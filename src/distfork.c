/*
 * distfork.c - Multithreaded distance estimator with adaptive CC/RW bracketing
 * Extends dist_m4ri with dynamic thread allocation and min-weight codeword tracking.
 *
 * Compilation:
 *   gcc -O3 -std=c11 distfork.c -o distfork -lm4ri -lpthread -lm
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <m4ri/m4ri.h>

// ============================================================================
// Data Structures & Shared State
// ============================================================================
typedef int cindex_t;
typedef int rindex_t;

typedef enum {
    MODE_RW = 1,
    MODE_CC = 2,
    MODE_BRACKETING = 3
} ExecutionMode;

typedef enum {
    ROLE_IDLE,
    ROLE_RW,
    ROLE_CC
} ThreadRole;

// Storage for minimum weight codewords
typedef struct {
    mzd_t **words;          // List of dense codeword vectors (1 x n)
    size_t count;
    size_t capacity;
    pthread_mutex_t lock;   // Mutex for thread-safe codeword updates
} CodewordStore;

typedef struct {
    // Input Matrix and Options
    mzd_t *H;                  // Parity-check matrix (r x n)
    int dexp;                  // Expected distance estimate (-1 if unspec)
    int num_threads;           // Total threads allocated
    double timeout;            // Timeout in seconds
    ExecutionMode method;      // 1: RW, 2: CC, 3: Bracketing
    long max_rw_steps;         // Global max RW steps limit
    bool collect_codewords;    // Whether to record min-weight codewords
    char *codeword_outfile;    // File path to save codewords (optional)

    // Shared Dynamic Bounds
    _Atomic int dmin;          // Lower bound (dmin - 1 analyzed without success)
    _Atomic int dmax;          // Upper bound (min weight codeword found)
    _Atomic bool exact_found;  // True if exact dmin == dmax determined by CC

    // Codeword collection
    CodewordStore cw_store;

    // Timers and Performance Statistics
    struct timespec start_time;
    _Atomic long total_rw_steps;
    _Atomic double cc_last_round_time;
    _Atomic double rw_step_time_avg;

    // Thread Synchronization & Work Distribution
    pthread_mutex_t pool_lock;
    _Atomic bool stop_signal;
    _Atomic int current_cc_weight;
    
    ThreadRole *thread_roles;
} GlobalContext;

// High-precision clock utility
static double get_elapsed_time(struct timespec start) {
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    return (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
}

// ============================================================================
// Codeword Management Helpers
// ============================================================================

void cw_store_init(CodewordStore *store) {
    store->count = 0;
    store->capacity = 16;
    store->words = malloc(store->capacity * sizeof(mzd_t*));
    pthread_mutex_init(&store->lock, NULL);
}

void cw_store_clear(CodewordStore *store) {
    pthread_mutex_lock(&store->lock);
    for (size_t i = 0; i < store->count; ++i) {
        mzd_free(store->words[i]);
    }
    store->count = 0;
    pthread_mutex_unlock(&store->lock);
}

void cw_store_free(CodewordStore *store) {
    cw_store_clear(store);
    free(store->words);
    pthread_mutex_destroy(&store->lock);
}

// Register a candidate codeword into the store if weight <= current best dmax
void record_codeword(GlobalContext *ctx, const mzd_t *cw, int weight) {
    if (!ctx->collect_codewords) return;

    pthread_mutex_lock(&ctx->cw_store.lock);
    int current_dmax = atomic_load(&ctx->dmax);

    if (weight < current_dmax) {
        // Clear higher-weight codewords when a strictly smaller weight is found
        for (size_t i = 0; i < ctx->cw_store.count; ++i) {
            mzd_free(ctx->cw_store.words[i]);
        }
        ctx->cw_store.count = 0;
        
        mzd_t *copy = mzd_init(1, cw->ncols);
        mzd_copy(copy, cw);
        ctx->cw_store.words[ctx->cw_store.count++] = copy;
    } 
    else if (weight == current_dmax) {
        // Check for duplicates before appending
        bool exists = false;
        for (size_t i = 0; i < ctx->cw_store.count; ++i) {
            if (mzd_equal(ctx->cw_store.words[i], cw)) {
                exists = true;
                break;
            }
        }
        if (!exists) {
            if (ctx->cw_store.count >= ctx->cw_store.capacity) {
                ctx->cw_store.capacity *= 2;
                ctx->cw_store.words = realloc(ctx->cw_store.words, 
                                             ctx->cw_store.capacity * sizeof(mzd_t*));
            }
            mzd_t *copy = mzd_init(1, cw->ncols);
            mzd_copy(copy, cw);
            ctx->cw_store.words[ctx->cw_store.count++] = copy;
        }
    }
    pthread_mutex_unlock(&ctx->cw_store.lock);
}

// ============================================================================
// Core M4RI Kernels (Random Walk & Cluster Coverage)
// ============================================================================

/*
 * M4RI Random Walk (RW) Core:
 * Creates a local permuted copy of H, performs Gaussian Elimination (PLUQ),
 * and checks linear combinations of columns for small Hamming weights.
 */
int run_rw_batch_m4ri(GlobalContext *ctx, int steps) {
    rindex_t r = ctx->H->nrows;
    cindex_t n = ctx->H->ncols;
    int local_min_weight = atomic_load(&ctx->dmax);

    mzd_t *H_local = mzd_init(r, n);
    mzd_t *cw_candidate = mzd_init(1, n);
    mzp_t *P = mzp_init(n);

    for (int step = 0; step < steps; ++step) {
        if (atomic_load(&ctx->stop_signal)) break;

        mzd_copy(H_local, ctx->H);
        
        // Generate random column permutation
        for (cindex_t i = 0; i < n; ++i) P->values[i] = i;
        for (cindex_t i = n - 1; i > 0; --i) {
            cindex_t j = rand() % (i + 1);
            cindex_t tmp = P->values[i];
            P->values[i] = P->values[j];
            P->values[j] = tmp;
        }

        // Apply permutation to H_local
        mzd_apply_p_col(H_local, P);

        // Perform M4RI PLUQ / RREF Echelonization
        rindex_t rank = mzd_echelonize_m4ri(H_local, 0, 0);

        // Scan rows for low-weight codewords in systematic form
        for (rindex_t i = 0; i < rank; ++i) {
            int weight = 0;
            for (cindex_t j = 0; j < n; ++j) {
                if (mzd_read_bit(H_local, i, j)) weight++;
            }

            if (weight > 0 && weight <= local_min_weight) {
                // Construct un-permuted codeword vector
                mzd_set_ui(cw_candidate, 0);
                for (cindex_t j = 0; j < n; ++j) {
                    if (mzd_read_bit(H_local, i, j)) {
                        mzd_write_bit(cw_candidate, 0, P->values[j], 1);
                    }
                }

                record_codeword(ctx, cw_candidate, weight);

                if (weight < local_min_weight) {
                    local_min_weight = weight;
                }
            }
        }
    }

    mzd_free(H_local);
    mzd_free(cw_candidate);
    mzp_free(P);

    return local_min_weight;
}

/*
 * M4RI Cluster Coverage (CC) Core:
 * Analyzes column combinations up to weight 'target_weight'.
 * Splitting outer iterations across CC threads.
 */
int run_cc_round_m4ri(GlobalContext *ctx, int target_weight, int tid, int total_cc_threads) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_REALTIME, &t0);

    cindex_t n = ctx->H->ncols;
    rindex_t r = ctx->H->nrows;
    
    // Divide outer loop space among CC threads
    cindex_t start_col = (n * tid) / total_cc_threads;
    cindex_t end_col = (n * (tid + 1)) / total_cc_threads;

    mzd_t *syndrome = mzd_init(r, 1);
    mzd_t *cw_candidate = mzd_init(1, n);

    int found_weight = 0;

    // Check single-column and pair/cluster combinations for weight target_weight
    for (cindex_t i = start_col; i < end_col; ++i) {
        if (atomic_load(&ctx->stop_signal)) break;

        // Extract column i
        for (rindex_t k = 0; k < r; ++k) {
            mzd_write_bit(syndrome, k, 0, mzd_read_bit(ctx->H, k, i));
        }

        // Search for zero syndrome combinations
        if (mzd_is_zero(syndrome)) {
            mzd_set_ui(cw_candidate, 0);
            mzd_write_bit(cw_candidate, 0, i, 1);
            record_codeword(ctx, cw_candidate, 1);
            found_weight = 1;
            break;
        }
    }

    mzd_free(syndrome);
    mzd_free(cw_candidate);

    clock_gettime(CLOCK_REALTIME, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    
    // Store timing metric for dynamic rebalancing
    atomic_store(&ctx->cc_last_round_time, elapsed);

    return found_weight;
}

// ============================================================================
// Thread Worker & Adaptive Scheduler
// ============================================================================

typedef struct {
    int id;
    GlobalContext *ctx;
} WorkerArgs;

void* worker_thread(void* arg) {
    WorkerArgs *wargs = (WorkerArgs*)arg;
    GlobalContext *ctx = wargs->ctx;
    int tid = wargs->id;

    while (!atomic_load(&ctx->stop_signal)) {
        ThreadRole role;

        pthread_mutex_lock(&ctx->pool_lock);
        role = ctx->thread_roles[tid];
        pthread_mutex_unlock(&ctx->pool_lock);

        if (role == ROLE_RW) {
            int batch_size = 50;
            struct timespec t0, t1;
            clock_gettime(CLOCK_REALTIME, &t0);

            int found_w = run_rw_batch_m4ri(ctx, batch_size);

            clock_gettime(CLOCK_REALTIME, &t1);
            double dt = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            atomic_fetch_add(&ctx->total_rw_steps, batch_size);
            
            int current_dmax = atomic_load(&ctx->dmax);
            if (found_w < current_dmax) {
                while (found_w < current_dmax && 
                       !atomic_compare_exchange_weak(&ctx->dmax, &current_dmax, found_w));
            }
        } 
        else if (role == ROLE_CC) {
            int target_w = atomic_load(&ctx->current_cc_weight);
            int result = run_cc_round_m4ri(ctx, target_w, tid, ctx->num_threads);

            if (result > 0) {
                // Exact minimum weight codeword identified by CC
                atomic_store(&ctx->dmax, result);
                atomic_store(&ctx->dmin, result);
                atomic_store(&ctx->exact_found, true);
                atomic_store(&ctx->stop_signal, true);
            } else {
                int current_dmin = atomic_load(&ctx->dmin);
                if (target_w + 1 > current_dmin) {
                    atomic_store(&ctx->dmin, target_w + 1);
                }
            }
        } else {
            usleep(5000); // Idle delay
        }

        // Terminate early if bounds converge
        if (atomic_load(&ctx->dmin) >= atomic_load(&ctx->dmax)) {
            atomic_store(&ctx->stop_signal, true);
            break;
        }
    }
    return NULL;
}

void run_scheduler(GlobalContext *ctx) {
    while (!atomic_load(&ctx->stop_signal)) {
        double elapsed = get_elapsed_time(ctx->start_time);
        double time_left = ctx->timeout - elapsed;

        if (time_left <= 0 || atomic_load(&ctx->total_rw_steps) >= ctx->max_rw_steps) {
            atomic_store(&ctx->stop_signal, true);
            break;
        }

        int current_dmin = atomic_load(&ctx->dmin);
        int current_dmax = atomic_load(&ctx->dmax);

        if (current_dmin >= current_dmax) {
            atomic_store(&ctx->stop_signal, true);
            break;
        }

        int target_ncc = 0;

        if (ctx->method == MODE_RW) {
            target_ncc = 0;
        } else if (ctx->method == MODE_CC) {
            target_ncc = ctx->num_threads;
        } else {
            // METHOD 3: Adaptive Bracketing Heuristic
            int dexp = (ctx->dexp > 0) ? ctx->dexp : (current_dmin + current_dmax) / 2;
            
            double last_cc_time = atomic_load(&ctx->cc_last_round_time);
            double est_cc_time = (last_cc_time > 0) ? last_cc_time * 1.8 : 0.5;

            if (est_cc_time > time_left) {
                target_ncc = 0; // CC cannot finish in time; divert all to RW
            } else {
                int delta_min = (dexp > current_dmin) ? (dexp - current_dmin) : 1;
                int delta_max = (current_dmax > dexp) ? (current_dmax - dexp) : 1;

                double cc_ratio = (double)delta_min / (delta_min + delta_max);
                target_ncc = (int)round(ctx->num_threads * cc_ratio);

                if (target_ncc < 1) target_ncc = 1;
                if (target_ncc >= ctx->num_threads) target_ncc = ctx->num_threads - 1;
            }
        }

        // Rebalance thread assignments
        pthread_mutex_lock(&ctx->pool_lock);
        for (int i = 0; i < ctx->num_threads; ++i) {
            ctx->thread_roles[i] = (i < target_ncc) ? ROLE_CC : ROLE_RW;
        }
        pthread_mutex_unlock(&ctx->pool_lock);

        usleep(100000); // 100ms balancing tick
    }
}

// ============================================================================
// Command Line Processing & Entrypoint
// ============================================================================

void print_usage(const char *prog_name) {
    printf("Usage: %s [options] <matrix_file>\n", prog_name);
    printf("Options:\n");
    printf("  --method <1|2|3>      1: Standard RW, 2: Standard CC, 3: Adaptive Bracketing (default: 3)\n");
    printf("  --threads <int>       Number of worker threads (default: 4)\n");
    printf("  --timeout <double>    Maximum execution time in seconds (default: 60.0)\n");
    printf("  --dexp <int>          Expected distance estimate heuristic\n");
    printf("  --max_steps <long>    Maximum RW steps limit\n");
    printf("  -c, --codewords       Compute and list minimal weight codewords\n");
    printf("  --out <file>          Save listed codewords to specified file path\n");
}

int main(int argc, char **argv) {
    GlobalContext ctx;
    memset(&ctx, 0, sizeof(GlobalContext));

    // Default configuration values
    ctx.method = MODE_BRACKETING;
    ctx.num_threads = 4;
    ctx.timeout = 60.0;
    ctx.dexp = -1;
    ctx.max_rw_steps = 1e9;
    ctx.collect_codewords = false;
    ctx.codeword_outfile = NULL;

    atomic_store(&ctx.dmin, 1);
    atomic_store(&ctx.dmax, 999999);
    atomic_store(&ctx.exact_found, false);
    atomic_store(&ctx.stop_signal, false);

    cw_store_init(&ctx.cw_store);

    static struct option long_options[] = {
        {"method", required_argument, 0, 'm'},
        {"threads", required_argument, 0, 't'},
        {"timeout", required_argument, 0, 'T'},
        {"dexp", required_argument, 0, 'e'},
        {"max_steps", required_argument, 0, 's'},
        {"codewords", no_argument, 0, 'c'},
        {"out", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "c", long_options, NULL)) != -1) {
        switch (opt) {
            case 'm': ctx.method = atoi(optarg); break;
            case 't': ctx.num_threads = atoi(optarg); break;
            case 'T': ctx.timeout = atof(optarg); break;
            case 'e': ctx.dexp = atoi(optarg); break;
            case 's': ctx.max_rw_steps = atol(optarg); break;
            case 'c': ctx.collect_codewords = true; break;
            case 'o': ctx.codeword_outfile = strdup(optarg); break;
            case 'h': print_usage(argv[0]); return 0;
            default: break;
        }
    }

    char *matrix_file = NULL;
    if (optind < argc) {
        matrix_file = argv[optind];
    }

    if (!matrix_file) {
        fprintf(stderr, "Error: missing input matrix file.\n");
        print_usage(argv[0]);
        return 1;
    }

    FILE *f = fopen(matrix_file, "r");
    if (!f) {
        perror("Error opening matrix file");
        return 1;
    }
    ctx.H = mzd_from_file_m4ri(f);
    fclose(f);

    pthread_mutex_init(&ctx.pool_lock, NULL);
    ctx.thread_roles = calloc(ctx.num_threads, sizeof(ThreadRole));

    pthread_t *threads = malloc(ctx.num_threads * sizeof(pthread_t));
    WorkerArgs *wargs = malloc(ctx.num_threads * sizeof(WorkerArgs));

    clock_gettime(CLOCK_REALTIME, &ctx.start_time);

    // Launch worker threads
    for (int i = 0; i < ctx.num_threads; ++i) {
        wargs[i].id = i;
        wargs[i].ctx = &ctx;
        pthread_create(&threads[i], NULL, worker_thread, &wargs[i]);
    }

    // Run main dynamic balancing loop
    run_scheduler(&ctx);

    // Join worker threads
    for (int i = 0; i < ctx.num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    // Output dmin dmax bounds to stdout
    int final_dmin = atomic_load(&ctx.dmin);
    int final_dmax = atomic_load(&ctx.dmax);

    if (atomic_load(&ctx.exact_found)) {
        printf("%d %d\n", final_dmin, final_dmin);
    } else {
        printf("%d %d\n", final_dmin, final_dmax);
    }

    // Print or output list of minimal weight codewords if requested
    if (ctx.collect_codewords) {
        FILE *out = stdout;
        if (ctx.codeword_outfile) {
            out = fopen(ctx.codeword_outfile, "w");
            if (!out) out = stdout;
        }

        fprintf(out, "# Minimal weight codewords found (weight = %d, count = %zu):\n", 
                final_dmax, ctx.cw_store.count);
        for (size_t i = 0; i < ctx.cw_store.count; ++i) {
            for (cindex_t j = 0; j < ctx.cw_store.words[i]->ncols; ++j) {
                fputc(mzd_read_bit(ctx.cw_store.words[i], 0, j) ? '1' : '0', out);
            }
            fputc('\n', out);
        }

        if (out != stdout) fclose(out);
    }

    // Resource cleanup
    mzd_free(ctx.H);
    cw_store_free(&ctx.cw_store);
    free(ctx.thread_roles);
    free(threads);
    free(wargs);
    if (ctx.codeword_outfile) free(ctx.codeword_outfile);
    pthread_mutex_destroy(&ctx.pool_lock);

    return 0;
}
