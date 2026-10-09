/** ********************************************************************** 
 * @file dist_m4ri.c
 * @brief Multithreaded distance calculation and bracketing (dist_m4ri)
 * 
 * The program implements multithreaded distance calculation:
 * - method=1: Multithreaded Random Window (RW) algorithm (upper bound)
 * - method=2: Multithreaded Connected Cluster (CC) algorithm (lower bound / exact)
 * - method=3: Bracketing mode (artillery fork / вилка) dynamically balancing
 *             CC and RW threads based on distance estimate (dexp/dest),
 *             current bounds [dmin, dmax], RW step count, timeout,
 *             and scaling characteristics.
 *
 * Output to stdout: "dmin dmax rw_steps"
 * where dmin-1 is the maximum cluster size analyzed without success by CC,
 * dmin=dmax if CC actually found a min-weight codeword of this size,
 * dmax is the smallest-weight codeword found by RW,
 * and rw_steps is the number of completed RW steps (0 if CC found a min-weight codeword
 * or if RW did not run in method=2).
 * NOTE: This 3-number output format is incompatible with legacy single-threaded dist_m4ri_old.
 *
 * Expert CC options (see `--morehelp`) restrict the CC search:
 * - noscan=1 (method=2): a single CC round at w=wmax.  Unless the supplied dmin equals wmax,
 *   lower weights are not scanned, so a codeword found only sets dmax and dmin is not raised.
 * - start=a,b,c: clusters are grown only from the listed columns (without the usual
 *   restriction to larger column indices); dmin is valid only under a code symmetry which
 *   maps every minimum-weight codeword to one containing a listed column.
 * - cbeg/cend: split runs; dmin only covers codewords whose smallest column is in [cbeg,cend].
 * A CC codeword gives dmin=dmax only if found in a round where all lower weights have been
 * analyzed (flag `cc_exact`).
 *
 * Threads: every worker can run CC rounds, while RW runs only on workers 0..rw_threads-1.  Thread
 * throttling (unless nothrottle=1) limits the RW workers by the memory of their dense matrices and by
 * the step count; in method=3 this limits only the RW share, and CC rounds can use all workers.
 * A worker without RW work (RW steps all claimed, RW stopped by min_hits) joins the current CC round.
 *
 * Timing model (method=2, 3): the CC work of a round is measured as the total CC thread time, and the work of the
 * next round is extrapolated with the growth factor of the last two rounds.  With a timeout, a CC round predicted not
 * to finish in time on all threads is not started in method=3 while RW runs (once RW has ended, the run ends); in
 * method=2 (and method=3 with steps=0), it is started anyway, as it can still find a codeword of weight w = dmin.  In
 * method=3, the RW step time is measured continuously; the round w = dmax-1 (which certifies dmin = dmax) gets all
 * threads, and while RW has not found any codeword, CC rounds w > dexp run only on the threads which cannot run RW.
 * The coordinator never ends CC while RW still runs.
 *
 * All debugging messages and confinement profile are sent to stderr.
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>
 ************************************************************************/

#define _GNU_SOURCE
#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <limits.h>
#include <m4ri/m4ri.h>

#include "mmio.h"
#include "uthash.h"
#include "util_hash.h"
#include "util_m4ri.h"
#include "util_io.h"
#include "dist_m4ri.h"
#include "dist_cc.h"

/* Mutex protecting M4RI's internal non-thread-safe MMC memory cache */
static pthread_mutex_t m4ri_mem_mutex = PTHREAD_MUTEX_INITIALIZER;

static inline mzd_t *safe_mzd_from_csr(mzd_t *dst, const csr_t *p) {
  pthread_mutex_lock(&m4ri_mem_mutex);
  mzd_t *res = mzd_from_csr(dst, p);
  pthread_mutex_unlock(&m4ri_mem_mutex);
  return res;
}

static inline mzd_t *safe_mzd_init(rci_t r, rci_t c) {
  pthread_mutex_lock(&m4ri_mem_mutex);
  mzd_t *res = mzd_init(r, c);
  pthread_mutex_unlock(&m4ri_mem_mutex);
  return res;
}

static inline void safe_mzd_free(mzd_t *M) {
  if (!M) return;
  pthread_mutex_lock(&m4ri_mem_mutex);
  mzd_free(M);
  pthread_mutex_unlock(&m4ri_mem_mutex);
}

static inline mzp_t *safe_mzp_init(rci_t length) {
  pthread_mutex_lock(&m4ri_mem_mutex);
  mzp_t *res = mzp_init(length);
  pthread_mutex_unlock(&m4ri_mem_mutex);
  return res;
}

static inline void safe_mzp_free(mzp_t *P) {
  if (!P) return;
  pthread_mutex_lock(&m4ri_mem_mutex);
  mzp_free(P);
  pthread_mutex_unlock(&m4ri_mem_mutex);
}

static inline double get_time_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static inline uint64_t splitmix64(uint64_t *state) {
  uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static inline int rand_uniform_thread(int max, uint64_t *state) {
  if (max <= 1) return 0;
  return (int)(splitmix64(state) % (uint64_t)max);
}

static inline mzp_t * mzp_rand_thread(mzp_t *q, rci_t length, uint64_t *state) {
  if (q == NULL) return NULL;
  for (int i = 0; i <= (int)length - 2; i++) {
    q->values[i] = i + rand_uniform_thread(length - i, state);
  }
  for (int i = length - 1; i < (int)q->length; i++) {
    q->values[i] = i;
  }
  return q;
}

typedef struct {
  params_t *p;
  int num_threads;
  int rw_threads; /* workers 0..rw_threads-1 may run RW (dense RW matrices); all workers may run CC */
  double timeout;
  double start_time;
  int dexp;

  /* Distance bounds & stop flags (cache-line isolated for read-mostly access) */
  _Alignas(64) atomic_int dmin; /* dmin-1 is max cluster size analyzed without success */
  atomic_int dmax;              /* smallest weight codeword found (0 if none) */
  atomic_int cc_found_weight;   /* smallest weight of a codeword found by CC (0 if none) */
  atomic_bool cc_exact;         /* cc_found_weight is certified exact (found in a round w=dmin) */
  atomic_bool stop_flag;        /* signals all threads to terminate */
  atomic_bool rw_stop_flag;     /* signals RW workers to stop (in method 3) */

  /* RW state (cache-line isolated from CC and read-mostly bounds) */
  _Alignas(64) long total_rw_steps;
  atomic_long rw_steps_started;
  atomic_long rw_steps_completed;
  atomic_llong rw_busy_ns;           /* thread time in ns of the timed RW steps, for the method-3 timing model */
  atomic_long rw_timed_steps;        /* number of RW steps timed in `rw_busy_ns` */

  /* CC state for current weight (cache-line isolated) */
  _Alignas(64) atomic_llong cc_next; /* weight of the current CC round and next start index, see `cc_pack()` */
  csr_t *mHT_cc;
  int max_col_W;
  atomic_int cc_active_workers;
  atomic_int cc_target_workers;
  atomic_int cc_round_active;
  atomic_llong cc_busy_ns;           /* total CC thread time in ns (all rounds), for the CC timing model */

  /* Codeword synchronization */
  pthread_mutex_t cw_mutex;

  /* Subspace sketching (ksub) state */
  mzd_t *N_global;
  int nu;
  pthread_rwlock_t basis_rwlock;
  atomic_long next_refresh_step;

  /* Timing stats: CC work of each completed round in thread-seconds (total CC thread time) */
  double cc_time_per_weight[MAX_W];

  /* Thread handles */
  pthread_t *threads;
} distfork_ctx_t;

typedef struct {
  distfork_ctx_t *ctx;
  int tid;
  int min_swei[MAX_W];
} worker_arg_t;

/* A CC round at weight w hands out start indices i = beg..end one at a time.  The weight and the next start index
 * are packed into the single atomic word `cc_next`, so that a worker always processes a claimed start index with
 * the weight of the round it belongs to, even if it joins just as the coordinator publishes the next round. */
#define CC_IDX_BITS 32
static inline long long cc_pack(const int w, const int idx) {
  return ((long long)w << CC_IDX_BITS) | (long long)(unsigned int)idx;
}
static inline int cc_unpack_w(const long long v) { return (int)(v >> CC_IDX_BITS); }
static inline int cc_unpack_idx(const long long v) { return (int)(v & 0xffffffffLL); }

/* First start index of a CC round: index into the expert `start` list, or the first start column (cbeg) */
static inline int cc_round_beg(const params_t * const p) {
  return (p->start_num > 0) ? 0 : ((p->cbeg >= 0) ? p->cbeg : 0);
}

/* Last start index of a CC round at weight w: a cluster started at column i only uses columns > i, so it
 * starts at column n-w at the latest (expert `start` list: all listed columns, clusters are unlimited) */
static inline int cc_round_end(const params_t * const p, const int w) {
  if (p->start_num > 0) return p->start_num - 1;
  const int nvar = p->spaH->cols;
  return (p->cend >= 0) ? minint(p->cend, nvar - w) : (nvar - w);
}

/* Recursive CC worker function (interruptible) */
static int start_CC_recurs_mt(one_vec_t *err, one_vec_t *urr, one_vec_t * const syn[],
                              const int w_limit, const int max_col_wt,
                              const csr_t * const mH, const csr_t * const mHT,
                              worker_arg_t *warg) {
  distfork_ctx_t *ctx = warg->ctx;
  if (atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed)) {
    return 0;
  }
  params_t * const p = ctx->p;
  const int w = err->wei;
  int current_limit = w_limit;
  int cur_dmax = atomic_load_explicit(&ctx->dmax, memory_order_relaxed);
  if (cur_dmax > 0 && p->dW >= 0) {
    current_limit = minint(w_limit, cur_dmax + p->dW);
  }
  if (w >= current_limit) {
    return 0;
  }

  const one_vec_t * const syn_w = syn[w];
  const int syn_w_wei = syn_w->wei;
  const int row = syn_w->vec[0];
  const csr_t * const mL = p->spaL;
  /* Clusters are grown only to columns larger than the start column, except with the expert `start`
   * list, where clusters are unlimited (columns already in `err` are skipped via `one_ordered_search`). */
  const int col_min = (p->start_num > 0) ? -1 : urr->vec[0];

  /* Leaf level: w + 1 == current_limit */
  if (w + 1 == current_limit) {
    int max_leaf_swei = 0;
    if (p->smax > 0 && current_limit < MAX_W) {
      int cur_min = warg->min_swei[current_limit] - 1;
      max_leaf_swei = (cur_min < p->smax) ? (cur_min > 0 ? cur_min : 0) : p->smax;
    }

    for (int i1 = mH->p[row]; i1 < mH->p[row + 1]; i1++) {
      const int col = mH->i[i1];
      if (col <= col_min) continue;

      const int p_beg = mHT->p[col];
      const int col_wt = mHT->p[col + 1] - p_beg;
      if (abs(syn_w_wei - col_wt) > max_leaf_swei) continue;

      int swei;
      if (max_leaf_swei == 0) {
        /* Only swei == 0 is of interest: check if mHT[col] == syn[w] */
        if (mHT->i[p_beg] != row ||
            mHT->i[p_beg + col_wt - 1] != syn_w->vec[col_wt - 1]) {
          continue;
        }
        if (memcmp(syn_w->vec, &mHT->i[p_beg], (size_t)col_wt * sizeof(int)) != 0) {
          continue;
        }
        if (one_ordered_search(err, col) != -1) continue;
        swei = 0;
      } else {
        if (one_ordered_search(err, col) != -1) continue;
        syn[w + 1]->wei = 0;
        swei = one_csr_row_combine(syn[w + 1], syn_w, mHT, col);
        if (swei > 0) {
          if (swei <= max_leaf_swei) {
            warg->min_swei[w + 1] = swei;
            int cur_min = swei - 1;
            max_leaf_swei = (cur_min < p->smax) ? (cur_min > 0 ? cur_min : 0) : p->smax;
          }
          continue;
        }
      }

      /* swei == 0: insert col into err to verify against mL and record codeword */
      int pos = one_ordered_ins(err, col);
      int nz = (!mL) || sparse_syndrome_non_zero(mL, err->wei, err->vec);
      if (nz) {
        bool stop = false;
        pthread_mutex_lock(&ctx->cw_mutex);
        p->codewords = codeword_add_maybe(p, err->vec, err->wei);
        int cur_d = atomic_load(&ctx->dmax);
        if (p->min_w < cur_d || cur_d == 0) {
          atomic_store(&ctx->dmax, p->min_w);
        }
        int cur_cc_found = atomic_load(&ctx->cc_found_weight);
        if (cur_cc_found == 0 || err->wei < cur_cc_found) {
          atomic_store(&ctx->cc_found_weight, err->wei);
        }
        if (!p->outC && p->maxC == 0) {
          atomic_store(&ctx->stop_flag, true);
          stop = true;
        }
        if (p->maxC && p->num_cws >= p->maxC) {
          atomic_store(&ctx->stop_flag, true);
          stop = true;
        }
        pthread_mutex_unlock(&ctx->cw_mutex);

        if (stop) {
          one_ordered_pos_del(err, col, pos);
          return 1;
        }
      }
      one_ordered_pos_del(err, col, pos);
    }
    return 0;
  }

  /* Internal level: w + 1 < current_limit */
  const int rem = current_limit - (w + 1);
  int max_s_needed = 0;
  if (p->smax > 0) {
    if (p->noscan || p->dmin > 1) {
      max_s_needed = p->smax;
    } else if (current_limit < MAX_W) {
      int cur_min = warg->min_swei[current_limit] - 1;
      max_s_needed = (cur_min < p->smax) ? (cur_min > 0 ? cur_min : 0) : p->smax;
    }
  }
  const int max_reach = rem * max_col_wt + max_s_needed;

  for (int i1 = mH->p[row]; i1 < mH->p[row + 1]; i1++) {
    const int col = mH->i[i1];
    if (col > col_min) {
      const int col_wt = mHT->p[col + 1] - mHT->p[col];
      if (syn_w_wei - col_wt > max_reach) continue;

      int pos = one_ordered_search(err, col);
      if (pos == -1) {
        syn[w + 1]->wei = 0;
        int swei = one_csr_row_combine(syn[w + 1], syn_w, mHT, col);

        if (p->smax && swei > 0 && swei <= p->smax && (w + 1 < MAX_W)) {
          if (swei < warg->min_swei[w + 1]) {
            warg->min_swei[w + 1] = swei;
          }
        }

        if (swei > 0 && swei <= max_reach) {
          urr->vec[w] = col;
          urr->wei++;
          pos = one_ordered_ins(err, col);
          int result = start_CC_recurs_mt(err, urr, syn, w_limit, max_col_wt,
                                          mH, mHT, warg);
          urr->wei--;
          one_ordered_pos_del(err, col, pos);
          if (result == 1) {
            return 1;
          }
        }
      }
    }
  }
  return 0;
}

/* Add the thread time `*acc_t` of `*acc_n` completed RW steps to the totals of the RW timing model, and reset them */
static inline void rw_time_publish(distfork_ctx_t *ctx, double *acc_t, int *acc_n) {
  if (*acc_n > 0) {
    atomic_fetch_add_explicit(&ctx->rw_busy_ns, (long long)(*acc_t * 1e9), memory_order_relaxed);
    atomic_fetch_add_explicit(&ctx->rw_timed_steps, *acc_n, memory_order_relaxed);
  }
  *acc_t = 0.0;
  *acc_n = 0;
}

/* Time of a completed RW step since `*t_last`: published at least every 2 ms (after every step if steps are slow),
 * and right away for the first step of a batch (`first`), so that the coordinator gets the step time early */
static inline void rw_time_step(distfork_ctx_t *ctx, double *t_last, double *acc_t, int *acc_n, const bool first) {
  const double t_now = get_time_sec();
  *acc_t += t_now - *t_last;
  *t_last = t_now;
  (*acc_n)++;
  if (first || *acc_t >= 0.002) rw_time_publish(ctx, acc_t, acc_n);
}

/* Run RW batch; returns the number of completed steps */
static int run_rw_steps(distfork_ctx_t *ctx, int n_steps,
                        mzd_t *mH, mzd_t *mHT, rci_t *ee,
                        mzp_t *perm, mzp_t *pivs, word *piv_mask,
                        int *eff_nrows_ptr,
                        int *visited_cols, int *visited_checks, int *col_queue,
                        int *visit_marker, uint64_t *rng_state, int tid) {
  params_t * const p = ctx->p;
  const csr_t * const spaL0 = p->spaL;
  const int nvar = p->spaH->cols;
  const int classical = p->classical;
  const int kwin = p->kwin;
  const int win_mode = p->win_mode;
  int eff_nrows = *eff_nrows_ptr;
  int n_done = 0;
  double t_last = get_time_sec(), acc_t = 0.0;
  int acc_n = 0;

  for (int step = 0; step < n_steps; step++) {
    if (atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed) ||
        atomic_load_explicit(&ctx->rw_stop_flag, memory_order_relaxed)) {
      break;
    }

    int eff_kwin = kwin;
    int eff_win_mode = win_mode;
    if (eff_kwin == 0 && nvar >= 500 && (step & 1) == 1) {
      eff_kwin = minint(512, (nvar * 3) / 4);
      eff_win_mode = 1;
    }

    if (eff_kwin > 0 && eff_kwin < nvar) {
      int seed_col = rand_uniform_thread(nvar, rng_state);
      localized_window_perm(perm, nvar, seed_col, eff_kwin, eff_win_mode,
                            p->spaH, ctx->mHT_cc,
                            visited_cols, visited_checks, col_queue,
                            visit_marker, rng_state);
    } else {
      pivs = mzp_rand_thread(pivs, nvar, rng_state);
      mzp_set_ui(perm, 1);
      perm = perm_p_trans(perm, pivs, 0);
    }

    memset(piv_mask, 0, mH->width * sizeof(word));
    int rank = 0;
    bool interrupted = false;
    for (int i = 0; i < nvar && rank < eff_nrows; i++) {
      /* a step on a large matrix can take many seconds: check the stop flags also every 16 columns */
      if ((i & 15) == 15 && (atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed) ||
                             atomic_load_explicit(&ctx->rw_stop_flag, memory_order_relaxed))) {
        interrupted = true;
        break;
      }
      int col = perm->values[i];
      if (gauss_one_rows(mH, col, rank, eff_nrows)) {
        pivs->values[rank++] = col;
        piv_mask[col >> 6] |= (word)1 << (col & 63);
      }
    }
    /* abandon the step: `eff_nrows` is kept, as rows rank..eff_nrows-1 are not processed (same row space) */
    if (interrupted) break;
    eff_nrows = rank;
    *eff_nrows_ptr = eff_nrows;

    mzd_transpose(mHT, mH);

    const int active_width = (rank + 63) >> 6;
    for (int col = 0; col < nvar; col++) {
      if ((piv_mask[col >> 6] >> (col & 63)) & 1) continue;
      int cnt = 0;
      ee[cnt++] = col;
      int limit = nvar + 1;
      int cur_dmax = atomic_load_explicit(&ctx->dmax, memory_order_relaxed);
      if (cur_dmax > 0) {
        if ((p->outC || p->maxC || p->dW > 0 || p->min_hits > 0 || (p->debug & 1)) &&
            p->dW >= 0) {
          limit = minint(limit, cur_dmax + p->dW + 1);
        } else {
          limit = minint(limit, cur_dmax);
        }
      }

      word *rawrow = mzd_row(mHT, col);
      rci_t j = -1;
      while (cnt < limit) {
        j = nextelement(rawrow, active_width, j);
        if (j == -1 || j >= rank) break;
        ee[cnt++] = pivs->values[j++];
      }

      if (cnt < limit) {
        rci_quick_sort(ee, cnt);
        int nz = classical ? 1 : sparse_syndrome_non_zero(spaL0, cnt, ee);
        if (nz) {
          pthread_mutex_lock(&ctx->cw_mutex);
          p->codewords = codeword_add_maybe(p, ee, cnt);
          if (cnt < p->min_w) p->min_w = cnt;
          int best = p->min_w;
          int old_dmax = atomic_load(&ctx->dmax);
          if (old_dmax == 0 || best < old_dmax) {
            atomic_store(&ctx->dmax, best);
            if (p->debug & 16) {
              int num_rw = (ctx->p->method == 1) ? ctx->num_threads
                           : minint(ctx->rw_threads, ctx->num_threads - atomic_load(&ctx->cc_target_workers));
              if (num_rw < 1) num_rw = 1;
              fprintf(stderr,
                      "# [thread %d] RW found new upper bound cw of weight %d "
                      "(using %d RW threads)\n", tid, best, num_rw);
            }
            int cur_dmin = atomic_load(&ctx->dmin);
            if (cur_dmin > 0 && best <= cur_dmin) {
              atomic_store(&ctx->stop_flag, true);
            }
          }
          if (p->wmin > 0 && best <= p->wmin) {
            atomic_store(&ctx->stop_flag, true);
          }
          if (p->maxC && p->num_cws >= p->maxC) {
            atomic_store(&ctx->stop_flag, true);
          }
          if (check_min_hits_convergence(p)) {
            int req = (p->cov_cws > 0) ? minint(p->cov_cws, 5) : 1;
            if ((p->debug & 1) && !atomic_load(&ctx->stop_flag) &&
                !atomic_load(&ctx->rw_stop_flag)) {
              int c_min = 0, c_max = 0;
              double c_avg = 0.0, c_std = 0.0;
              compute_min_w_hit_stats(p, &c_min, &c_max, &c_avg, &c_std);
              fprintf(stderr,
                      "# RW convergence reached: >= %d min-weight (w=%d) cws "
                      "(cov=%d/%d) each hit >= %d times "
                      "(hits min=%d, max=%d, avg=%.2f, stdev=%.2f)\n",
                      req, best, p->min_w_cov, p->min_w_cws, p->min_hits,
                      c_min, c_max, c_avg, c_std);
            }
            if (p->method == 1) {
              atomic_store(&ctx->stop_flag, true);
            } else {
              atomic_store(&ctx->rw_stop_flag, true);
            }
          }
          pthread_mutex_unlock(&ctx->cw_mutex);
        }
      }
    }
    atomic_fetch_add(&ctx->rw_steps_completed, 1);
    n_done++;
    rw_time_step(ctx, &t_last, &acc_t, &acc_n, n_done == 1);
  }
  rw_time_publish(ctx, &acc_t, &acc_n);
  return n_done;
}

/* Run compact subspace RW batch (ksub > 0); returns the number of completed steps */
static int run_rw_steps_ksub(distfork_ctx_t *ctx, int n_steps,
                             mzd_t *M_sub, rci_t *ee, mzp_t *perm, mzp_t *pivs,
                             int *visited_cols, int *visited_checks, int *col_queue,
                             int *visit_marker, uint64_t *rng_state, int tid) {
  params_t * const p = ctx->p;
  const csr_t * const spaL0 = p->spaL;
  const int nvar = p->spaH->cols;
  const int classical = p->classical;
  const int nu = ctx->nu;
  const int ksub = M_sub->nrows;
  const int kwin = p->kwin;
  const int win_mode = p->win_mode;

  if (nu <= 0 || ksub <= 0 || !ctx->N_global) { /* not reached: no subspace to sample, count the steps as done */
    atomic_fetch_add(&ctx->rw_steps_completed, n_steps);
    return 0;
  }

  int n_done = 0;
  double t_last = get_time_sec(), acc_t = 0.0;
  int acc_n = 0;
  pthread_rwlock_rdlock(&ctx->basis_rwlock);
  const mzd_t * const N = ctx->N_global;

  for (int step = 0; step < n_steps; step++) {
    if (atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed) ||
        atomic_load_explicit(&ctx->rw_stop_flag, memory_order_relaxed)) {
      break;
    }

    /* 1. Generate column permutation (localized window or uniform) */
    int eff_kwin = kwin;
    int eff_win_mode = win_mode;
    if (eff_kwin == 0 && nvar >= 500 && (step & 1) == 1) {
      eff_kwin = minint(512, (nvar * 3) / 4);
      eff_win_mode = 1;
    }

    if (eff_kwin > 0 && eff_kwin < nvar) {
      int seed_col = rand_uniform_thread(nvar, rng_state);
      localized_window_perm(perm, nvar, seed_col, eff_kwin, eff_win_mode,
                            p->spaH, ctx->mHT_cc,
                            visited_cols, visited_checks, col_queue,
                            visit_marker, rng_state);
    } else {
      pivs = mzp_rand_thread(pivs, nvar, rng_state);
      mzp_set_ui(perm, 1);
      perm = perm_p_trans(perm, pivs, 0);
    }

    /* 2. Sample ksub rows from N into M_sub (with optional window overlap preference) */
    const int marker = *visit_marker;
    for (int i = 0; i < ksub; i++) {
      int r = rand_uniform_thread(nu, rng_state);
      if (eff_kwin > 0 && eff_kwin < nvar && visited_cols) {
        for (int attempt = 0; attempt < 4; attempt++) {
          const word *raw_n = mzd_row_cons(N, r);
          int j_bit = nextelement(raw_n, N->width, -1);
          int overlaps = 0;
          while (j_bit >= 0 && j_bit < nvar) {
            if (visited_cols[j_bit] == marker) {
              overlaps = 1;
              break;
            }
            j_bit = nextelement(raw_n, N->width, j_bit + 1);
          }
          if (overlaps) break;
          r = rand_uniform_thread(nu, rng_state);
        }
      }
      mzd_copy_row(M_sub, i, N, r);
    }

    /* 3. Echelonize M_sub directly in L1/L2 cache using permuted column order */
    int rank = 0;
    for (int i = 0; i < nvar && rank < ksub; i++) {
      int col = perm->values[i];
      if (gauss_one(M_sub, col, rank)) {
        rank++;
      }
    }

    /* 4. Extract candidate codewords directly from reduced rows (already in sorted order) */
    for (int ir = 0; ir < rank; ir++) {
      int cnt = 0;
      int limit = nvar + 1;
      int cur_dmax = atomic_load_explicit(&ctx->dmax, memory_order_relaxed);
      if (cur_dmax > 0) {
        if ((p->outC || p->maxC || p->dW > 0 || p->min_hits > 0 || (p->debug & 1)) &&
            p->dW >= 0) {
          limit = minint(limit, cur_dmax + p->dW + 1);
        } else {
          limit = minint(limit, cur_dmax);
        }
      }

      word *rawrow = mzd_row(M_sub, ir);
      rci_t j = -1;
      while (cnt < limit) {
        j = nextelement(rawrow, M_sub->width, j);
        if (j == -1 || j >= nvar) break;
        ee[cnt++] = j++;
      }

      if (cnt > 0 && cnt < limit) {
        int nz = classical ? 1 : sparse_syndrome_non_zero(spaL0, cnt, ee);
        if (nz) {
          pthread_mutex_lock(&ctx->cw_mutex);
          p->codewords = codeword_add_maybe(p, ee, cnt);
          if (cnt < p->min_w) p->min_w = cnt;
          int best = p->min_w;
          int old_dmax = atomic_load(&ctx->dmax);
          if (old_dmax == 0 || best < old_dmax) {
            atomic_store(&ctx->dmax, best);
            if (p->debug & 16) {
              int num_rw = (ctx->p->method == 1) ? ctx->num_threads
                           : minint(ctx->rw_threads, ctx->num_threads - atomic_load(&ctx->cc_target_workers));
              if (num_rw < 1) num_rw = 1;
              fprintf(stderr,
                      "# [thread %d] RW (ksub=%d) found new upper bound cw of weight %d "
                      "(using %d RW threads)\n", tid, ksub, best, num_rw);
            }
            int cur_dmin = atomic_load(&ctx->dmin);
            if (cur_dmin > 0 && best <= cur_dmin) {
              atomic_store(&ctx->stop_flag, true);
            }
          }
          if (p->wmin > 0 && best <= p->wmin) {
            atomic_store(&ctx->stop_flag, true);
          }
          if (p->maxC && p->num_cws >= p->maxC) {
            atomic_store(&ctx->stop_flag, true);
          }
          if (check_min_hits_convergence(p)) {
            int req = (p->cov_cws > 0) ? minint(p->cov_cws, 5) : 1;
            if ((p->debug & 1) && !atomic_load(&ctx->stop_flag) &&
                !atomic_load(&ctx->rw_stop_flag)) {
              int c_min = 0, c_max = 0;
              double c_avg = 0.0, c_std = 0.0;
              compute_min_w_hit_stats(p, &c_min, &c_max, &c_avg, &c_std);
              fprintf(stderr,
                      "# RW convergence reached: >= %d min-weight (w=%d) cws "
                      "(cov=%d/%d) each hit >= %d times "
                      "(hits min=%d, max=%d, avg=%.2f, stdev=%.2f)\n",
                      req, best, p->min_w_cov, p->min_w_cws, p->min_hits,
                      c_min, c_max, c_avg, c_std);
            }
            if (p->method == 1) {
              atomic_store(&ctx->stop_flag, true);
            } else {
              atomic_store(&ctx->rw_stop_flag, true);
            }
          }
          pthread_mutex_unlock(&ctx->cw_mutex);
        }
      }
    }
    atomic_fetch_add(&ctx->rw_steps_completed, 1);
    n_done++;
    rw_time_step(ctx, &t_last, &acc_t, &acc_n, n_done == 1);
  }
  rw_time_publish(ctx, &acc_t, &acc_n);

  pthread_rwlock_unlock(&ctx->basis_rwlock);

  /* 5. Periodic adaptive basis refresh if enabled */
  if (p->refresh > 0 && !atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed) &&
      !atomic_load_explicit(&ctx->rw_stop_flag, memory_order_relaxed)) {
    long done = atomic_load_explicit(&ctx->rw_steps_completed, memory_order_relaxed);
    long next_ref = atomic_load_explicit(&ctx->next_refresh_step, memory_order_relaxed);
    if (done >= next_ref && next_ref > 0) {
      if (atomic_compare_exchange_strong(&ctx->next_refresh_step, &next_ref,
                                         done + p->refresh)) {
        int cw_wt = 0;
        pthread_mutex_lock(&ctx->cw_mutex);
        if (p->codewords && p->min_w < nvar) {
          cw_vec_t *cw, *tmp;
          HASH_ITER(hh, p->codewords, cw, tmp) {
            if (cw->weight == p->min_w) {
              cw_wt = cw->weight;
              for (int k = 0; k < cw_wt; k++) ee[k] = cw->arr[k];
              break;
            }
          }
        }
        pthread_mutex_unlock(&ctx->cw_mutex);

        pthread_rwlock_wrlock(&ctx->basis_rwlock);
        pthread_mutex_lock(&m4ri_mem_mutex);
        refresh_nullspace_basis(ctx->N_global, cw_wt > 0 ? ee : NULL, cw_wt, rng_state);
        pthread_mutex_unlock(&m4ri_mem_mutex);
        pthread_rwlock_unlock(&ctx->basis_rwlock);
      }
    }
  }
  return n_done;
}

/* Worker thread main loop */
static void *worker_thread_func(void *arg) {
  worker_arg_t *warg = (worker_arg_t *)arg;
  distfork_ctx_t *ctx = warg->ctx;
  int tid = warg->tid;
  const int nvar = ctx->p->spaH->cols;
  /* only workers 0..rw_threads-1 run RW (method 1 or 3) and allocate dense RW matrices; all workers may run CC */
  const bool enable_rw = ((ctx->p->method & 1) != 0) && (tid < ctx->rw_threads);
  const bool use_ksub = enable_rw && (ctx->p->ksub > 0);
  const int ksub_eff = (use_ksub && ctx->nu > 0) ? minint(ctx->p->ksub, ctx->nu) : 0;

  /* Initialize min_swei for this thread */
  for (int i = 0; i < MAX_W; i++) {
    warg->min_swei[i] = ctx->p->spaH->rows + 1;
  }

  /* Thread-local RW matrices (allocated safely only if RW is enabled) */
  mzd_t *mH = NULL;
  mzd_t *mHT_rw = NULL;
  mzd_t *M_sub = NULL;
  rci_t *ee = NULL;
  mzp_t *perm = NULL;
  mzp_t *pivs = NULL;
  word *piv_mask = NULL;
  int eff_nrows = ctx->p->spaH->rows;
  int *visited_cols = NULL;
  int *visited_checks = NULL;
  int *col_queue = NULL;
  int visit_marker = 0;
  uint64_t rng_state = (uint64_t)ctx->p->seed
                       + (uint64_t)tid * 0x9e3779b97f4a7c15ULL
                       + 0x517cc1b727220a95ULL;

  if (enable_rw) {
    ee = malloc((nvar + 2) * sizeof(rci_t));
    perm = safe_mzp_init(nvar);
    pivs = safe_mzp_init(nvar);
    if (ctx->p->kwin > 0 || nvar >= 500) {
      visited_cols = calloc(nvar, sizeof(int));
      visited_checks = calloc(ctx->p->spaH->rows, sizeof(int));
      col_queue = calloc(nvar, sizeof(int));
    }
    if (use_ksub) {
      if (ksub_eff > 0) {
        M_sub = safe_mzd_init(ksub_eff, nvar);
      }
    } else {
      mH = safe_mzd_from_csr(NULL, ctx->p->spaH);
      mHT_rw = safe_mzd_init(nvar, ctx->p->spaH->rows);
      piv_mask = calloc(mH->width, sizeof(word));
    }
  }

  /* Thread-local CC memory */
  const int wmax_alloc = MAX_W - 1;
  one_vec_t *err = calloc(
      1, sizeof(one_vec_t) + sizeof(int) * (wmax_alloc + 2)
  );
  one_vec_t *urr = calloc(
      1, sizeof(one_vec_t) + sizeof(int) * (wmax_alloc + 2)
  );
  one_vec_t **syn = calloc(wmax_alloc + 3, sizeof(one_vec_t *));
  for (int i = 0; i <= wmax_alloc + 2; i++) {
    syn[i] = calloc(
        1, sizeof(one_vec_t) + sizeof(int) * (ctx->p->spaH->rows + 1)
    );
  }

  while (!atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed)) {
    if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    bool did_work = false;

    /* 1. Try to take CC work if CC is active (method 2 or 3) */
    if (ctx->p->method >= 2 && atomic_load_explicit(&ctx->cc_round_active, memory_order_acquire)) {
      int active = atomic_load_explicit(&ctx->cc_active_workers, memory_order_relaxed);
      int target = atomic_load_explicit(&ctx->cc_target_workers, memory_order_relaxed);
      /* A worker without RW work joins any CC round: method 2, a CC-only worker (tid >= rw_threads), RW stopped by
       * `min_hits`, or all RW steps already claimed by other workers (it would otherwise idle until the round ends) */
      if (!enable_rw || atomic_load_explicit(&ctx->rw_stop_flag, memory_order_relaxed) ||
          atomic_load_explicit(&ctx->rw_steps_started, memory_order_relaxed) >= ctx->total_rw_steps) {
        target = ctx->num_threads;
      }
      const long long next = atomic_load_explicit(&ctx->cc_next, memory_order_relaxed);
      if (active < target && cc_unpack_idx(next) <= cc_round_end(ctx->p, cc_unpack_w(next))) {
        if (atomic_compare_exchange_weak(&ctx->cc_active_workers, &active, active + 1)) {
          const double t_join = get_time_sec();
          while (!atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed)) {
            if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
              atomic_store(&ctx->stop_flag, true);
              break;
            }
            /* claim the next start index together with the weight of its round */
            const long long claim = atomic_fetch_add_explicit(&ctx->cc_next, 1, memory_order_relaxed);
            const int w = cc_unpack_w(claim);
            const int idx = cc_unpack_idx(claim);
            if (idx > cc_round_end(ctx->p, w)) break;
            /* with the expert `start` list, `idx` enumerates the listed start columns */
            const int col = (ctx->p->start_num > 0) ? ctx->p->start_list[idx] : idx;

            err->vec[0] = urr->vec[0] = col;
            err->wei = urr->wei = 1;
            syn[1]->wei = 0;
            int swei = one_csr_row_combine(syn[1], syn[0], ctx->mHT_cc, col);

            if (ctx->p->smax && swei > 0 && swei <= ctx->p->smax) {
              if (swei < warg->min_swei[1]) {
                warg->min_swei[1] = swei;
              }
            }

            if (w > 1) {
              if (swei > 0 && swei <= (w - 1) * ctx->max_col_W + ctx->p->smax) {
                start_CC_recurs_mt(err, urr, syn, w, ctx->max_col_W,
                                   ctx->p->spaH, ctx->mHT_cc, warg);
              }
            } else {
              if (!swei) {
                int nz = (!ctx->p->spaL) ||
                         sparse_syndrome_non_zero(ctx->p->spaL, 1, err->vec);
                if (nz) {
                  pthread_mutex_lock(&ctx->cw_mutex);
                  ctx->p->codewords = codeword_add_maybe(ctx->p, err->vec, 1);
                  atomic_store(&ctx->cc_found_weight, 1);
                  atomic_store(&ctx->cc_exact, true); /* weight 1 is always the exact distance */
                  atomic_store(&ctx->dmin, 1);
                  atomic_store(&ctx->dmax, 1);
                  atomic_store(&ctx->stop_flag, true);
                  pthread_mutex_unlock(&ctx->cw_mutex);
                }
              }
            }
            err->wei = urr->wei = 0;
          }
          /* CC thread time (added before leaving the round, so the coordinator sees it once the round is over) */
          atomic_fetch_add_explicit(&ctx->cc_busy_ns, (long long)((get_time_sec() - t_join) * 1e9),
                                    memory_order_relaxed);
          atomic_fetch_sub_explicit(&ctx->cc_active_workers, 1, memory_order_release);
          did_work = true;
          continue;
        }
      }
    }

    /* 2. Try to take RW work if RW is active (method 1 or 3) */
    if (enable_rw && !atomic_load(&ctx->stop_flag) && !atomic_load(&ctx->rw_stop_flag)) {
      long cur_s = atomic_load(&ctx->rw_steps_started);
      if (cur_s < ctx->total_rw_steps) {
        long batch_size;
        if (ctx->p->chunk_size > 0) {
          batch_size = ctx->p->chunk_size;
        } else {
          int nvar = ctx->p->nvar;
          if (nvar < 500) {
            if (ctx->total_rw_steps >= 50000) batch_size = 500;
            else if (ctx->total_rw_steps >= 1000) batch_size = 250;
            else batch_size = 50;
          } else if (nvar < 5000) {
            if (ctx->total_rw_steps >= 10000) batch_size = 100;
            else batch_size = 50;
          } else {
            /* Large matrices (e.g. n >= 5000): keep chunk bounded to ~2-3s */
            if (ctx->total_rw_steps >= 10000) batch_size = 50;
            else batch_size = 25;
          }
        }
        /* Prevent thread starvation: ensure chunk size doesn't monopolize steps across the RW threads */
        if (ctx->rw_threads > 1 && ctx->total_rw_steps > 0) {
          long max_chunk = (ctx->total_rw_steps + ctx->rw_threads - 1) / ctx->rw_threads;
          if (max_chunk >= 1 && batch_size > max_chunk) {
            batch_size = max_chunk;
          }
        }
        long target_s = cur_s + batch_size;
        if (target_s > ctx->total_rw_steps) target_s = ctx->total_rw_steps;
        if (atomic_compare_exchange_weak(&ctx->rw_steps_started, &cur_s, target_s)) {
          int n_steps = (int)(target_s - cur_s);
          if (use_ksub) {
            if (ksub_eff > 0) {
              run_rw_steps_ksub(ctx, n_steps, M_sub, ee, perm, pivs,
                                visited_cols, visited_checks, col_queue,
                                &visit_marker, &rng_state, tid);
            } else {
              atomic_fetch_add(&ctx->rw_steps_completed, n_steps);
            }
          } else {
            run_rw_steps(ctx, n_steps, mH, mHT_rw, ee, perm, pivs,
                         piv_mask, &eff_nrows, visited_cols, visited_checks,
                         col_queue, &visit_marker, &rng_state, tid);
          }
          did_work = true;
          continue;
        }
      }
    }

    if (!did_work) {
      usleep(100);
    }
  }

  if (enable_rw) {
    free(piv_mask);
    safe_mzp_free(perm);
    safe_mzp_free(pivs);
    free(ee);
    free(visited_cols);
    free(visited_checks);
    free(col_queue);
    safe_mzd_free(M_sub);
    safe_mzd_free(mHT_rw);
    safe_mzd_free(mH);
  }

  for (int i = 0; i <= wmax_alloc + 2; i++) free(syn[i]);
  free(syn);
  free(err);
  free(urr);

  return NULL;
}

/* Method 1 coordinator */
static void run_method1_coordinator(distfork_ctx_t *ctx) {
  if (ctx->p->debug & 2) {
    fprintf(stderr, "# running method=1 (multithreaded RW) with %d threads, total steps=%ld\n",
            ctx->num_threads, ctx->total_rw_steps);
  }

  while (!atomic_load(&ctx->stop_flag)) {
    if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }
    if (atomic_load(&ctx->rw_steps_completed) >= ctx->total_rw_steps) {
      break;
    }
    usleep(1000);
  }
}

/* Estimated work of the CC round at weight w in thread-seconds (the unit of `cc_time_per_weight`, i.e., wall time
 * times the number of CC threads): the measured work W(w-1) of the previous round times the growth factor
 * W(w-1)/W(w-2), clamped to [2, 10] (4 if W(w-2) is not known), also returned in `*growth`.  Returns 0 if W(w-1)
 * has not been measured (the first round at a supplied dmin > 1). */
static double cc_work_estimate(const distfork_ctx_t * const ctx, const int w, double * const growth) {
  *growth = 4.0;
  if (w <= 1) return 1e-4; /* one syndrome per column */
  if (w > MAX_W) return 0.0;
  const double prev = ctx->cc_time_per_weight[w - 1];
  if (!(prev > 0.0)) return 0.0;
  if (w >= 3 && ctx->cc_time_per_weight[w - 2] > 1e-4) {
    const double g = prev / ctx->cc_time_per_weight[w - 2];
    *growth = (g < 2.0) ? 2.0 : ((g > 10.0) ? 10.0 : g);
  }
  return ((prev > 1e-4) ? prev : 1e-4) * (*growth);
}

/* Method 2 coordinator */
static void run_method2_coordinator(distfork_ctx_t *ctx) {
  const int wmax = ctx->p->wmax;
  const int w_start = ctx->p->noscan ? wmax : (ctx->p->dmin > 1 ? ctx->p->dmin : 1);
  /* with outC or maxC, CC does not stop at the first codeword but collects all codewords of a round */
  const bool collecting = (ctx->p->outC != NULL) || (ctx->p->maxC > 0);
  /* with outC and dW > 0, extra CC rounds export the codewords of weight up to d + dW */
  const int extra_w = (ctx->p->outC && ctx->p->dW > 0) ? ctx->p->dW : 0;

  if (ctx->p->debug & 2) {
    fprintf(stderr, "# running method=2 (multithreaded CC) with %d threads, w_start=%d wmax=%d\n",
            ctx->num_threads, w_start, wmax);
  }

  int w_limit = wmax;
  for (int w = w_start; ; w++) {
    if (atomic_load(&ctx->stop_flag)) break;
    /* The current upper bound dmax (given, read from finC, or found by CC) limits the rounds: since the CC kernel
     * caps the cluster weight at dmax + dW, a round w > dmax + extra_w would only enumerate the same clusters
     * again.  Unless codewords are collected, the round w = dmax is not needed once dmin = dmax. */
    const int cur_dmax = atomic_load(&ctx->dmax);
    if (cur_dmax > 0) {
      if (!collecting && atomic_load(&ctx->dmin) >= cur_dmax) {
        if (ctx->p->debug & 1) {
          fprintf(stderr, "# bounds coincide: dmin = dmax = %d (CC round w=%d not needed)\n", cur_dmax, w);
        }
        break;
      }
      w_limit = minint(w_limit, cur_dmax + extra_w);
    }
    if (w > w_limit) break;

    double now = get_time_sec();
    double remaining_time = ctx->timeout - (now - ctx->start_time);
    if (ctx->timeout > 0.0 && remaining_time <= 0.0) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    /* With a timeout, a round predicted not to finish in time is started anyway: it can still end early with a
     * codeword of weight w, which gives the exact distance (all lower weights have been analyzed) */
    if (ctx->timeout > 0.0 && (ctx->p->debug & 1)) {
      double growth;
      const double t_cc_est = cc_work_estimate(ctx, w, &growth) / ctx->num_threads; /* wall time */
      if (t_cc_est > remaining_time) {
        fprintf(stderr, "# CC for w=%d (est %.2fs) may not finish in the remaining time %.2fs, searching anyway "
                "(dmin=%d)\n", w, t_cc_est, remaining_time, atomic_load(&ctx->dmin));
      }
    }

    const int beg = cc_round_beg(ctx->p);
    const int end = cc_round_end(ctx->p, w);
    /* all weights below w analyzed (or a supplied dmin): a codeword of weight w found is exact */
    const bool certified = (w == atomic_load(&ctx->dmin));

    atomic_store(&ctx->cc_target_workers, ctx->num_threads);
    const long long busy_ns0 = atomic_load(&ctx->cc_busy_ns);
    atomic_store(&ctx->cc_next, cc_pack(w, beg));
    atomic_store(&ctx->cc_round_active, 1);

    double cc_start = get_time_sec();

    if (ctx->p->debug & 2) {
      const char *note = (certified || atomic_load(&ctx->cc_exact)) ? "" : " (lower weights not scanned)";
      if (ctx->p->start_num > 0) {
        fprintf(stderr, "# searching w=%d with %d CC threads, start list (%d columns)%s\n",
                w, ctx->num_threads, ctx->p->start_num, note);
      } else {
        fprintf(stderr, "# searching w=%d with %d CC threads, columns [%d, %d]%s\n",
                w, ctx->num_threads, beg, end, note);
      }
    }

    bool round_completed = false;
    while (!atomic_load(&ctx->stop_flag)) {
      if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
        atomic_store(&ctx->stop_flag, true);
        break;
      }
      if (cc_unpack_idx(atomic_load(&ctx->cc_next)) > end && atomic_load(&ctx->cc_active_workers) == 0) {
        round_completed = true;
        break;
      }
      usleep(100);
    }

    atomic_store(&ctx->cc_round_active, 0);

    double cc_dur = get_time_sec() - cc_start;
    if (w < MAX_W) { /* CC work of the round in thread-seconds, as measured by the workers */
      ctx->cc_time_per_weight[w] = 1e-9 * (double)(atomic_load(&ctx->cc_busy_ns) - busy_ns0);
    }

    int cw_found = atomic_load(&ctx->cc_found_weight);
    if (cw_found > 0 && certified && cw_found == w && !atomic_load(&ctx->cc_exact)) {
      atomic_store(&ctx->cc_exact, true);
    }
    const bool exact = (cw_found > 0) && atomic_load(&ctx->cc_exact);
    if (exact) {
      atomic_store(&ctx->dmin, cw_found);
      atomic_store(&ctx->dmax, cw_found);
      /* the distance is known: no rounds beyond w = cw_found + extra_w */
      w_limit = minint(w_limit, cw_found + extra_w);
      const bool more_rounds = round_completed && (w < w_limit);
      if (ctx->p->debug & 1) {
        if (w == cw_found && more_rounds) {
          fprintf(stderr,
                  "# CC round w=%d finished in %.3fs (%d CC threads): found min-weight codewords "
                  "(dmin=%d, continuing up to w=%d for dW=%d, total %lld cws)\n",
                  w, cc_dur, ctx->num_threads, cw_found, w_limit, ctx->p->dW, ctx->p->num_cws);
        } else if (w == cw_found) {
          fprintf(stderr, "# CC found min-weight codeword: d=%d (using %d CC threads, total %lld cws)\n",
                  cw_found, ctx->num_threads, ctx->p->num_cws);
        } else if (round_completed) {
          fprintf(stderr,
                  "# CC round w=%d finished in %.3fs (%d CC threads): extra dW round completed "
                  "(dmin=%d, total %lld cws)\n",
                  w, cc_dur, ctx->num_threads, cw_found, ctx->p->num_cws);
        }
      }
      if (!more_rounds) {
        atomic_store(&ctx->stop_flag, true);
        break;
      }
    } else if (cw_found > 0) {
      /* noscan=1 without dmin=wmax: lower weights not scanned, the codeword only gives an upper bound */
      if (ctx->p->debug & 1) {
        fprintf(stderr,
                "# CC round w=%d finished in %.3fs (%d CC threads): found codeword of weight %d -> dmax=%d "
                "(dmin=%d not certified: lower weights not scanned)\n",
                w, cc_dur, ctx->num_threads, cw_found, atomic_load(&ctx->dmax), atomic_load(&ctx->dmin));
      }
      atomic_store(&ctx->stop_flag, true);
      break;
    } else {
      if (!round_completed) {
        break;
      }
      if (certified) {
        /* Weight w analyzed without success */
        atomic_store(&ctx->dmin, w + 1);
        if (ctx->p->debug & 1) {
          fprintf(stderr, "# CC w=%d completed in %.3fs (%d CC threads): no codewords found -> dmin=%d\n",
                  w, cc_dur, ctx->num_threads, w + 1);
        }
      } else if (ctx->p->debug & 1) {
        fprintf(stderr, "# CC w=%d completed in %.3fs (%d CC threads): no codewords found "
                "(dmin=%d not raised: lower weights not scanned)\n",
                w, cc_dur, ctx->num_threads, atomic_load(&ctx->dmin));
      }
    }
  }
}

/* Method 3: RW still runs (RW threads, RW steps not all completed, RW not stopped by min_hits) */
static inline bool m3_rw_running(distfork_ctx_t * const ctx) {
  return ctx->rw_threads > 0 && !atomic_load(&ctx->rw_stop_flag) &&
         atomic_load(&ctx->rw_steps_completed) < ctx->total_rw_steps;
}

/* Method 3: RW thread time per step, measured continuously by the RW workers (before the first step is completed, a
 * step takes at least the time elapsed since the start) */
static double m3_rw_step_time(distfork_ctx_t * const ctx) {
  const long n = atomic_load(&ctx->rw_timed_steps);
  if (n > 0) return 1e-9 * (double)atomic_load(&ctx->rw_busy_ns) / (double)n;
  const double elapsed = get_time_sec() - ctx->start_time;
  return (elapsed > 5e-5) ? elapsed : 5e-5;
}

/* Method 3: the largest CC weight needed: dmax-1 once an upper bound dmax is known (dmax+dW when exporting codewords
 * with outC), otherwise wmax (or n); at most MAX_W-2.  (`dexp` only affects the thread split, see
 * m3_plan_cc_threads().) */
static int m3_cc_target_w(const distfork_ctx_t * const ctx, const int cur_dmin, const int cur_dmax) {
  const params_t * const p = ctx->p;
  int target;
  if (cur_dmax > 0) {
    if (p->outC && (p->dW > 0 || cur_dmin >= cur_dmax)) {
      target = cur_dmax + (p->dW > 0 ? p->dW : 0);
    } else {
      target = cur_dmax - 1;
    }
  } else {
    target = (p->wmax > 0) ? p->wmax : p->spaH->cols;
  }
  const int max_allowed_w = (p->wmax > 0) ? minint(p->wmax, MAX_W - 2) : (MAX_W - 2);
  return minint(target, max_allowed_w);
}

#define M3_DEADLINE_FRACTION 0.75 /* with a timeout, a round gets enough CC threads to finish in 3/4 of the time left */
#define M3_OVERRUN_FACTOR 1.5     /* a round taking 1.5 times (+20 ms) longer than predicted gets all threads */

/* Method 3: number of CC threads for the round at weight w with the estimated work `est` (thread-seconds, 0 if not
 * known) and the growth factor `growth` (see cc_work_estimate()); the other threads run RW.  Returns 0 if CC is
 * paused (`dexp`, below). */
static int m3_plan_cc_threads(distfork_ctx_t * const ctx, const int w, const double est, const double growth,
                              const double remaining_time) {
  const int nthr = ctx->num_threads;
  if (nthr <= 1) return 1;
  /* RW finished, or all RW steps claimed: the RW workers join the round anyway */
  if (!m3_rw_running(ctx) || atomic_load(&ctx->rw_steps_started) >= ctx->total_rw_steps) return nthr;
  const int cur_dmin = atomic_load(&ctx->dmin);
  const int cur_dmax = atomic_load(&ctx->dmax);
  /* The round w = dmax-1 certifies dmin = dmax (later rounds only collect codewords).  RW cannot improve this result,
   * it could only find a codeword of weight w before CC does: all threads run CC. */
  if (cur_dmax > 0 && (w >= cur_dmax - 1 || cur_dmin >= cur_dmax)) return nthr;
  /* `dexp` hint: as long as RW has not found any codeword, CC rounds w > dexp run only on the threads which cannot
   * run RW (none: CC pauses) */
  if (cur_dmax == 0 && ctx->dexp > 0 && w > ctx->dexp) return nthr - ctx->rw_threads;

  int n_cc;
  if (!(est > 0.0)) {
    n_cc = nthr / 2; /* the first round at a supplied dmin: its work is not known */
  } else if (est < 0.005) {
    n_cc = (nthr >= 4) ? 2 : 1;
  } else {
    /* Split by the remaining work in thread-seconds: the CC rounds w, w+1, w+2 up to the horizon (dmax-1, or dexp
     * before RW has found a codeword), and the remaining RW steps (at most 2000 once dmax is known), but not more RW
     * work than the RW threads can do before the timeout */
    int horizon = m3_cc_target_w(ctx, cur_dmin, cur_dmax);
    if (cur_dmax == 0 && ctx->dexp > 0 && ctx->dexp < horizon) horizon = ctx->dexp;
    double t_cc = est, t_k = est;
    for (int k = w + 1; k <= horizon && k <= w + 2; k++) {
      t_k *= growth;
      t_cc += t_k;
    }
    const long steps_rem = ctx->total_rw_steps - atomic_load(&ctx->rw_steps_completed);
    const long eff_steps = (cur_dmax > 0 && steps_rem > 2000) ? 2000 : steps_rem;
    double t_rw = (double)eff_steps * m3_rw_step_time(ctx);
    if (ctx->timeout > 0.0 && t_rw > remaining_time * ctx->rw_threads) t_rw = remaining_time * ctx->rw_threads;
    n_cc = (int)round((double)nthr * t_cc / (t_cc + t_rw));
    if (n_cc < 1) n_cc = 1;
    if (n_cc > nthr - 1) n_cc = nthr - 1; /* keep at least one RW thread */
  }
  /* timeout: enough CC threads to finish the round in time */
  if (ctx->timeout > 0.0 && est > 0.0 && remaining_time > 0.0) {
    const double n_min = ceil(est / (M3_DEADLINE_FRACTION * remaining_time));
    if (n_min > n_cc) n_cc = (n_min < nthr) ? (int)n_min : nthr;
  }
  /* only `rw_threads` workers can run RW (memory / steps limits): the other workers always run CC */
  if (nthr - n_cc > ctx->rw_threads) n_cc = nthr - ctx->rw_threads;
  return n_cc;
}

/* Method 3: wait while RW runs, until the upper bound dmax differs from `dmax0`, RW ends, or the run stops */
static void m3_wait_rw(distfork_ctx_t * const ctx, const int dmax0) {
  while (!atomic_load(&ctx->stop_flag) && m3_rw_running(ctx) && atomic_load(&ctx->dmax) == dmax0) {
    if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }
    usleep(1000);
  }
}

/* Method 3 coordinator */
static void run_method3_coordinator(distfork_ctx_t *ctx) {
  const int nvar = ctx->p->spaH->cols;
  const int nthr = ctx->num_threads;
  int w = ctx->p->noscan ? ctx->p->wmax : (ctx->p->dmin > 1 ? ctx->p->dmin : 1);

  if (ctx->p->debug & 2) {
    fprintf(stderr, "# running method=3 (bracketing mode) with %d threads, timeout=%.1fs, dexp=%d\n",
            ctx->num_threads, ctx->timeout, ctx->dexp);
  }

  int init_dmax = atomic_load(&ctx->dmax);
  int init_dmin = atomic_load(&ctx->dmin);
  if (init_dmax > 0 && init_dmin >= init_dmax && !ctx->p->outC) {
    atomic_store(&ctx->stop_flag, true);
    return;
  }

  /* No RW probe: the RW step time is measured continuously by the RW workers (m3_rw_step_time()).  CC is never ended
   * while RW still runs: if no CC round can be started (CC done up to wmax, a round predicted to exceed the timeout,
   * or CC paused by `dexp`), the coordinator waits for RW and re-plans when dmax changes or RW ends. */
  int wait_msg_w = 0; /* weight w for which a waiting message was printed (once per round) */
  while (!atomic_load(&ctx->stop_flag)) {
    double now = get_time_sec();
    double remaining_time = (ctx->timeout > 0.0) ? (ctx->timeout - (now - ctx->start_time)) : 1e9;
    if (ctx->timeout > 0.0 && remaining_time <= 0.0) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    int cur_dmax = atomic_load(&ctx->dmax);
    int cur_dmin = atomic_load(&ctx->dmin);
    const bool rw_on = m3_rw_running(ctx);

    /* Target cluster size for CC */
    const int target_cc_w = m3_cc_target_w(ctx, cur_dmin, cur_dmax);

    if (cur_dmax > 0 && cur_dmin >= cur_dmax && w > target_cc_w) {
      /* Bracketing converged and all requested dW rounds completed */
      atomic_store(&ctx->dmin, cur_dmax);
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    if (w > target_cc_w) {
      /* CC done up to wmax: only RW can still lower dmax */
      if (!rw_on) break;
      if ((ctx->p->debug & 2) && wait_msg_w != w) {
        wait_msg_w = w;
        fprintf(stderr, "# CC done up to w=%d, waiting for RW (bounds [%d, %d])\n", target_cc_w, cur_dmin, cur_dmax);
      }
      m3_wait_rw(ctx, cur_dmax);
      continue;
    }

    /* Estimated CC work of the round (thread-seconds).  A round predicted not to finish in time on all threads is not
     * started while RW runs; once RW has ended, the run ends.  Without RW (steps=0), as in method 2, the round is
     * started anyway: it can still end early with a codeword of weight w (exact if w = dmin). */
    double growth;
    const double t_cc_est = cc_work_estimate(ctx, w, &growth);
    if (ctx->timeout > 0.0 && t_cc_est / nthr > remaining_time) {
      if (rw_on) {
        if ((ctx->p->debug & 2) && wait_msg_w != w) {
          wait_msg_w = w;
          fprintf(stderr, "# CC for w=%d (est %.2fs) exceeds remaining timeout %.2fs, devoting %d threads to RW\n",
                  w, t_cc_est / nthr, remaining_time, ctx->rw_threads);
        }
        m3_wait_rw(ctx, cur_dmax);
        continue;
      }
      if (ctx->rw_threads > 0) {
        if (ctx->p->debug & 1) {
          fprintf(stderr, "# CC for w=%d (est %.2fs) exceeds remaining timeout %.2fs, RW ended: terminating early "
                  "(dmin=%d)\n", w, t_cc_est / nthr, remaining_time, cur_dmin);
        }
        atomic_store(&ctx->stop_flag, true);
        break;
      }
      if ((ctx->p->debug & 1) && wait_msg_w != w) {
        wait_msg_w = w;
        fprintf(stderr, "# CC for w=%d (est %.2fs) may not finish in the remaining time %.2fs, searching anyway "
                "(dmin=%d)\n", w, t_cc_est / nthr, remaining_time, cur_dmin);
      }
    }

    /* Calculate thread balancing */
    int n_cc = m3_plan_cc_threads(ctx, w, t_cc_est, growth, remaining_time);
    if (n_cc <= 0) {
      if ((ctx->p->debug & 2) && wait_msg_w != w) {
        wait_msg_w = w;
        fprintf(stderr, "# CC paused before w=%d > dexp=%d until RW finds a codeword (or ends)\n", w, ctx->dexp);
      }
      m3_wait_rw(ctx, cur_dmax);
      continue;
    }
    const int n_cc0 = n_cc;
    int n_rw = nthr - n_cc;
    const long steps_rem = rw_on ? ctx->total_rw_steps - atomic_load(&ctx->rw_steps_completed) : 0;

    const int beg = cc_round_beg(ctx->p);
    const int end = cc_round_end(ctx->p, w);
    /* all weights below w analyzed (or a supplied dmin): a codeword of weight w found is exact */
    const bool certified = (w == atomic_load(&ctx->dmin));

    atomic_store(&ctx->cc_target_workers, n_cc);
    const long long busy_ns0 = atomic_load(&ctx->cc_busy_ns);
    atomic_store(&ctx->cc_next, cc_pack(w, beg));
    atomic_store(&ctx->cc_round_active, 1);

    if (ctx->p->debug & 2) {
      char est_str[64] = "unknown";
      char rw_str[64] = "";
      if (t_cc_est > 0.0) snprintf(est_str, sizeof(est_str), "%.3g thread-s", t_cc_est);
      if (rw_on) snprintf(rw_str, sizeof(rw_str), ", RW step %.3g s", m3_rw_step_time(ctx));
      fprintf(stderr,
              "# CC round w=%d started: %d CC threads, %d RW threads "
              "(bounds [%d, %d], rem_rw=%ld, rem_time=%.2fs, est %s%s)\n",
              w, n_cc, n_rw, cur_dmin, cur_dmax, steps_rem, remaining_time, est_str, rw_str);
    }

    double cc_start = get_time_sec();
    bool round_completed = false;
    int dmax_seen = cur_dmax;
    unsigned int poll = 0;

    while (!atomic_load(&ctx->stop_flag)) {
      const double t = get_time_sec();
      if (ctx->timeout > 0.0 && (t - ctx->start_time >= ctx->timeout)) {
        atomic_store(&ctx->stop_flag, true);
        break;
      }
      if (cc_unpack_idx(atomic_load(&ctx->cc_next)) > end && atomic_load(&ctx->cc_active_workers) == 0) {
        round_completed = true;
        break;
      }
      /* Re-plan during the round (only adding CC threads, workers do not leave a round): when RW finds a new upper
       * bound (e.g., this round now certifies dmin = dmax), and when the round takes much longer than predicted */
      if ((++poll & 15) == 0 && n_cc < nthr) {
        int n_new = n_cc;
        const char *why = "";
        const int d_now = atomic_load(&ctx->dmax);
        if (d_now != dmax_seen) {
          dmax_seen = d_now;
          const double rem = (ctx->timeout > 0.0) ? ctx->timeout - (t - ctx->start_time) : 1e9;
          n_new = m3_plan_cc_threads(ctx, w, t_cc_est, growth, rem);
          why = "new upper bound";
        }
        const bool paused = (d_now == 0 && ctx->dexp > 0 && w > ctx->dexp && m3_rw_running(ctx));
        if (!paused && t_cc_est > 0.0 && t - cc_start > M3_OVERRUN_FACTOR * t_cc_est / n_cc0 + 0.02) {
          n_new = nthr;
          why = "longer than predicted";
        }
        if (n_new > n_cc) {
          if (ctx->p->debug & 2) {
            fprintf(stderr, "# CC round w=%d: %d -> %d CC threads (%s, bounds [%d, %d], %.3fs into the round)\n",
                    w, n_cc, n_new, why, atomic_load(&ctx->dmin), d_now, t - cc_start);
          }
          n_cc = n_new;
          n_rw = nthr - n_cc;
          atomic_store(&ctx->cc_target_workers, n_cc);
        }
      }
      usleep(100);
    }

    atomic_store(&ctx->cc_round_active, 0);

    double cc_dur = get_time_sec() - cc_start;
    if (w < MAX_W) {
      /* CC work of the round in thread-seconds, as measured by the workers: the number of CC threads may exceed
       * n_cc (CC-only workers, and RW workers that ran out of RW steps join the round) */
      ctx->cc_time_per_weight[w] = 1e-9 * (double)(atomic_load(&ctx->cc_busy_ns) - busy_ns0);
    }

    int cw_found = atomic_load(&ctx->cc_found_weight);
    if (cw_found > 0 && certified && cw_found == w && !atomic_load(&ctx->cc_exact)) {
      atomic_store(&ctx->cc_exact, true);
    }
    const bool exact = (cw_found > 0) && atomic_load(&ctx->cc_exact);
    if (exact) {
      atomic_store(&ctx->dmin, cw_found);
      atomic_store(&ctx->dmax, cw_found);

      int max_w_lim = minint(ctx->p->wmax > 0 ? ctx->p->wmax : nvar, cw_found + ctx->p->dW);
      if (ctx->p->outC && ctx->p->dW > 0 && w < max_w_lim) {
        if (ctx->p->debug & 1) {
          if (w == cw_found) {
            fprintf(stderr,
                    "# CC round w=%d finished in %.3fs (%d CC threads, %d RW threads): found codewords "
                    "(dmin=%d, continuing up to w=%d for dW=%d, total %lld cws)\n",
                    w, cc_dur, n_cc, n_rw, cw_found, max_w_lim, ctx->p->dW, ctx->p->num_cws);
          } else if (round_completed) {
            fprintf(stderr,
                    "# CC round w=%d finished in %.3fs (%d CC threads, %d RW threads): "
                    "extra dW round completed (dmin=%d, total %lld cws)\n",
                    w, cc_dur, n_cc, n_rw, cw_found, ctx->p->num_cws);
          }
        }
      } else {
        if (ctx->p->debug & 1) {
          if (w > cw_found) {
            if (round_completed) {
              fprintf(stderr,
                      "# CC round w=%d finished in %.3fs (%d CC threads, %d RW threads): "
                      "extra dW round completed (dmin=%d, total %lld cws)\n",
                      w, cc_dur, n_cc, n_rw, cw_found, ctx->p->num_cws);
            }
          } else {
            fprintf(stderr, "# CC found min-weight codeword: d=%d (using %d CC threads, total %lld cws)\n",
                    cw_found, n_cc, ctx->p->num_cws);
          }
        }
        atomic_store(&ctx->stop_flag, true);
        break;
      }
    } else if (cur_dmax > 0 && cur_dmin >= cur_dmax) {
      /* Extra dW round completed */
      if (round_completed && (ctx->p->debug & 1)) {
        fprintf(stderr,
                "# CC round w=%d finished in %.3fs (%d CC threads, %d RW threads): "
                "extra dW round completed (dmin=%d, total %lld cws)\n",
                w, cc_dur, n_cc, n_rw, cur_dmin, ctx->p->num_cws);
      }
      if (!round_completed) {
        break;
      }
    } else {
      if (!round_completed) {
        break;
      }

      /* Weight w analyzed without success */
      int new_dmin = w + 1;
      atomic_store(&ctx->dmin, new_dmin);
      if (ctx->p->debug & 1) {
        fprintf(stderr, "# CC round w=%d finished in %.3fs (%d CC threads, %d RW threads): no codewords -> dmin=%d\n",
                w, cc_dur, n_cc, n_rw, new_dmin);
      }

      cur_dmax = atomic_load(&ctx->dmax);
      if (cur_dmax > 0 && new_dmin >= cur_dmax) {
        atomic_store(&ctx->dmin, cur_dmax);
        if (ctx->p->outC && ctx->p->dW > 0 && cur_dmax + ctx->p->dW > cur_dmax) {
          if (ctx->p->debug & 1) {
            fprintf(stderr, "# bracketing bounds coincide: dmin = dmax = %d (continuing up to w=%d for dW=%d)\n",
                    cur_dmax, cur_dmax + ctx->p->dW, ctx->p->dW);
          }
        } else {
          atomic_store(&ctx->stop_flag, true);
          if (ctx->p->debug & 1) {
            fprintf(stderr, "# bracketing bounds coincide: dmin = dmax = %d\n", cur_dmax);
          }
          break;
        }
      }
    }

    w++;
  }
}

/* Maximum number of RW workers for thread throttling (all workers in method 1, the RW share in method 3):
 * - dense memory: a full-matrix RW worker holds dense copies of H (m x n) and of its transpose (n x m), a worker
 *   with ksub > 0 only the sampled subspace matrix (ksub x n).  Above 15 MB per worker, the total is kept at
 *   ~1.5 GB (2 to 32 workers) to prevent DRAM bus and L3 cache thrashing;
 * - small step counts: at most ceil(steps/10) RW workers. */
static int rw_thread_cap(const params_t * const p) {
  int cap = INT_MAX;
  const int nvar = p->nvar;
  const int nrows = p->spaH ? p->spaH->rows : 0;
  if (nrows > 0 && nvar > 0) {
    const size_t words_n = ((size_t)nvar + 63) / 64;
    size_t dense_bytes_per_thread;
    if (p->ksub > 0) {
      dense_bytes_per_thread = (size_t)minint(p->ksub, nrows) * words_n * 8;
    } else {
      dense_bytes_per_thread = (size_t)nrows * words_n * 8 + (size_t)nvar * (((size_t)nrows + 63) / 64) * 8;
    }
    if (dense_bytes_per_thread > 15ULL * 1024 * 1024) {
      int max_mem_threads = (int)((1536ULL * 1024 * 1024) / dense_bytes_per_thread);
      if (max_mem_threads < 2) max_mem_threads = 2;
      if (max_mem_threads > 32) max_mem_threads = 32;
      cap = max_mem_threads;
    }
  }
  if (p->steps > 0) {
    const int max_step_threads = (int)(((long)p->steps + 9) / 10);
    cap = minint(cap, max_step_threads);
  }
  return cap;
}

int main(int argc, char **argv) {
  params_t * const p = &prm;

  var_init(argc, argv, p);

  if (p->finC) {
    nzlist_read(p->finC, p);
  }

  /* Check ksub feasibility before thread throttling (require m >= nu >= n - m) */
  if ((p->method & 1) && p->ksub > 0 && p->spaH) {
    int nrows = p->spaH->rows;
    int nvar = p->spaH->cols;
    if (2 * nrows < nvar) {
      if (p->debug & 1) {
        fprintf(stderr,
                "# Warning: ksub=%d requested, but m=%d < n-m=%d (<= nu); "
                "falling back to full-matrix RW (ksub=0)\n",
                p->ksub, nrows, nvar - nrows);
      }
      p->ksub = 0;
    }
  }

  /* Determine number of threads */
  int requested_threads = p->threads;
  int num_threads = p->threads;
  if (num_threads <= 0) {
    long nprocs = sysconf(_SC_NPROCESSORS_ONLN);
    num_threads = (nprocs > 0) ? (int)nprocs : 4;
    if (num_threads > 64) num_threads = 64;
  }
  /* Workers which may run RW (each allocates dense RW matrices): all workers in method 1, none in method 2.
   * In method 3, the RW share may be limited further while CC rounds can use all workers. */
  int rw_threads = 0;

  if (!p->nothrottle) {
    /* Thread throttling for small codes and large-memory matrices */
    int nvar = p->nvar;
    int nrows = p->spaH ? p->spaH->rows : 0;
    unsigned long long n_elements = (unsigned long long)nrows * (unsigned long long)nvar;

    /* 1. Small-code throttling: avoid thread spawn/join overhead on tiny workloads */
    double rw_work = (p->method & 1)
                     ? ((double)n_elements * (double)(p->steps > 0 ? p->steps : 1))
                     : 0.0;
    if (nvar < 60 || (n_elements < 20000ULL && rw_work < 5.0e7)) {
      if (num_threads > 4) num_threads = 4;
    } else if ((nvar < 150 || n_elements < 100000ULL) && rw_work < 2.0e8) {
      if (num_threads > 16) num_threads = 16;
    }

    /* 2. Large-matrix memory throttling and 3. steps throttling: these limit only the RW workers, i.e., all
     *    workers in method 1 and the RW share in method 3; CC (method 2, CC rounds in method 3) is not limited */
    if (p->method & 1) {
      const int rw_cap = rw_thread_cap(p);
      if (p->method == 1 && num_threads > rw_cap) num_threads = rw_cap;
      rw_threads = minint(num_threads, rw_cap);
    }

    if (requested_threads > 0 && num_threads != requested_threads && (p->debug & 2)) {
      fprintf(stderr, "# Note: throttled threads from %d to %d (n=%d, r=%d)\n",
              requested_threads, num_threads, nvar, nrows);
    }
  } else {
    if (p->method & 1) rw_threads = num_threads;
    if (p->debug & 2) {
      fprintf(stderr, "# Thread throttling disabled (nothrottle=1); running with %d threads\n", num_threads);
    }
  }
  if (p->method == 3 && p->steps == 0) {
    rw_threads = 0; /* no RW steps: pure CC through the bracketing coordinator, no RW matrices needed */
  }
  if (p->method == 3 && rw_threads < num_threads && (p->debug & 2)) {
    if (rw_threads == 0) {
      fprintf(stderr, "# Note: steps=0, no RW; CC rounds use all %d threads\n", num_threads);
    } else {
      fprintf(stderr, "# Note: RW limited to %d of %d threads (memory, steps); CC rounds can use all threads\n",
              rw_threads, num_threads);
    }
  }

  double timeout = (p->timeout > 0.0) ? p->timeout : 0.0;

  distfork_ctx_t ctx;
  memset(&ctx, 0, sizeof(ctx));
  ctx.p = p;
  ctx.num_threads = num_threads;
  ctx.rw_threads = rw_threads;
  ctx.timeout = timeout;
  ctx.start_time = get_time_sec();
  ctx.dexp = p->dexp;
  ctx.total_rw_steps = (p->steps >= 0) ? p->steps : 1;

  /* Initialize dmin and dmax */
  atomic_init(&ctx.dmin, p->dmin > 1 ? p->dmin : 1);
  int init_dmax = 0;
  if (p->dmax > 0) {
    init_dmax = p->dmax;
  }
  if (p->min_w != INT_MAX) {
    if (init_dmax == 0 || p->min_w < init_dmax) {
      init_dmax = p->min_w;
    }
  }
  atomic_init(&ctx.dmax, init_dmax);

  if (init_dmax > 0 && p->wmin > 0 && init_dmax <= p->wmin && !p->outC) {
    if (p->debug & 2) {
      fprintf(stderr, "# early termination due to wmin=%d (known dmax=%d <= wmin)\n", p->wmin, init_dmax);
    }
    printf("%d %d 0\n", p->dmin > 1 ? p->dmin : 1, init_dmax);
    var_kill(p);
    return 0;
  }

  if (p->method == 3 && init_dmax > 0 && p->dmin > 1 && p->dmin >= init_dmax && !p->outC) {
    if (p->debug & 2) {
      fprintf(stderr, "# running method=3 (bracketing mode) with %d threads, timeout=%.1fs, dexp=%d\n",
              num_threads, timeout, p->dexp);
    }
    printf("%d %d 0\n", p->dmin, init_dmax);
    var_kill(p);
    return 0;
  }

  atomic_init(&ctx.cc_found_weight, 0);
  atomic_init(&ctx.cc_exact, false);
  atomic_init(&ctx.stop_flag, false);
  atomic_init(&ctx.rw_stop_flag, false);
  atomic_init(&ctx.rw_steps_started, 0);
  atomic_init(&ctx.rw_steps_completed, 0);
  atomic_init(&ctx.rw_busy_ns, 0);
  atomic_init(&ctx.rw_timed_steps, 0);
  atomic_init(&ctx.cc_next, cc_pack(0, 0));
  atomic_init(&ctx.cc_active_workers, 0);
  atomic_init(&ctx.cc_target_workers, 0);
  atomic_init(&ctx.cc_round_active, 0);
  atomic_init(&ctx.cc_busy_ns, 0);
  atomic_init(&ctx.next_refresh_step, p->refresh > 0 ? p->refresh : 0);

  pthread_mutex_init(&ctx.cw_mutex, NULL);
  pthread_rwlock_init(&ctx.basis_rwlock, NULL);

  ctx.mHT_cc = csr_transpose(NULL, p->spaH);
  ctx.max_col_W = csr_max_row_wght(ctx.mHT_cc);

  if ((p->method & 1) && p->ksub > 0 && ctx.rw_threads > 0) {
    double t_ker0 = get_time_sec();
    ctx.N_global = mzd_nullspace(p->spaH);
    ctx.nu = ctx.N_global ? ctx.N_global->nrows : 0;
    if (p->spaH->rows < ctx.nu || ctx.nu <= 0) {
      if (p->debug & 1) {
        fprintf(stderr,
                "# Warning: ksub=%d requested, but m=%d < nu=%d; "
                "falling back to full-matrix RW (ksub=0)\n",
                p->ksub, p->spaH->rows, ctx.nu);
      }
      safe_mzd_free(ctx.N_global);
      ctx.N_global = NULL;
      ctx.nu = 0;
      p->ksub = 0;
      if (!p->nothrottle) { /* full-matrix RW workers need more memory: re-apply the RW thread limit */
        const int rw_cap = rw_thread_cap(p);
        if (p->method == 1 && num_threads > rw_cap) {
          num_threads = rw_cap;
          ctx.num_threads = rw_cap;
        }
        if (ctx.rw_threads > rw_cap) ctx.rw_threads = rw_cap;
        if (p->debug & 2) {
          fprintf(stderr, "# Note: %d threads, %d of them for RW (full-matrix RW)\n", num_threads, ctx.rw_threads);
        }
      }
    } else if (p->debug & 2) {
      fprintf(stderr,
              "# computed ker(H) in %.4fs: nu=%d, n=%d, ksub=%d (eff=%d)\n",
              get_time_sec() - t_ker0, ctx.nu, p->spaH->cols, p->ksub,
              minint(p->ksub, ctx.nu));
    }
  }

  /* Allocate and launch worker threads */
  ctx.threads = malloc(num_threads * sizeof(pthread_t));
  worker_arg_t *args = malloc(num_threads * sizeof(worker_arg_t));

  for (int i = 0; i < num_threads; i++) {
    args[i].ctx = &ctx;
    args[i].tid = i;
    for (int k = 0; k < MAX_W; k++) {
      args[i].min_swei[k] = p->spaH->rows + 1;
    }
    pthread_create(&ctx.threads[i], NULL, worker_thread_func, &args[i]);
  }

  if (p->method == 1) {
    run_method1_coordinator(&ctx);
  } else if (p->method == 2) {
    run_method2_coordinator(&ctx);
  } else if (p->method == 3) {
    run_method3_coordinator(&ctx);
  } else {
    ERROR("invalid method %d\n", p->method);
  }

  /* Signal stop and wait for all workers */
  atomic_store(&ctx.stop_flag, true);
  for (int i = 0; i < num_threads; i++) {
    pthread_join(ctx.threads[i], NULL);
  }

  int final_dmin = atomic_load(&ctx.dmin);
  int final_dmax = atomic_load(&ctx.dmax);
  int cc_found = atomic_load(&ctx.cc_found_weight);
  /* A CC codeword gives the exact distance only if found in a round where all lower weights had
   * been analyzed (not with noscan=1 unless dmin=wmax is supplied); see `cc_exact`. */
  const bool cc_exact = (cc_found > 0) && atomic_load(&ctx.cc_exact) &&
                        (final_dmax == 0 || cc_found <= final_dmax);

  if (cc_exact) {
    final_dmin = cc_found;
    final_dmax = cc_found;
  } else if (final_dmax > 0 && final_dmin >= final_dmax) {
    final_dmin = final_dmax;
  }

  if (p->wmin > 0 && final_dmax > 0 && final_dmax <= p->wmin) {
    fprintf(stderr, "# early termination due to wmin=%d (cw of weight %d <= wmin found)\n", p->wmin, final_dmax);
  }

  /* Confinement profile output (if smax > 0 and CC was run) */
  if (p->smax && p->method >= 2) {
    int max_w_analyzed = (final_dmin > 1) ? (final_dmin - 1) : ((p->wmax > 0) ? p->wmax : 0);
    if (cc_exact) max_w_analyzed = cc_found;
    if (max_w_analyzed > 0) {
      int global_swei[MAX_W];
      for (int i = 0; i < MAX_W; i++) global_swei[i] = p->spaH->rows + 1;
      for (int t = 0; t < num_threads; t++) {
        for (int i = 1; i <= max_w_analyzed; i++) {
          if (args[t].min_swei[i] < global_swei[i]) {
            global_swei[i] = args[t].min_swei[i];
          }
        }
      }
      int skipped = 0;
      if (p->debug & 1) {
        for (int i = 1; i <= max_w_analyzed; i++) {
          if (global_swei[i] <= p->spaH->rows) {
            fprintf(stderr, "# w=%d min non-zero syndrome weight %d\n", i, global_swei[i]);
          } else {
            skipped = 1;
          }
        }
      } else {
        fprintf(stderr, "# confinement: ");
        for (int i = 1; i <= max_w_analyzed; i++) {
          if (global_swei[i] <= p->spaH->rows) {
            fprintf(stderr, "%d%s", global_swei[i], i < max_w_analyzed ? "," : "");
          } else {
            skipped = 1;
            fprintf(stderr, "?%s", i < max_w_analyzed ? "," : "");
          }
        }
        fprintf(stderr, "\n");
      }
      if (skipped) {
        fprintf(stderr,
                "# Note: Some weights were skipped in confinement profile. "
                "Try increasing smax (current: %d)\n", p->smax);
      }
    }
  }

  long reported_rw_steps = 0;
  if (p->method != 2 && !cc_exact) {
    reported_rw_steps = atomic_load(&ctx.rw_steps_completed);
  }

  if (p->debug & 1) {
    print_codeword_stats(stderr, p);
  }

  /* Output to stdout: dmin dmax rw_steps */
  printf("%d %d %ld\n", final_dmin, final_dmax, reported_rw_steps);
  fflush(stdout);

  /* Codeword export */
  if (p->outC) {
    char comment[256];
    sprintf(comment, "generated by dist_m4ri");
    nzlist_write(p->outC, comment, p);
  }

  if (p->debug & 32) {
    cw_vec_t *cw;
    for (cw = p->codewords; cw != NULL; cw = (cw_vec_t *)(cw->hh.next)) {
      fprintf(stderr, "# cw: [ ");
      for (int i = 0; i < cw->weight; i++) fprintf(stderr, "%d ", 1 + cw->arr[i]);
      fprintf(stderr, "] cnt=%d\n", cw->cnt);
    }
  }

  /* Cleanup */
  safe_mzd_free(ctx.N_global);
  csr_free(ctx.mHT_cc);
  free(ctx.threads);
  free(args);
  pthread_rwlock_destroy(&ctx.basis_rwlock);
  pthread_mutex_destroy(&ctx.cw_mutex);

  var_kill(p);

  return 0;
}
