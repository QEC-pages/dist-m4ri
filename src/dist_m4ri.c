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
  double timeout;
  double start_time;
  int dexp;

  /* Distance bounds & stop flags (cache-line isolated for read-mostly access) */
  _Alignas(64) atomic_int dmin; /* dmin-1 is max cluster size analyzed without success */
  atomic_int dmax;              /* smallest weight codeword found (0 if none) */
  atomic_int cc_found_weight;   /* weight of codeword if CC found exact */
  atomic_bool stop_flag;        /* signals all threads to terminate */
  atomic_bool rw_stop_flag;     /* signals RW workers to stop (in method 3) */

  /* RW state (cache-line isolated from CC and read-mostly bounds) */
  _Alignas(64) long total_rw_steps;
  atomic_long rw_steps_started;
  atomic_long rw_steps_completed;

  /* CC state for current weight (cache-line isolated) */
  _Alignas(64) atomic_int cc_weight;
  atomic_int cc_col_next;
  int cc_col_beg;
  int cc_col_end;
  csr_t *mHT_cc;
  int max_col_W;
  atomic_int cc_active_workers;
  atomic_int cc_target_workers;
  atomic_int cc_round_active;

  /* Codeword synchronization */
  pthread_mutex_t cw_mutex;

  /* Subspace sketching (ksub) state */
  mzd_t *N_global;
  int nu;
  pthread_rwlock_t basis_rwlock;
  atomic_long next_refresh_step;

  /* Timing stats */
  double cc_time_per_weight[MAX_W];
  double avg_rw_step_time;

  /* Thread handles */
  pthread_t *threads;
} distfork_ctx_t;

typedef struct {
  distfork_ctx_t *ctx;
  int tid;
  int min_swei[MAX_W];
} worker_arg_t;

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
  const int col_min = urr->vec[0];

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

/* Run RW batch */
static void run_rw_steps(distfork_ctx_t *ctx, int n_steps,
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
    for (int i = 0; i < nvar && rank < eff_nrows; i++) {
      int col = perm->values[i];
      if (gauss_one_rows(mH, col, rank, eff_nrows)) {
        pivs->values[rank++] = col;
        piv_mask[col >> 6] |= (word)1 << (col & 63);
      }
    }
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
        if ((p->outC || p->maxC || p->dW > 0 || p->min_hits > 0) && p->dW >= 0) {
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
                           : (ctx->num_threads - atomic_load(&ctx->cc_target_workers));
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
              fprintf(stderr,
                      "# RW convergence reached: >= %d min-weight (w=%d) cws "
                      "(cov=%d/%d) each hit >= %d times\n",
                      req, best, p->min_w_cov, p->min_w_cws, p->min_hits);
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
  }
}

/* Run compact subspace RW batch (ksub > 0) */
static void run_rw_steps_ksub(distfork_ctx_t *ctx, int n_steps,
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

  if (nu <= 0 || ksub <= 0 || !ctx->N_global) return;

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
        if ((p->outC || p->maxC || p->dW > 0 || p->min_hits > 0) && p->dW >= 0) {
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
                           : (ctx->num_threads - atomic_load(&ctx->cc_target_workers));
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
              fprintf(stderr,
                      "# RW convergence reached: >= %d min-weight (w=%d) cws "
                      "(cov=%d/%d) each hit >= %d times\n",
                      req, best, p->min_w_cov, p->min_w_cws, p->min_hits);
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
  }

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
}

/* Worker thread main loop */
static void *worker_thread_func(void *arg) {
  worker_arg_t *warg = (worker_arg_t *)arg;
  distfork_ctx_t *ctx = warg->ctx;
  int tid = warg->tid;
  const int nvar = ctx->p->spaH->cols;
  const bool enable_rw = (ctx->p->method & 1) != 0;
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
      if (atomic_load_explicit(&ctx->rw_stop_flag, memory_order_relaxed)) {
        target = ctx->num_threads;
      }
      if (active < target &&
          atomic_load_explicit(&ctx->cc_col_next, memory_order_relaxed) <= ctx->cc_col_end) {
        if (atomic_compare_exchange_weak(&ctx->cc_active_workers, &active, active + 1)) {
          int w = atomic_load_explicit(&ctx->cc_weight, memory_order_relaxed);
          int end = ctx->cc_col_end;
          while (!atomic_load_explicit(&ctx->stop_flag, memory_order_relaxed)) {
            if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
              atomic_store(&ctx->stop_flag, true);
              break;
            }
            int col = atomic_fetch_add_explicit(&ctx->cc_col_next, 1, memory_order_relaxed);
            if (col > end) break;

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
                  atomic_store(&ctx->dmin, 1);
                  atomic_store(&ctx->dmax, 1);
                  atomic_store(&ctx->stop_flag, true);
                  pthread_mutex_unlock(&ctx->cw_mutex);
                }
              }
            }
            err->wei = urr->wei = 0;
          }
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
        /* Prevent thread starvation: ensure chunk size doesn't monopolize steps across threads */
        if (ctx->num_threads > 1 && ctx->total_rw_steps > 0) {
          long max_chunk = (ctx->total_rw_steps + ctx->num_threads - 1) / ctx->num_threads;
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

/* Method 2 coordinator */
static void run_method2_coordinator(distfork_ctx_t *ctx) {
  const int nvar = ctx->p->spaH->cols;
  const int wmax = ctx->p->wmax;
  const int w_start = ctx->p->noscan ? wmax : (ctx->p->dmin > 1 ? ctx->p->dmin : 1);

  if (ctx->p->debug & 2) {
    fprintf(stderr, "# running method=2 (multithreaded CC) with %d threads, w_start=%d wmax=%d\n",
            ctx->num_threads, w_start, wmax);
  }

  int w_limit = wmax;
  if (ctx->p->dmax > 0) {
    if (ctx->p->outC && ctx->p->dW > 0) {
      w_limit = minint(wmax > 0 ? wmax : ctx->p->dmax + ctx->p->dW, ctx->p->dmax + ctx->p->dW);
    } else {
      w_limit = minint(wmax > 0 ? wmax : ctx->p->dmax, ctx->p->dmax);
    }
  }

  for (int w = w_start; w <= w_limit; w++) {
    if (atomic_load(&ctx->stop_flag)) break;
    double now = get_time_sec();
    double remaining_time = ctx->timeout - (now - ctx->start_time);
    if (ctx->timeout > 0.0 && remaining_time <= 0.0) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    /* Estimate CC time for weight w if timeout > 0 */
    if (ctx->timeout > 0.0 && w > w_start) {
      double prev = ctx->cc_time_per_weight[w - 1];
      if (prev <= 0.0001) prev = 0.001;
      double growth = 4.0;
      if (w >= 3 && ctx->cc_time_per_weight[w - 2] > 0.0001) {
        growth = ctx->cc_time_per_weight[w - 1] / ctx->cc_time_per_weight[w - 2];
        if (growth < 2.0) growth = 2.0;
        if (growth > 10.0) growth = 10.0;
      }
      double t_cc_est = prev * growth;
      if ((t_cc_est / ctx->num_threads) > remaining_time * 1.5) {
        if (ctx->p->debug & 1) {
          fprintf(stderr, "# CC for w=%d (est %.2fs) exceeds remaining timeout %.2fs, terminating early (dmin=%d)\n",
                  w, t_cc_est / ctx->num_threads, remaining_time, atomic_load(&ctx->dmin));
        }
        atomic_store(&ctx->stop_flag, true);
        break;
      }
    }

    int beg = (ctx->p->cbeg >= 0) ? ctx->p->cbeg : 0;
    int end = (ctx->p->cend >= 0) ? minint(ctx->p->cend, nvar - w) : (nvar - w);

    atomic_store(&ctx->cc_weight, w);
    ctx->cc_col_beg = beg;
    ctx->cc_col_end = end;
    atomic_store(&ctx->cc_col_next, beg);
    atomic_store(&ctx->cc_target_workers, ctx->num_threads);
    atomic_store(&ctx->cc_round_active, 1);

    double cc_start = get_time_sec();

    if (ctx->p->debug & 2) {
      fprintf(stderr, "# searching w=%d with %d CC threads, columns [%d, %d]\n",
              w, ctx->num_threads, beg, end);
    }

    bool round_completed = false;
    while (!atomic_load(&ctx->stop_flag)) {
      if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
        atomic_store(&ctx->stop_flag, true);
        break;
      }
      if (atomic_load(&ctx->cc_col_next) > end && atomic_load(&ctx->cc_active_workers) == 0) {
        round_completed = true;
        break;
      }
      usleep(100);
    }

    atomic_store(&ctx->cc_round_active, 0);

    double cc_dur = get_time_sec() - cc_start;
    if (w < MAX_W) {
      ctx->cc_time_per_weight[w] = cc_dur;
    }

    int cw_found = atomic_load(&ctx->cc_found_weight);
    if (cw_found > 0) {
      atomic_store(&ctx->dmin, cw_found);
      atomic_store(&ctx->dmax, cw_found);

      if (ctx->p->outC && ctx->p->dW > 0 && w < minint(wmax, cw_found + ctx->p->dW)) {
        w_limit = minint(wmax, cw_found + ctx->p->dW);
        if (ctx->p->debug & 1) {
          if (w == cw_found) {
            fprintf(stderr,
                    "# CC round w=%d finished in %.3fs (%d CC threads): found min-weight codewords "
                    "(dmin=%d, continuing up to w=%d for dW=%d, total %lld cws)\n",
                    w, cc_dur, ctx->num_threads, cw_found, w_limit, ctx->p->dW, ctx->p->num_cws);
          } else if (round_completed) {
            fprintf(stderr,
                    "# CC round w=%d finished in %.3fs (%d CC threads): extra dW round completed "
                    "(dmin=%d, total %lld cws)\n",
                    w, cc_dur, ctx->num_threads, cw_found, ctx->p->num_cws);
          }
        }
      } else {
        if (ctx->p->debug & 1) {
          if (w > cw_found) {
            if (round_completed) {
              fprintf(stderr,
                      "# CC round w=%d finished in %.3fs (%d CC threads): extra dW round completed "
                      "(dmin=%d, total %lld cws)\n",
                      w, cc_dur, ctx->num_threads, cw_found, ctx->p->num_cws);
            }
          } else {
            fprintf(stderr, "# CC found min-weight codeword: d=%d (using %d CC threads, total %lld cws)\n",
                    cw_found, ctx->num_threads, ctx->p->num_cws);
          }
        }
        if (w >= w_limit || !round_completed) {
          atomic_store(&ctx->stop_flag, true);
          break;
        }
      }
    } else {
      if (!round_completed) {
        break;
      }
      /* Weight w analyzed without success */
      atomic_store(&ctx->dmin, w + 1);
      if (ctx->p->debug & 1) {
        fprintf(stderr, "# CC w=%d completed in %.3fs (%d CC threads): no codewords found -> dmin=%d\n",
                w, cc_dur, ctx->num_threads, w + 1);
      }
    }
  }
}

/* Method 3 coordinator */
static void run_method3_coordinator(distfork_ctx_t *ctx) {
  const int nvar = ctx->p->spaH->cols;
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

  /* Initial RW probe to measure average step time */
  double t_rw_start = get_time_sec();
  usleep(2000);
  double t_rw_dur = get_time_sec() - t_rw_start;
  long initial_steps = atomic_load(&ctx->rw_steps_completed);
  if (initial_steps > 0) {
    ctx->avg_rw_step_time = (t_rw_dur * (double)ctx->num_threads) / (double)initial_steps;
  } else {
    ctx->avg_rw_step_time = 0.00005;
  }

  while (!atomic_load(&ctx->stop_flag)) {
    double now = get_time_sec();
    double remaining_time = (ctx->timeout > 0.0) ? (ctx->timeout - (now - ctx->start_time)) : 1e9;
    if (ctx->timeout > 0.0 && remaining_time <= 0.0) {
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    int cur_dmax = atomic_load(&ctx->dmax);
    int cur_dmin = atomic_load(&ctx->dmin);

    /* Target cluster size for CC */
    int target_cc_w;
    if (cur_dmax > 0) {
      if (ctx->p->outC && (ctx->p->dW > 0 || cur_dmin >= cur_dmax)) {
        target_cc_w = cur_dmax + (ctx->p->dW > 0 ? ctx->p->dW : 0);
      } else {
        target_cc_w = cur_dmax - 1;
      }
    } else if (ctx->dexp > 0) {
      target_cc_w = ctx->dexp;
    } else if (ctx->p->wmax > 0) {
      target_cc_w = ctx->p->wmax;
    } else {
      target_cc_w = nvar;
    }

    int max_allowed_w = (ctx->p->wmax > 0)
                        ? minint(ctx->p->wmax, MAX_W - 2) : (MAX_W - 2);
    if (target_cc_w > max_allowed_w) {
      target_cc_w = max_allowed_w;
    }

    if (cur_dmax > 0 && cur_dmin >= cur_dmax && w > target_cc_w) {
      /* Bracketing converged and all requested dW rounds completed */
      atomic_store(&ctx->dmin, cur_dmax);
      atomic_store(&ctx->stop_flag, true);
      break;
    }

    if (w > target_cc_w) {
      /* Let remaining RW steps finish */
      while (!atomic_load(&ctx->stop_flag) && !atomic_load(&ctx->rw_stop_flag)) {
        if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) break;
        if (atomic_load(&ctx->rw_steps_completed) >= ctx->total_rw_steps) break;
        usleep(1000);
      }
      break;
    }

    /* Estimate CC time for weight w */
    double t_cc_est;
    if (w == 1) {
      t_cc_est = 0.0001;
    } else if (w == 2) {
      t_cc_est = 0.001;
    } else {
      double prev = ctx->cc_time_per_weight[w - 1];
      if (prev <= 0.0001) prev = 0.001;
      double growth = 4.0;
      if (w >= 3 && ctx->cc_time_per_weight[w - 2] > 0.0001) {
        growth = ctx->cc_time_per_weight[w - 1] / ctx->cc_time_per_weight[w - 2];
        if (growth < 2.0) growth = 2.0;
        if (growth > 10.0) growth = 10.0;
      }
      t_cc_est = prev * growth;
    }

    /* Check if CC for weight w can finish within timeout */
    if (ctx->timeout > 0.0 && (t_cc_est / ctx->num_threads > remaining_time * 1.5)) {
      if (ctx->p->debug & 2) {
        fprintf(stderr, "# CC for w=%d (est %.2fs) exceeds remaining timeout %.2fs, devoting %d threads to RW\n",
                w, t_cc_est / ctx->num_threads, remaining_time, ctx->num_threads);
      }
      while (!atomic_load(&ctx->stop_flag) && !atomic_load(&ctx->rw_stop_flag)) {
        if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) break;
        if (atomic_load(&ctx->rw_steps_completed) >= ctx->total_rw_steps) break;
        usleep(1000);
      }
      break;
    }

    /* Calculate thread balancing */
    long steps_done = atomic_load(&ctx->rw_steps_completed);
    long steps_rem = (!atomic_load(&ctx->rw_stop_flag) && ctx->total_rw_steps > steps_done)
                     ? (ctx->total_rw_steps - steps_done) : 0;

    int n_cc;
    if (ctx->num_threads == 1) {
      n_cc = 1;
    } else if (steps_rem == 0) {
      n_cc = ctx->num_threads;
    } else if (t_cc_est < 0.005) {
      n_cc = (ctx->num_threads >= 4) ? 2 : 1;
    } else {
      long eff_steps_rem = (cur_dmax > 0 && steps_rem > 2000) ? 2000 : steps_rem;
      double t_rw_total_1t = (double)eff_steps_rem * ctx->avg_rw_step_time;
      double t_cc_total_1t = t_cc_est;
      double est_accum = t_cc_est;
      for (int k = w + 1; k <= target_cc_w && k <= w + 2; k++) {
        est_accum *= 4.0;
        t_cc_total_1t += est_accum;
      }
      double ratio = t_cc_total_1t / (t_cc_total_1t + t_rw_total_1t);
      n_cc = (int)round((double)ctx->num_threads * ratio);
      if (n_cc < 1) n_cc = 1;
      if (n_cc >= ctx->num_threads && steps_rem > 0) n_cc = ctx->num_threads - 1;
    }

    int n_rw = ctx->num_threads - n_cc;

    int beg = (ctx->p->cbeg >= 0) ? ctx->p->cbeg : 0;
    int end = (ctx->p->cend >= 0) ? minint(ctx->p->cend, nvar - w) : (nvar - w);

    atomic_store(&ctx->cc_weight, w);
    ctx->cc_col_beg = beg;
    ctx->cc_col_end = end;
    atomic_store(&ctx->cc_col_next, beg);
    atomic_store(&ctx->cc_target_workers, n_cc);
    atomic_store(&ctx->cc_round_active, 1);

    if (ctx->p->debug & 2) {
      fprintf(stderr,
              "# CC round w=%d started: %d CC threads, %d RW threads "
              "(bounds [%d, %d], rem_rw=%ld, rem_time=%.2fs)\n",
              w, n_cc, n_rw, cur_dmin, cur_dmax, steps_rem, remaining_time);
    }

    double cc_start = get_time_sec();
    bool round_completed = false;

    while (!atomic_load(&ctx->stop_flag)) {
      if (ctx->timeout > 0.0 && (get_time_sec() - ctx->start_time >= ctx->timeout)) {
        atomic_store(&ctx->stop_flag, true);
        break;
      }
      if (atomic_load(&ctx->cc_col_next) > end && atomic_load(&ctx->cc_active_workers) == 0) {
        round_completed = true;
        break;
      }
      usleep(100);
    }

    atomic_store(&ctx->cc_round_active, 0);

    double cc_dur = get_time_sec() - cc_start;
    if (w < MAX_W) {
      ctx->cc_time_per_weight[w] = cc_dur * (double)n_cc;
    }

    int cw_found = atomic_load(&ctx->cc_found_weight);
    if (cw_found > 0) {
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

  if (!p->nothrottle) {
    /* Thread throttling for small codes and large-memory matrices */
    int nvar = p->nvar;
    int nrows = p->spaH ? p->spaH->rows : 0;
    unsigned long long n_elements = (unsigned long long)nrows * (unsigned long long)nvar;

    /* 1. Small-code throttling: avoid thread spawn/join and lock overhead */
    if (nvar < 100 || n_elements < 100000ULL) {
      if (num_threads > 4) num_threads = 4;
    } else if (nvar < 300 || n_elements < 500000ULL) {
      if (num_threads > 16) num_threads = 16;
    }

    /* 2. Large-matrix memory throttling: prevent DRAM bus & L3 cache thrashing */
    int eff_rows = (p->ksub > 0 && p->ksub < nrows) ? p->ksub : nrows;
    if (eff_rows > 0 && nvar > 0) {
      size_t words_H = ((size_t)nvar + 63) / 64;
      size_t bytes_H = (size_t)eff_rows * words_H * 8;
      size_t words_HT = ((size_t)eff_rows + 63) / 64;
      size_t bytes_HT = (size_t)nvar * words_HT * 8;
      size_t dense_bytes_per_thread = bytes_H + bytes_HT;

      /* If dense memory per thread exceeds 15 MB, cap total memory at ~1.5 GB */
      if (dense_bytes_per_thread > 15ULL * 1024 * 1024) {
        int max_mem_threads = (int)((1536ULL * 1024 * 1024) / dense_bytes_per_thread);
        if (max_mem_threads < 2) max_mem_threads = 2;
        if (max_mem_threads > 32) max_mem_threads = 32;
        if (num_threads > max_mem_threads) num_threads = max_mem_threads;
      }
    }

    /* 3. Steps throttling for RW: don't spawn more threads than needed for small step counts */
    if ((p->method & 1) && p->steps > 0) {
      int max_step_threads = (p->steps + 9) / 10;
      if (max_step_threads < 1) max_step_threads = 1;
      if (num_threads > max_step_threads) num_threads = max_step_threads;
    }

    if (requested_threads > 0 && num_threads != requested_threads && (p->debug & 2)) {
      fprintf(stderr, "# Note: throttled threads from %d to %d (n=%d, r=%d)\n",
              requested_threads, num_threads, nvar, nrows);
    }
  } else if (p->debug & 2) {
    fprintf(stderr, "# Thread throttling disabled (nothrottle=1); running with %d threads\n", num_threads);
  }

  double timeout = (p->timeout > 0.0) ? p->timeout : 0.0;

  distfork_ctx_t ctx;
  memset(&ctx, 0, sizeof(ctx));
  ctx.p = p;
  ctx.num_threads = num_threads;
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
  atomic_init(&ctx.stop_flag, false);
  atomic_init(&ctx.rw_stop_flag, false);
  atomic_init(&ctx.rw_steps_started, 0);
  atomic_init(&ctx.rw_steps_completed, 0);
  atomic_init(&ctx.cc_weight, 1);
  atomic_init(&ctx.cc_col_next, 0);
  atomic_init(&ctx.cc_active_workers, 0);
  atomic_init(&ctx.cc_target_workers, 0);
  atomic_init(&ctx.cc_round_active, 0);
  atomic_init(&ctx.next_refresh_step, p->refresh > 0 ? p->refresh : 0);

  pthread_mutex_init(&ctx.cw_mutex, NULL);
  pthread_rwlock_init(&ctx.basis_rwlock, NULL);

  ctx.mHT_cc = csr_transpose(NULL, p->spaH);
  ctx.max_col_W = csr_max_row_wght(ctx.mHT_cc);

  if ((p->method & 1) && p->ksub > 0) {
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

  if (cc_found > 0) {
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
    if (cc_found > 0) max_w_analyzed = cc_found;
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
  if (p->method != 2 && cc_found == 0) {
    reported_rw_steps = atomic_load(&ctx.rw_steps_completed);
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
