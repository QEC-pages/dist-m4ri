#ifndef UTIL_IO_H
/************************************************************************ 
 * qLDPC code input utility routines for distance/decoder package               
 *                                                                      
 * currently: CSS only 
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 ************************************************************************/
#define UTIL_IO_H

#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <limits.h>
#include <m4ri/m4ri.h>

#include "mmio.h"
#include "uthash.h"
#include "util_hash.h"
#include "util_m4ri.h"

#define _maybe_unused __attribute__((unused))

/* Debug output bitmap (`debug=[int]`, see MORE_HELP), sorted by usability (lower bits are more informative).
 * Diagnostic output goes to stderr (the legacy dist_m4ri_old also prints to stdout).  Multiple `debug=` arguments
 * are OR-combined, the first one replaces the default DBG_DEFAULT.  Bits 512, 1024, and 8192..32768 are reserved;
 * bits 65536 and above are reserved for the Python wrapper dist_m4ri.py (ignored by the binaries).
 * (util_hash.h uses the literal value of DBG_LEGACY, as it may be included first.) */
#define DBG_SUMMARY   1    /**< input summary, warnings, stop reason, codeword and hit statistics, export */
#define DBG_PROGRESS  2    /**< run plan, finished CC rounds, new upper bounds, bounds changes, RW convergence */
#define DBG_STATUS    4    /**< periodic status (after 1, 2, 4, ... s, then every 60 s); legacy: every 1000 RW steps */
#define DBG_TIMING    8    /**< thread allocation and timing: throttling, CC round planning and work, waits */
#define DBG_PARAMS    16   /**< code parameters rank(H), rank(L) or rank(G), k (dense elimination at startup) */
#define DBG_CODEWORDS 32   /**< support of each new lightest codeword (RW, CC, finC) */
#define DBG_ARGS      64   /**< command-line arguments, rng seed, files read */
#define DBG_MATRICES  128  /**< dump the matrices H, G, L (if n < 150 unless DBG_LARGE) */
#define DBG_CWDUMP    256  /**< dump all codewords in hash with their hit counts at the end */
#define DBG_LARGE     2048 /**< no size cutoff for the matrices and codeword supports */
#define DBG_LEGACY    4096 /**< dist_m4ri_old only: CC recursion and confinement hash traces */
#define DBG_DEFAULT   (DBG_SUMMARY | DBG_PROGRESS)
#define DBG_CW_MAX    25   /**< without DBG_LARGE, at most this many positions of a codeword support are printed */

//static const int max_row_wt=10; 

#define MAX_W 100 
struct CW_VEC_T;
typedef struct CW_VEC_T cw_vec_t;
typedef struct{
  int debug; /* debug output bitmap, see DBG_SUMMARY etc. (default DBG_DEFAULT) */
  int classical; /* 1 for a classical code, i.e., no `G=Hz` matrix*/
  int css; /* 1: css, 0: non-css -- currently not supported */
  int method; /* bitmap. 1: random window; 2: cluster; 3: both */
  int steps; /* how many RW decoding steps */
  int smax; /** max syndrome weight of interest for `confinement`
		calculation.  When `smax=0` (default), do not
		calculate confinement or use hashing storage.
	     */
  int wmax; /** max cluster size to try for `CC`; */
  int dmin; /** known lower bound on distance (w starts from dmin in CC) */
  int dmax; /** known upper bound on distance (RW ignores codewords of weight >= dmax unless collecting or for
                `min_hits`) */
  int wmin; /** min distance below which we are not interested: once a codeword of weight w <= wmin is found, the
                run ends immediately (not when collecting codewords with outC or maxC) */
  int noscan; /** 1: start CC directly with wmax (no scan over w).  Expert option: the
                  result is a certified lower bound only if `dmin=wmax` is supplied */
  int seed;/* rng seed, set=0 for automatic */
  int dist; /* target distance of the code */
  int dist_max; /* distance actually checked */
  int dist_min; /* distance actually checked */
  int max_row_wgt_H; /* needed for C */
  int max_col_wgt_H; /* needed ? */
  //! int max_row_wt;  /* WARNING: this is defined in `util_io.h` as `static const int` */
  int swei[MAX_W]; /** minimum syndrome weight for each error weight */
  int *start_list; /** sorted list of distinct CC start columns from `start=a,b,c` (NULL: not set).
                       Clusters grown from a listed column are not limited to larger column indices */
  int start_num;   /** number of entries in `start_list` (0: not set, all columns are used) */
  int start_scan;  /** 1: a cluster grown from a listed start column c only uses the columns > c, as in the scan of
                       all columns (orbit exclusion when the listed columns are the smallest ones of the orbits) */
  int cbeg; /** first CC start column for split runs (-1: column 0) */
  int cend; /** last CC start column for split runs (-1: column n-1) */
  //  int linear; /* not supported */
  int n0;  /* code length, =nvar for css, (nvar/2) for non-css */
  int nvar; /* actual n = matrix size */
  int nchk; /* actual k = number of codewords */
  long long int maxC;
  int dW;
  char *finC;
  char *outC;
  int irred; /* 1 (default): only irreducible codewords are recorded, i.e., not the union of a codeword and a disjoint
              * nonzero zero-syndrome vector (CC codewords and RW pairs are filtered, finC codewords are replaced by an
              * irreducible nontrivial part); 0: expert, all codewords as before 0.14.8 */
  cw_vec_t *codewords; /* hash of codewords: the collection window (weight up to min_w + dW with outC or maxC)
                          and the representative set for `min_hits`, see codeword_add_maybe() */
  long long int num_cws; /* number of codewords in hash */
  int min_w;             /* minimum weight of the codewords in hash (INT_MAX: none) */
  int cw_max_w;          /* maximum weight of the codewords in hash (0: none) */
  long long int *cw_cnt_w;  /* number of codewords of each weight in hash (size nvar+2, allocated on first use) */
  long long int *cw_hits_w; /* total hits of the codewords of each weight in hash (same size) */
  char *fdem;
  double pmin;
  int merge;       /* 1: merge the equal columns of a DEM (same detectors and observables), see dem_merge_columns() */
  int nvar_orig;   /* number of columns before merging (and dropping) the columns of a DEM (0: not changed) */
  int *col_orig;   /* column -> its original column (the first one of its class; NULL: not changed) */
  int *col_class;  /* original column -> column (-1: dropped, see dem_drop_detectors(); NULL: not changed) */
  int nvar_merged; /* number of columns after merging the equal columns of a DEM (0: not merged) */
  char *fdrop;     /* file with the detectors to drop from a DEM if this is exact, see dem_drop_detectors() */
  int drop_mode;   /* drop=0: keep all detectors, drop=1 or drop=auto (default): find the detectors to drop
                      (dem_drop_auto()), 2: drop=FILE (fdrop) */
  int drop_set;    /* 1: drop= given explicitly */
  int drop_none;   /* 1: drop=auto found no detector which can be dropped exactly */
  int *drop_list;  /* the detectors dropped (0-based, increasing; ndrop of them), for outD */
  char *outD;      /* file to write the detectors dropped from a DEM to (outD=FILE) */
  int ndrop;       /* number of detectors dropped (0: none, or not verified) */
  int ndrop_rows;  /* number of detectors (rows of H) before dropping */
  int ndrop_cols;  /* number of columns which touch the dropped detectors (replaced by their twins) */
  int drop_A;      /* 1: the twins have other observables (L(twin) = L + A H with A != 0) */
  int sort_rows;   /* 1: sort the rows of H by increasing weight (stable) for CC, see csr_sort_rows_by_weight();
                      -1 (default): 1 for matrices, 0 for a DEM */
  int rows_sorted; /* 1: sort_rows changed the order of the rows of H */
  int cc_prune;    /* the CC search with smax = 0: 0 the walk of the first unsatisfied row (cc_walk(), before 0.15.0),
                      1 (default) the branch and bound (cc_bnb()), 2 the branch and bound with the cover bounds */
  char *finH;
  char *finG;
  char *finL;
  char *fin;
  csr_t *spaH;
  csr_t *spaG;
  csr_t *spaL;
  colmask_t *maskL; /* column masks of L for fast checks L c != 0 (NULL for a classical code), see var_init() */
  int threads; /* number of threads to use (0 for auto) */
  int dexp;    /* expected distance value (0 for auto/none) */
  int dstop;   /* stop target for the lower bound: the run ends once CC certifies dmin >= dstop (0: off); unlike dmax,
                  not an upper bound, and not reported as dmax */
  double timeout; /* timeout in seconds (default 60.0, 0 for infinite) */
  int nothrottle; /* 1: disable automatic thread throttling */
  int chunk_size; /* RW chunk/batch size (0 for auto) */
  int ksub;       /* RW subspace dimension sampled from ker(H) (0 for full H); experimental, should not be used */
  int kwin;       /* RW localized window size W (0 for uniform permutation) */
  int win_mode;   /* RW window mode: 0 = Tanner BFS, 1 = index proximity */
  int win_hyb;    /* window size of the hybrid RW steps (kwin=0, n >= 500: an index window on every other step); 0:
                     256 columns for a DEM, 512 for matrices (see rw_hybrid_window()) */
  int pairs;      /* RW pairs of non-pivot columns (p=2) agreeing on `pairs` random pivot rows (0: off; -1: auto, see
                     rw_pairs_bits()) */
  int round_D;    /* the rounds of a DEM (debug&1, see dem_round_estimate()): detectors per round (0: none found) */
  int round_L;    /* columns per round, n D / m (0: none found) */
  double round_score; /* the score of the estimate (accepted if >= 0.5; -1: not estimated) */
  int min_hits;   /* RW stopping criterion: average hits <n> per lowest-weight cw (0 = off) */
  int cov_cws;    /* RW stopping criterion: size of the representative set of lowest-weight cws (default 100) */
  int refresh;    /* RW steps interval for adaptive basis refresh (0 = off, auto 5000 if ksub>0) */
} params_t;

static inline int minint(const int a, const int b) { return (a < b) ? a : b; }
// #define MININT(a,b) do{ int t1=(a); int t2=(b); t1<t2? t1 :t2; } while(0)

/* The original (1-based in output: + 1) column of the column j: different only after merging the equal columns of a
 * DEM (see dem_merge_columns()) */
static inline int col_orig_index(const params_t * const p, const int j) { return p->col_orig ? p->col_orig[j] : j; }

extern params_t prm;
/**
 * @brief Initialize parameters and load matrices from command line arguments.
 * 
 * Parses command line arguments, sets up the parameter structure,
 * loads matrices from specified files (Matrix Market or DEM), 
 * constructs logical matrices if needed, and performs consistency checks.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 * @param p Pointer to the params_t structure to initialize.
 */
void var_init(int argc, char **argv, params_t * const p);

/**
 * @brief Clean up and free memory allocated in the params_t structure.
 * 
 * Frees sparse matrices (spaH, spaG, spaL) and codeword lists.
 *
 * @param p Pointer to the params_t structure to clean up.
 */
void var_kill(params_t * const p);

/**
 * @brief Read a Detector Error Model (DEM) file and construct H and L matrices.
 * 
 * Parses a DEM file (e.g. from Stim), filters error events based on pmin,
 * and builds the corresponding sparse check matrix H and logical matrix L.
 *
 * @param fnam Path to the DEM file.
 * @param p_spaH Pointer to store the constructed sparse check matrix H.
 * @param p_spaL Pointer to store the constructed sparse logical matrix L.
 * @param pmin Minimum error probability threshold to keep an error event.
 * @param debug Debug print level bitmap.
 */
void read_dem_file(char *fnam, csr_t **p_spaH, csr_t **p_spaL, double pmin, int debug);

/**
 * @brief Read codewords from a .nz list file and add them to the codeword hash.
 * 
 * Reads the file, verifies that each codeword satisfies the code requirements
 * (orthogonal to H, not orthogonal to L for quantum codes), and adds valid ones
 * to the hash table in params_t.
 *
 * @param fnam Path to the .nz file.
 * @param p Pointer to the params_t structure containing the code matrices and hash.
 * @return Number of valid codewords successfully read and added, or -1 on error.
 */
long long int nzlist_read(const char fnam[], params_t *p);

/**
 * @brief Write the collected codewords from the hash table to a .nz file.
 * 
 * Exports the codewords in the collection window, i.e., of weight up to min_w + dW (min_w: the
 * minimum weight found), in NZLIST format.  Heavier codewords kept in the hash only for the
 * `min_hits` statistic are not exported.
 *
 * @param fnam Path to the output .nz file.
 * @param comment An optional comment string to include in the file header.
 * @param p Pointer to the params_t structure containing the codeword hash.
 * @return Number of codewords written, or -1 on error.
 */
long long int nzlist_write(const char fnam[], const char comment[], params_t *p);

/**
 * @brief Add a candidate codeword to the hash table, or count one more hit of a known codeword.
 * 
 * A codeword already in the hash gets one more hit (`cnt`).  A new codeword is added if:
 * - its weight is below the minimum weight min_w (always, also with maxC; min_w is updated);
 * - with outC or maxC (collecting), its weight is at most min_w + max(dW,0), unless maxC
 *   codewords of such weights are collected (see codeword_maxc_reached());
 * - without outC and maxC, its weight equals min_w, unless cov_cws > 0 such codewords are in hash;
 * - with min_hits > 0 and cov_cws > 0, it is heavier but helps to keep a representative set of
 *   cov_cws lowest-weight codewords: the hash holds fewer than cov_cws codewords, or the codeword
 *   is lighter than the heaviest codeword in hash.
 * Heavier codewords are then gradually removed (whole weight classes, starting from the heaviest)
 * as long as the remaining ones keep the collection window and at least cov_cws codewords.
 *
 * @param p Pointer to the params_t structure.
 * @param arr Array of indices representing the support of the codeword (sorted).
 * @param weight Weight of the codeword (length of arr).
 * @param src Where the codeword was found, for the hit counts of each kind of RW step (CW_SRC_RW_UNIFORM,
 *            CW_SRC_RW_WINDOW; CW_SRC_OTHER for CC and finC), see print_rw_step_stats().
 * @return The codeword hash (`p->codewords`).
 */
cw_vec_t * codeword_add_maybe(params_t * const p, const int arr[], int weight, int src);

/** @brief The source of a codeword added to the hash, see codeword_add_maybe(). */
enum { CW_SRC_OTHER = 0, CW_SRC_RW_UNIFORM = 1, CW_SRC_RW_WINDOW = 2 };

/**
 * @brief DBG_SUMMARY: the codewords of the lowest three weights in hash found by the RW steps with a uniform random
 *        permutation and by those with a window (the distinct codewords and their hits), e.g., to see whether some
 *        codewords are found only by the uniform steps (not local in the column order).
 */
void print_rw_step_stats(FILE *stream, const params_t * const p);

/**
 * @brief Window size of the hybrid RW steps (kwin=0, n >= 500: an index window on every other step): win_hyb if
 *        given, otherwise 256 columns for a DEM and 512 for matrices; at most 3n/4.
 *
 * The window columns come first in the elimination order, so that they are pivots unless they depend on earlier
 * window columns.  A window of W consecutive DEM columns (errors close in time) touches only so many detectors: the
 * windows of 192-256 columns find the minimum-weight codewords of circuit DEMs 2-10 times faster than 512, and windows
 * of one round of columns are much slower (most of their columns cannot be pivots).
 */
static inline int rw_hybrid_window(const params_t * const p, const int nvar) {
  const int w = (p->win_hyb > 0) ? p->win_hyb : (p->fdem ? 256 : 512);
  const int cap = (nvar * 3) / 4;
  return (w < cap) ? w : cap;
}

/** @brief The largest number of random pivot rows of the buckets of the RW pairs (`pairs`, a 64-bit key). */
#define RW_PAIRS_MAX_BITS 64

/** @brief The smallest chance to keep a codeword of the weight of the upper bound in the buckets of the RW pairs with
 *         `pairs=-1` (auto), see rw_pairs_bits(). */
#define RW_PAIRS_KEEP 0.8

/** @brief The elimination of a RW step (rank r, n columns) takes about r^2 n / (64 RW_PAIRS_WORK) word operations, the
 *         scale of the work of its RW pairs with `pairs=-1` (auto), see rw_pairs_min_bits().  (Measured with l = 1:
 *         the pairs took 0.5 and 2.5 times the time of the elimination for c1920H and a random LDPC code with
 *         n = 3000, against 0.67 and 2.0 times from this estimate.) */
#define RW_PAIRS_WORK 16.0

/** @brief The number of pairs of a RW step allowed with `pairs=-1` (auto) whatever the work of its elimination (small
 *         matrices), see rw_pairs_min_bits(). */
#define RW_PAIRS_FREE 1048576.0

/**
 * @brief The chance C(r - l, w - 2) / C(r, w - 2) = C(r - w + 2, l) / C(r, l) that a codeword of weight w with exactly
 *        two non-pivot positions avoids l random pivot rows (rank r), see rw_pairs_bits().
 */
static inline double rw_pairs_keep(const int rank, const int w, const int l) {
  double q = 1.0;
  for (int i = 0; i < w - 2; i++) {
    if (rank - l - i <= 0) return 0.0;
    q *= (double)(rank - l - i) / (double)(rank - i);
  }
  return q;
}

/**
 * @brief The fewest random pivot rows l of the buckets of the RW pairs with `pairs=-1` (auto), for a matrix with n
 *        columns and rank r: the n - r non-pivot columns give about (n - r)^2 / 2^(l+1) pairs per RW step, each
 *        compared in at most r/64 words, which should be at most about the work of the elimination of the step,
 *        r n / RW_PAIRS_WORK pairs (or RW_PAIRS_FREE pairs, for small matrices).
 *
 * This matters while the upper bound w is far above the distance, early in a run: the chance to keep a codeword of
 * weight w (rw_pairs_keep()) then asks for l = 1, i.e., all pairs in two buckets, so that a step on a circuit DEM
 * (n - r about 10 r) would check millions (mitten_270, n = 8910) or billions (bb_756, n = 129276) of pairs; with this
 * bound, l >= 6 for these DEMs.  Codes with n - r < r (e.g., c1920H) keep l = 1.
 */
static inline int rw_pairs_min_bits(const int nvar, const int rank) {
  const double k = (double)(nvar - rank);
  double budget = (double)rank * (double)nvar / RW_PAIRS_WORK;
  if (budget < RW_PAIRS_FREE) budget = RW_PAIRS_FREE;
  int l = 1;
  double pairs = k * k / 4.0; /* (n - r)^2 / 2^(l+1) */
  while ((l < RW_PAIRS_MAX_BITS) && (pairs > budget)) {
    l++;
    pairs /= 2.0;
  }
  return l;
}

/**
 * @brief The number l of random pivot rows of the buckets of the RW pairs (p=2) for a matrix with n columns and rank
 *        r, with the upper bound w on the distance (0: no pairs).
 *
 * A RW step with pairs also finds the codewords with exactly two non-pivot columns a, b: the sum of the codewords of
 * a and b (p=1), whose pivot part is the sum of the columns a and b of the reduced matrix.  Only the pairs which agree
 * on l random pivot rows are checked (the non-pivot columns are bucketed by these l bits, as in Stern's algorithm),
 * which keeps a codeword of weight w with the chance rw_pairs_keep().  `pairs` if given (at most r); auto (-1): about
 * one non-pivot column per bucket, l = ceil(log2(n - r)), but fewer rows (at least one) while the chance to keep a
 * codeword of weight w is below RW_PAIRS_KEEP, and not fewer than rw_pairs_min_bits() (the work of the pairs of a
 * step).  No pairs without an upper bound (w = 0).
 */
static inline int rw_pairs_bits(const params_t * const p, const int nvar, const int rank, const int w) {
  const int k = nvar - rank;
  if ((p->pairs == 0) || (k < 2) || (rank < 1) || (w < 2)) return 0;
  int l = p->pairs;
  if (l < 0) {
    l = 1;
    while ((l < RW_PAIRS_MAX_BITS) && (((uint64_t)1 << l) < (uint64_t)k)) l++;
    const int lmin = rw_pairs_min_bits(nvar, rank);
    while ((l > 1) && (l > lmin) && (rw_pairs_keep(rank, w, l) < RW_PAIRS_KEEP)) l--;
  }
  return (l < rank) ? l : rank;
}

/**
 * @brief Number of codewords in the collection window (weight up to min_w + dW with outC or maxC,
 *        otherwise weight min_w): the codewords exported with outC and counted for maxC.
 */
long long int codeword_window_count(const params_t * const p);

/** @brief Returns 1 if maxC > 0 and maxC codewords of the collection window are collected. */
int codeword_maxc_reached(const params_t * const p);

/**
 * @brief Weight limit (exclusive) of RW codewords needed for the `min_hits` statistic.
 *
 * Returns 0 if min_hits <= 0, the weight min_w + 1 if only the minimum-weight codewords are used
 * (cov_cws <= 0), INT_MAX (any weight) while the representative set is not filled, and otherwise
 * the maximum weight in hash plus one.
 */
int codeword_feed_limit(const params_t * const p);

/**
 * @brief Weight limit (exclusive) for RW candidate codewords (the same with any `debug` value).
 *
 * Without an upper bound (cur_dmax <= 0), any weight is of interest.  Otherwise, codewords of
 * weight >= cur_dmax are of no interest, except when collecting (outC or maxC with dW >= 0:
 * weight up to cur_dmax + dW) and for the `min_hits` statistic (`feed`, see codeword_feed_limit()).
 *
 * @param p Pointer to the params_t structure.
 * @param cur_dmax Current upper bound on the distance (0 if none).
 * @param feed Weight limit for the `min_hits` statistic, see codeword_feed_limit().
 * @return Weight limit, at most nvar + 1.
 */
static inline int codeword_rw_limit(const params_t * const p, const int cur_dmax, const int feed) {
  const int lim_max = p->nvar + 1;
  if (cur_dmax <= 0) return lim_max;
  long long int lim = cur_dmax;
  if ((p->outC != NULL || p->maxC > 0) && p->dW >= 0) lim = (long long int)cur_dmax + p->dW + 1;
  if (feed > lim) lim = feed;
  return (lim < lim_max) ? (int)lim : lim_max;
}

/** @brief Hit statistics for the `min_hits` stopping criterion, see codeword_hit_stats(). */
typedef struct {
  long long int min_cws;  /**< number of codewords of the minimum weight min_w */
  long long int min_hits; /**< their total number of hits */
  long long int set_cws;  /**< number of codewords in the representative set (weights w_lo..w_hi) */
  long long int set_hits; /**< their total number of hits */
  int w_lo;               /**< smallest weight in the representative set (min_w) */
  int w_hi;               /**< largest weight in the representative set */
  double avg;             /**< <n>: the smaller of the average hits per codeword of weight min_w and
                               of the representative set (0 if no codewords) */
} cw_hit_stats_t;

/**
 * @brief Compute the hit statistics of the `min_hits` stopping criterion.
 *
 * The representative set consists of the lowest weight classes (starting with min_w) in hash
 * with at least cov_cws codewords in total (all classes if fewer; only min_w if cov_cws <= 0 or
 * min_hits <= 0).  As in QDistRnd, <n> estimates the average number of times a codeword is found;
 * a lighter codeword is missed with probability about exp(-<n>), assuming lighter codewords are
 * found at least as often.  Heavier codewords make the estimate more conservative.
 *
 * @param p Pointer to the params_t structure.
 * @param st Output statistics.
 */
void codeword_hit_stats(const params_t * const p, cw_hit_stats_t * const st);

/**
 * @brief Check whether the QDistRnd-style hit count stopping condition is met.
 *
 * Returns 1 if p->min_hits > 0 and the average number of hits <n> (see codeword_hit_stats())
 * reaches p->min_hits, both for the codewords of weight p->min_w and for the representative set
 * of up to p->cov_cws lowest-weight codewords.  A single codeword (e.g., for a classical code
 * with k=1) suffices.  Otherwise returns 0.
 *
 * @param p Pointer to the params_t structure.
 * @return 1 if convergence criterion is met, 0 otherwise.
 */
int check_min_hits_convergence(const params_t * const p);

/**
 * @brief Check the uniformity of the hit counts within each weight class of the representative set.
 *
 * Codewords of the same weight found with the same probability per RW step have (zero-truncated)
 * Poisson-distributed hit counts.  For each weight class of the representative set with at least
 * 2 codewords, the variance of the hit counts is compared with the Poisson variance.  A warning is
 * printed if it is significantly larger (p-value < 0.001), the spread (standard deviation) of the
 * hit rates is at least their mean, and, for hit rates with this relative spread (gamma
 * distribution), the probability to miss a lighter codeword is at least 3 times larger than
 * exp(-mu) with equal hit rates (mu: average hits): some codewords are found much more often than
 * others of the same weight, and the estimate exp(-<n>) is optimistic.  Takes a single pass over
 * the hash.
 *
 * @param stream Output stream for the warnings (typically stderr).
 * @param p Pointer to the params_t structure.
 * @return Number of weight classes with non-uniform hit counts.
 */
int codeword_hit_check(FILE *stream, const params_t * const p);

/**
 * @brief Information-set estimate of the probability that a RW step finds a given codeword of weight w.
 *
 * Model: uniform random information sets.  A RW step reduces H (rank r) with a random column order;
 * a codeword is found if exactly one of its w positions is among the k = n - r non-pivot columns:
 * P1(w) = w C(n-w, k-1) / C(n, k).  Codewords of weight w > r + 1 are never found (P1 = 0).  With the
 * RW pairs (`pairs_l` > 0, see rw_pairs_bits()), also if exactly two of its positions are non-pivot
 * columns and its other w - 2 positions avoid the pivot rows of the buckets:
 * P2(w) = C(w, 2) C(n-w, k-2) / C(n, k) * C(r - w + 2, l) / C(r, l), and the probability is P1 + P2.
 *
 * @param n Block length (number of columns of H).
 * @param rank Rank r of H.
 * @param w Codeword weight.
 * @param pairs_l The number l of random pivot rows of the buckets of the RW pairs (0: no pairs).
 * @return P1(w) + P2(w).
 */
double rw_infoset_find_prob(const int n, const int rank, const int w, const int pairs_l);

/**
 * @brief Print the information-set estimate of the probability that RW missed a codeword lighter than
 *        the upper bound found.
 *
 * Among the weights w_lo..w_hi (not excluded by CC), the weight w with the smallest P1(w) (see
 * rw_infoset_find_prob()) is used: a single codeword of this weight is missed in `steps` RW steps
 * with uniform random permutations with probability (1 - P1)^steps.  Also prints the total number
 * of RW steps for a 1% miss probability (with the fraction steps / steps_total of uniform steps, as
 * in this run; localized windows, used in about half of the steps for n >= 500, are not counted)
 * and, if `t_step` > 0, the corresponding time.
 *
 * @param stream Output stream (typically stderr).
 * @param n Block length (number of columns of H).
 * @param rank Rank of H.
 * @param w_lo Smallest weight to consider (dmin, at least 1).
 * @param w_hi Largest weight to consider (dmax - 1).
 * @param steps Number of completed RW steps with uniform random permutations.
 * @param steps_total Number of all completed RW steps.
 * @param t_step Wall time per RW step in seconds (0: unknown).
 * @param pairs_l The number l of random pivot rows of the buckets of the RW pairs (0: no pairs).
 */
void print_rw_infoset_estimate(FILE *stream, const int n, const int rank, const int w_lo, const int w_hi,
                               const long steps, const long steps_total, const double t_step, const int pairs_l);

/**
 * @brief Compute hit count statistics (min, max, avg, stdev) for minimum-weight codewords.
 *
 * @param p Pointer to the params_t structure.
 * @param min_cnt Output minimum hit count.
 * @param max_cnt Output maximum hit count.
 * @param avg_cnt Output average hit count.
 * @param stdev_cnt Output standard deviation of hit counts.
 */
void compute_min_w_hit_stats(const params_t * const p, int *min_cnt, int *max_cnt,
                             double *avg_cnt, double *stdev_cnt);

/**
 * @brief Print accumulated codeword and hit statistics to the given stream.
 *
 * One line for the minimum weight (with the `min_hits` statistic), and one line for each heavier weight class in
 * hash; without DBG_LARGE in `p->debug`, the weight classes after the 4 lightest heavier ones share one line.
 *
 * @param stream Output stream (typically stderr).
 * @param p Pointer to the params_t structure.
 */
void print_codeword_stats(FILE *stream, const params_t * const p);

/**
 * @brief Print the support of a codeword (1-based column indices) in one line, for DBG_CODEWORDS.
 *
 * Prints `prefix`, the weight, and the first DBG_CW_MAX positions (all of them with DBG_LARGE).
 *
 * @param stream Output stream (typically stderr).
 * @param prefix Line prefix, e.g., "# RW codeword".
 * @param arr Sorted 0-based column indices.
 * @param weight Number of indices.
 * @param debug Debug bitmap (DBG_LARGE: no size cutoff).
 */
void print_codeword_support(FILE *stream, const char *prefix, const int arr[], const int weight, const int debug);

/**
 * @brief Print a binary CSR matrix row by row ('1': nonzero, '.': zero), for DBG_MATRICES.
 *
 * Nothing is printed without DBG_MATRICES in `debug`; for n >= 150 columns, only the size is printed unless
 * DBG_LARGE is set.
 *
 * @param stream Output stream (typically stderr).
 * @param M The matrix.
 * @param name Name of the matrix, e.g., "H".
 * @param debug Debug bitmap.
 */
void csr_dump(FILE *stream, const csr_t * const M, const char name[], const int debug);

#define DIST_M4RI_VERSION "0.15.1"

/**
 * @brief Print short help message listing all allowed parameters to stderr.
 *
 * @param prog Program name (argv[0]).
 */
void print_short_help(const char *prog);

#define SHORT_HELP \
  "%s (version %s): calculate distance of a classical or quantum CSS code\n" \
  "Usage: %s [method=1|2|3] [parameter=value ...]\n\n" \
  "Allowed parameters:\n" \
  "  method, finH, finG, finL, fin, fdem, pmin, merge, drop, outD, classical, css,\n" \
  "  steps, wmin, wmax, dmin, dmax, dstop, dexp (dest), smax, start, start_scan,\n" \
  "  sort_rows, cc_prune, cbeg, cend, noscan, threads, timeout, nothrottle,\n" \
  "  chunk_size (batch), ksub, kwin (win), win_mode, win_hyb, pairs, min_hits,\n" \
  "  cov_cws, refresh, finC, outC, maxC, dW, irred, seed, debug\n\n" \
  "Help options:\n" \
  "  -h, --help    : display help for commonly used parameters (fits 80 rows; --help TOPIC: as --morehelp TOPIC)\n" \
  "  --morehelp    : display full help for all available parameters (or for TOPIC: --morehelp TOPIC)\n"

#define USAGE \
  "%s (version %s): calculate distance of a classical or quantum CSS code\n" \
  "Usage: %s [method=1|2|3] [parameter=value ...]\n\n" \
  "Calculation method:\n" \
  "  method=[int]       1: Random Window (RW) algorithm (upper bound)\n" \
  "                     2: Connected Cluster (CC) algorithm (lower bound / exact)\n" \
  "                     3: Bracketing mode (concurrent RW and CC) (default: 3)\n\n" \
  "Input matrices (Matrix Market .mmx/.mtx format or Stim DEM):\n" \
  "  finH=[file]        Parity check matrix H (classical) or Hx (CSS quantum)\n" \
  "  finG=[file]        Hz check matrix (quantum CSS code only)\n" \
  "  finL=[file]        Lx logical operator matrix (quantum CSS code only)\n" \
  "                     Note: Either L=Lx or G=Hz is required for quantum CSS codes\n" \
  "  fin=[str]          Base name for CSS matrices (loads ${fin}X.mtx, ${fin}Z.mtx)\n" \
  "  fdem=[file]        Stim detector error model (DEM) file\n" \
  "  pmin=[float]       Minimum error probability threshold to keep for DEM (0.0)\n" \
  "  merge=[0|1]        Merge equal DEM columns (exact; codewords use the first one) (1)\n" \
  "  drop=[auto|file]   DEM detectors to drop if exact: auto=1 (default), 0 (none), or a file\n" \
  "  outD=[file]        Write the DEM detectors dropped (0-based) to file\n" \
  "  classical=[0|1]    1: classical code (Hx only), 0: quantum CSS (auto-detected)\n\n" \
  "Distance bounds and guidance:\n" \
  "  dmin=[int]         Certified lower bound on distance (CC starts from dmin) (1)\n" \
  "  dmax=[int]         Known upper bound on distance (CC up to w=dmax-1 only) (0)\n" \
  "  dstop=[int]        Stop once CC certifies dmin>=dstop (not an upper bound) (0)\n" \
  "  dexp=[int]         Expected distance for method=3 thread allocation (alias: dest) (0)\n\n" \
  "Search limits and stopping criteria:\n" \
  "  steps=[int]        Maximum RW decoding steps / information sets (100000)\n" \
  "  wmax=[int]         Maximum cluster weight to search in CC (0=until bound/timeout)\n" \
  "  wmin=[int]         Stop immediately if cw with weight <= wmin is found (1)\n" \
  "  min_hits=[int]     Stop RW when lowest-wt cws are hit min_hits times on average (5)\n" \
  "  timeout=[sec]      Execution timeout in seconds, 0 for infinite (60.0)\n\n" \
  "Multithreading and RW optimization:\n" \
  "  threads=[int]      Max worker threads to use (0: all logical CPUs) (0)\n" \
  "  ksub=[int]         EXPERIMENTAL, do not use: subspace dimension sampled from ker(H) for RW (0)\n" \
  "  kwin=[int]         Localized column permutation window size W (0: auto/hybrid) (0)\n\n" \
  "Codeword collection:\n" \
  "  outC=[file]        Export found min-weight codewords (up to +dW) to file (.nz format)\n" \
  "  finC=[file]        Import initial codewords from file (.nz format)\n" \
  "  maxC=[int]         Maximum number of codewords to collect (0 for unlimited) (0)\n" \
  "  dW=[int]           Collect codewords up to the minimum weight found + dW (0)\n\n" \
  "Extra parameters (see --morehelp for details):\n" \
  "  smax=[int] (0)         Max syndrome weight for confinement profile (0 to disable)\n" \
  "  noscan=[0|1] (0)       Expert, method 2: CC at w=wmax only (no lower bound!)\n" \
  "  start=[list] (-1)      Expert: CC from listed columns only, e.g. start=0,48\n" \
  "  start_scan=[0|1] (0)   Expert, with start: clusters from column c use only columns > c\n" \
  "  cbeg/cend=[int] (-1)   Expert: CC start column range, for split runs\n" \
  "  irred=[0|1] (1)        Expert: 0 records also reducible codewords (cw + disjoint 0-syndrome vector)\n" \
  "  sort_rows=[0|1] (1)    Sort the rows of H by weight for CC (default 0 with fdem)\n" \
  "  cc_prune=[0|1|2] (1)   CC search: 0 first-row walk, 1 branch and bound, 2 also cover bounds\n" \
  "  nothrottle=[0|1] (0)   Disable thread throttling (also --no-throttle)\n" \
  "  chunk_size=[int] (0)   RW batch chunk size (0: auto, alias: batch)\n" \
  "  win_mode=[0|1] (0)     Window mode: 0=Tanner BFS, 1=index proximity\n" \
  "  win_hyb=[int] (0)      Hybrid RW window size (0: 256 for a DEM, 512 otherwise)\n" \
  "  pairs=[int] (-1)       RW pairs of non-pivot columns (p=2; -1: auto, 0: off)\n" \
  "  cov_cws=[int] (100)    Number of lowest-wt cws used for the min_hits statistic\n" \
  "  refresh=[int] (0)      RW steps between ker(H) basis refreshes (experimental ksub only; 5000 if ksub>0)\n" \
  "  seed=[int] (0)         RNG seed [0 for time(NULL)]\n" \
  "  debug=[int] (3)        Debug bitmask (0: silent, 1: summary, 2: progress, ...)\n\n" \
  "Help options:\n" \
  "  -h, --help         Display this help message (--help TOPIC: the same as --morehelp TOPIC)\n" \
  "  --morehelp         Display full help for all parameters (--morehelp TOPIC: only TOPIC, e.g., maxC or windows)\n" \
  "  --version          Display program version\n"

#define MORE_HELP \
  "%s (version %s): calculate distance of a classical or quantum CSS code\n" \
  "Usage: %s [method=1|2|3] [parameter=value ...]\n\n" \
  "Calculation method:\n" \
  "  method=[int]       Bitmap / identifier for calculation method (default: 3):\n" \
  "                     1: Random Window (RW) algorithm (upper bound)\n" \
  "                        Finds an upper bound on distance by testing random\n" \
  "                        information sets. Fast for finding small codewords.\n" \
  "                     2: Connected Cluster (CC) algorithm (lower bound / exact)\n" \
  "                        Exhaustive cluster search finding a lower bound or exact\n" \
  "                        minimum distance. Guaranteed to find the true code\n" \
  "                        distance if run to completion.\n" \
  "                     3: Bracketing mode (concurrent RW and CC)\n" \
  "                        Dynamically balances CC (lower bound) and RW (upper\n" \
  "                        bound) worker threads based on measured CC and RW times,\n" \
  "                        current bounds [dmin, dmax], remaining timeout, and the\n" \
  "                        distance estimate (dexp).\n\n" \
  "Input matrices and code specification:\n" \
  "  finH=[file]        Parity check matrix H (for classical codes) or Hx (for\n" \
  "                     quantum CSS codes) in Matrix Market (.mmx / .mtx) format.\n" \
  "  finG=[file]        Generator / Hz matrix for quantum CSS codes in Matrix Market\n" \
  "                     format. Used to verify orthogonality and construct logicals.\n" \
  "  finL=[file]        Logical operator matrix Lx for quantum CSS codes in Matrix\n" \
  "                     Market format.\n" \
  "                     Note: For a quantum CSS code, either finL (Lx) or finG (Hz)\n" \
  "                     must be provided.\n" \
  "  fin=[str]          Base name for the matrix files (default: ../examples/try).\n" \
  "                     Automatically looks for \"${fin}X.mtx\" as finH and\n" \
  "                     \"${fin}Z.mtx\" as finG.  A later finH, finG, or finL\n" \
  "                     replaces fin (with a warning); fin after them is an error.\n" \
  "  fdem=[file]        Detector Error Model file (.dem) generated by Stim.\n" \
  "                     Automatically constructs H and L matrices. Cannot be\n" \
  "                     combined with finH, finG, finL, or fin.\n" \
  "  pmin=[float]       Minimum error probability threshold for DEM parsing (0.0).\n" \
  "                     Error mechanisms with probability < pmin are ignored.\n" \
  "                     Only valid when fdem is specified. The column indices\n" \
  "                     (start, finC, outC) count only the error mechanisms kept:\n" \
  "                     they change with pmin.\n" \
  "  merge=[0|1]        Merge the equal columns of a DEM (error mechanisms with the\n" \
  "                     same detectors and observables) when reading it (default:\n" \
  "                     1; ignored without fdem). Exact: dropping two equal columns\n" \
  "                     of a codeword gives a lighter one, so the distance does not\n" \
  "                     change. Column indices (outC, finC, start, cbeg, cend, and\n" \
  "                     the printed codewords) are those of the DEM; the codewords\n" \
  "                     found use the first one of each set of equal columns (use\n" \
  "                     merge=0 to enumerate the codewords with each of them, or\n" \
  "                     expand=1 of dist_m4ri.py to expand the list found).\n" \
  "  drop=[auto|file]   Detectors of a DEM to drop if this is verified to be exact\n" \
  "                     (default: auto; also drop=1; drop=0 keeps all detectors;\n" \
  "                     only with fdem; not with classical=1): every column\n" \
  "                     touching them has a twin, a column with the same other\n" \
  "                     detectors which touches none of them (no column if it has\n" \
  "                     no other detector), with observables L + A H for one\n" \
  "                     matrix A.  Then the twins replace the columns touching the\n" \
  "                     dropped detectors in any codeword, without increasing its\n" \
  "                     weight, and the search runs on the other detectors and the\n" \
  "                     columns which touch none of them: every codeword found\n" \
  "                     (outC) is a codeword of the full DEM (using only these\n" \
  "                     columns), and the distance is exact.  finC codewords are\n" \
  "                     mapped to the twins.\n" \
  "                     auto: a search for the largest such set (with A = 0) in the\n" \
  "                     DEM, e.g., the minority-basis detectors of a CSS circuit\n" \
  "                     (with the usual noise): a column which contains no other\n" \
  "                     column with the same observables is dropped all or nothing,\n" \
  "                     and the detectors of the other columns follow by propagation\n" \
  "                     (see README.md).\n" \
  "                     file: the detectors listed (0-based, k for Dk, separated by\n" \
  "                     white space or commas; lines starting with '#' or '%%' are\n" \
  "                     comments; write ./1 for a file named 1); if dropping them is\n" \
  "                     not exact, a note is printed and the full DEM is used.\n" \
  "  outD=[file]        Write the detectors dropped from the DEM (0-based, in the\n" \
  "                     format of a drop file; only the comment line if none) to\n" \
  "                     file.  Only valid with fdem.\n" \
  "  classical=[0|1]    Code type override:\n" \
  "                     1: Treat as a classical linear code (Hx only; ignores\n" \
  "                        or discards logical matrix L).\n" \
  "                     0: Treat as a quantum CSS code (requires L or G matrix).\n" \
  "                     Default: 1 if only finH is given; 0 if finG, finL, fin,\n" \
  "                     or fdem is provided.\n" \
  "  css=[int]          Reserved for future non-CSS quantum code support (1).\n\n" \
  "Distance bounds and search guidance:\n" \
  "  dmin=[int]         Known certified lower bound on distance (default: 1).\n" \
  "                     In CC (method 2/3), cluster search begins at w = dmin.\n" \
  "  dmax=[int]         Known upper bound on distance (default: 0).\n" \
  "                     CC (method 2/3) ends once dmin = dmax (with outC, after the\n" \
  "                     rounds w = dmax..dmax+dW).  In RW (method 1/3), codewords\n" \
  "                     of weight >= dmax are ignored unless collecting codewords\n" \
  "                     or needed for min_hits.\n" \
  "  dstop=[int]        Stop target for the lower bound (default: 0 = off): the run\n" \
  "                     (method 2/3) ends once CC has certified dmin >= dstop, that\n" \
  "                     is, that there is no codeword of weight < dstop (with outC,\n" \
  "                     after the rounds w = dstop..dstop+dW), unless a lighter\n" \
  "                     codeword is found (then the search goes on as usual).\n" \
  "                     Unlike dmax, dstop is not an upper bound and is not\n" \
  "                     reported as dmax: the output is 'dmin w steps', where w is\n" \
  "                     the weight of the lightest codeword found (0: none).  E.g.,\n" \
  "                     dstop=U wmin=U-1 decides whether d >= U, and for a CSS code\n" \
  "                     with d = min(dX, dZ) <= U, dstop=U is all that is needed.\n" \
  "                     With method=1, it only applies to a supplied dmin >= dstop.\n" \
  "  dexp=[int]         Expected code distance (alias: dest) (default: 0).\n" \
  "                     Hint for method=3 (bracketing): as long as RW has found no\n" \
  "                     codeword, CC rounds at w > dexp run only on threads which\n" \
  "                     cannot run RW (CC pauses if RW can use all threads); CC\n" \
  "                     resumes once RW finds a codeword or ends.\n\n" \
  "Search limits and Connected Cluster (CC) options:\n" \
  "  steps=[int]        Maximum number of RW decoding steps / information sets\n" \
  "                     (default: 100000). Ignored in method=2.\n" \
  "  wmin=[int]         Minimum distance threshold (default: 1).\n" \
  "                     If a codeword of weight w <= wmin is found, execution\n" \
  "                     terminates immediately (not when collecting codewords with\n" \
  "                     outC or maxC). Useful for screening codes.\n" \
  "  wmax=[int]         Maximum cluster weight to analyze in CC (default: 0).\n" \
  "                     In method=2, CC terminates after checking weight wmax.\n" \
  "                     0 means continue until a codeword is found, the bounds\n" \
  "                     meet, or timeout (method=2 then needs timeout>0 or dmax>0).\n" \
  "  smax=[int]         Maximum syndrome weight for confinement profile (default: 0).\n" \
  "                     When smax > 0, tracks minimum syndrome weights for each\n" \
  "                     error weight. When smax=0, confinement is not computed.\n" \
  "  sort_rows=[0|1]    Sort the rows of H by increasing weight (stable) after\n" \
  "                     reading the input (default: 1 for matrices, finH or fin;\n" \
  "                     0 for a DEM, fdem).  With cc_prune=0, CC extends a cluster\n" \
  "                     by the columns of the first row of its syndrome (in row\n" \
  "                     order), so that light rows first reduce the CC work for\n" \
  "                     matrices with very different row weights, e.g., 1.7-4x\n" \
  "                     with the hook columns of an extended CSS distance, while\n" \
  "                     the detectors of a DEM are in time order, which is better\n" \
  "                     (sorting them made CC 1.3-4x slower).  With cc_prune=1 or\n" \
  "                     2, the order of the rows only breaks ties.  The distance\n" \
  "                     and RW do not depend on the order of the rows; CC finds the\n" \
  "                     same codewords (outC), possibly in another order (with\n" \
  "                     cc_prune=0, also with other hit counts).\n" \
  "  cc_prune=[0|1|2]   The CC search with smax=0 (default: 1).  0: the walk of\n" \
  "                     the versions before 0.15.0: a cluster is extended by the\n" \
  "                     columns of the first unsatisfied row of its syndrome, and\n" \
  "                     a branch is pruned if its syndrome weight exceeds what the\n" \
  "                     columns left can cancel.  1: branch and bound: it branches\n" \
  "                     on the unsatisfied row with the fewest available columns,\n" \
  "                     the columns of the earlier branches are forbidden in the\n" \
  "                     later ones (each codeword is found once: one CC hit), and\n" \
  "                     a branch is pruned by lower bounds on the number of columns\n" \
  "                     needed (an unsatisfied row without available columns, or a\n" \
  "                     packing of unsatisfied rows with disjoint sets of available\n" \
  "                     columns).  2: also the coverage and fractional cover bounds\n" \
  "                     (fewer nodes, more work per node).  All values find the\n" \
  "                     same codewords and distances (doc/cc_completeness.tex);\n" \
  "                     with the expert start list (and start_scan=0), cc_prune=1\n" \
  "                     and 2 forbid the columns listed before the start column.\n" \
  "  noscan=[0|1]       Expert option: 1: run CC only at weight w = wmax, skipping\n" \
  "                     the scan over weights w < wmax (default: 0). Only valid for\n" \
  "                     method=2. WARNING: lower weights are not checked, so a\n" \
  "                     codeword found only sets dmax, and dmin is not raised above\n" \
  "                     the supplied dmin. Certified only if dmin=wmax is given.\n" \
  "  start=[list]       Expert option: comma-separated list of CC start columns,\n" \
  "                     e.g., start=0,48,96 (default: -1 = all columns). Clusters\n" \
  "                     grown from a listed column are NOT limited to larger column\n" \
  "                     indices (unless start_scan=1). Use, e.g., one column per\n" \
  "                     block of a quasi-cyclic code. WARNING: dmin is valid only\n" \
  "                     if a code symmetry maps every minimum-weight codeword to one\n" \
  "                     containing a listed column. Cannot be combined with\n" \
  "                     cbeg/cend.\n" \
  "  start_scan=[0|1]   Expert option, with start (default: 0).  1: a cluster grown\n" \
  "                     from a listed column c only uses the columns > c, as in the\n" \
  "                     scan of all columns (see cbeg).  WARNING: dmin is valid\n" \
  "                     only if a code symmetry maps every minimum-weight codeword\n" \
  "                     to one whose smallest column is listed, e.g., if the listed\n" \
  "                     columns are the smallest columns of the orbits of a group of\n" \
  "                     code automorphisms, one per orbit (start=0,l for the two\n" \
  "                     blocks of a two-block group-algebra code of length 2l over\n" \
  "                     an abelian group).  The later start columns then need less\n" \
  "                     work, the least if each orbit is a range of consecutive\n" \
  "                     columns.  Not with fdem.\n" \
  "  cbeg=[int]         Expert option for split runs: first CC start column\n" \
  "                     (default: -1 = column 0). A cluster started at column i\n" \
  "                     only uses columns > i, so every codeword is found starting\n" \
  "                     from its smallest column.\n" \
  "  cend=[int]         Expert option for split runs: last CC start column\n" \
  "                     (default: -1 = column n-1). WARNING: with cbeg/cend, dmin\n" \
  "                     only covers codewords whose smallest column is in\n" \
  "                     [cbeg,cend]; combine runs covering all columns 0..n-1 by\n" \
  "                     taking the minimum of dmin (and of dmax) over the runs.\n\n" \
  "Codeword collection and export:\n" \
  "  outC=[file]        Export found codewords of weight up to min_w + dW (min_w:\n" \
  "                     minimum weight found) to file in .nz list format.  In\n" \
  "                     method=2/3, once dmin = dmax = d, CC rounds w = d..d+dW\n" \
  "                     enumerate all such irreducible codewords (see irred; unless\n" \
  "                     the timeout is hit; for a DEM, up to the choice of equal\n" \
  "                     columns: see merge).\n" \
  "  finC=[file]        Import initial candidate codewords from file in .nz format.\n" \
  "                     A reducible codeword is replaced by an irreducible part with\n" \
  "                     a nonzero syndrome of L (see irred).\n" \
  "  maxC=[int]         Maximum number of codewords to collect (default: 0).\n" \
  "                     0 means collect all valid codewords found up to wmax or stop\n" \
  "                     flag. If maxC > 0, halts when maxC codewords are collected.\n" \
  "  dW=[int]           Extra weight window above minimum distance (default: 0).\n" \
  "                     Collects codewords with weight up to min_w + dW.\n" \
  "  irred=[0|1]        Expert option: 1 (default) records only irreducible\n" \
  "                     codewords: a codeword which is the union of a codeword and a\n" \
  "                     disjoint nonzero vector with zero syndrome (e.g., a trivial\n" \
  "                     weight-3 vector of a circuit DEM, from weight d+3 on) is not\n" \
  "                     recorded.  Minimum-weight codewords are always irreducible.\n" \
  "                     Tested for the CC codewords (the rank of their columns of\n" \
  "                     H) and for the RW pairs (the sum of two disjoint vectors\n" \
  "                     with zero syndrome is skipped); a finC codeword is replaced\n" \
  "                     by an irreducible part with a nonzero syndrome of L.\n" \
  "                     0: all codewords are recorded, as before version 0.14.8.\n" \
  "                     (dist_m4ri_old: finC only.)\n\n" \
  "Execution, multithreading, and timing:\n" \
  "  threads=[int]      Maximum number of worker threads to use (default: 0 = all\n" \
  "                     logical CPUs available, i.e., the CPU affinity mask, e.g.,\n" \
  "                     with taskset). Subject to throttling for small codes or\n" \
  "                     large matrices unless nothrottle=1 is set.\n" \
  "  timeout=[sec]      Execution timeout in seconds (default: 60.0, 0 = infinite).\n" \
  "                     In method=3, it guides the CC vs RW thread balance, and a\n" \
  "                     CC round predicted not to finish in time is not started\n" \
  "                     while RW runs (in method=2, and with steps=0, it is started\n" \
  "                     anyway: it may still find a codeword of weight w=dmin).\n" \
  "  nothrottle=[0|1]   Disable automatic thread throttling (default: 0).\n" \
  "                     Aliases: --no-throttle, -no-throttle, nothrottle.\n" \
  "                     By default, threads are throttled for very small codes or\n" \
  "                     large matrices to avoid cache and memory bus thrashing.\n" \
  "                     The limits for small codes (by the RW work) apply in methods\n" \
  "                     1 and 3; method=2 starts 4 threads, and the others once a CC\n" \
  "                     round takes longer (> 2 ms). The limits for large matrices\n" \
  "                     (above 15 MB of dense RW memory per thread: at most 32 RW\n" \
  "                     threads, within a quarter of the physical memory) and for\n" \
  "                     small step counts (ceil(steps/10)) apply to RW threads only;\n" \
  "                     in method=3, CC rounds can use all threads.\n" \
  "  chunk_size=[int]   Maximum RW batch per claim of a worker (default: 0 = adaptive,\n" \
  "                     25 to 500 steps by n and steps). Near the end of the steps,\n" \
  "                     batches shrink to ceil(remaining/(4*RW threads)) steps, so\n" \
  "                     that the RW threads finish together. Alias: batch=[int].\n" \
  "  ksub=[int]         EXPERIMENTAL, should not be used (a warning is printed):\n" \
  "                     subspace dimension sampled from ker(H) in RW (default: 0 =\n" \
  "                     original full-matrix RW). E.g. ksub=32 or 64 keeps working\n" \
  "                     matrices inside L1/L2 cache. Automatically falls back to\n" \
  "                     ksub=0 (with a warning) if m < nu = dim(ker(H)).  Each step\n" \
  "                     reduces the span of ksub distinct rows of a fixed basis of\n" \
  "                     ker(H), so that a few codewords are found much more often\n" \
  "                     than others, and min_hits may end RW early.\n" \
  "  kwin=[int]         Localized column permutation window size W, used in every RW\n" \
  "                     step (alias: win=[int]).  Default: 0 = automatic hybrid for\n" \
  "                     n>=500: every other step uses an index window (win_mode=1)\n" \
  "                     of win_hyb columns, the other steps uniform random\n" \
  "                     permutations (also for codewords which are not local).\n" \
  "  win_mode=[0|1]     Window construction mode when kwin > 0 (default: 0):\n" \
  "                     0: Tanner graph BFS neighbors around random seed column.\n" \
  "                     1: Contiguous index proximity window around seed column.\n" \
  "  win_hyb=[int]      Window size of the hybrid RW steps (kwin=0, n>=500; at most\n" \
  "                     3n/4).  Default: 0 = 256 columns for a DEM, 512 for\n" \
  "                     matrices.  The window columns come first in the elimination\n" \
  "                     order; in a DEM, they are errors close in time, which touch\n" \
  "                     few detectors: windows of 192-256 columns find circuit\n" \
  "                     codewords 2-10 times faster than 512, and windows of one\n" \
  "                     round of columns much slower.  With debug&1, the binary\n" \
  "                     prints the rounds of a DEM (estimated from the time\n" \
  "                     translation of its errors: D detectors and L columns per\n" \
  "                     round) and the codewords found by the uniform and by the\n" \
  "                     window steps.\n" \
  "  pairs=[int]        RW pairs of non-pivot columns (p=2, as in Lee-Brickell and\n" \
  "                     Stern): a RW step also finds the codewords with exactly two\n" \
  "                     non-pivot columns, checking only the pairs of columns which\n" \
  "                     agree on pairs=l random pivot rows (at most 64; 0 = off).\n" \
  "                     Default: -1 = auto, about one non-pivot column per bucket,\n" \
  "                     with fewer rows while a codeword of the weight of the upper\n" \
  "                     bound would be kept with probability below 0.8, but enough\n" \
  "                     rows to keep the pairs of a step within about the work of\n" \
  "                     its elimination (e.g., l >= 6 for large circuit DEMs while\n" \
  "                     the upper bound is far above the distance).  The pairs\n" \
  "                     find the minimum-weight codewords of circuit DEMs 7 to 25\n" \
  "                     times faster, and lighter codewords of codes with a large\n" \
  "                     distance.  Not with ksub.  With debug&1, the binary prints\n" \
  "                     the pairs checked per step and the codewords found by pairs\n" \
  "                     and by single columns.\n" \
  "  min_hits=[int]     QDistRnd-style RW stopping criterion (default: 5, 0 = off).\n" \
  "                     Stops RW when the average number of times <n> that RW has\n" \
  "                     found each codeword reaches min_hits, both for the codewords\n" \
  "                     of the minimum weight found and for a representative set of\n" \
  "                     cov_cws lowest-weight codewords (heavier codewords are kept\n" \
  "                     until enough lighter ones are found).  A single codeword\n" \
  "                     suffices.  A lighter codeword is missed with probability\n" \
  "                     about exp(-<n>), assuming that lighter codewords are found\n" \
  "                     at least as often.  With debug&1, a warning is printed at\n" \
  "                     the end if the hit counts of codewords of the same weight\n" \
  "                     are strongly non-uniform (then exp(-<n>) is optimistic),\n" \
  "                     together with an information-set estimate (uniform random\n" \
  "                     information sets) of the probability to miss a lighter\n" \
  "                     codeword and of the RW steps needed for 1%%.\n" \
  "  cov_cws=[int]      Size of the representative set of lowest-weight codewords\n" \
  "                     for min_hits (default: 100; 0 = only the minimum weight).\n" \
  "                     Without outC and maxC, at most cov_cws minimum-weight\n" \
  "                     codewords are tracked (all if cov_cws=0).  With fewer\n" \
  "                     minimum-weight codewords, a smaller cov_cws may stop RW\n" \
  "                     earlier.\n" \
  "  refresh=[int]      RW steps interval for adaptive ker(H) basis refresh via\n" \
  "                     low-weight codeword exchange and re-echelonization, only\n" \
  "                     with the experimental ksub > 0 (default: 5000 when\n" \
  "                     ksub > 0, 0 = off).\n" \
  "  seed=[int]         Random number generator seed (default: 0; seed<=0 is\n" \
  "                     replaced by time(NULL) + 10*pid - 1000*seed).\n\n" \
  "Debug output bitmap (debug=[int], default: 3):\n" \
  "  Diagnostic output goes to stderr (stdout only has the result line); lower\n" \
  "  bits are more informative.  The value is decimal or hexadecimal (0x...).\n" \
  "  Multiple debug arguments are OR-combined, where the first one replaces the\n" \
  "  default: e.g., debug=0 alone is silent, while debug=4 and debug=0 debug=4\n" \
  "  both give only the periodic status.\n" \
  "    0    : Silent: errors and warnings on the validity of the result (expert\n" \
  "           options, experimental ksub, invalid finC codewords) only, and the\n" \
  "           confinement profile with smax > 0\n" \
  "    1    : Summary (default): input matrices (and the rounds of a DEM, RW\n" \
  "           windows), warnings, stop reason, codeword and hit statistics,\n" \
  "           RW information sets, pairs, and step types, export\n" \
  "    2    : Progress (default): run plan (threads, steps, timeout, seed),\n" \
  "           finished CC rounds, new upper bounds, RW convergence\n" \
  "    4    : Periodic status (after 1, 2, 4, ... s, then every 60 s): bounds,\n" \
  "           RW steps and rate, <n>, progress of the current CC round\n" \
  "    8    : Thread allocation and timing: throttling, CC round planning and\n" \
  "           re-planning, waits for RW, predicted vs measured CC work, and the\n" \
  "           nodes of each CC round (calls of the recursion)\n" \
  "    16   : Code parameters rank(H), rank(L) or rank(G), and k (dense\n" \
  "           elimination at startup, slow for large matrices)\n" \
  "    32   : Support (1-based) of each new lightest codeword (RW, CC, finC)\n" \
  "    64   : Command-line arguments, rng seed, and files read\n" \
  "    128  : Matrices H, G, L (if n < 150, unless 2048 is set)\n" \
  "    256  : All codewords in hash with their hit counts (at the end)\n" \
  "    2048 : No size cutoff for matrices and codeword supports\n" \
  "    4096 : dist_m4ri_old only: CC recursion and confinement hash traces\n" \
  "  Bits 512, 1024, and 8192 to 32768 are reserved.  Bits 65536 and above are\n" \
  "  reserved for the Python wrapper dist_m4ri.py (ignored by the binary).\n" \
  "  E.g., debug=7 adds the periodic status to the default output, and debug=11\n" \
  "  the thread allocation and timing.\n\n" \
  "Output format (stdout):\n" \
  "  Standard output produces a single line with three space-separated integers:\n" \
  "    dmin dmax rw_steps\n" \
  "  where:\n" \
  "    dmin-1   : Maximum cluster weight analyzed by CC without finding any codeword.\n" \
  "               If CC finds an exact minimum-weight codeword, dmin = dmax = weight.\n" \
  "               With the expert options noscan, start, or cbeg/cend, see the\n" \
  "               warnings above for the validity of dmin.\n" \
  "    dmax     : Smallest weight of a codeword found (by RW or CC, or read with\n" \
  "               finC), or the supplied dmax if smaller (0 if none).\n" \
  "    rw_steps : Number of completed RW steps (0 if CC found exact or method=2).\n\n" \
  "Option groups (--help GROUP prints all entries of GROUP):\n" \
  "  windows            kwin, win_mode, win_hyb\n" \
  "  stopping, stop     steps, wmin, timeout, min_hits, cov_cws\n" \
  "  bounds             dmin, dmax, dstop, dexp, wmax\n" \
  "  dem                fdem, pmin, merge, drop, outD\n" \
  "  matrices, matrix   finH, finG, finL, fin, classical\n" \
  "  codewords          outC, finC, maxC, dW, irred\n" \
  "  threading          threads, nothrottle, chunk_size\n" \
  "  expert             noscan, start, start_scan, cbeg, cend\n" \
  "  experimental       ksub, refresh\n\n" \
  "Help options:\n" \
  "  -h, --help         Display help for commonly used parameters (fits 80 rows);\n" \
  "                     --help TOPIC is the same as --morehelp TOPIC.\n" \
  "  --morehelp         Display this full help message listing all parameters;\n" \
  "                     --morehelp TOPIC (also --morehelp=TOPIC) displays only the\n" \
  "                     entries of TOPIC, e.g., --morehelp maxC: one or more names\n" \
  "                     separated by commas (e.g., maxC,dW), case-insensitive,\n" \
  "                     without '-', '_', or '=VALUE', including aliases such as\n" \
  "                     dest or batch (and debug for the debug output bitmap); an\n" \
  "                     unknown name gives the closest names on stderr, and the\n" \
  "                     exit code 1.  A group name (see Option groups, also as\n" \
  "                     singular or plural, e.g., window) gives all entries of the\n" \
  "                     group, and a parameter of a group also the other entries\n" \
  "                     of its group (see also).\n" \
  "  --version          Display program version\n"

#define BRIEF_HELP \
  "try \"%s -h\" for help"

#endif /* UTIL_IO_H */
