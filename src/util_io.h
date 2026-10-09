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

//static const int max_row_wt=10; 

#define MAX_W 100 
struct CW_VEC_T;
typedef struct CW_VEC_T cw_vec_t;
typedef struct{
  int debug; /* debug information */ 
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
  int wmin; /** min distance below which we are not interested 
		if w <= wmin found in RW, terminate immediately 
		start clusters with `wmin` for `CC`
	     */
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
  cw_vec_t *codewords; /* hash of codewords: the collection window (weight up to min_w + dW with outC or maxC)
                          and the representative set for `min_hits`, see codeword_add_maybe() */
  long long int num_cws; /* number of codewords in hash */
  int min_w;             /* minimum weight of the codewords in hash (INT_MAX: none) */
  int cw_max_w;          /* maximum weight of the codewords in hash (0: none) */
  long long int *cw_cnt_w;  /* number of codewords of each weight in hash (size nvar+2, allocated on first use) */
  long long int *cw_hits_w; /* total hits of the codewords of each weight in hash (same size) */
  char *fdem;
  double pmin;
  char *finH;
  char *finG;
  char *finL;
  char *fin;
  csr_t *spaH;
  csr_t *spaG;
  csr_t *spaL;
  int threads; /* number of threads to use (0 for auto) */
  int dexp;    /* expected distance value (0 for auto/none) */
  double timeout; /* timeout in seconds (default 60.0, 0 for infinite) */
  int nothrottle; /* 1: disable automatic thread throttling */
  int chunk_size; /* RW chunk/batch size (0 for auto) */
  int ksub;       /* RW subspace dimension sampled from ker(H) (0 for full H) */
  int kwin;       /* RW localized window size W (0 for uniform permutation) */
  int win_mode;   /* RW window mode: 0 = Tanner BFS, 1 = index proximity */
  int min_hits;   /* RW stopping criterion: average hits <n> per lowest-weight cw (0 = off) */
  int cov_cws;    /* RW stopping criterion: size of the representative set of lowest-weight cws (default 100) */
  int refresh;    /* RW steps interval for adaptive basis refresh (0 = off, auto 5000 if ksub>0) */
} params_t;

static inline int minint(const int a, const int b) { return (a < b) ? a : b; }
// #define MININT(a,b) do{ int t1=(a); int t2=(b); t1<t2? t1 :t2; } while(0)

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
 * @return The codeword hash (`p->codewords`).
 */
cw_vec_t * codeword_add_maybe(params_t * const p, const int arr[], int weight);

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
 * P1(w) = w C(n-w, k-1) / C(n, k).  Codewords of weight w > r + 1 are never found (P1 = 0).
 *
 * @param n Block length (number of columns of H).
 * @param rank Rank r of H.
 * @param w Codeword weight.
 * @return P1(w).
 */
double rw_infoset_find_prob(const int n, const int rank, const int w);

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
 */
void print_rw_infoset_estimate(FILE *stream, const int n, const int rank, const int w_lo, const int w_hi,
                               const long steps, const long steps_total, const double t_step);

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
 * @param stream Output stream (typically stderr).
 * @param p Pointer to the params_t structure.
 */
void print_codeword_stats(FILE *stream, const params_t * const p);

#define DIST_M4RI_VERSION "0.10.2"

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
  "  method, finH, finG, finL, fin, fdem, pmin, classical, css,\n" \
  "  steps, wmin, wmax, dmin, dmax, dexp (dest), smax, start, cbeg,\n" \
  "  cend, noscan, threads, timeout, nothrottle, chunk_size (batch),\n" \
  "  ksub, kwin (win), win_mode, min_hits, cov_cws, refresh,\n" \
  "  finC, outC, maxC, dW, seed, debug\n\n" \
  "Help options:\n" \
  "  -h, --help    : display help for commonly used parameters (fits 80 rows)\n" \
  "  --morehelp    : display full help for all available parameters\n"

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
  "  classical=[0|1]    1: classical code (Hx only), 0: quantum CSS (auto-detected)\n\n" \
  "Distance bounds and guidance:\n" \
  "  dmin=[int]         Certified lower bound on distance (CC starts from dmin) (1)\n" \
  "  dmax=[int]         Known upper bound on distance (RW ignores cw wt >= dmax) (0)\n" \
  "  dexp=[int]         Expected distance for method=3 thread allocation (alias: dest) (0)\n\n" \
  "Search limits and stopping criteria:\n" \
  "  steps=[int]        Maximum RW decoding steps / information sets (100000)\n" \
  "  wmax=[int]         Maximum cluster weight to search in CC (0=until bound/timeout)\n" \
  "  wmin=[int]         Stop immediately if cw with weight <= wmin is found (1)\n" \
  "  min_hits=[int]     Stop RW when lowest-wt cws are hit min_hits times on average (5)\n" \
  "  timeout=[sec]      Execution timeout in seconds, 0 for infinite (60.0)\n\n" \
  "Multithreading and RW optimization:\n" \
  "  threads=[int]      Max worker threads to use (0: auto CPU count) (0)\n" \
  "  ksub=[int]         Subspace dimension sampled from ker(H) for RW (0: full H) (0)\n" \
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
  "  cbeg/cend=[int] (-1)   Expert: CC start column range, for split runs\n" \
  "  nothrottle=[0|1] (0)   Disable thread throttling (also --no-throttle)\n" \
  "  chunk_size=[int] (0)   RW batch chunk size (0: auto, alias: batch)\n" \
  "  win_mode=[0|1] (0)     Window mode: 0=Tanner BFS, 1=index proximity\n" \
  "  cov_cws=[int] (100)    Number of lowest-wt cws used for the min_hits statistic\n" \
  "  refresh=[int] (0)      RW steps between adaptive ker(H) basis refreshes (auto 5000 if ksub>0)\n" \
  "  seed=[int] (0)         RNG seed [0 for time(NULL)]\n" \
  "  debug=[int] (3)        Debug bitmask (0: silent, 1: general, 2: verbose, ...)\n\n" \
  "Help options:\n" \
  "  -h, --help         Display this help message (commonly used parameters)\n" \
  "  --morehelp         Display full help with all parameter descriptions\n" \
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
  "  fin=[str]          Base name for matrix input files (default: \"try\").\n" \
  "                     Automatically looks for \"${fin}X.mtx\" as finH and\n" \
  "                     \"${fin}Z.mtx\" as finG.\n" \
  "  fdem=[file]        Detector Error Model file (.dem) generated by Stim.\n" \
  "                     Automatically constructs H and L matrices. Cannot be\n" \
  "                     combined with finH, finG, finL, or fin.\n" \
  "  pmin=[float]       Minimum error probability threshold for DEM parsing (0.0).\n" \
  "                     Error mechanisms with probability < pmin are ignored.\n" \
  "                     Only valid when fdem is specified.\n" \
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
  "                     In RW (method 1/3), codewords of weight >= dmax are ignored\n" \
  "                     unless collecting codewords or needed for min_hits.\n" \
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
  "                     0 means continue until codeword found, bounds meet, or timeout.\n" \
  "  smax=[int]         Maximum syndrome weight for confinement profile (default: 0).\n" \
  "                     When smax > 0, tracks minimum syndrome weights for each\n" \
  "                     error weight. When smax=0, confinement is not computed.\n" \
  "  noscan=[0|1]       Expert option: 1: run CC only at weight w = wmax, skipping\n" \
  "                     the scan over weights w < wmax (default: 0). Only valid for\n" \
  "                     method=2. WARNING: lower weights are not checked, so a\n" \
  "                     codeword found only sets dmax, and dmin is not raised above\n" \
  "                     the supplied dmin. Certified only if dmin=wmax is given.\n" \
  "  start=[list]       Expert option: comma-separated list of CC start columns,\n" \
  "                     e.g., start=0,48,96 (default: -1 = all columns). Clusters\n" \
  "                     grown from a listed column are NOT limited to larger column\n" \
  "                     indices. Use, e.g., one column per block of a quasi-cyclic\n" \
  "                     code. WARNING: dmin is valid only if a code symmetry maps\n" \
  "                     every minimum-weight codeword to one containing a listed\n" \
  "                     column. Cannot be combined with cbeg/cend.\n" \
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
  "                     enumerate all such codewords (unless the timeout is hit).\n" \
  "  finC=[file]        Import initial candidate codewords from file in .nz format.\n" \
  "  maxC=[int]         Maximum number of codewords to collect (default: 0).\n" \
  "                     0 means collect all valid codewords found up to wmax or stop\n" \
  "                     flag. If maxC > 0, halts when maxC codewords are collected.\n" \
  "  dW=[int]           Extra weight window above minimum distance (default: 0).\n" \
  "                     Collects codewords with weight up to min_w + dW.\n\n" \
  "Execution, multithreading, and timing:\n" \
  "  threads=[int]      Maximum number of worker threads to use (default: 0 =\n" \
  "                     hardware concurrency, capped at 64). Subject to throttling\n" \
  "                     for small codes or large matrices unless nothrottle=1 is set.\n" \
  "  timeout=[sec]      Execution timeout in seconds (default: 60.0, 0 = infinite).\n" \
  "                     In method=3, it guides the CC vs RW thread balance, and a\n" \
  "                     CC round predicted not to finish in time is not started\n" \
  "                     while RW runs (in method=2, and with steps=0, it is started\n" \
  "                     anyway: it may still find a codeword of weight w=dmin).\n" \
  "  nothrottle=[0|1]   Disable automatic thread throttling (default: 0).\n" \
  "                     Aliases: --no-throttle, -no-throttle, nothrottle.\n" \
  "                     By default, threads are throttled for very small codes or\n" \
  "                     large matrices to avoid cache and memory bus thrashing.\n" \
  "                     The limits for large matrices (dense RW memory) and for\n" \
  "                     small step counts (ceil(steps/10)) apply to RW threads only;\n" \
  "                     in method=3, CC rounds can use all threads.\n" \
  "  chunk_size=[int]   RW batch chunk size per worker (default: 0 = adaptive).\n" \
  "                     Alias: batch=[int].\n" \
  "  ksub=[int]         Subspace dimension sampled from ker(H) in RW (default: 0 =\n" \
  "                     original full-matrix RW). E.g. ksub=32 or 64 keeps working\n" \
  "                     matrices inside L1/L2 cache. Automatically falls back to\n" \
  "                     ksub=0 (with a warning) if m < nu = dim(ker(H)).\n" \
  "  kwin=[int]         Localized column permutation window size W (default: 0 =\n" \
  "                     automatic hybrid window/uniform for n>=500; alias: win=[int]).\n" \
  "  win_mode=[0|1]     Window construction mode when kwin > 0 (default: 0):\n" \
  "                     0: Tanner graph BFS neighbors around random seed column.\n" \
  "                     1: Contiguous index proximity window around seed column.\n" \
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
  "                     low-weight codeword exchange and re-echelonization\n" \
  "                     (default: 5000 when ksub > 0, 0 = off).\n" \
  "  seed=[int]         Random number generator seed (default: 0 = initialize from\n" \
  "                     current time).\n\n" \
  "Debug output bitmap (debug=[int], default: 3):\n" \
  "  The debug parameter accepts a bitmask controlling diagnostic output to stderr:\n" \
  "    0    : Clear entire debug bitmap (completely silent execution)\n" \
  "    1    : General progress and round summary information (on by default)\n" \
  "    2    : Detailed thread allocation and timing information (on by default)\n" \
  "    4    : Command-line argument parsing diagnostics\n" \
  "    8    : Progress reports every 1000 RW steps\n" \
  "    16   : Output new minimum-weight codewords as they are found\n" \
  "    32   : Dump matrices and full codeword lists\n" \
  "    64   : Debug confinement hash table updates (swei changes)\n" \
  "    128  : Debug duplicate syndromes in confinement hash (debug build)\n" \
  "    2048 : Allow large matrix and vector output (bypasses size cutoff)\n" \
  "  Multiple debug arguments are XOR-combined (except debug=0 which clears all).\n" \
  "  Tip: Place debug=0 as the first argument to silence all diagnostic output.\n\n" \
  "Output format (stdout):\n" \
  "  Standard output produces a single line with three space-separated integers:\n" \
  "    dmin dmax rw_steps\n" \
  "  where:\n" \
  "    dmin-1   : Maximum cluster weight analyzed by CC without finding any codeword.\n" \
  "               If CC finds an exact minimum-weight codeword, dmin = dmax = weight.\n" \
  "               With the expert options noscan, start, or cbeg/cend, see the\n" \
  "               warnings above for the validity of dmin.\n" \
  "    dmax     : Smallest weight of any codeword found (0 if no codeword found).\n" \
  "    rw_steps : Number of completed RW steps (0 if CC found exact or method=2).\n\n" \
  "Help options:\n" \
  "  -h, --help         Display help for commonly used parameters (fits 80 rows)\n" \
  "  --morehelp         Display this full help message listing all parameters\n" \
  "  --version          Display program version\n"

#define BRIEF_HELP \
  "try \"%s -h\" for help"

#endif /* UTIL_IO_H */
