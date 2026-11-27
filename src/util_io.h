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
#include <stdatomic.h>
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
#define DBG_LARGE     2048 /**< no size cutoff for the matrices, codeword supports, and per-weight cw statistics */
#define DBG_LEGACY    4096 /**< dist_m4ri_old only: CC recursion and confinement hash traces */
#define DBG_DEFAULT   (DBG_SUMMARY | DBG_PROGRESS)
#define DBG_CW_MAX    25   /**< without DBG_LARGE, at most this many positions of a codeword support are printed */

//static const int max_row_wt=10; 

#define MAX_W 100 
struct CW_VEC_T;
typedef struct CW_VEC_T cw_vec_t;

/* The stages of the RW search selected by `rw_opt=[int]` (a bitmask; -1: automatic, 0: none, as in version 0.16.1),
 * see rw_triplet.h and doc/rw_triplet_degeneracies.tex (Sec. 9, Stages 1-5 and the dispatch rules D1-D5; the automatic
 * mask follows D1-D5, except that Stage 2b is explicit only and the row slabs of D4 are automatic only for a DEM with a
 * single round or without round structure, see rw_slab_auto() and rw_stages_plan()) */
#define RW_OPT_T3     0x01 /**< Stage 1 (T_3 enumeration) and 2a (contraction, hoisted L syndromes, Stern filter) */
#define RW_OPT_GEN    0x02 /**< Stage 2b: shortening by the rows of G (finG) and the realized-generator filter */
#define RW_OPT_SLAB   0x04 /**< Stage 3: row-slab windows of a DEM with MIS push-back order */
#define RW_OPT_DUAL   0x08 /**< Stage 4: sparse-dual m x m inverse elimination */
#define RW_OPT_CANON  0x10 /**< Stage 5: canonical classes, coalesced hits, Good-Turing / Chao2 stopping (the Harris-FKG
                            *   score is printed with any stage, on the dense full-matrix steps with pairs) */
#define RW_OPT_EXPAND 0x20 /**< the expansion of the collected codewords by T_3 before the export (outC) */
#define RW_OPT_ALL    0x3f /**< all stages */

/** @brief Stage 5 (RW_OPT_CANON): the statistics of the canonical classes of the minimum weight w in the uniform RW
 *         steps, see rw_canon_tab_t in rw_triplet.h (each class counted at most once per step). */
typedef struct {
  int active;           /**< 1: the canonical stopping rule replaces the `min_hits` rule of codeword_hit_stats() */
  int full;             /**< 1: too many classes for the table (the rule of codeword_hit_stats() is used) */
  int w;                /**< the weight of the classes (the current min_w; INT_MAX: none yet) */
  long long classes;    /**< C_obs: the canonical classes of weight w seen in the uniform steps */
  long long raw;        /**< the distinct raw codewords of weight w seen in the uniform steps */
  long long n_unif;     /**< N_unif: the sum of the uniform-step incidences of the classes */
  long long f1, f2;     /**< Q1, Q2: the classes with exactly one and two uniform-step incidences */
} rw_canon_stats_t;

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
  int ortho; /* 1 (default): with finG, check that the rows of H and G are orthogonal (H G^T = 0); 0: expert, not
              * checked (the caller guarantees it) */
  cw_vec_t *codewords; /* hash of codewords: the collection window (weight up to min_w + dW with outC or maxC)
                          and the representative set for `min_hits`, see codeword_add_maybe() */
  long long int num_cws; /* number of codewords in hash */
  int min_w;             /* minimum weight of the codewords in hash (INT_MAX: none) */
  int cw_max_w;          /* maximum weight of the codewords in hash (0: none) */
  long long int *cw_cnt_w;  /* number of codewords of each weight in hash (size nvar+2, allocated on first use) */
  long long int *cw_hits_w; /* total hits of the codewords of each weight in hash (same size) */
  long long int *cw_nhits_w; /* their hits counted for the `min_hits` statistic <n> (the sources in nhits_src; same
                                size, see codeword_hit_stats()) */
  long long int *cw_ncnt_w; /* the codewords of each weight in hash with such a hit (same size) */
  int nhits_src;            /* the sources of the hits counted for <n>: a bitmask of 1 << CW_SRC_* (default: all; set by
                               dist_m4ri once the RW stages are resolved, see rw_nhits_sources() in dist_m4ri.c) */
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
  int drop_none;   /* 1: drop=auto found no detector which can be dropped exactly; 2: found none, but its search was
                      not exhaustive (a group not verified, or the work limit reached) */
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
  int k_zero;       /* k = 0 (no codeword at all): 1 yes (no search, the result is 0 0 0), 0 no (or not checked), -1 not
                       known yet (the rank check of k_zero_rank() follows), see var_init() */
  char k0_reason[64]; /* why k = 0, for the note of k_zero_note() */
  int threads; /* number of threads to use (0 for auto) */
  int dexp;    /* expected distance value (0 for auto/none) */
  int dstop;   /* stop target for the lower bound: the run ends once CC certifies dmin >= dstop (0: off); unlike dmax,
                  not an upper bound, and not reported as dmax */
  double timeout; /* timeout in seconds (default 60.0, 0 for infinite) */
  int nothrottle; /* 1: disable automatic thread throttling */
  int chunk_size; /* RW chunk/batch size (0 for auto) */
  int ksub;       /* RW subspace dimension sampled from ker(H) (0 for full H); experimental, should not be used */
  int kwin;       /* RW localized window size W, used in every step (0: uniform random permutations, and for n >= 500
                     an index window of win_hyb columns on every other step, see rw_hybrid_window(); < 0: no
                     windows) */
  int win_mode;   /* RW window mode: 0 = Tanner BFS, 1 = index proximity */
  int win_hyb;    /* window size of the hybrid RW steps (kwin=0, n >= 500: an index window on every other step); 0:
                     256 columns for a DEM, 512 for matrices (see rw_hybrid_window()) */
  int pairs;      /* RW pairs of non-pivot columns (p=2) agreeing on `pairs` random pivot rows (0: off; -1: auto, see
                     rw_pairs_bits()) */
  int rw_opt;     /* the RW stages (bitmask of RW_OPT_T3 etc.): -1 (default) automatic, 0 none (as in version 0.16.1);
                     ignored by dist_m4ri_old */
  int rw_active;  /* the RW stages in use, resolved from rw_opt by dist_m4ri at the start of the run (0 before) */
  int rw_planned; /* 1 once rw_active is resolved (dist_m4ri, see rw_hybrid_on()) */
  int slab_rows;  /* study option, the row slabs of rw_opt 4: detectors per slab (0, the default: m/2) */
  int slab_bnd;   /* study option: the probability in percent of each of the two boundary slabs (default 50: only the
                     boundary slabs, the two detector layers of a single-round DEM) */
  int slab_var;   /* study option: the random height reduction of a slab in percent (0, the default: none) */
  int slab_every; /* study option: a slab on every slab_every-th RW step (default 2: every other step; at least 2, never
                     slabs only) */
  char *steplog;  /* study option steplog=FILE: the log of the RW steps and of their codewords (NULL: none), see
                     steplog_t in dist_m4ri.c */
  int dem_rounds; /* the rounds T of the DEM if known (fdem; 0, the default: unknown): 1 adds the row slabs of rw_opt 4
                     to the automatic rw_opt, 2 or more excludes them, see rw_slab_auto() */
  rw_canon_stats_t canon; /* Stage 5 (RW_OPT_CANON): the statistics of the canonical classes (rw_canon_converged()) */
  int round_D;    /* the rounds of a DEM (see dem_round_estimate(); with debug&1, or for rw_slab_auto()): detectors per
                     round (0: none found) */
  int round_L;    /* columns per round, n D / m (0: none found) */
  double round_score; /* the score of the estimate (accepted if >= 0.5; -1: not estimated) */
  int min_hits;   /* RW stopping criterion: average hits <n> per lowest-weight cw (0 = off) */
  int cov_cws;    /* RW stopping criterion: size of the representative set of lowest-weight cws (default 100) */
  int refresh;    /* RW steps interval for adaptive basis refresh (0 = off, auto 5000 if ksub>0) */
  int stamp;      /* opt-in time stamps on progress and stop lines: 1 appends ' [t=S.SSSs]' (0: off, default) */
  int json;       /* opt-in JSON summary on stderr: 1 writes '# RESULT {json}' as its last line (0: off, default) */
  char *json_file;/* opt-in JSON summary file (json=FILE; NULL: none); may be combined with json=1 */
  int m_full;     /* rows of H as read (after pmin, before DEM merge and drop), for the JSON summary */
  int n_full;     /* columns of H as read (after pmin, before DEM merge and drop), for the JSON summary */
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
 * constructs logical matrices if needed, and performs consistency checks.  If k = 0, i.e., there is no codeword
 * (rank H = n for a classical code; for a quantum code or a DEM, every row of L is in the row space of H, also if L
 * has no rows, e.g., a DEM with no L left after the pmin cut, merge, and drop), and this is known at once (L without
 * rows, also L computed from H and G, or L = 0), it sets p->k_zero = 1 and returns before the warnings for the search
 * options: the binaries then print the note of k_zero_note() and the result 0 0 0.  If k is not known at once, it sets
 * p->k_zero = -1: the binaries then run k_zero_rank() (dist_m4ri in its own thread during the search).
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
 * @brief The note on stderr (at any debug level) for k = 0 (p->k_zero, see var_init()): "# Note: k=0: no logical
 * codeword (<reason>); the result is <result>".
 *
 * @param p Pointer to the params_t structure (p->k0_reason).
 * @param result The result line of the binary ("0 0 0" for dist_m4ri, "0" for dist_m4ri_old).
 */
void k_zero_note(const params_t * const p, const char * const result);

/**
 * @brief DBG_SUMMARY: the line of the hybrid RW windows of the input summary (see rw_hybrid_on()); dist_m4ri prints it
 * after resolving an automatic rw_opt to 0.
 */
void print_window_summary(FILE *stream, const params_t * const p);

/**
 * @brief The rank check of k = 0 (when not known at once, p->k_zero = -1 after var_init()): the row echelon form of a
 * dense copy of H (a classical code: rank H = n), then the rows of L reduced by it (every row of L in the row space of
 * H).  Thread-safe: it reads only H and L, writes only `reason`, and allocates with malloc() only (M4RI's memory cache
 * is not thread-safe); the dense copy has m*n bits (see the help of the output format).
 *
 * @param H The matrix H (at least one row and one column).
 * @param L The matrix L, NULL for a classical code.
 * @param reason Why k = 0 (written only if k = 0), for p->k0_reason.
 * @param len The size of `reason`.
 * @param stop If not NULL: the check ends early (returns -1) once *stop is nonzero.
 * @return 1 if k = 0, 0 if k > 0, -1 if stopped.
 */
int k_zero_rank(const csr_t * const H, const csr_t * const L, char * const reason, const size_t len,
                atomic_int * const stop);

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
 *            CW_SRC_RW_WINDOW, CW_SRC_RW_SLAB; CW_SRC_OTHER for CC and finC), see print_rw_step_stats(); the hits of
 *            the sources in p->nhits_src are counted for the `min_hits` statistic <n> (see codeword_hit_stats()).
 * @return The codeword hash (`p->codewords`).
 */
cw_vec_t * codeword_add_maybe(params_t * const p, const int arr[], int weight, int src);

/**
 * @brief The hits of a codeword in the hash (see codeword_add_maybe()), 0 if it is not there (for steplog=FILE).
 *
 * @param p Pointer to the params_t structure.
 * @param arr The support of the codeword (sorted).
 * @param weight Its weight (length of arr).
 * @return The hit count `cnt` of the codeword in hash, or 0.
 */
int codeword_count(const params_t * const p, const int arr[], const int weight);

/** @brief The source of a codeword added to the hash, see codeword_add_maybe(): CC, finC, and the other sources; the RW
 *         steps with a uniform random permutation, with a window (kwin, the hybrid index windows), and the row slabs of
 *         rw_opt 4. */
enum { CW_SRC_OTHER = 0, CW_SRC_RW_UNIFORM = 1, CW_SRC_RW_WINDOW = 2, CW_SRC_RW_SLAB = 3 };

/** @brief All sources (p->nhits_src: all hits count for <n>, as in version 0.16.1). */
#define CW_SRC_ALL ((1 << CW_SRC_OTHER) | (1 << CW_SRC_RW_UNIFORM) | (1 << CW_SRC_RW_WINDOW) | (1 << CW_SRC_RW_SLAB))

/**
 * @brief DBG_SUMMARY: the codewords of the lowest three weights in hash found by the RW steps with a uniform random
 *        permutation, by those with a window (index windows, kwin), and by the row slabs of rw_opt 4: the distinct
 *        codewords and their hits, e.g., to see whether some codewords are found only by the uniform steps (not
 *        local in the column order).  Only the codewords in hash are counted: without outC and maxC, at most cov_cws
 *        codewords of the minimum weight (the first ones found) and heavier ones only for the representative set of
 *        min_hits (see codeword_add_maybe()); "both" counts the codewords found by the uniform steps and by the window
 *        or slab steps, and "of N cws" all codewords of the weight in hash (also those found only by CC or read with
 *        finC).
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

/**
 * @brief Whether the hybrid RW steps (kwin=0, n >= 500) use index windows on every other step: always when the RW
 *        stages in use are none (rw_opt=0, as in version 0.16.1, or an automatic rw_opt which resolves to 0, e.g., for
 *        a classical code, steps <= 2, or ksub > 0); otherwise only for a DEM (columns in circuit order), or with
 *        win_hyb > 0 (Rule D3 of doc/rw_triplet_degeneracies.tex: windows do not help the matrices of codes and
 *        extended matrices, whose column order is not local).  Keyed on the resolved mask p->rw_active once
 *        p->rw_planned is set by dist_m4ri (before: on the argument rw_opt, an automatic -1 not yet known).
 */
static inline int rw_hybrid_on(const params_t * const p, const int nvar) {
  const int mask = p->rw_planned ? p->rw_active : p->rw_opt;
  return (p->kwin == 0) && (nvar >= 500) && ((mask == 0) || (p->fdem != NULL) || (p->win_hyb > 0));
}

/** @brief The most detectors of a DEM without round structure with the row slabs in the automatic rw_opt, see
 *         rw_slab_auto() (the default slab height m/2 is then at most 320, as in the row-slab study). */
#define RW_SLAB_AUTO_MMAX 640

/**
 * @brief Whether the automatic rw_opt (-1) has the row slabs of rw_opt 4 (Stage 3, Rule D4) for a DEM (fdem; with
 *        kwin=0 and the rows not sorted, see rw_stages_plan() in dist_m4ri.c): if the DEM is known to have a single
 *        round (dem_rounds=1, e.g., passed by dist_m4ri.py for a Stim circuit with one round), or if its rounds are not
 *        known (dem_rounds=0) and it has no round structure and at most RW_SLAB_AUTO_MMAX detectors.  A DEM has round
 *        structure if the estimate of dem_round_estimate() finds a time translation of its errors by D detectors, 2 <=
 *        D <= m/3 (at least three detector layers; its score at least 1/2): a single-round DEM, whose detectors form
 *        two layers (the round and the final readout), never has one, nor do the folded schedules whose rounds
 *        alternate between two forms, or a DEM whose detectors are not in time order.  On such a DEM, the boundary
 *        slabs (the first and the last m/2 detectors, see slab_bnd) are the two detector layers of a single round:
 *        they reach the codewords confined to one layer, which the uniform steps and the index windows hardly find
 *        (the row-slab study: the weight-16 codewords of the double-gross DEM with one round).  With dem_rounds >= 2,
 *        no automatic slabs (on deep DEMs, the slabs found fewer codewords than the index windows).
 */
static inline int rw_slab_auto(const params_t * const p) {
  if (p->dem_rounds == 1) return 1;
  return (p->dem_rounds == 0) && (p->spaH != NULL) && (p->spaH->rows <= RW_SLAB_AUTO_MMAX) &&
         (p->round_score >= 0.0) && (p->round_D == 0);
}

/**
 * @brief Whether the row slabs of rw_opt 4 apply to the input (Stage 3, Rule D4, see rw_stages_plan() in dist_m4ri.c):
 *        a DEM (fdem) with its detectors in time order (rows not sorted, sort_rows=0), kwin=0 (the slabs replace the
 *        index windows of the hybrid steps), at least 4 detectors, and the stage explicit (`expl`: rw_opt with bit 4)
 *        or in the automatic mask (rw_slab_auto()).
 */
static inline int rw_slab_input(const params_t * const p, const int expl) {
  return (p->fdem != NULL) && (p->kwin == 0) && !p->rows_sorted && (p->spaH != NULL) && (p->spaH->rows >= 4) &&
         (expl || rw_slab_auto(p));
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
  long long int min_cws;  /**< number of codewords of the minimum weight min_w with a counted hit */
  long long int min_hits; /**< their total number of counted hits */
  long long int set_cws;  /**< number of codewords in the representative set (weights w_lo..w_hi) with a counted hit */
  long long int set_hits; /**< their total number of counted hits */
  long long int all_cws;  /**< number of all codewords in the representative set (also those without a counted hit) */
  long long int excl_hits; /**< the hits of the representative set not counted (sources not in p->nhits_src) */
  int w_lo;               /**< smallest weight in the representative set (min_w) */
  int w_hi;               /**< largest weight in the representative set */
  double avg;             /**< <n>: the smaller of the average counted hits per codeword (with a counted hit) of weight
                               min_w and of the representative set (0 if no codeword of weight min_w has one) */
} cw_hit_stats_t;

/**
 * @brief Compute the hit statistics of the `min_hits` stopping criterion.
 *
 * The representative set consists of the lowest weight classes (starting with min_w) in hash
 * with at least cov_cws codewords in total (all classes if fewer; only min_w if cov_cws <= 0 or
 * min_hits <= 0).  As in QDistRnd, <n> estimates the average number of times a codeword is found;
 * a lighter codeword is missed with probability about exp(-<n>), assuming lighter codewords are
 * found at least as often.  Heavier codewords make the estimate more conservative.  Only the hits of the sources in
 * p->nhits_src count, over the codewords with such a hit: with a nonzero resolved rw_opt, not the hits of the row
 * slabs and of the window steps, which find the same codewords over and over (unless every RW step has a window,
 * kwin > 0), see rw_nhits_sources() in dist_m4ri.c.
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
 * with k=1) suffices.  Otherwise returns 0.  With the canonical classes of rw_opt (RW_OPT_CANON, p->canon.active),
 * the rule of rw_canon_converged() is used instead (see min_hits_rule()).
 *
 * @param p Pointer to the params_t structure.
 * @return 1 if convergence criterion is met, 0 otherwise.
 */
int check_min_hits_convergence(const params_t * const p);

/** @brief The rule of the `min_hits` stop, see min_hits_rule(). */
enum { MIN_HITS_OFF = 0, MIN_HITS_VECTOR = 1, MIN_HITS_CANON = 2 };

/**
 * @brief The rule of the `min_hits` stop of check_min_hits_convergence(): MIN_HITS_OFF (min_hits <= 0),
 *        MIN_HITS_CANON (the canonical classes of rw_opt bit 16, rw_canon_converged(), unless their table is full), or
 *        MIN_HITS_VECTOR (the average hits <n> of codeword_hit_stats()).  With MIN_HITS_CANON, only the codewords of
 *        the minimum weight found in the uniform RW steps count: while there is none (e.g., when only the window steps
 *        find codewords of the minimum weight), the stop cannot fire (the status line and the summary say so).
 */
int min_hits_rule(const params_t * const p);

/**
 * @brief The canonical stopping rule of rw_opt (RW_OPT_CANON, Sec. 7.3 of doc/rw_triplet_degeneracies.tex): with the
 *        statistics p->canon of the canonical classes of the minimum weight in the uniform RW steps (each class
 *        counted at most once per step), returns 1 if N_unif >= max(15, 3 min_hits), N_unif / C_obs >= min_hits, the
 *        Good-Turing estimate of the unseen probability mass p_GT = Q1 / N_unif <= min(0.02, 3 exp(-min_hits)), and
 *        the unseen fraction C0 / (C_obs + C0) <= 0.05 of the bias-corrected Chao2 estimate C0 = Q1 (Q1 - 1) /
 *        (2 (Q2 + 1)) of the unseen classes.
 */
int rw_canon_converged(const params_t * const p);

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
 * the hash.  Only the hits counted for <n> (see codeword_hit_stats()) and the codewords with such a hit are tested.
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
 * Among the weights w_lo..w_hi (not excluded by CC), the weight w with the smallest probability P(w) = P1(w) + P2(w)
 * (see rw_infoset_find_prob(); P2 = 0 without pairs) is used: a single codeword of this weight is missed in `steps`
 * RW steps with uniform random permutations with probability (1 - P)^steps.  Also prints the total number of RW steps
 * for a 1% miss probability (with the fraction steps / steps_total of uniform steps, as in this run; the window
 * steps, i.e., the index windows of kwin or win_hyb and the row slabs of rw_opt 4, are not counted) and, if `t_step`
 * > 0, the corresponding time.
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

#define DIST_M4RI_VERSION "0.17.5"

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
  "  method, finH, finG, finL, fin, fdem, pmin, merge, drop, outD, dem_rounds,\n" \
  "  classical, css, steps, wmin, wmax, dmin, dmax, dstop, dexp (dest), smax,\n" \
  "  start, start_scan, sort_rows, cc_prune, cbeg, cend, noscan, threads, timeout,\n" \
  "  nothrottle, chunk_size (batch), ksub, kwin (win), win_mode, win_hyb, pairs,\n" \
  "  rw_opt, min_hits (max_hits), cov_cws, refresh, finC, outC, maxC, dW, irred,\n" \
  "  ortho, seed, debug, stamp, json, steplog, slab_rows, slab_bnd, slab_var,\n" \
  "  slab_every\n\n" \
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
  "  ortho=[0|1] (1)        Expert: 0 skips the check H G^T = 0 of finG (valid only if it holds)\n" \
  "  sort_rows=[0|1] (1)    Sort the rows of H by weight for CC (default 0 with fdem)\n" \
  "  cc_prune=[0|1|2] (1)   CC search: 0 first-row walk, 1 branch and bound, 2 also cover bounds\n" \
  "  nothrottle=[0|1] (0)   Disable thread throttling (also --no-throttle)\n" \
  "  chunk_size=[int] (0)   RW batch chunk size (0: auto, alias: batch)\n" \
  "  win_mode=[0|1] (0)     Window mode: 0=Tanner BFS, 1=index proximity\n" \
  "  win_hyb=[int] (0)      Hybrid RW window size (0: 256 for a DEM, 512 otherwise)\n" \
  "  pairs=[int] (-1)       RW pairs of non-pivot columns (p=2; -1: auto, 0: off)\n" \
  "  rw_opt=[int] (-1)      RW triplet stages: -1 auto, 0 off (as 0.16.1), or a bitmask 1..63\n" \
  "  slab_rows=[int] (0)    Study, row slabs of rw_opt 4: detectors per slab (0: m/2); also\n" \
  "                         slab_bnd, slab_var, slab_every (see --morehelp study)\n" \
  "  steplog=[file]         Study: log the RW steps and the codewords found in them (see --morehelp study)\n" \
  "  cov_cws=[int] (100)    Number of lowest-wt cws used for the min_hits statistic\n" \
  "  refresh=[int] (0)      RW steps between ker(H) basis refreshes (experimental ksub only; 5000 if ksub>0)\n" \
  "  seed=[int] (0)         RNG seed [0 for time(NULL)]\n" \
  "  stamp=[0|1] (0)        Append elapsed search time ' [t=S.SSSs]' to progress and stop lines\n" \
  "  json=[0|1|file] (0)    Emit JSON summary (dm4ri-result/1): 1 '# RESULT' on stderr, or file\n" \
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
  "                     format. Used to verify orthogonality (see ortho) and\n" \
  "                     construct logicals.\n" \
  "  ortho=[0|1]        Expert option: 1 (default) checks that the rows of H and G\n" \
  "                     (finG) are orthogonal, H G^T = 0, and stops with an error\n" \
  "                     otherwise; 0 skips the check (a WARNING is printed), e.g.,\n" \
  "                     to screen many codes which are orthogonal by construction\n" \
  "                     (the check takes O(nnz(G) times the column weight of H)\n" \
  "                     operations).  The result is valid only if H G^T = 0.\n" \
  "                     Ignored (with a WARNING) without finG.\n" \
  "  finL=[file]        Logical operator matrix Lx for quantum CSS codes in Matrix\n" \
  "                     Market format.\n" \
  "                     Note: For a quantum CSS code, either finL (Lx) or finG (Hz)\n" \
  "                     must be provided.\n" \
  "  fin=[str]          Base name for the matrix files (default: ../examples/try).\n" \
  "                     Automatically looks for \"${fin}X.mtx\" as finH and\n" \
  "                     \"${fin}Z.mtx\" as finG.  A later finH, finG, or finL\n" \
  "                     cancels fin (with a warning; then no file is taken from\n" \
  "                     fin, so that finH must be given); fin after them is an\n" \
  "                     error.\n" \
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
  "  dem_rounds=[int]   The rounds T of the DEM, if known (default: 0 = unknown;\n" \
  "                     only valid with fdem): 1, a single round of syndrome\n" \
  "                     measurements and the final readout, adds the row slabs of\n" \
  "                     rw_opt 4 to the automatic rw_opt; 2 or more: no automatic\n" \
  "                     row slabs; 0: the automatic rw_opt has them for a DEM\n" \
  "                     without round structure (see win_hyb; never for more than\n" \
  "                     640 detectors).  dist_m4ri.py passes the rounds of a Stim\n" \
  "                     circuit (the largest number of measurements of a qubit).\n" \
  "  classical=[0|1]    Code type override:\n" \
  "                     1: Treat as a classical linear code (H only: finG,\n" \
  "                        finL, or fin is an error, and the observables L of a\n" \
  "                        DEM are discarded).\n" \
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
  "                     meet, dmin >= dstop, or timeout (method=2 then needs\n" \
  "                     timeout>0, dmax>0, or dstop>0), up to the weight 99 at\n" \
  "                     most (MAX_W - 1, see util_io.h; also the largest wmax).\n" \
  "  smax=[int]         Maximum syndrome weight for confinement profile (default: 0).\n" \
  "                     When smax > 0, tracks minimum syndrome weights for each\n" \
  "                     error weight. When smax=0, confinement is not computed.\n" \
  "                     Disabled (with a WARNING) with noscan=1 or dmin > 1, as the\n" \
  "                     small cluster weights are then skipped.\n" \
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
  "                     cbeg/cend.  With drop (fdem), a listed column touching a\n" \
  "                     dropped detector is replaced by its twin (or removed if it\n" \
  "                     has none, see drop).\n" \
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
  "                     a nonzero syndrome of L (see irred).  If finC is the same\n" \
  "                     file name as outC, a missing file is not an error.  A\n" \
  "                     column index below 1 is an error; codewords with a column\n" \
  "                     beyond n, not orthogonal to H, or orthogonal to L are\n" \
  "                     skipped (with a warning), and so are those heavier than\n" \
  "                     wmax > 0 (with method=2 and wmax=0, heavier than 99).\n" \
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
  "                     In method=3, adaptive batches take at most about 0.5 s of\n" \
  "                     RW steps (by the measured step time), and while RW steps\n" \
  "                     remain to be claimed, an RW worker leaves its batch between\n" \
  "                     two steps for a CC round which needs more threads (the rest\n" \
  "                     of the batch runs after the round).\n" \
  "  ksub=[int]         EXPERIMENTAL, should not be used (a warning is printed):\n" \
  "                     subspace dimension sampled from ker(H) in RW (default: 0 =\n" \
  "                     original full-matrix RW). E.g. ksub=32 or 64 keeps working\n" \
  "                     matrices inside L1/L2 cache. Automatically falls back to\n" \
  "                     ksub=0 if m < nu = dim(ker(H)) (with a warning if debug&1).\n" \
  "                     Each step reduces the span of ksub distinct rows of a fixed\n" \
  "                     basis of ker(H), so that a few codewords are found much\n" \
  "                     more often than others, and min_hits may end RW early.\n" \
  "  kwin=[int]         Localized column permutation window size W, used in every RW\n" \
  "                     step (alias: win=[int]).  Default: 0 = automatic hybrid for\n" \
  "                     n>=500: every other step uses an index window (win_mode=1)\n" \
  "                     of win_hyb columns, the other steps uniform random\n" \
  "                     permutations (also for codewords which are not local);\n" \
  "                     with a nonzero resolved rw_opt (the run plan), only for a\n" \
  "                     DEM or with win_hyb > 0, and with the row slabs of rw_opt\n" \
  "                     4 (e.g., automatic for a single round), a slab instead.\n" \
  "                     kwin=-1 (any negative value): no windows, uniform random\n" \
  "                     permutations in every step.\n" \
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
  "                     translation of its errors by D detectors, 2 <= D <= m/3,\n" \
  "                     i.e., at least three detector layers: D detectors and L\n" \
  "                     columns per round; a single round, a folded schedule,\n" \
  "                     whose rounds alternate between two forms, or detectors\n" \
  "                     not in time order give no round structure, which the\n" \
  "                     automatic row slabs of rw_opt 4 use, see dem_rounds) and\n" \
  "                     the codewords found by the uniform steps, by the window\n" \
  "                     steps, and by the row slabs (line '# RW step types': the\n" \
  "                     codewords in hash, i.e., without outC and maxC at most\n" \
  "                     cov_cws of the minimum weight, the first found; 'both':\n" \
  "                     found by the uniform and by the other steps).\n" \
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
  "                     find the minimum-weight codewords of circuit DEMs 7 to 26\n" \
  "                     times faster, and lighter codewords of codes with a large\n" \
  "                     distance.  Not with ksub.  With debug&1, the binary prints\n" \
  "                     the pairs checked per step and the codewords found by pairs\n" \
  "                     and by single columns.\n" \
  "  rw_opt=[int]       RW stages with the weight-3 degeneracies T3 (triplets of\n" \
  "                     columns with zero syndromes of H and L) and the known\n" \
  "                     generators (doc/rw_triplet_degeneracies.tex).  Default: -1\n" \
  "                     = automatic (the stages below which apply to the input);\n" \
  "                     0: none (the RW search of version 0.16.1); or the sum of\n" \
  "                     the stages to use (1 to 63, also as 0x..):\n" \
  "                     1: enumerate T3 (after merge and drop) and contract the\n" \
  "                        candidates (2 -> 1, 3 -> 0 columns): the single-column\n" \
  "                        candidates up to the weight min(dmax+8, 64), the pairs\n" \
  "                        up to min(dmax+1, 64) (+dW when collecting); in the\n" \
  "                        pairs, skip the columns whose circuit (the column and\n" \
  "                        its pivot columns) has a zero L syndrome and at most 3\n" \
  "                        columns (the Stern filter: loss-free for the distance,\n" \
  "                        i.e., the minimum weight, with the contraction, hence\n" \
  "                        only while dmax+8 <= 64, +dW when collecting; when\n" \
  "                        collecting codewords with outC or maxC, only with 32:\n" \
  "                        32 restores the codewords the filter removes through\n" \
  "                        2 <-> 2 exchanges, 1 -> 2 splits, and 1 <-> 1 swaps of\n" \
  "                        identical columns, e.g., in finH/finL exported from a\n" \
  "                        DEM without merging, or with merge=0); a\n" \
  "                        zero column of H with a nonzero column of L gives d=1\n" \
  "                        at once (automatic: quantum codes, steps > 2);\n" \
  "                     2: with finG, also shorten the candidates by the rows of\n" \
  "                        G (steepest descent, before the contraction of 1; the\n" \
  "                        single-column candidates up to min(dmax+w_G-2, 64),\n" \
  "                        with 1 up to min(dmax+max(8, w_G-2), 64)) and skip the\n" \
  "                        columns whose circuits are rows of G in the pairs (only\n" \
  "                        while that weight is at most 64; not when collecting;\n" \
  "                        automatic: off);\n" \
  "                     4: row-slab windows of a DEM in time order (kwin=0, not\n" \
  "                        with sort_rows=1): every other step (the others use\n" \
  "                        uniform random permutations), the columns inside the\n" \
  "                        first or the last m/2 detectors, i.e., the two\n" \
  "                        detector layers of a single round (see the study\n" \
  "                        options slab_rows, slab_bnd, slab_var, slab_every), a\n" \
  "                        random independent set of the triplet graph first\n" \
  "                        (automatic: a DEM with a single round, dem_rounds=1, or\n" \
  "                        without round structure, see win_hyb, and with at most\n" \
  "                        640 detectors: the slabs find the codewords confined to\n" \
  "                        one layer, e.g., weight 16 of the double-gross DEM with\n" \
  "                        one round, which the other steps hardly find; not on\n" \
  "                        deep circuit DEMs, where the index windows with 1 found\n" \
  "                        more minimum-weight codewords per second; the slabs\n" \
  "                        find the same codewords over and over, and min_hits\n" \
  "                        does not count their hits);\n" \
  "                     8: sparse-dual elimination: an m x m inverse instead of\n" \
  "                        the m x n matrix (the same information sets, less\n" \
  "                        memory per RW thread; quantum codes with up to 64 rows\n" \
  "                        of L; automatic: m x n above 512 KiB, and n >= 20 m or\n" \
  "                        m >= 1400; no '# RW FKG score');\n" \
  "                     16: canonical classes: the codewords of the minimum weight\n" \
  "                        found in the uniform RW steps are counted by their\n" \
  "                        canonical form (2 <-> 2 exchanges across two triplets\n" \
  "                        with a common column, or half-overlap exchanges with\n" \
  "                        the even rows of G), each at most once per step, and\n" \
  "                        min_hits stops RW by the Good-Turing / Chao2 rule (see\n" \
  "                        min_hits; only with min_hits > 0, not with kwin > 0;\n" \
  "                        automatic: quantum codes); the tables hold up to 2^24\n" \
  "                        classes and codewords (up to about 640 MiB; then the\n" \
  "                        rule of <n> is used); the codewords of a RW step are\n" \
  "                        recorded at its end (lighter ones at once), each at\n" \
  "                        most once (also with an active 1 or 2);\n" \
  "                     32: with outC or maxC, add the codewords obtained from the\n" \
  "                        codewords found by T3 before the export: 2 <-> 2\n" \
  "                        exchanges and 1 -> 2 splits, up to the weight min_w+dW;\n" \
  "                        and by 1 <-> 1 swaps of the identical columns of H and\n" \
  "                        L, also without T3 (see 1, irred, and maxC; automatic:\n" \
  "                        quantum codes with outC or maxC).\n" \
  "                     When the resolved rw_opt is not 0 (see the run plan),\n" \
  "                     matrices (finH) use uniform random permutations in all\n" \
  "                     RW steps unless kwin > 0 or win_hyb > 0 (the hybrid\n" \
  "                     windows are for DEMs, see kwin); an automatic rw_opt which\n" \
  "                     resolves to 0 (rw_opt=0x0 in the run plan, e.g., for a\n" \
  "                     classical code) keeps the hybrid windows of rw_opt=0.\n" \
  "                     Only with full-matrix RW (ksub=0) and\n" \
  "                     steps > 0 (automatic: steps > 2).  With debug&1, the\n" \
  "                     binary prints the stages used and their statistics\n" \
  "                     (lines '# RW triplets', '# RW canonical', and '# RW FKG\n" \
  "                     score': the mean Harris-FKG bound -sum ln(1 + |S_a|) on\n" \
  "                     ln Pr(J) of the dense full-matrix steps with pairs, with\n" \
  "                     any stage but not with 8 or on the slab steps, and\n" \
  "                     -ln C(n, r)).  dist_m4ri_old ignores rw_opt.\n" \
  "  min_hits=[int]     QDistRnd-style RW stopping criterion (default: 5, 0 = off;\n" \
  "                     alias: max_hits).  Stops RW when the average number of hits\n" \
  "                     <n> per codeword (over the codewords hit) reaches min_hits,\n" \
  "                     both for the codewords of the minimum weight found and for\n" \
  "                     a representative set of cov_cws lowest-weight codewords\n" \
  "                     (heavier codewords are kept until enough lighter ones are\n" \
  "                     found).  A single codeword suffices.  The hits counted:\n" \
  "                     with a resolved rw_opt of 0 (rw_opt=0, or e.g. a classical\n" \
  "                     code), all, as in version 0.16.1; otherwise (see the run\n" \
  "                     plan) those of the uniform RW steps, of CC, of finC, and of\n" \
  "                     the windows of kwin > 0 (then no step is uniform), but not\n" \
  "                     those of the row slabs of rw_opt 4 and of the hybrid index\n" \
  "                     windows (kwin=0), which find the same codewords over and\n" \
  "                     over (they are shown apart, see '# RW step types').  While\n" \
  "                     no codeword of the minimum weight has a counted hit, the\n" \
  "                     min_hits stop is inactive (the status line and the summary\n" \
  "                     say so).  A lighter codeword is missed with probability\n" \
  "                     about exp(-<n>), assuming that lighter codewords are found\n" \
  "                     at least as often.  With debug&1, a warning is printed at\n" \
  "                     the end if the hit counts of codewords of the same weight\n" \
  "                     are strongly non-uniform (then exp(-<n>) is optimistic),\n" \
  "                     together with an information-set estimate (uniform random\n" \
  "                     information sets) of the probability to miss a lighter\n" \
  "                     codeword and of the RW steps needed for 1%%.\n" \
  "                     With min_hits=0, RW records only the codewords lighter than\n" \
  "                     the upper bound (unless collecting codewords with outC or\n" \
  "                     maxC), so that the repeated hits are not counted.\n" \
  "                     With the canonical classes of rw_opt (16, automatic for\n" \
  "                     quantum codes, kwin <= 0), RW stops instead once the\n" \
  "                     canonical classes of the minimum weight found in the\n" \
  "                     uniform steps (C of them, Q1 and Q2 found in exactly one\n" \
  "                     and two steps, N step incidences in total) give N >=\n" \
  "                     max(15, 3 min_hits), N/C >= min_hits, the Good-Turing\n" \
  "                     unseen mass Q1/N <= min(0.02, 3 exp(-min_hits)), and the\n" \
  "                     Chao2 unseen fraction C0/(C+C0) <= 0.05 with C0 = Q1(Q1-1)\n" \
  "                     /(2(Q2+1)).  Only the uniform steps count: while none of\n" \
  "                     them has found a codeword of the minimum weight (e.g., on\n" \
  "                     deep circuit DEMs, where the window steps find them), the\n" \
  "                     min_hits stop is inactive and RW runs until the steps or\n" \
  "                     the timeout end (the status line and the summary say so;\n" \
  "                     the window steps are shown for information only, see\n" \
  "                     '# RW step types').  Without 16, <n> above.\n" \
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
  "Study options (the statistics of the RW step types, e.g., the row slabs):\n" \
  "  slab_rows=[int]    The detectors per row slab of rw_opt 4 (default: 0 = m/2;\n" \
  "                     at most m).  The slab options are ignored (with a warning)\n" \
  "                     in a run without row slabs.\n" \
  "  slab_bnd=[int]     The probability in percent of each of the two boundary\n" \
  "                     slabs (the first and the last slab_rows detectors; default:\n" \
  "                     50, i.e., only the boundary slabs; at most 50); the other\n" \
  "                     slabs start at a uniform random detector.\n" \
  "  slab_var=[int]     A random reduction of the slab height by up to slab_var\n" \
  "                     percent: the height of each slab is uniform in [slab_rows\n" \
  "                     (100 - slab_var)/100, slab_rows], at least 2 (default: 0, a\n" \
  "                     fixed height; at most 99).\n" \
  "  slab_every=[int]   A slab on every slab_every-th RW step (default: 2, every\n" \
  "                     other step; at least 2: never slabs only, as the boundary\n" \
  "                     slabs miss the codewords across the two layers); the\n" \
  "                     other steps use uniform random permutations.\n" \
  "  steplog=[file]     Write a line for each full-matrix RW step (ksub=0) and for\n" \
  "                     each codeword recorded in it to file (the statistics of the\n" \
  "                     step types; the file is written once the search starts):\n" \
  "                     'S tid step type a b ns rank ok ncw t_us dt_us': the\n" \
  "                     thread, the global step, its type u (uniform), i (index\n" \
  "                     window of the columns a..b-1 mod n), b (Tanner BFS window:\n" \
  "                     seed column a, b columns), or s (row slab: detectors\n" \
  "                     a..b-1), its ns columns, the rank of its elimination, ok=0\n" \
  "                     if abandoned at a stop, the ncw codeword lines which\n" \
  "                     follow, the time of its end since the start and its\n" \
  "                     duration (microseconds); then\n" \
  "                     'C w key ckey flag o [lsig rlo rhi c_1 ... c_w]': the\n" \
  "                     weight, the hash of the support (the same in all runs on\n" \
  "                     the same matrix), that of its canonical form (see rw_opt\n" \
  "                     16; 0 without the triplets T3 and G), n (new in the\n" \
  "                     codeword hash), o (one more hit of a codeword in hash), or\n" \
  "                     x (not kept), its origin o: 1 a single column (p=1), 2 a\n" \
  "                     pair (p=2), 5 and 6 their shortened forms (see rw_opt 1),\n" \
  "                     and on the first line of a key in the file its L\n" \
  "                     syndrome (hex; a hash with more than 64 rows of L), its\n" \
  "                     smallest and largest detector (row), and its columns\n" \
  "                     (0-based, after merge and drop).  The lines of a step are\n" \
  "                     written together at its end by its thread; a header gives\n" \
  "                     the parameters, and '# end:' the result.  The codewords are\n" \
  "                     logged as RW records them: with min_hits=0, only the\n" \
  "                     lighter ones (min_hits=1000000000 cov_cws=0: every hit of\n" \
  "                     the minimum weight).  Default: none (no log); not with\n" \
  "                     ksub > 0; dist_m4ri_old ignores it.\n\n" \
  "Machine-readable output (opt-in; since version 0.17.4):\n" \
  "  stamp=[0|1]        Append the elapsed search time ' [t=S.SSSs]' (seconds\n" \
  "                     since the start of the search, the clock of timeout and\n" \
  "                     of '# stopped after') to every progress line (debug bit\n" \
  "                     2) and the stop line (debug bit 1).  Default: 0 (off;\n" \
  "                     default output unchanged).  Needs debug bits 1 or 2 for\n" \
  "                     those lines to appear (e.g., debug=3).  dist_m4ri_old\n" \
  "                     ignores it.\n" \
  "  json=[0|1|file]    Emit a single-line JSON summary (schema 'dm4ri-result/1')\n" \
  "                     on every normal exit (including timeout, k=0, and\n" \
  "                     no-search stops; nothing on errors with exit code 255 or\n" \
  "                     signals):\n" \
  "                     0: off (default);\n" \
  "                     1: write '# RESULT {json}' as the last line of stderr\n" \
  "                        (at any debug level);\n" \
  "                     file: write the JSON object (one line) to file (opened\n" \
  "                        after input loading, empty if killed; write ./1 or\n" \
  "                        ./0 for a file named 1 or 0).\n" \
  "                     Repeating json= adds destinations (e.g., json=1\n" \
  "                     json=out.json writes both); json=0 clears both.  Fields:\n" \
  "                     schema, version, dmin, dmax, rw_steps, k0, k0_reason,\n" \
  "                     search, stop {code, text}, t_input (process start to\n" \
  "                     search start), t_stop, t_ub, t_lb (search clock),\n" \
  "                     timeout, seed (resolved), method, threads {requested,\n" \
  "                     cpus, used, spawned, rw, cc}, rw_opt (resolved),\n" \
  "                     rw_opt_requested, rw_steps_done, size {full, searched},\n" \
  "                     cc_rounds [{w, thread_s, wall_s, nodes, threads,\n" \
  "                     complete}], events [[t, dmin, dmax], ...], witness\n" \
  "                     {available, weight, source, base, support (all 1-based\n" \
  "                     original columns)}.  Callers probing an older binary can\n" \
  "                     check 'dist_m4ri --help json' (exit 0 vs 1) or catch\n" \
  "                     exit code 255 with 'unrecognized parameter' on stderr.\n" \
  "                     dist_m4ri_old ignores it.\n\n" \
  "Debug output bitmap (debug=[int], default: 3):\n" \
  "  Diagnostic output goes to stderr (stdout only has the result line); lower\n" \
  "  bits are more informative.  The value is decimal or hexadecimal (0x...).\n" \
  "  Multiple debug arguments are OR-combined, where the first one replaces the\n" \
  "  default: e.g., debug=0 alone is silent, while debug=4 and debug=0 debug=4\n" \
  "  both give only the periodic status.\n" \
  "    0    : Silent: errors, warnings on the validity of the result (expert\n" \
  "           options, experimental ksub, invalid finC codewords) and on ignored\n" \
  "           or overridden parameters, and the confinement profile with smax > 0\n" \
  "    1    : Summary (default): input matrices (and the rounds of a DEM, RW\n" \
  "           windows), warnings, stop reason, codeword and hit statistics,\n" \
  "           RW information sets, pairs, step types, and stages (rw_opt),\n" \
  "           export\n" \
  "    2    : Progress (default): run plan (threads, steps, timeout, seed),\n" \
  "           finished CC rounds, new upper bounds, RW convergence\n" \
  "    4    : Periodic status (after 1, 2, 4, ... s, then every 60 s): bounds,\n" \
  "           RW steps and rate, the min_hits statistic (<n>, or the\n" \
  "           canonical classes of rw_opt 16), progress of the current CC round\n" \
  "    8    : Thread allocation and timing: throttling, CC round planning and\n" \
  "           re-planning, waits for RW, predicted vs measured CC work, the\n" \
  "           nodes of each CC round (calls of the recursion), and in method=3\n" \
  "           the RW workers which left their RW batch to join it\n" \
  "    16   : Code parameters rank(H), rank(L) or rank(G), and k (dense\n" \
  "           elimination at startup, slow for large matrices)\n" \
  "    32   : Support (1-based) of each new lightest codeword (RW, CC, finC)\n" \
  "    64   : Command-line arguments, rng seed, and files read\n" \
  "    128  : Matrices H, G, L (if n < 150, unless 2048 is set)\n" \
  "    256  : All codewords in hash with their hit counts (at the end)\n" \
  "    2048 : No size cutoff for matrices, codeword supports, and per-weight\n" \
  "           codeword statistics\n" \
  "    4096 : dist_m4ri_old only: CC recursion and confinement hash traces\n" \
  "  Bits 512, 1024, and 8192 to 32768 are reserved.  Bits 65536 and above are\n" \
  "  reserved for the Python wrapper dist_m4ri.py (ignored by the binary).\n" \
  "  E.g., debug=7 adds the periodic status to the default output, and debug=11\n" \
  "  the thread allocation and timing.\n\n" \
  "Output format (stdout):\n" \
  "  Standard output produces a single line with three space-separated integers:\n" \
  "    dmin dmax rw_steps\n" \
  "  where:\n" \
  "    dmin-1   : The largest weight w certified to have no codeword of weight\n" \
  "               <= w: CC found none from the supplied dmin (default 1) up to w\n" \
  "               (odd weights are skipped if all codewords have even weight).\n" \
  "               If CC finds an exact minimum-weight codeword, dmin = dmax = weight.\n" \
  "               With the expert options noscan, start, or cbeg/cend, see the\n" \
  "               warnings above for the validity of dmin.\n" \
  "    dmax     : Smallest weight of a codeword found (by RW or CC, or read with\n" \
  "               finC), or the supplied dmax if smaller (0 if none).\n" \
  "    rw_steps : Number of completed RW steps (0 if CC found exact or method=2).\n" \
  "  If k = 0, i.e., there is no codeword (rank H = n for a classical code; every\n" \
  "  row of L in the row space of H, e.g., a DEM with no L left after the pmin\n" \
  "  cut), the result is 0 0 0, with a note on stderr.  Unless k is known at once\n" \
  "  (e.g., L without rows; with finG, k is the number of rows of the computed\n" \
  "  L), one extra thread (also with threads=1) checks the rank on a dense copy\n" \
  "  of H during the search (not if m*n > 2^33, i.e., over 1 GiB): the run ends\n" \
  "  as soon as it finds k = 0, while a codeword found proves k > 0; a search\n" \
  "  without a codeword waits for it.  dist_m4ri_old checks before the search.\n" \
  "  The legacy single-threaded dist_m4ri_old (the same options) prints one number\n" \
  "  and also its own progress and codeword lines to stdout, and ignores threads,\n" \
  "  timeout, nothrottle, chunk_size, dexp, dstop, cc_prune, pairs, rw_opt,\n" \
  "  dem_rounds, stamp, json, and the study options (with a warning).\n\n" \
  "Option groups (--help GROUP prints all entries of GROUP):\n" \
  "  windows            kwin, win_mode, win_hyb\n" \
  "  stopping, stop     steps, wmin, timeout, min_hits, cov_cws\n" \
  "  bounds             dmin, dmax, dstop, dexp, wmax\n" \
  "  dem                fdem, pmin, merge, drop, outD, dem_rounds\n" \
  "  matrices, matrix   finH, finG, finL, fin, classical, ortho\n" \
  "  codewords          outC, finC, maxC, dW, irred\n" \
  "  threading          threads, nothrottle, chunk_size\n" \
  "  expert             noscan, start, start_scan, cbeg, cend\n" \
  "  experimental       ksub, refresh\n" \
  "  study              slab_rows, slab_bnd, slab_var, slab_every, steplog\n" \
  "  output             stamp, json\n\n" \
  "Help options:\n" \
  "  -h, --help         Display help for commonly used parameters (fits 80 rows);\n" \
  "                     --help TOPIC is the same as --morehelp TOPIC.\n" \
  "  --morehelp         Display this full help message listing all parameters;\n" \
  "                     --morehelp TOPIC (also --morehelp=TOPIC) displays only the\n" \
  "                     entries of TOPIC, e.g., --morehelp maxC: one or more names\n" \
  "                     separated by commas (e.g., maxC,dW), case-insensitive,\n" \
  "                     without '-' or '_' (and '=VALUE' in --help=TOPIC and\n" \
  "                     --morehelp=TOPIC: a separate argument with '=' is not a\n" \
  "                     TOPIC), including aliases such as dest or batch (and debug\n" \
  "                     for the debug output bitmap); an unknown name gives the\n" \
  "                     closest names on stderr, and the exit code 1.  A group\n" \
  "                     name (see Option groups, also as singular or plural, e.g.,\n" \
  "                     window) gives all entries of the group, and a parameter of\n" \
  "                     a group also the other entries of its group (see also).\n" \
  "  --version          Display program version\n"

#define BRIEF_HELP \
  "try \"%s -h\" for help"

#endif /* UTIL_IO_H */
