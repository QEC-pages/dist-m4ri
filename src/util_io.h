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
  int dmax; /** known upper bound on distance (RW ignores codewords of weight >= dmax unless collecting) */
  int wmin; /** min distance below which we are not interested 
		if w <= wmin found in RW, terminate immediately 
		start clusters with `wmin` for `CC`
	     */
  int noscan; /** 1: start CC directly with wmax (no scan over w) */
  int seed;/* rng seed, set=0 for automatic */
  int dist; /* target distance of the code */
  int dist_max; /* distance actually checked */
  int dist_min; /* distance actually checked */
  int max_row_wgt_H; /* needed for C */
  int max_col_wgt_H; /* needed ? */
  //! int max_row_wt;  /* WARNING: this is defined in `util_io.h` as `static const int` */
  int swei[MAX_W]; /** minimum syndrome weight for each error weight */
  int start;
  int cbeg;
  int cend;
  //  int linear; /* not supported */
  int n0;  /* code length, =nvar for css, (nvar/2) for non-css */
  int nvar; /* actual n = matrix size */
  int nchk; /* actual k = number of codewords */
  long long int maxC;
  int dW;
  char *finC;
  char *outC;
  cw_vec_t *codewords;
  long long int num_cws;
  int min_w;
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
 * @brief Write the found codewords from the hash table to a .nz file.
 * 
 * Exports all codewords currently stored in the hash table to a file in NZLIST format.
 *
 * @param fnam Path to the output .nz file.
 * @param comment An optional comment string to include in the file header.
 * @param p Pointer to the params_t structure containing the codeword hash.
 * @return Number of codewords written, or -1 on error.
 */
long long int nzlist_write(const char fnam[], const char comment[], params_t *p);

/**
 * @brief Add a candidate codeword to the hash table if it meets weight limits.
 * 
 * Compares the candidate codeword weight with the current minimum weight and dW limit.
 * If it is within the limits, it is added to the hash. If a new strictly smaller minimum
 * weight is found, it updates the global minimum weight and prunes heavier codewords
 * from the hash.
 *
 * @param p Pointer to the params_t structure.
 * @param arr Array of indices representing the support of the codeword.
 * @param weight Weight of the codeword (length of arr).
 * @return Pointer to the added/existing codeword structure, or NULL if not added.
 */
cw_vec_t * codeword_add_maybe(params_t * const p, const int arr[], int weight);

#define DIST_M4RI_VERSION "0.9.0"

/**
 * @brief Print short help message listing all allowed parameters to stderr.
 *
 * @param prog Program name (argv[0]).
 */
void print_short_help(const char *prog);

#define SHORT_HELP \
  "%s (version %s): calculate distance of a classical or quantum CSS code\n" \
  "Usage: %s method=[1|2|3] [parameter=value ...]\n\n" \
  "Allowed parameters:\n" \
  "  method, finH, finG, finL, fin, fdem, pmin, classical, css,\n" \
  "  steps, wmin, wmax, dmin, dmax, dexp (dest), smax, start, cbeg,\n" \
  "  cend, noscan, threads, timeout, nothrottle, chunk_size (batch),\n" \
  "  finC, outC, maxC, dW, seed, debug\n\n" \
  "Help options:\n" \
  "  -h, --help    : display help for commonly used parameters (fits 80 rows)\n" \
  "  --morehelp    : display full help for all available parameters\n"

#define USAGE \
  "%s (version %s): calculate distance of a classical or quantum CSS code\n" \
  "Usage: %s method=[1|2|3] [parameter=value ...]\n\n" \
  "Required parameter:\n" \
  "  method=[int]       1: Random Window (RW) algorithm (upper bound)\n" \
  "                     2: Connected Cluster (CC) algorithm (lower bound / exact)\n" \
  "                     3: Bracketing mode (concurrent RW and CC)\n\n" \
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
  "  steps=[int]        Maximum RW decoding steps / information sets (1000)\n" \
  "  wmax=[int]         Maximum cluster weight to search in CC (0=until bound/timeout)\n" \
  "  wmin=[int]         Stop immediately if cw with weight <= wmin is found (1)\n" \
  "  timeout=[sec]      Execution timeout in seconds, 0 for infinite (60.0)\n\n" \
  "Multithreading:\n" \
  "  threads=[int]      Max worker threads to use (0: auto CPU count) (0)\n" \
  "                     (subject to throttling unless nothrottle=1)\n\n" \
  "Codeword collection:\n" \
  "  outC=[file]        Export found minimum-weight codewords to file (.nz format)\n" \
  "  finC=[file]        Import initial codewords from file (.nz format)\n" \
  "  maxC=[int]         Maximum number of codewords to collect (0 for unlimited) (0)\n" \
  "  dW=[int]           Collect codewords up to weight dmin + dW (default: 0)\n\n" \
  "Extra parameters (see --morehelp for details):\n" \
  "  smax=[int] (5)         Max syndrome weight for confinement profile (0 to disable)\n" \
  "  noscan=[0|1] (0)       CC method 2: start directly at wmax, skip scanning w<wmax\n" \
  "  start/cbeg/cend=[int]  Limit CC search to specific column(s) (-1: all)\n" \
  "  nothrottle=[0|1] (0)   Disable thread throttling (also --no-throttle)\n" \
  "  chunk_size=[int] (0)   RW batch chunk size (0: auto, alias: batch)\n" \
  "  seed=[int] (0)         RNG seed [0 for time(NULL)]\n" \
  "  debug=[int] (3)        Debug bitmask (0: silent, 1: general, 2: verbose, ...)\n" \
  "  css=[int] (1)          Reserved for future use\n\n" \
  "Help options:\n" \
  "  -h, --help         Display this help message (commonly used parameters)\n" \
  "  --morehelp         Display full help with all parameter descriptions\n" \
  "  --version          Display program version\n"

#define MORE_HELP \
  "%s (version %s): calculate distance of a classical or quantum CSS code\n" \
  "Usage: %s method=[1|2|3] [parameter=value ...]\n\n" \
  "Required parameter:\n" \
  "  method=[int]       Bitmap / identifier for calculation method (no default):\n" \
  "                     1: Random Window (RW) algorithm (upper bound)\n" \
  "                        Finds an upper bound on distance by testing random\n" \
  "                        information sets. Fast for finding small codewords.\n" \
  "                     2: Connected Cluster (CC) algorithm (lower bound / exact)\n" \
  "                        Exhaustive cluster search finding a lower bound or exact\n" \
  "                        minimum distance. Guaranteed to find the true code\n" \
  "                        distance if run to completion.\n" \
  "                     3: Bracketing mode (concurrent RW and CC)\n" \
  "                        Dynamically balances CC (lower bound) and RW (upper\n" \
  "                        bound) worker threads based on distance estimate (dexp),\n" \
  "                        current bounds [dmin, dmax], and remaining timeout.\n\n" \
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
  "                     unless collecting codewords.\n" \
  "  dexp=[int]         Expected code distance (alias: dest) (default: 0).\n" \
  "                     Used in method=3 (bracketing) to balance worker threads\n" \
  "                     between CC and RW and estimate search feasibility.\n\n" \
  "Search limits and Connected Cluster (CC) options:\n" \
  "  steps=[int]        Maximum number of RW decoding steps / information sets\n" \
  "                     (default: 1000). Ignored in method=2.\n" \
  "  wmin=[int]         Minimum distance threshold (default: 1).\n" \
  "                     If a codeword of weight w <= wmin is found, execution\n" \
  "                     terminates immediately. Useful for screening codes.\n" \
  "  wmax=[int]         Maximum cluster weight to analyze in CC (default: 0).\n" \
  "                     In method=2, CC terminates after checking weight wmax.\n" \
  "                     0 means continue until codeword found, bounds meet, or timeout.\n" \
  "  smax=[int]         Maximum syndrome weight for confinement profile (default: 5).\n" \
  "                     When smax > 0, tracks minimum syndrome weights for each\n" \
  "                     error weight. Set smax=0 to disable confinement calculation.\n" \
  "  noscan=[0|1]       1: Start CC directly at weight wmax, skipping scan over\n" \
  "                     weights w < wmax (default: 0). Only valid for method=2.\n" \
  "  start=[int]        Restrict CC search to start column index (default: -1).\n" \
  "                     Equivalent to setting cbeg=start cend=start.\n" \
  "  cbeg=[int]         Beginning column index for CC search (default: -1, start at 0).\n" \
  "  cend=[int]         Ending column index for CC search (default: -1, end at n-1).\n\n" \
  "Codeword collection and export:\n" \
  "  outC=[file]        Export found codewords to file in .nz list format.\n" \
  "  finC=[file]        Import initial candidate codewords from file in .nz format.\n" \
  "  maxC=[int]         Maximum number of codewords to collect (default: 0).\n" \
  "                     0 means collect all valid codewords found up to wmax or stop\n" \
  "                     flag. If maxC > 0, halts when maxC codewords are collected.\n" \
  "  dW=[int]           Extra weight window above minimum distance (default: 0).\n" \
  "                     Collects codewords with weight up to min_w + dW.\n\n" \
  "Execution, multithreading, and timing:\n" \
  "  threads=[int]      Maximum number of worker threads to use (default: 0 =\n" \
  "                     hardware concurrency). Subject to throttling for small\n" \
  "                     codes or large matrices unless nothrottle=1 is set.\n" \
  "  timeout=[sec]      Execution timeout in seconds (default: 60.0, 0 = infinite).\n" \
  "                     In method=3, dynamically guides CC vs RW thread balance.\n" \
  "  nothrottle=[0|1]   Disable automatic thread throttling (default: 0).\n" \
  "                     Aliases: --no-throttle, -no-throttle, nothrottle.\n" \
  "                     By default, threads are throttled for very small codes or\n" \
  "                     large matrices to avoid cache and memory bus thrashing.\n" \
  "  chunk_size=[int]   RW batch chunk size per worker (default: 0 = adaptive).\n" \
  "                     Alias: batch=[int].\n" \
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
  "    dmax     : Smallest weight of any codeword found (0 if no codeword found).\n" \
  "    rw_steps : Number of completed RW steps (0 if CC found exact or method=2).\n\n" \
  "Help options:\n" \
  "  -h, --help         Display help for commonly used parameters (fits 80 rows)\n" \
  "  --morehelp         Display this full help message listing all parameters\n" \
  "  --version          Display program version\n"

#define BRIEF_HELP \
  "try \"%s -h\" for help"

#endif /* UTIL_IO_H */
