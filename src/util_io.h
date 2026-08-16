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
  double timeout; /* timeout in seconds (default 60.0) */
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

#define USAGE								\
  "%s: distance of a classical or quantum CSS code\n"			\
  "\tusage: %s parameter=value [...]\n\n"				\
  "   Required parameter:\n"						\
  "\tmethod=[int]: bitmap for method used (no default): \n" \
  "\n"									\
  "\t\t1: random window (RW) algorithm. Options:\n"			\
  "\t\t   steps=[int]: how many information sets to use (1000)\n"		\
  "\t\t   wmin=[int]:  minimum distance of interest (1)\n"		\
  "\t\t\t immediately stop and return '-w' on a cw of weight w<=wmin\n" \
  "\t\t\t use this option to quickly scan over a large number of codes\n" \
  "\n"									\
  "\t\t2: connected cluster (CC) algorithm.  Options:\n"		\
  "\t\t   wmax=[int]:  maximum cluster weight to construct, inclusive (0)\n" \
  "\t\t\t must be non-zero for CC only, otherwise use upper bound from RW\n" \
  "\t\t   smax=[int]:  maximum syndrome weight of interest, inclusive (5)\n" \
  "\t\t\t must be non-zero to calculate confinement profile\n"          \
  "\t\t   start=[int]: use only this position to start (equiv. to cbeg=cend=start) (-1)\n" \
  "\t\t   cbeg=[int]:  start column to begin CC search (-1)\n"		\
  "\t\t   cend=[int]:  end column to limit CC search (-1)\n"		\
  "\t\t   noscan=[int]: start CC directly with wmax (0)\n" \
  "\t\t3: bracketing mode (balanced concurrent RW and CC)\n" \
  "\n"									\
  "   Multithreading parameters (distfork):\n"				\
  "\tthreads=[int]: number of threads to use (0 for auto) (0)\n"	\
  "\tdexp=[int]:    expected distance value (alias: dest) (0)\n"	\
  "\ttimeout=[sec]: timeout in seconds (60)\n"				\
  "\n"									\
  "   Distance bounds parameters:\n"					\
  "\tdmin=[int]:    known lower bound on distance, inclusive (w starts from dmin in CC) (1)\n" \
  "\tdmax=[int]:    known upper bound on distance, inclusive (RW ignores codewords of weight >= dmax) (0)\n" \
  "\n"									\
  "   General parameters:\n"						\
  "\tfdem=[str]: detector error model (DEM) file from stim (NULL)\n" \
  "\tpmin=[float]: minimum error probability to keep for DEM (0.0)\n" \
  "\tfinH=[str]: parity check matrix Hx (NULL)\n"			\
  "\tfinG=[str]: matrix Hz (quantum CSS code only) (NULL)\n"		\
  "\tfinL=[str]: matrix Lx (quantum CSS code only) (NULL)\n"		\
  "\t\t Either L=Lx or G=Hz matrix is required for a quantum CSS code\n" \
  "\tfin=[str]:  base name for input files (\"try\")\n"			\
  "\t\t set finH->\"${fin}X.mtx\"  finG->\"${fin}Z.mtx\"\n"		\
  "\tcss=[int]:  reserved for future use (1)\n"				\
  "\tseed=[int]: rng seed [use 0 for time(NULL)] (0)\n"			\
  "\tdebug=[int]:\t bitmap for aux information to output (3)\n"		\
  "\t\t0: clear the entire debug bitmap to 0.\n"			\
  "\t\t1: output misc general info (on by default)\n"			\
  "\t\t2: output more general info (on by default)\n"			\
  "\t\t4: debug command line arguments parsing\n"			\
  "\t\t8: output progress reports every 1000 steps\n"			\
  "\t\t16: output new min-weight codewords found (cut large vectors)\n"	\
  "\t\t32: output matrices (unless n is large)\n"			\
  "\t\t64: debug confinement hash updates (swei changes)\n" \
  "\t\t128: debug duplicate syndromes in confinement hash (debug build only)\n" \
  "\t\t256: reserved\n"							\
  "\t\t512: reserved\n"							\
  "\t\t1024: reserved\n"						\
  "\t\t2048: allow big matrix / large vector output\n"			\
  "\t\t   see the source code for more options\n"			\
  "\t  Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t  Use debug=0 as the 1st argument to suppress all debug messages.\n"\
  "   -h gives this help (also '--help')\n"

#define BRIEF_HELP				\
  "try \"%s -h\" for help"	       

#endif /* UTIL_IO_H */
