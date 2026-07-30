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
#include <m4ri/m4ri.h>

#include "mmio.h"
#include "uthash.h"
#include "util_hash.h"
#include "util_m4ri.h"

#define _maybe_unused __attribute__((unused))

//static const int max_row_wt=10; 

#define MAX_W 100 
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
  int dmax; /** only look for errors of weight < dmax for `RW` warning: do not
              set this up unless you know what you are doing!  This should only
              be an actual upper bound on the distance from a previous
              calculation!  Setting it will speed things up a little, but may
              make results invalid! */
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
  //  int linear; /* not supported */
  int n0;  /* code length, =nvar for css, (nvar/2) for non-css */
  int nvar; /* actual n = matrix size */
  int nchk; /* actual k = number of codewords */
  int maxC;
  char *fdem;
  double pmin;
  char *finH;
  char *finG;
  char *finL;
  char *fin;
  csr_t *spaH;
  csr_t *spaG;
  csr_t *spaL;
} params_t;

static inline int minint(const int a, const int b) { return (a < b) ? a : b; }
// #define MININT(a,b) do{ int t1=(a); int t2=(b); t1<t2? t1 :t2; } while(0)

extern params_t prm;
void var_init(int argc, char **argv, params_t * const p);
void var_kill(params_t * const p);
void read_dem_file(char *fnam, csr_t **p_spaH, csr_t **p_spaL, double pmin, int debug);

#define USAGE								\
  "%s: distance of a classical or quantum CSS code\n"			\
  "\tusage: %s parameter=value [...]\n\n"				\
  "   Required parameter:\n"						\
  "\tmethod=[int]: bitmap for method used (no default): \n" \
  "\n"									\
  "\t\t1: random window (RW) algorithm. Options:\n"			\
  "\t\t   steps=[int]: how many information sets to use (1)\n"		\
  "\t\t   wmin=[int]:  minimum distance of interest (1)\n"		\
  "\t\t\t immediately stop and return '-w' on a cw of weight w<=wmin\n" \
  "\t\t\t use this option to quickly scan over a large number of codes\n" \
  "\t\t   dmax=[int]:  if non-zero, ignore vectors of this and larger wgt (0)\n" \
  "\t\t\t this option accelerates the search somewhat\n" \
  "\n"									\
  "\t\t2: connected cluster (CC) algorithm.  Options:\n"		\
  "\t\t   wmax=[int]:  maximum cluster weight to construct, inclusive (0)\n" \
  "\t\t\t must be non-zero for CC only, otherwise use upper bound from RW\n" \
  "\t\t   smax=[int]:  maximum syndrome weight of interest, inclusive (20)\n" \
  "\t\t\t must be non-zero to calculate confinement profile\n"          \
  "\t\t   start=[int]: use only this position to start (-1)\n"		\
  "\t\t   noscan=[int]: start CC directly with wmax (0)\n" \
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
  "\t\t64: reserved\n"							\
  "\t\t128: reserved\n"							\
  "\t\t256: print out neighbor lists\n"					\
  "\t\t512: print out vectors/syndrome weights during recursion\n"	\
  "\t\t1024: print piv/skip_pivs/reserved\n"						\
  "\t\t2048: allow big matrix / large vector output\n"			\
  "\t\t   see the source code for more options\n"			\
  "\t  Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t  Use debug=0 as the 1st argument to suppress all debug messages.\n"\
  "   -h gives this help (also '--help')\n"

#define BRIEF_HELP				\
  "try \"%s -h\" for help"	       

#endif /* UTIL_IO_H */
