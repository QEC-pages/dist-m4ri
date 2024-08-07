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
#include "util_m4ri.h"

#define _maybe_unused __attribute__((unused))

//static const int max_row_wt=10; 

typedef struct{
  int debug; /* debug information */ 
  int classical; /* 1 for a classical code, i.e., no `G=Hz` matrix*/
  int css; /* 1: css, 0: non-css -- currently not supported */
  int method; /* bitmap. 1: random window; 2: cluster; 3: both */
  int steps; /* how many RW decoding steps */
  int wmax; /* max cluster size to try */
  int wmin; /* min distance below which we are not interested at all */
  int seed;/* rng seed, set=0 for automatic */
  int dist; /* target distance of the code */
  int dist_max; /* distance actually checked */
  int dist_min; /* distance actually checked */
  int max_row_wgt_H; /* needed for C */
  int max_col_wgt_H; /* needed ? */
  //! int max_row_wt;  /* WARNING: this is defined in `util_io.h` as `static const int` */
  int start;
  //  int linear; /* not supported */
  int n0;  /* code length, =nvar for css, (nvar/2) for non-css */
  int nvar; /* actual n = matrix size */
  int nchk; /* actual k = number of codewords */
  int maxC;
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

#define USAGE								\
  "%s: distance of a classical or quantum CSS code\n"			\
  "\tusage: %s parameter=value [...]\n\n"				\
  "   Required parameter:\n"						\
  "\tmethod=[int]: bitmap for method used (required, default 0: none): \n" \
  "\t\t1: random window (RW) algorithm. Options:\n"			\
  "\t\t   steps=[int]: how many information sets to use (1)\n"		\
  "\t\t   wmin=[int]:  minimum distance of interest (1)\n"		\
  "\t\t2: connected cluster (CC) algorithm.  Options:\n"		\
  "\t\t   wmax=[int]:  maximum cluster weight (5) \n"			\
  "\t\t   start=[int]: use only this position to start (-1)\n\n"	\
  "   General parameters:\n"						\
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
