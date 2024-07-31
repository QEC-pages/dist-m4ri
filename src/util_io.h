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
  int max_row_wgt_G; /* needed for C */
  //! int maxrow;  /* WARNING: this is defined in `dist_m4ri.h` as `static const int` */
  int start;
  //  int linear; /* not supported */
  int n0;  /* code length, =nvar for css, (nvar/2) for non-css */
  int nvar; /* actual n = matrix size */
  int nchk; /* actual k = number of codewords */
  int swait;
  int maxC;
  char *finH;
  char *finG;
  char *finL;
  char *fin;
  csr_t *spaH;
  csr_t *spaG;
  csr_t *spaL;
} params_t; 

extern params_t prm;
void var_init(int argc, char **argv, params_t * const p);
void var_kill(params_t * const p);

#define USAGE								\
  "%s: calculate the minumum distance of a q-LDPC code\n"		\
  "\tusage: %s [arguments [...]]\n"					\
  "Supported parameters:\n"						\
  "\tdebug=[int]:\t bitmap for aux information (3)\n"			\
  "\tfin=[string]: base name for input files (\"try\")\n"		\
  "\t\t finH->\"${try}X.mtx\"  finG->\"${try}X.mtx\"\n"			\
  "\tfinH=[str]: parity check matrix Hx (NULL)\n"			\
  "\tfinG=[str]: matrix Hz or NULL for classical code (NULL)\n"		\
  "\tfinL=[str]: matrix Lx or NULL for classical code (NULL)\n"		\
  "\t\t Either L=Lx or G=Hz matrix is required for a quantum CSS code\n" \
  "\tcss=1: this is a CSS code (the only supported one) (1)\n"		\
  "\tseed=[int]: rng seed  [0 for time(NULL)]\n"			\
  "\tmethod=[int]: bitmap for method used: \n"				\
  "\t\t1: random window (RW) algorithm\n"				\
  "\t\t2: cluster (C) algorithm\n"					\
  "\tsteps=[int]: how many RW decoding cycles to use (1)\n"		\
  "\twmax=[int]: max cluster weight in C (5) \n"			\
  "\twmin=[int]: min distance of interest in RW (1)\n"			\
  "\t-h or --help gives this help\n"



#endif /* UTIL_IO_H */
