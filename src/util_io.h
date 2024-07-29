/* 	$Id:  $	 */
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

#define USAGE_IO						\
  "Supported input parameters:\n"				\
  "\tdebug=1: bitmap for aux information\n"			\
  "\tstring: fin=\"try\": base name for input files\n"		\
  "\tstring: finP=\"tryX.mtx\": check matrix\n"			\
  "\tstring: finG=\"tryZ.mtx\": dual generator matrix\n"	\
  "\tswitch=1: switch X and Z generators\n"			\
  "\tcss=1: this is a CSS code (the only supported one)\n"

typedef struct{
  int css; /* 1: css, 0: non-css -- currently not supported */
  //  int linear; /* not supported */
  int debug; /* debug information */ 
  //  char *finP, *finG; /* generators */
  int n0;  /* code length, =n for css, (n/2) for non-css */
  int n; /* actual n = matrix size */
} params_io_t; 

extern params_io_t pio;

int local_io_init(int argc, char **argv, csr_t **spaP, csr_t **spaG);
/* initialized variables needed and the sparse matrices P and G;
 * return the position of the first argument which failed to interpret
 * (to be treated by a subsequent routine specific to the program) 
 */

#endif /* UTIL_IO_H */
