
/** ********************************************************************** 
 * @brief distance of a classical or quantum CSS code
 * 
 * The program implements two methods:
 * 1. Random information set (random window) algorithm (upper bound).  
 *    This works with any code (LDPC or not).
 * (2) depth-first codeword enumeration (connected cluster) algorithm
 * (Lower bound or actual distance if a codeword is found.)  
 * 
 * A. Dumer, A. A. Kovalev, and L. P. Pryadko "Distance verification..."
 * in IEEE Trans. Inf. Th., vol. 63, p. 4675 (2017). 
 * doi: 10.1109/TIT.2017.2690381
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>, Weilei Zeng
 ************************************************************************/
// #include <m4ri/config.h>
#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <time.h>
#include <m4ri/m4ri.h>

#include "mmio.h"
#include "util_m4ri.h"
#include "util_io.h"
#include "dist_m4ri.h"

#ifdef STANDALONE

int main(int argc, char **argv){
  params_t * const p = &prm;

  var_init(argc,argv,p);

  printf("# test functions here\n");

  printf("# test finished\n");
  var_kill(p);
  
  return 0;
}

#endif /* STANDALONE */
