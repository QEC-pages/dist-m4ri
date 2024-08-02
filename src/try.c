
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
// #include "util.h"

typedef struct ONE_VEC_T{
  int wei; /** current weight */
  int max; /** allocated */
  int vec[0];
} one_vec_t;

// v1[:] = v0[:] + mat[row,:]
static inline int csr_row_combine(one_vec_t * const v1, const one_vec_t * const v0,
				  const csr_t * const mat, const int row){
#ifndef NDEBUG
  if ((!v1) || (!v0) || (!mat))
    ERROR("all arguments must be allocated: v1=%p v0=%p mat=%p\n",v1,v0,mat);
  if(v1 == v0)
    ERROR("the two vectors should not be the same !");
  if((row<0) || (row >= mat->rows) ||
     (v1->max < mat->cols) ||
     (v0->max < mat->cols))
    ERROR("this should not happen\n");
#endif
  int iM, i0=0, i1=0; /** iterators */
  for (iM = mat->p[row]; iM < mat->p[row+1]; iM++){
    const int ic = mat->i[iM];
    while((i0 < v0->wei) && (v0->vec[i0] < ic))
      v1->vec[i1++] = v0->vec[i0++];
    if(i0 >= v0->wei)
      break;
    if(v0->vec[i0]==ic)
      i0++; /** `1+1=0` just skip this position */
  }
  for (                      ; iM < mat->p[row+1]; iM++){
    const int ic = mat->i[iM];
    v1->vec[i1++] = ic;
  }
  v1->wei = i1;
  return i1; /** weith of the out vector */
}

static int start_CC_rec(const int w, const int wmax, int err[], one_vec_t * const syn[],
			const csr_t * const mH, const csr_t * const mHT, const mzd_t * const L){  
  int row = syn[w-1]->vec[0]; /** row with the first non-zero syndrome bit */
  int col0 = err[w-1]; /** position already set */
  for(int i1 = mH->p[row]; i1 < mH->p[row+1]; i1++){
    const int col = mH->i[i1];
    if(!bsearch(&col, err->vec, err->wei, sizeof(int), cmp_rci_t)){
      /** TODO: make a copy of `vec` with ordered positions */
    if(col0 != col)
  }
}


//! rewrite of the cluster method function using only sparse matrices
//! try recursive version first 
int do_CC_dist(const csr_t * const mH, const csr_t * const mHT, const csr_t * mL,
	   const int wmax, const int start, const int debug){

  const int nchk = mH->cols, nvar = mH->rows;
  int lev;
  if((start<0) || (start>=nvar))
    ERROR("invalid start=%d for H[%d,%d]\n",start, nchk, nvar);

  one_vec_t *err = calloc(1, sizeof(one_vec_t)+sizeof(int)*wmax);  
  one_vec_t **syn = calloc(wmax+1, sizeof(one_vec_t *));
  if((!syn) || (!err))
    ERROR("memory allocation");
  err->max = wmax;
  for(int i=0; i <= wmax; i++){
    syn[i]=calloc(1, sizeof(one_vec_t)+sizeof(int)*mH->rows);
    if(!syn[i]) ERROR("i=%d memory allocation",i);
    syn[i]->max = mH->rows;
  }
  
  for(int w=2; w <= wmax; w++){ /* cluster weight */
    int beg = 0, end = nvar - 1;
    if (start >= 0)
      beg = end = start;
    
    for(int i = beg; i <= end; i++){ /* start column position */
      err->vec[0]=i;
      err->wei=1;
      lev=1;
      int syn_wt = csr_row_combine(syn[1], syn[0], mHT, i);
      if(!syn_wt){
	/** verify the vector */
	int nz = sparse_syndrome_non_zero(mL, err->wei, err->vec);
	/** TODO: check the syndrome vector with NDEBUG ??? */
	if(nz)
	  goto all_done;
      }
      else{
	/** go up */
      }
      syn_wt = csr_row_combine(syn[1], syn[0], mHT, i);
      assert(syn_wt==0);
    }
  }

  
 all_done:
  for(int i=0; i<= wmax; i++)
    free(syn[i]);
  free(syn);
  free(err);
  return lev;

}


