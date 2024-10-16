
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
//#include "dist_m4ri.h"
// #include "util.h"

typedef struct ONE_VEC_T{
  int wei; /** current weight */
  int max; /** allocated */
  int vec[0];
} one_vec_t;

/** @brief print entire `one_vec_t` structure by pointer */
void one_vec_print(const one_vec_t * const pvec){
  printf(" w=%d [ ",pvec->wei);
  for(int i=0; i < pvec->wei; i++)
    printf("%d ",pvec->vec[i]);
  printf("]\n");
}


// v1[:] = v0[:] + mat[row,:]
static inline int one_csr_row_combine(one_vec_t * const v1, const one_vec_t * const v0,
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
    else
      v1->vec[i1++] = ic;
  }
  if(i0 >= v0->wei) /** remaining `mat[row]` entries */
    for (                      ; iM < mat->p[row+1]; iM++){
      const int ic = mat->i[iM];
      v1->vec[i1++] = ic;
    }
  else /** remaining `v0` entries */
    while(i0 < v0->wei)
      v1->vec[i1++] = v0->vec[i0++];
  
  v1->wei = i1;
  return i1; /** weith of the out vector */
}

/** @brief insert `j` (originally absent) into ordered array, return position */
static inline int one_ordered_ins(one_vec_t * const err, const int j){
  int pos=err->wei-1;
  while(j < err->vec[pos]){
    err->vec[pos+1] = err->vec[pos];
    pos--;
  }
#ifndef NDEBUG  
  if (j == err->vec[pos]) 
    ERROR("Unexpected! vec[%d]=%d is already present!",pos,j);
#endif   
  err->vec[pos+1]=j;
#ifndef NDEBUG
  if(err->wei >= err->max)
    ERROR("unexpected!"); /** before increment */
  for(int i=0; i < err->wei; i++)
    if(err->vec[i] >= err->vec[i+1]){
      printf("check ordering at i=%d! ",i); /** before increment */
      one_vec_print(err);
      ERROR("unexpected");
    }
#endif   
  err->wei ++;
  return pos+1;
}


/** @brief find `val` in ordered array 
 *  @return position of `val` if found, -1 otherwise 
*/
static inline int one_ordered_search(one_vec_t * const err, const int val){
  /** binary search for pos of `val` */
  int bot=0, top=err->wei , mid=0;
#ifndef NDEBUG  
  if (!top)
    return -1;
#endif   
  while(top - bot > 1){
    mid = (top+bot) >> 1;
#ifndef NDEBUG  
  if (mid>=err->wei)
    ERROR("this should not happen");
#endif   
    if (err->vec[mid] <= val)
      bot = mid;
    else
      top = mid;
    //    printf("bot=%d mid=%d top=%d err->vec[mid]=%d val=%d\n",bot,mid,top,err->vec[mid],val);
  }
  if ( err->vec[bot] == val)
    return bot;
  else
    return -1;
} 

/** @brief delete `val` (if originally present) from ordered array 
 *  @return 1 of `val` was found, 0 otherwise 
*/
static inline int one_ordered_find_del(one_vec_t * const err, const int val){
  /** binary search for pos of `j` */
  int bot=0, top=err->wei, mid=0;
#ifndef NDEBUG  
  if (!top)
    return 0;
#endif   
  while(top - bot > 1){
    mid = (top+bot) >> 1;
    if (err->vec[mid] <= val)
      bot = mid;
    else
      top = mid;
  }
  if ( err->vec[mid] != val)
    return 0;
    //    ERROR("unexpected! value val=%d not found",val);
  
  for(int i=mid; i < err->wei; i++)
      err->vec[i] = err->vec[i+1];
  err->wei --;
  return 1;
}


/** @brief delete `val` in known position `pos` from ordered array */
static inline void one_ordered_pos_del(one_vec_t * const err, _maybe_unused const int val, const int pos){
#ifndef NDEBUG
  if ((pos<0) || (pos >= err->wei) || (err->wei == 0) || (err->vec[pos] != val))
    ERROR("this should not happen!");
#endif
  err->wei --; 
  for(int i=pos; i < err->wei; i++)
      err->vec[i] = err->vec[i+1];
}

/** @brief recursively construct codewords 
 * 
 * @param err error vector with sorted components 
 * @param urr unsorted vector so far 
 * @param syn array of syndrome vectors with sorted components (indexed by weight of error)
 * @param wmax max recursion level (max weight of an error to process)
 * @param max_col_wt maximum column weight (used to predict early termination)
 * @param mH matrix `H` (check matrix of the code or `Hx` for a CSS code)
 * @param mHT matrix `H` transposed
 * @param mL matrix `L=Lx` for a CSS code, or `NULL` for a classical binary code, used to check whether zero-syndrome error is trivial or not
 * @param p_swei minimum syndrome weight array
 * @param debug bitmap 
 */
int start_CC_recurs(one_vec_t *err, one_vec_t *urr, one_vec_t * const syn[],
		    const int wmax, const int max_col_wt, 
		    const csr_t * const mH, const csr_t * const mHT, const csr_t * const mL, int p_swei[], 
		    const int debug){
  const int w=err->wei;
  int row = syn[w]->vec[0]; /** row with the first non-zero syndrome bit */
#ifndef NDEBUG  
  if(debug&64){
    printf("starting CC recurs w=%d row=%d:\n urr: ",w,row);
    one_vec_print(urr);
    printf(" err: ");
    one_vec_print(err);
    for(int i=0; i <= w; i++){
      printf("i=%d ",i);
      one_vec_print(syn[i]);
    }
  }
#endif   
  const int col_min=urr->vec[0]; /** all valid positions should be to the right of here */
  for(int i1 = mH->p[row]; i1 < mH->p[row+1]; i1++){
    const int col = mH->i[i1];
    if(col > col_min){
      int pos = one_ordered_search(err, col);
      if(pos == -1){ /** not there */
	urr->vec[w] = col;
	urr->wei++;
#ifndef NDEBUG
	if(debug&64){
	  printf(" pos=%d urr: ",pos);
	  one_vec_print(urr);
	}
#endif 	
	pos = one_ordered_ins(err,col);
	int swei = one_csr_row_combine(syn[w+1],syn[w], mHT, col);
	if(p_swei[err->wei] > swei){
	  //#ifndef NDEBUG
	  if(debug&64){
	    printf("# swei[%d]=%d -> %d change\n# err: ",
		   err->wei,p_swei[err->wei],swei);
	    one_vec_print(err);
	    printf("# syn: ");
	    one_vec_print(syn[w+1]);
	  }
	  //#endif 	  
	  p_swei[err->wei]=swei;
	}
	int result = 0;
	if (err->wei < wmax){
	  if (swei){ /** go up */
	    if(swei <= (wmax - err->wei)*max_col_wt){ /** reachable goal? */
	      result = start_CC_recurs(err,urr,syn,wmax,max_col_wt,
				       mH,mHT,mL,p_swei,debug);
	      if(result == 1)
		return 1;
	    }
	  }
	  // swei == 0 means it is a degenerate vector
	}
	else{ // wei == wmax
	  if(!swei){
	    if((!mL) ||  /** classical code */
	       (sparse_syndrome_non_zero(mL, err->wei, err->vec))){
	      if(debug&32){
		printf("swei=%d *** success ***\n",swei);
		one_vec_print(syn[w+1]);
	      }
	      return 1; /** success, just get out fast */
	    }
	  }
	}
	urr->wei--;
	one_ordered_pos_del(err,col,pos);
	if(debug&32){
	  printf("xerr: ");
	  one_vec_print(err);
	}
      }
    }
  }
  if(debug&32)
    printf("exiting CC recurs\n\n");
  return 0; /** nothing found */
}

//! rewrite of the cluster method function using only sparse matrices
//! try recursive version first 
int do_CC_dist(const csr_t * const mH, const csr_t * mL,
	       const int wmax, const int start, int p_swei[], const int debug){

  const int nchk = mH->rows, nvar = mH->cols;
  if((start<-1) || (start>=nvar))
    ERROR("invalid start=%d for H[%d,%d]\n",start, nchk, nvar);

  csr_t * const mHT = csr_transpose(NULL,mH);
  if(debug&32){
    if((mHT->cols<150)||(debug&2048))
      csr_print(mHT,"HT");
  }
  int max_col_W = csr_max_row_wght(mHT);
  
  one_vec_t *err = calloc(1, sizeof(one_vec_t)+sizeof(int)*wmax);  
  one_vec_t *urr = calloc(1, sizeof(one_vec_t)+sizeof(int)*wmax);  
  one_vec_t **syn = calloc(wmax+1, sizeof(one_vec_t *));
  if((!syn) || (!err) || (!urr))
    ERROR("memory allocation");
  err->max = wmax;
  for(int i=0; i <= wmax; i++){
    syn[i]=calloc(1, sizeof(one_vec_t)+sizeof(int)*mH->rows);
    if(!syn[i]) ERROR("i=%d memory allocation",i);
    syn[i]->max = mH->rows;    
  }
  int result = 0;
  for(int w=1; w <= wmax; w++){ /* cluster weight */
    int beg = 0, end = nvar - wmax ;
    if (start >= 0)
      beg = end = start;
    if(debug&2)
      printf("# recursively searching for w=%d codewords wmax=%d beg=%d end=%d\n",w,wmax,beg,end);
    for(int i = beg; i <= end; i++){ /* start column position */
      /** prepare the 1st error vector and the syndrome */
      err->vec[0] = urr->vec[0] = i;
      err->wei = urr->wei = 1;
      int swei = one_csr_row_combine(syn[1], syn[0], mHT, i);
      if(p_swei[1] > swei)
	p_swei[1]=swei;
      if (1<w){
	if (swei){ /** go up */
	  result = start_CC_recurs(err,urr,syn,w,max_col_W,mH,mHT,mL,p_swei,debug);
	  if(result == 1)
	    break;
	}
      }
      else{ // w==1
	if(!swei){	/** verify the vector */
	  if((!mL) ||  /** classical code */
	     (sparse_syndrome_non_zero(mL, err->wei, err->vec))){
	    result = 1; /** success */
	    break;
	  }
	}
      }
      err->wei = urr->wei = 0;
    }
    if(result == 1)
      break;
  }
  
  if(result==1){
    result = err->wei; /** codeword weight found */
    if(debug&16){
      printf("# wmax=%d found cw of weight %d: [",wmax,result);
      int max = ((result<50) || (debug&2048)) ? result : 50 ;
      for(int i=0; i< max; i++)
	printf("%d%s",err->vec[i], i+1!=max?" ": (result==max ? "]\n" : "...]\n"));
    }
  }
  else
    result = -wmax; /** not found a codeword up to wmax */

  if(debug&2){
    for(int i=1;i<=wmax; i++)
      if(p_swei[i] <= mH->rows)
	printf("# w=%d min syndrome weight %d\n",i,p_swei[i]);
  }
  
  for(int i=0; i<= wmax; i++)
    free(syn[i]);
  free(syn);
  free(err);
  free(urr);
  csr_free(mHT);
  return result;
}


