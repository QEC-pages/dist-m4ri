
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
#include "uthash.h"
#include "util_hash.h"
#include "util_m4ri.h"
#include "util_io.h"
#include "dist_cc.h"

/** @brief print entire `one_vec_t` structure by pointer */
void one_vec_print(const one_vec_t * const pvec){
  printf(" w=%d [ ",pvec->wei);
  for(int i=0; i < pvec->wei; i++)
    printf("%d ",pvec->vec[i]);
  printf("]\n");
}

two_vec_t *errors=NULL;

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
		    const int w_limit, const int max_col_wt, 
		    const csr_t * const mH, const csr_t * const mHT,
                    params_t * const p){
  const int w=err->wei;
  assert(urr->wei == err->wei);
  int row = syn[w]->vec[0]; /** row with the first non-zero syndrome bit */
  const csr_t * const mL = p->spaL;
  int *p_swei = p->swei;
  const int smax = p->smax;
  const int debug = p->debug;
#ifndef NDEBUG  
  if(debug&64){
    printf("starting CC recurs w=%d row=%d:\n urr: ",w,row);
    one_vec_print(urr);
    printf(" err: ");
    one_vec_print(err);
    for(int i=0; i <= w; i++){
      printf("# syn[i=%d] ",i);
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
	syn[w+1]->wei=0;
	int swei = one_csr_row_combine(syn[w+1],syn[w], mHT, col);
	if(debug&64){
	  printf(" syn: ");
	  one_vec_print(syn[w+1]);
	}
	int result = 0;
        int current_limit = w_limit;
        if (p->min_w != INT_MAX && p->dW >= 0) {
          current_limit = minint(w_limit, p->min_w + p->dW);
        }
	if (err->wei < current_limit){
	  if (swei){ /** go up */
//	    if(swei <= (w_limit - err->wei)*max_col_wt){ /** reachable goal? */
	      result = start_CC_recurs(err,urr,syn,w_limit,max_col_wt,
				       mH,mHT,p);
	      if(result == 1)
		return 1;
//	    }
	  }
	  // swei == 0 means it is a degenerate vector
	  // do not go up in this case 
	}
	else{ // wei >= current_limit
	  assert(err->wei == current_limit);
	  if(!swei){
	    if((!mL) ||  /** classical code */
	       (sparse_syndrome_non_zero(mL, err->wei, err->vec))){
	      if(debug&32){
		printf("swei=%d *** success ***\n",swei);
		one_vec_print(err);
		one_vec_print(syn[w+1]);
	      }
              p->codewords = codeword_add_maybe(p, err->vec, err->wei);
              if (p->maxC && p->num_cws >= p->maxC) {
                return 1;
              }
              if (!p->outC && p->maxC == 0) {
                return 1;
              }
	    }
	  }
	  else if(swei <= smax){/** update p_swei if not in hash yet */
	    if(debug&64){
	      printf("# try adding to the hash ewei=%d swei=%d smax=%d p_swei[%d]=%d\n",err->wei,swei,smax,w+1,p_swei[w+1]);
	    }
	    //	    
	    errors = hash_add_maybe(syn[w+1],err,errors, p_swei, debug);
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
//! p_swei[]: min syndrome weight distribution to return (`confinement`).
int do_CC_dist(params_t * const p){
  const csr_t * const mH = p->spaH;
  const csr_t * const mL = p->spaL;
  const int wmax = p->wmax;
  const int noscan = p->noscan;
  int *p_swei = p->swei;
  const int smax = p->smax;
  const int debug = p->debug;

  const int nvar = mH->cols;

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
  //  err->max = wmax;
  for(int i=0; i <= wmax; i++){
    syn[i]=calloc(1, sizeof(one_vec_t)+sizeof(int)*mH->rows);
    if(!syn[i]) ERROR("i=%d memory allocation",i);
    //    syn[i]->max = mH->rows;    
  }
  int result = 0;
  const int w_start = noscan ? wmax : (p->dmin > 1 ? p->dmin : 1);
  int w_limit_dynamic = wmax;
  if (p->dmax > 0) {
    if (p->outC && p->dW > 0) {
      w_limit_dynamic = minint(wmax > 0 ? wmax : p->dmax + p->dW, p->dmax + p->dW);
    } else {
      w_limit_dynamic = minint(wmax > 0 ? wmax : p->dmax, p->dmax);
    }
  }
  for(int w=w_start; w <= w_limit_dynamic; w++){ /* cluster weight */
    if (p->min_w != INT_MAX && p->dW >= 0) {
      w_limit_dynamic = minint(wmax > 0 ? wmax : p->min_w + p->dW, p->min_w + p->dW);
    }
    if (w > w_limit_dynamic) {
      break;
    }
    int beg = (p->cbeg >= 0) ? p->cbeg : 0;
    int end = (p->cend >= 0) ? minint(p->cend, nvar - w) : nvar - w;
    if(debug&2)
      printf("# recursively searching for w=%d codewords wmax=%d beg=%d end=%d\n",w,wmax,beg,end);
    for(int i = beg; i <= end; i++){ /* start column position */
      /** prepare the 1st error vector and the syndrome */
      err->vec[0] = urr->vec[0] = i;
      err->wei = urr->wei = 1;
      syn[1]->wei=0;
      int swei = one_csr_row_combine(syn[1], syn[0], mHT, i);
      if (w>1){
	if (swei){ /** go up */
	  result = start_CC_recurs(err,urr,syn,w,max_col_W,mH,mHT,p);
	  if(result == 1)
	    break;
	}
      }
      else{ // w==1
	if(!swei){	/** verify the vector */
	  if((!mL) ||  /** classical code */
	     (sparse_syndrome_non_zero(mL, err->wei, err->vec))){
            p->codewords = codeword_add_maybe(p, err->vec, err->wei);
            if (p->maxC && p->num_cws >= p->maxC) {
              result = 1;
              break;
            }
            if (!p->outC && p->maxC == 0) {
              result = 1;
              break;
            }
	  }
	}
	else if(swei <= smax) /** update p_swei if not in hash yet */
	  errors = hash_add_maybe(syn[1],err,errors,p_swei, debug);
      }
      err->wei = urr->wei = 0;
    }
    if(result == 1)
      break;

    if (p->min_w > w && w < w_limit_dynamic) {
      printf("-%d\n", w);
      fflush(stdout);
    } else if (p->min_w <= w && p->dW >= 0 && w <= w_limit_dynamic) {
      if (debug & 1) {
        fprintf(stderr, "# CC round w=%d finished: searched with dW=%d (min_w=%d, total %lld codewords)\n",
                w, p->dW, p->min_w, p->num_cws);
      }
    }
  }
  
  if(result==1){
    result = err->wei; /** codeword weight found */
#ifndef NDEBUG
    assert(  0 == sparse_syndrome_non_zero(mH, err->wei, err->vec)); 
    if(mL) assert(sparse_syndrome_non_zero(mL, err->wei, err->vec)); 
#endif     
    if(debug&16){
      printf("# wmax=%d found cw of weight %d: [",wmax,result);
      int max = ((result<50) || (debug&2048)) ? result : 50 ;
      for(int i=0; i< max; i++)
	printf("%d%s",err->vec[i], i+1!=max?" ": (result==max ? "]\n" : "...]\n"));
    }
  }
  else {
    if (p->min_w <= wmax) {
      result = p->min_w;
    } else {
      result = -wmax; /** not found a codeword up to wmax */
    }
  }

  if(smax){
    int skipped = 0;
    if(debug){
      for(int i=1;i<=wmax; i++) {
        if(p_swei[i] <= mH->rows) {
          printf("# w=%d min non-zero syndrome weight %d\n",i,p_swei[i]);
        } else {
          skipped = 1;
        }
      }
    }
    else{
      printf("# confinement: ");
      for(int i=1;i<=wmax; i++) {
        if(p_swei[i] <= mH->rows) {
          printf("%d%s",p_swei[i],i<wmax?",":"");
        } else {
          skipped = 1;
          printf("?%s",i<wmax?",":"");
        }
      }
      printf("\n");
    }
    if (skipped) {
      printf("# Note: Some weights were skipped in confinement profile. Try increasing smax (current: %d)\n", smax);
    }
  }
  
  for(int i=0; i<= wmax; i++)
    free(syn[i]);
  free(syn);
  free(err);
  free(urr);
  csr_free(mHT);

  /** prescribed way to clean the hashing table */
  two_vec_t *terr, *tmp;
  HASH_ITER(hh, errors, terr, tmp) {
    //    two_vec_print(terr);
    HASH_DEL(errors, terr);
    free(terr);
  }
  
  return result;
}


