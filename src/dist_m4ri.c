
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


/** @brief Random Information Set search for small-E logical operators.
4 *
 * @param dW weight increment from the minimum found
 * @param p pointer to global parameters structure
 * @param classical set to `1` for classical code (do not use `L` matrix), `0` otherwise  
 * @return minimum `weight` of a CW found (or `-weigt` if early termination condition is reached). Or `0` if no codewords wit `w<wmax` have been found.
 */
int do_RW_dist(const csr_t * const spaH0, const csr_t * const spaL0,
	       const int steps, const int wmin, const int wmax,
	       const int classical, const int debug){
  /** whether to verify logical ops as a vector or individually */
  const int nvar = spaH0->cols;
  if(((!classical)&&(spaL0==NULL)) ||
     ((classical)&&(spaL0!=NULL))){
    printf("L0 %s NULL classical=%d\n",spaL0==NULL ? "=" : "!=", classical);	   
    ERROR("L0 should be non-NULL only for classical code!\n");
  }

  int minW = wmax > 0 ? wmax : nvar+1;

  if(debug&2)
    printf("# running do_RW_dist() with steps=%d wmin=%d wmax=%d classical=%d nvar=%d\n",
	   steps, wmin, wmax, classical, nvar);
  
  mzd_t * mH = mzd_from_csr(NULL, spaH0);
  mzd_t *mHT = NULL;
  /** actual `vector` in sparse form */
  rci_t *ee = malloc(nvar*sizeof(rci_t)); 
  
  if((!mH) || (!ee))
    ERROR("memory allocation failed!\n");

  /** 1. Construct random column permutation P */
  mzp_t * perm=mzp_init(nvar); /** identity column permutation */
  mzp_t * pivs=mzp_init(nvar); /** list of pivot columns */
  mzp_t * pivs_srtd=mzp_init(nvar); /** list of pivot columns */
  mzp_t * skip_pivs=mzp_init(nvar); /** list of pivot columns */
  if((!pivs) || (!perm))
    ERROR("memory allocation failed!\n");

  for (int ii=0; ii< steps; ii++){
    pivs=mzp_rand(pivs); /** random pivots LAPAC-style */
    mzp_set_ui(perm,1);
    perm=perm_p_trans(perm,pivs,0); /**< corresponding permutation */

    /** full row echelon form of `H` (gauss) using the order in `perm` */
    int rank=0;
    for(int i=0; i< nvar; i++){
      int col=perm->values[i];
      int ret=gauss_one(mH, col, rank);
      if(ret)
        pivs->values[rank++]=col;
    }

    /** construct skip-pivot permutation */
    pivs_srtd = mzp_copy(pivs_srtd,pivs);
    qsort(pivs_srtd->values, rank, sizeof(pivs->values[0]), cmp_rci_t);
    int end=-1, num=0;
    for(int i=0; i < rank; i++){
      int beg = end + 1;
      end = pivs_srtd->values[i];
      for(int j = beg; j < end; j++)
	skip_pivs->values[num++] = j;
    }
    for(int j = end + 1 ; j < nvar; j++)
      skip_pivs->values[num++] = j;
    
#ifndef NDEBUG
    if (num + rank != nvar)
      ERROR("mismatch: rank=%d and num=%d do not add to nvar=%d\n",rank,num,nvar);
#endif
    skip_pivs->length = num;

#ifndef NEW
# define NEW 1
#endif 
#if (NEW==1)
    /** it is a bit faster to transpose `mH` first. */
    mHT = mzd_transpose(mHT,mH);
#endif     
    /** calculate sparse version of each vector (list of positions)
     *  `p`    `p``p`               # pivot columns marked with `p`   
     *  [1  a1        b1 ] ->  [a1  1  a2 a3 0 ]
     *  [   a2  1     b2 ]     [b1  0  b2 b3 1 ]
     *  [   a3     1  b3 ]
     */
    int k = nvar - rank;
    for (int ir=0; ir< k; ir++){ /** each row in the dual matrix */
      int cnt=0; /** how many non-zero elements */
      const int col = ee[cnt++] = skip_pivs->values[ir];
#if (NEW==0) /** older version going over columns of `H` */
      for(int ix=0; ix<rank; ix++){
        if(mzd_read_bit(mH,ix,col))
          ee[cnt++] = pivs->values[ix];
	if (cnt >= minW) /** `cw` of no interest */
	  break;
      }
#elif (NEW==2) /** 
		   function `mzd_find_pivot()` walks over columns
		   one-by-one to find a non-zero bit.  Returns `1` if
		   a non-zero bit was found.
		   WARNING: this is the slowest option!!!
	       */
      rci_t ic=col, ix=0;
      while (ix < rank){
	int res = mzd_find_pivot(mH, ix, col, &ix, &ic);
	if((res)&&(ic==col)){
	  ee[cnt++] = pivs->values[ix++];
	  //	  printf("cnt=%d j=%d\n",cnt,ix); 
	  if (cnt >= minW) /** `cw` of no interest */
	    break;
	}
	else
	  break;
      }
#else /** NEW==1, use transposed `H` -- the `fastest` version of the code*/
      word * rawrow = mHT->rows[col];  
      rci_t j=-1;
      while(cnt < minW){/** `cw` of no interest */
	j=nextelement(rawrow,mHT->width,j);
	if(j==-1) // empty line after simplification
	  break; 
	ee[cnt++] = pivs->values[j++];
      }
#endif /* NEW */              
      if (cnt < minW){
	/** sort the column indices */
	qsort(ee, cnt, sizeof(rci_t), cmp_rci_t);
#ifndef NDEBUG
	/** expensive: verify orthogonality */
	if(sparse_syndrome_non_zero(spaH0, cnt, ee)){
	  printf("# cw of weight %d: [",cnt);
	  for(int i=0; i<cnt;i++)
	    printf("%d%s",1+ee[i],i+1==cnt?" ":"]\n");
	  ERROR("this should not happen: cw not orthogonal to H");
	}
#endif /* NDEBUG */
      
	/** verify logical operator */
	int nz;
	if (classical)
	  nz=1; /** no need to verify */
	else
	  nz = sparse_syndrome_non_zero(spaL0, cnt, ee);	
	if(nz){ /** we got non-trivial codeword! */
	  /** TODO: try local search to `lerr` (if 2 or larger) */
	  /** at this point we have `cnt` codeword indices in `ee` */
	  if(debug&16){
	    printf("# step=%d row=%d minW=%d found cw of W=%d: [",ii,ir,minW,cnt);
	    const int max = ((cnt<25) || (debug&2048)) ?  cnt : 25 ;
	    for(int i=0; i< max; i++)
	      printf("%d%s", 1+ee[i], i+1!=max?" ": (cnt==max ? "]\n" : "...]\n"));
	  }
	  minW=cnt;
	  if (minW <= wmin){ /** early termination condition */
	    minW = - minW;   /** this distance value is of little interest; */
	    goto alldone; /** stop right away */
	  }
	}
      }      
    } /** end of the dual matrix rows loop */
    if(debug&8){
      if(ii%1000==999)
	printf("# round=%d of %d minW=%d\n", ii+1, steps, minW);
    }
    
  }/** end of `steps` random window */

 alldone: /** early termination label */

  /** clean up */
  if(skip_pivs)
    mzp_free(skip_pivs);
  if(pivs_srtd)
    mzp_free(pivs_srtd);
  mzp_free(perm);
  mzp_free(pivs);
  free(ee);
  if(mHT)
    mzd_free(mHT);
  mzd_free(mH);
  
  return ((wmax>0) && (minW==wmax)) ? 0 :minW ;
}


#ifdef STANDALONE

int do_CC_dist(const csr_t * const mH, const csr_t * mL,
	       const int wmax, const int start, int p_swei[],
	       const int smax, const int debug);


int main(int argc, char **argv){
  params_t * const p = &prm;

  var_init(argc,argv,p);

  //  const int n=p->nvar;

  if (prm.method & 1){ /* RW method */
    
    prm.dist_max=do_RW_dist(p->spaH,p->spaL,p->steps, p->wmin, p->wmax, p->classical, p->debug);

    if (prm.debug&1){
      printf("### RW upper bound on the distance: %d\n",prm.dist_max);
    if(prm.dist_max <0)
      printf("### negative distance due to wmin=%d set (early termination)\n",prm.wmin);
    else if (prm.dist_max ==0)
            printf("### no vectors below wmax=%d found\n",prm.wmax);
    }
    prm.wmax=minint(prm.wmax, abs(prm.dist_max)-1);
  }

  if (prm.method & 2){ /* cluster method */
    int dmin=do_CC_dist(p->spaH,p->spaL,p->wmax,p->start,p->swei,p->smax, p->debug);

    if (dmin>0){ 
      if (prm.debug&1)
	printf("### Cluster (actual min-weight codeword found): dmin=%d\n",dmin);
      prm.dist_min = dmin; /* actual distance found */
      prm.dist_max = dmin;
    }
    else if (dmin<0){
      if (prm.debug&1)
	printf("### Cluster dmin=%d  (no codewords of weight up to %d)\n",dmin,-dmin);
      if (-dmin==abs(prm.dist_max)-1)
	prm.dist_min=abs(prm.dist_max); /* OK */
      else 
	prm.dist_min=-dmin;
    }
    else
      ERROR("unexpected dmin=0\n");

    if (prm.dist_min==abs(prm.dist_max)){
      if(p->method==3)
	printf("success  (two distance bounds coincide) d=%d\n",prm.dist_min);
      else
	printf("success  (found min-weight codeword) d=%d\n",prm.dist_min);
    }      
    else if (prm.dist_max>prm.dist_min)
      printf("distance in the interval (inclusive) %d to %d\n", prm.dist_min,prm.dist_max);
    else
      printf("cluster algorithm failed to find a codeword up to wmax=%d\n",-dmin);
  }
  else{ /* just RW */
    printf("RW algorithm upper bound for the distance d=%d\n", prm.dist_max);
  }
  var_kill(p);
  
  return 0;
}

#endif /* STANDALONE */
