
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
 *
 * @param dW weight increment from the minimum found
 * @param p pointer to global parameters structure
 * @param classical set to `1` for classical code (do not use `L` matrix), `0` otherwise  
 * @return minimum `weight` of a CW found (or `-weight` if early termination condition is reached),
 *         or `0` if no codewords with `w<wmax` have been found.
 */
int do_RW_dist(params_t * const p){
  const csr_t * const spaH0 = p->spaH;
  const csr_t * const spaL0 = p->spaL;
  const int steps = p->steps;
  const int wmin = p->wmin;
  const int wmax = p->wmax;
  const int classical = p->classical;
  const int debug = p->debug;
  /** whether to verify logical ops as a vector or individually */
  const int nvar = spaH0->cols;
  if(((!classical)&&(spaL0==NULL)) ||
     ((classical)&&(spaL0!=NULL))){
    printf("L0 %s NULL classical=%d\n",spaL0==NULL ? "=" : "!=", classical);	   
    ERROR("L0 should be non-NULL only for classical code!\n");
  }

  int minW = nvar + 1;
  if (p->dmax > 0 && p->dmax < minW) {
    minW = p->dmax;
  }
  if (p->min_w != INT_MAX && p->min_w < minW) {
    minW = p->min_w;
  }

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
  word * piv_mask=calloc(mH->width, sizeof(word)); /** bitmask of pivot columns */
  if((!pivs) || (!perm) || (!piv_mask))
    ERROR("memory allocation failed!\n");
  int eff_nrows = spaH0->rows;

  mzd_t *N_ker = NULL;
  mzd_t *M_sub = NULL;
  int nu = 0;
  int ksub_eff = 0;
  if (p->ksub > 0) {
    if (2 * spaH0->rows < nvar) {
      if (debug & 1) {
        fprintf(stderr,
                "# Warning: ksub=%d requested, but m=%d < n-m=%d (<= nu); "
                "falling back to full-matrix RW (ksub=0)\n",
                p->ksub, spaH0->rows, nvar - spaH0->rows);
      }
      p->ksub = 0;
    } else {
      N_ker = mzd_nullspace(spaH0);
      nu = N_ker ? N_ker->nrows : 0;
      if (spaH0->rows < nu || nu <= 0) {
        if (debug & 1) {
          fprintf(stderr,
                  "# Warning: ksub=%d requested, but m=%d < nu=%d; "
                  "falling back to full-matrix RW (ksub=0)\n",
                  p->ksub, spaH0->rows, nu);
        }
        if (N_ker) {
          mzd_free(N_ker);
          N_ker = NULL;
        }
        nu = 0;
        p->ksub = 0;
      } else {
        ksub_eff = minint(p->ksub, nu);
        if (ksub_eff > 0) {
          M_sub = mzd_init(ksub_eff, nvar);
        }
      }
    }
  }

  csr_t *mHT_csr = NULL;
  int *visited_cols = NULL;
  int *visited_checks = NULL;
  int *col_queue = NULL;
  int visit_marker = 0;
  uint64_t rng_state = (uint64_t)p->seed + 0x517cc1b727220a95ULL;
  if (p->kwin > 0 && p->kwin < nvar) {
    mHT_csr = csr_transpose(NULL, spaH0);
    visited_cols = calloc(nvar, sizeof(int));
    visited_checks = calloc(spaH0->rows, sizeof(int));
    col_queue = calloc(nvar, sizeof(int));
  }

  for (int ii=0; ii< steps; ii++){
    if (p->kwin > 0 && p->kwin < nvar) {
      int seed_col = rand_uniform(nvar);
      localized_window_perm(perm, nvar, seed_col, p->kwin, p->win_mode,
                            spaH0, mHT_csr, visited_cols, visited_checks,
                            col_queue, &visit_marker, &rng_state);
    } else {
      pivs=mzp_rand(pivs); /** random pivots LAPAC-style */
      mzp_set_ui(perm,1);
      perm=perm_p_trans(perm,pivs,0); /**< corresponding permutation */
    }

    if (p->ksub > 0) {
      if (ksub_eff <= 0) break;
      for (int i = 0; i < ksub_eff; i++) {
        int r = rand_uniform(nu);
        mzd_copy_row(M_sub, i, N_ker, r);
      }
      int rank = 0;
      for (int i = 0; i < nvar && rank < ksub_eff; i++) {
        int col = perm->values[i];
        if (gauss_one(M_sub, col, rank)) {
          rank++;
        }
      }
      for (int ir = 0; ir < rank; ir++) {
        int cnt = 0;
        int limit = nvar + 1;
        int cur_d = (minW <= nvar) ? minW : 0;
        if (cur_d > 0) {
          if ((p->outC || p->maxC || p->dW > 0 || p->min_hits > 0) && p->dW >= 0) {
            limit = minint(limit, cur_d + p->dW + 1);
          } else {
            limit = minint(limit, cur_d);
          }
        }
        word *rawrow = mzd_row(M_sub, ir);
        rci_t j = -1;
        while (cnt < limit) {
          j = nextelement(rawrow, M_sub->width, j);
          if (j == -1 || j >= nvar) break;
          ee[cnt++] = j++;
        }
        if (cnt > 0 && cnt < limit) {
          int nz = classical ? 1 : sparse_syndrome_non_zero(spaL0, cnt, ee);
          if (nz) {
            p->codewords = codeword_add_maybe(p, ee, cnt);
            if (cnt < minW) minW = cnt;
            if (p->maxC && p->num_cws >= p->maxC) goto alldone;
            if (check_min_hits_convergence(p)) goto alldone;
            if (minW <= wmin) {
              minW = -minW;
              goto alldone;
            }
          }
        }
      }
      if (p->refresh > 0 && (ii + 1) % p->refresh == 0) {
        int cw_wt = 0;
        if (p->codewords && p->min_w < nvar) {
          cw_vec_t *cw, *tmp;
          HASH_ITER(hh, p->codewords, cw, tmp) {
            if (cw->weight == p->min_w) {
              cw_wt = cw->weight;
              for (int k = 0; k < cw_wt; k++) ee[k] = cw->arr[k];
              break;
            }
          }
        }
        refresh_nullspace_basis(N_ker, cw_wt > 0 ? ee : NULL, cw_wt, &rng_state);
      }
      continue;
    }

    /** full row echelon form of `H` (gauss) using the order in `perm` */
    memset(piv_mask, 0, mH->width * sizeof(word));
    int rank = 0;
    for (int i = 0; i < nvar && rank < eff_nrows; i++) {
      int col = perm->values[i];
      if (gauss_one_rows(mH, col, rank, eff_nrows)) {
        pivs->values[rank++] = col;
        piv_mask[col >> 6] |= (word)1 << (col & 63);
      }
    }
    eff_nrows = rank;

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
    for (int col = 0, ir = 0; col < nvar; col++){ /** each row in the dual matrix */
      if ((piv_mask[col >> 6] >> (col & 63)) & 1) continue;
      ir++;
      int cnt=0; /** how many non-zero elements */
      ee[cnt++] = col;
      int limit = nvar + 1;
      int cur_d = (minW <= nvar) ? minW : 0;
      if (cur_d > 0) {
        if ((p->outC || p->maxC || p->dW > 0 || p->min_hits > 0) && p->dW >= 0) {
          limit = minint(limit, cur_d + p->dW + 1);
        } else {
          limit = minint(limit, cur_d);
        }
      }
#if (NEW==0) /** older version going over columns of `H` */
      for(int ix=0; ix<rank; ix++){
        if(mzd_read_bit(mH,ix,col))
          ee[cnt++] = pivs->values[ix];
	if (cnt >= limit) /** `cw` of no interest */
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
	  if (cnt >= limit) /** `cw` of no interest */
	    break;
	}
	else
	  break;
      }
#else /** NEW==1, use transposed `H` -- the `fastest` version of the code*/
      word * rawrow = mzd_row(mHT,col);  
      rci_t j=-1;
      const int active_width = (rank + 63) >> 6;
      while(cnt < limit){/** `cw` of no interest */
	j=nextelement(rawrow,active_width,j);
	if(j==-1 || j >= rank) // empty line after simplification
	  break; 
	ee[cnt++] = pivs->values[j++];
      }
#endif /* NEW */              
      if (cnt < limit){
	/** sort the column indices */
	rci_quick_sort(ee, cnt);
#ifndef NDEBUG
	/** expensive: verify orthogonality */
	if(sparse_syndrome_non_zero(spaH0, cnt, ee)){
	  printf("# cw of weight %d: [",cnt);
	  for(int i=0; i<cnt;i++)
	    printf("%d%s",ee[i],i+1==cnt?" ":"]\n");
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
          p->codewords = codeword_add_maybe(p, ee, cnt);
	  if(debug&16){
	    printf("# step=%d row=%d minW=%d found cw of W=%d: [",ii,ir,minW,cnt);
	    const int max = ((cnt<25) || (debug&2048)) ?  cnt : 25 ;
	    for(int i=0; i< max; i++)
	      printf("%d%s", ee[i], i+1!=max?" ": (cnt==max ? "]\n" : "...]\n"));
	  }
          if (cnt < minW) {
            minW = cnt;
          }
          if (p->maxC && p->num_cws >= p->maxC) {
            goto alldone;
          }
          if (check_min_hits_convergence(p)) {
            goto alldone;
          }
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
  if (M_sub) mzd_free(M_sub);
  if (N_ker) mzd_free(N_ker);
  if (mHT_csr) csr_free(mHT_csr);
  free(visited_cols);
  free(visited_checks);
  free(col_queue);
  free(piv_mask);
  mzp_free(perm);
  mzp_free(pivs);
  free(ee);
  if(mHT)
    mzd_free(mHT);
  mzd_free(mH);
  
  if (minW < 0) {
    return minW;
  }
  if (p->min_w == INT_MAX) {
    return 0;
  }
  if (wmax > 0 && p->min_w > wmax) {
    return 0;
  }
  return p->min_w;
}


#ifdef STANDALONE

int do_CC_dist(params_t * const p);


int main(int argc, char **argv){
  params_t * const p = &prm;

  var_init(argc,argv,p);

  if (p->finC) {
    nzlist_read(p->finC, p);
  }

  //  const int n=p->nvar;

  if (prm.method & 1){ /* RW method */
    
    prm.dist_max=do_RW_dist(p);

    if (prm.debug&1){
      if (prm.dist_max != 0)
        printf("### RW upper bound on the distance: %d\n",prm.dist_max);
      else
        printf("### RW: no upper bound found up to wmax = %d\n", prm.wmax);
      if(prm.dist_max <0)
        printf("### negative distance due to wmin=%d set (early termination)\n",prm.wmin);
      else if (prm.dist_max ==0)
        printf("### no codewords of weight <= %d found\n",prm.wmax);
    }
    if (prm.dist_max < 0) {
      if(prm.debug) {
        if (prm.method == 1) 
          printf("RW algorithm upper bound for the distance d=%d\n", prm.dist_max);
        else
          printf("RW algorithm upper bound for the distance d=%d "
                 "(early termination due to wmin=%d)\n", prm.dist_max, prm.wmin);
      }
      printf("%d\n",prm.dist_max);
      prm.dist_min = prm.dist_max;
      goto end_all;
    }
    if (prm.method==1){ /** just RW */
      if(prm.debug) {
        if (prm.dist_max > 0)
          printf("RW algorithm upper bound for the distance d=%d\n", prm.dist_max);
        else
          printf("RW algorithm: no upper bound found up to wmax = %d\n", prm.wmax);
      }
      printf("%d\n",prm.dist_max);
    }
    else{
      if (prm.wmax==0) {
        if (prm.dist_max == 0) {
          ERROR("RW found no codewords, cannot run CC with wmax=0. Please specify wmax > 0.\n");
        }
        prm.wmax=abs(prm.dist_max)-1;
      }
      else if (prm.dist_max != 0) {
        prm.wmax=minint(prm.wmax, abs(prm.dist_max)-1);   
      }
      if (prm.wmax == 0) {
        prm.dist_min = 1;
        prm.dist_max = 1;
        if (prm.debug & 1) {
          printf("success (distance is 1) d=1\n");
        }
        printf("1\n");
        goto end_all;
      }
    }
  }
  
  if (prm.method & 2){ /* cluster method */
    int dmin=do_CC_dist(p);

    if (dmin>0){ 
      if (prm.debug&1)
	printf("### Cluster (actual min-weight codeword found): d=%d\n",dmin);
      printf("%d\n",dmin);
      prm.dist_min = dmin; /* actual distance found */
      prm.dist_max = dmin;
      goto end_all;
    }
    else if (dmin<0){
      if (prm.debug&1)
	printf("### Cluster dmin=%d  (no codewords of weight up to %d)\n",dmin,-dmin);
      if (-dmin==abs(prm.dist_max)-1){
	prm.dist_min=abs(prm.dist_max); /* OK */
        if (prm.debug&1)
          printf("success (two distance bounds coincide) d=%d\n",prm.dist_min);
        printf("%d\n",prm.dist_min);       
        goto end_all;
      }
      else{
	prm.dist_min=-dmin;
        if(prm.debug){
          if (prm.dist_max>prm.dist_min)
            printf("# distance in the interval (inclusive) %d to %d\n", prm.dist_min,prm.dist_max);
          else
            printf("# cluster algorithm failed to find a codeword up to wmax=%d\n",-dmin);
        }
        printf("%d\n",dmin);
      }
    }
    else
      ERROR("unexpected dmin=0\n");
  }
 end_all:
    if (p->outC) {
      char comment[256];
      sprintf(comment, "generated by dist_m4ri");
      nzlist_write(p->outC, comment, p);
    }
    if (p->debug & 32) {
      cw_vec_t *cw;
      for(cw = p->codewords; cw != NULL; cw = (cw_vec_t *)(cw->hh.next)){
        printf("# cw: [ ");
        for(int i=0; i<cw->weight; i++) printf("%d ", 1 + cw->arr[i]);
        printf("] cnt=%d\n", cw->cnt);
      }
    }
    var_kill(p);
    
    return 0;
  }

#endif /* STANDALONE */
