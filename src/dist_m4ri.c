
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

/** 
 * @brief Add row=i of sparse matrix spaQ to row.
 * WARNING: no range check is done 
 */ 
static inline void addto_inline(mzd_t *row, const csr_t *spaQ, const int i){
#ifndef NDEBUG  
  if (i>=spaQ->rows)
    ERROR("addto: attempt to get row=%d of %d",i,spaQ->rows);
  if (row->ncols != spaQ->cols)
    ERROR("addto: column number mismatch");
  //  mzd_print(row);
#endif
  for (int j=spaQ->p[i]; j<spaQ->p[i+1]; j++)
    mzd_flip_bit(row, 0,spaQ->i[j]); /* flip that bit */
}

/** 
 * given current value bits and syndrome bits, see which bits in the
 * vicinity of the check z0 can be flipped to make it happy.
 * Return number of entries to process (length of nei) 
 *  input: z0 is the row of P to check 
 *  nei: neighborhood vector to construct 
 *  v: value bits 
 *  s: syndrome bits 
 *  P: LDPC parity check (with row weight <= max_row_wt )
 */

static inline int prep_neis(const int z0, int * nei, const mzd_t * v, _maybe_unused const mzd_t * s, const csr_t * P){ 
  int cnt=0, max=P->p[z0+1];
  for (int j=P->p[z0]; j<max; j++){
    if (!mzd_read_bit(v,0,P->i[j])) /* never flip a bit twice */
      nei[cnt++]=P->i[j]; 
  }
#ifndef NDEBUG
  //  if(prm.debug&256){
    printf("nei=[");
    for(int i=0; i< cnt; i++)
      printf(" %d ",nei[i]);
    printf("]\n");
    //  }
#endif 
  return cnt;
}

/** recursive function. 
 *  return: -1: fail; 0: success; 1: immediate termination 
 *  input: 
 * 	w=current recursion level
 * 	wmax = max recursion level 
 *      v,s: current value and syndrome bit vectors 
 * 	P, PT, G: matrices as in do_dist_clus
 */
static int start_rec(const int w, const int wmax, mzd_t * v, mzd_t * s,
		     const csr_t * const P, const csr_t * const PT, const mzd_t * const G, const int rankG){
  int res=0, all_zero=1;
  mzd_t *v0=NULL;
  word * rawrow = s->rows[0];  
  rci_t ns=v->ncols;
#ifndef NDEBUG  
  //  if(prm.debug & 512){
    printf("w=%d wgh(Pv)=%d  \nv=",w,syndrome_bit_count(v, P));
    mzd_print(v);
    printf("s=");  mzd_print(s);
    //  }
#endif
    
  for(rci_t i=-1 ; i+1< ns ; ){
    i=nextelement(rawrow,s->width,i+1);
    if (i<0)
      break; /* no more non-zero syndrome bits */
    else if(i<ns){ /* found a syndrome to fix */
      if(w>=wmax)
	return -1; /* failed to find a cluster */
      all_zero=0;
      int nei[max_row_wt];
      int neisize=prep_neis(i,nei,v,s,P);
      //      cout << " row="<< i<< " nei=" << nei<< endl;
      //      if(neisize==0) ERROR("neisize=0");  /* this is actually OK: just continue */
      for(int p=0; p< neisize; p++){
	int j=nei[p];
	mzd_write_bit(v,0,j,1); 
	addto_inline(s,PT,j); 
	res=start_rec(w+1,wmax,v,s,P,PT,G, rankG);
	if(res==1)
	  return 1; /* just get out fast */
	addto_inline(s,PT,j); 
	mzd_write_bit(v,0,j,0);  /* clean up */
	//	printf("s=");	mzd_print(s);
      }      
    }
    else break;
  }
  //  printf("done with w=%d all_zero=%d\n",w,all_zero);
  if(all_zero){ // syndrome vector was zero
    if(G){ /** quantum code */
    //    cout<< "found v="<< v<< endl;
    v0=mzd_copy(v0,v);
    int result=do_reduce(v0,G,rankG);
    if (result==-1)  
      return -1; /* go back down, the found row was trivial */
#ifndef NDEBUG    
    if (((int)mzd_weight(v)!=wmax)||(0!=syndrome_bit_count(v, P)))
      ERROR(" start_rec: something is wrong %zu != wmax=%d ",mzd_weight(v),wmax);
#endif
    }
    return 1; /* success */
  }
  return 0; /* keep going */
}

/** @brief lower bound on the minimum distance by cluster enumeration
 * 
 * WARNING: only intended for LDPC codes
 */
int do_dist_clus(const csr_t * const P, const mzd_t * const G, int debug, int wmax, int start, const int rankG){
  // input: wmax=max cluster weight,
  // start=initial position (-1 to scan all),
  // P:  LDPC parity check 
  // PT:  LDPC parity check *** TRANSPOSED ***
  // G: systematic form of the dual generator
  // these must be orthogonal (no check)

  int nc=P->cols, ns=P->rows;
  csr_t* PT=csr_transpose(NULL,P); 
  //  csr_out(spaG0T);

  mzd_t *v=mzd_init(1,nc), *s=mzd_init(1,ns);
  for(int w=2;w<=wmax;w++){ // cluster weight loop 
    if (debug&8)
      printf( "# starting w=%d\n", w);
    int beg=0, end=nc-w-1; 
    if (start>=0)
      beg=end=start; 
    for(int i=beg;i<=end;i++){ // starting bit loop 
      if ((debug&8)&&(w==wmax))
	printf( "# w=%d start=%d\n", w,i);
      mzd_set_ui(v,0);
      mzd_set_ui(s,0); 
      mzd_write_bit(v,0,i,1);      
      //      printf("v="); mzd_print(v);
      addto_inline(s,PT,i); //  s+=col[i];
      //      printf("s=");  mzd_print(s);
      int done=start_rec(1,w,v,s,P,PT,G,rankG);
      if(done==1){
	//	cout << "distance="<< weight(v)<< endl;
	//	cout << "prod="<< P.get_H()*v<< endl;
	csr_free(PT);
	mzd_free(v);
	return w;
      }
    }
  }
  csr_free(PT);  
  mzd_free(v);
  return -wmax; /* failed up to wmax */
}


/** @brief prepare an ordered pivot-skip list of length `n-rank` */
mzp_t * do_skip_pivs(const size_t rank, const mzp_t * const pivs){
  const rci_t n=pivs->length;
  rci_t j1=rank; /** position to insert the next result */
  mzp_t * ans = mzp_copy(NULL,pivs);
  qsort(ans->values, rank, sizeof(pivs->values[0]), cmp_rci_t);

  for(rci_t j=0; j<n; j++){
    if(!bsearch(&j, ans->values, rank, sizeof(ans->values[0]), cmp_rci_t)){
      ans->values[j1++]=j;
    }
  }
  assert(j1==n);

  int j=rank;
  for(size_t i=0; j<n; j++)
    ans->values[i++] = ans->values[j];
  ans->length = n-rank;

  if(prm.debug&8){/** in skip_pivs */
    printf("skip_pivs of len=%d: ",ans->length);
    for(int i=0; i< ans->length; i++)
      printf(" %d%s",ans->values[i],i+1 == ans->length ?"\n":"");
    printf("pivs of len=%d, rank=%zu: ",pivs->length, rank);
    for(size_t i=0; i< rank; i++)
      printf(" %d%s",pivs->values[i],i+1 == rank ?"\n":"");
  }
  return ans;
}


/** @brief Random Information Set search for small-E logical operators.
 *
 * @param dW weight increment from the minimum found
 * @param p pointer to global parameters structure
 * @param classical set to `1` for classical code (do not use `L` matrix), `0` otherwise  
 * @return minimum `weight` of a CW found (or `-weigt` if early termination condition is reached). 
 */
int do_RW_dist(const csr_t * const spaH0, const csr_t * const spaL0,
	       const int steps, const int wmin, const int classical, const int debug){
  /** whether to verify logical ops as a vector or individually */
  const int nvar = spaH0->cols;
  if(((!classical)&&(spaL0==NULL)) ||
     ((classical)&&(spaL0!=NULL))){
    printf("L0 %s NULL classical=%d\n",spaL0==NULL ? "=" : "!=", classical);	   
    ERROR("L0 should be non-NULL present only for classical code!\n");
  }

  int minW=nvar+1;

  if(debug&16)
    printf("running do_RW_dist() with steps=%d wmin=%d classical=%d nvar=%d\n",steps, wmin, classical, nvar);
  
  mzd_t * mH = mzd_from_csr(NULL, spaH0);
  //  mzd_t *mLt = NULL, *eemLt = NULL; //, *mL = NULL;
  rci_t *ee = malloc(nvar*sizeof(rci_t)); /** actual `vector` */
  
  if((!mH) || (!ee))
    ERROR("memory allocation failed!\n");
  //  if(p->debug & 16)  mzd_print(mH);
  /** 1. Construct random column permutation P */

  mzp_t * perm=mzp_init(nvar); /** identity column permutation */
  mzp_t * pivs=mzp_init(nvar); /** list of pivot columns */
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
    mzp_t * skip_pivs = do_skip_pivs(rank, pivs);

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
      for(int ix=0; ix<rank; ix++){
        if(mzd_read_bit(mH,ix,col))
          ee[cnt++] = pivs->values[ix];
	if (cnt >= minW)
	  break;
      }
      if(cnt < minW){
	/** sort the column indices */
	qsort(ee, cnt, sizeof(rci_t), cmp_rci_t);
      
	/** verify logical operator */
	int nz;
	if (classical)
	  nz=1; /** no need to verify */
	else
	  nz = sparse_syndrome_non_zero(spaL0, cnt, ee);	
	if(nz){ /** we got non-trivial codeword! */
	  /** TODO: try local search to `lerr` (if 2 or larger) */
	  /** calculate the energy and compare */
	  /** at this point we have `cnt` codeword indices in `ee`, and its `energ` */
	  //        if (cnt < minW){
	  minW=cnt;
	  if (minW <= wmin){ /** early termination condition */
	    minW = - minW; /** this distance value is of little interest; */
	  }
	}
      }
    } /** end of the dual matrix rows loop */
    if(debug & 16)
      printf(" round=%d of %d minW=%d\n", ii+1, steps, minW);
    
    mzp_free(skip_pivs);
    
  }/** end of `steps` random window */

  //alldone: /** early termination label */

  /** clean up */
  mzp_free(perm);
  mzp_free(pivs);
  free(ee);

  mzd_free(mH);
  
  return minW;
}


int do_dist_rnd(csr_t *spaG0, mzd_t *matP0, int debug,int steps, int wmin){
  rci_t n = spaG0-> cols;

  mzd_t *matG0perp=NULL;
  int weimin=n, rt=0;  
  csr_t *spaG1=NULL;
  mzp_t *piv1=mzp_init(n), *q1=mzp_init(n);

  for (int ii=0; ii < steps ; ii++){ /* random window decoding loop */
    mzp_rand(piv1); /* random permutation in terms of pivots */
    mzp_set_ui(q1,1); q1=perm_p_trans(q1,piv1,0);    /* corresponding permutation */
    spaG1=csr_apply_perm(spaG1,spaG0,q1);     /* G with permuted cols */
    matG0perp=mzd_generator_from_csr(matG0perp,spaG1);
    mzd_apply_p_right(matG0perp, piv1); // permute cols back to order in G0
    //    mzd_info(matG0perp,0);
    // generator matrix dual to G0
    if((debug & 1)&&(ii==0))
      printf("# rankG=%d\n",matG0perp->nrows);
    if (((debug &512)||(debug & 2048)) ){
      if (debug & 2048) printf("current G\n");
      mzd_t *matG0=mzd_from_csr(NULL,spaG0);
      if ((debug & 2048) && (n<=150)){
	mzd_print(matG0); 
	printf("current Gperp=\n");
	mzd_print(matG0perp);
	printf(" G*Gperp_T=\n");
      }
      mzd_t *prod1;
      //      mzd_info(prod1,0);
      //      mzd_info(matG0,0);
      //      mzd_info(matG0perp,0);
      mzd_t * matG0ptran = mzd_transpose(NULL,matG0perp);
      prod1=csr_mzd_mul(NULL,spaG0,matG0ptran,1);
      mzd_free(matG0ptran);
      if ((debug & 2048) && (n<=150))
	mzd_print(prod1);
      printf("weigt of G0*G0perp_T=%d\n",(int)mzd_weight(prod1));
      if ((debug & 2048) && (n<=150)){
	printf("current P=\n");
	mzd_print(matP0);
	printf(" G*P_T=\n");
      }
      mzd_free(prod1);
      prod1=mzd_init(matG0->nrows,matP0->nrows);
      prod1=_mzd_mul_naive(prod1,matG0,matP0,1);
      if ((debug & 2048) && (n<=150))
	mzd_print(prod1);
      printf("weigt of G*P^T=%d\n",(int)mzd_weight(prod1));
      mzd_free(matG0);
      mzd_free(prod1);
    }

    rt = matG0perp->nrows; 
    if (rt==matP0->nrows)
      ERROR("This is an empty code, distance infinite!");
    rci_t wei;
    mzd_t *row, *row0=NULL;
    int width = matG0perp -> width;
    if((debug & 1)&&(ii==0))
      printf("# n=%d width=%d\n",n,width);
    //    printf("##### here ########\n");
    for(int i=0;i<rt;i++){
      fflush(stdout);
      row=mzd_init_window(matG0perp,i,0,i+1,n);
      wei=mzd_weight(row); 
      if((wei<weimin)){ 
	//	if (debug & 1024)
	  row0=mzd_copy(row0,row);
	  //	else
	  //	  row0=row;
	rci_t j=do_reduce(row,matP0,matP0->nrows);
	if (j!=-1){ /* not an empty row */
	  if (debug & 1024){
	    int wei0=mzd_weight_naive(row0);
	    if (wei!=wei0){
	      printf("### naive wei=%d vs std_bitcount %d\n",wei0,wei);
	      mzd_print(row0);
	    }
	  }
	  if (debug & 3)
	    printf("# round=%d i=%d w=%d s=%d " "s0=%d "
		   "min=%d\n",ii,i,wei,
		   debug &2 ? syndrome_bit_count(row0,spaG0):0,syndrome_bit_count(row,spaG0),
		   weimin);	    	  
	  weimin=wei;
	}
	/*	else
	  printf("# wei=%d, syndrome should be zero: s0=%d \n" , wei,
		   syndrome_bit_count(row0,spaG0));	    	  
	*/
      }
      mzd_free(row);
      if(weimin<=wmin){ // no need to continue 
	ii=steps;
	weimin=-weimin;
	break;
      }

    }
    if (row0!=NULL){
      mzd_free(row0);
      row0=NULL;
    }
    if (weimin<0)
      break;
  }
  if (debug & 2) 
    printf("# n=%d k=%d weimin=%d\n",n,rt-matP0->nrows,weimin);

  return weimin;
}

#ifdef STANDALONE

int main(int argc, char **argv){
  params_t * const p = &prm;
  var_init(argc,argv,p);

  const int n=p->nvar;

  

  if (prm.method & 1){ /* RW method */
    if (p->debug & 1)
      printf("# starting RW method with wmin=%d steps=%d classical=%d\n",p->wmin,p->steps,p->classical);
#if 0 /** older version, may have bugs */    
    prm.dist_max=do_dist_rnd(spaH0,matG0,prm.debug,prm.steps,prm.wmin);
#else //! `new` distance-finding routine 
    prm.dist_max=do_RW_dist(p->spaH,p->spaL,p->steps, p->wmin, p->classical, p->debug);
#endif /* 0 */    
    if (prm.debug & 1){
      printf("### RW upper bound on the distance: %d\n",prm.dist_max);
    if(prm.dist_max <0)
      printf("### negative distance due to wmax=%d set (early termination)\n",prm.wmax);
    }
    prm.wmax=minint(prm.wmax, abs(prm.dist_max)-1);
  }

  if (prm.method & 2){ /* cluster method */
    /* convert G to standard form */
    mzp_t *piv0=mzp_init(n);  //  mzp_out(piv0);
    mzd_t *matG0=mzd_from_csr(NULL,p->spaG); 
    rci_t rankG=mzd_gauss_naive(matG0,piv0,1); 
    if(prm.debug & 1)
      printf("# rankG=%d\n",rankG);
    mzd_apply_p_right_trans(matG0,piv0);
    //    matG0->nrows=rankG;
      
    mzp_t *q0=perm_p_trans(NULL,piv0,1);    // permutation equiv to piv0 
    csr_t *spaH0=csr_apply_perm(NULL,p->spaH,q0); // permuted sparse H
    //  csr_t* spaG0=csr_apply_perm(NULL,p->spaG,q0); // permuted sparse G -- not needed here
    if(prm.debug & 2048){
      if ((n<=150))
	{
	  printf("matG0=\n");
	  mzd_print(matG0);
	    
	  printf("matP0=\n");
	  mzd_t *matP0=mzd_from_csr(NULL,spaH0);  
	  mzd_print(matP0);
	  mzd_free(matP0);
	}
      int wei=product_weight_csr_mzd(spaH0, matG0,1);    
      printf("weigt of H0*G0_T=%d\n",wei);
      if(wei>0)
	ERROR("expected zero weight product!");
    }
    mzp_free(piv0);
    mzp_free(q0);
    
    int dmin=do_dist_clus(spaH0,matG0,prm.debug,prm.wmax,prm.start,rankG);
    csr_free(spaH0);
    mzd_free(matG0);
    if (dmin>0){ 
      if (prm.debug & 1)
	printf("### Cluster (actual min-weight codeword found): dmin=%d\n",dmin);
      prm.dist_min = dmin; /* actual distance found */
      prm.dist_max = dmin;
    }
    else if (dmin<0){
      if (prm.debug & 1)
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
    //    if(prm.debug &1)
    printf("RW algorithm upper bound for the distance d=%d\n", prm.dist_max);
      
      //    }
  }
  var_kill(p);
  
  return 0;
}

#endif /* STANDALONE */
