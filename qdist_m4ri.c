/* 	$Id$	 */

/************************************************************************ 
 * calculate distance of a quantum LDPC code using random window algorithm
 * A. Dumer, A. A. Kovalev, and L. P. Pryadko "Distance verification..."
 * in IEEE Trans. Inf. Th., vol. 63, p. 4675 (2017). 
 * doi: 10.1109/TIT.2017.2690381
 *
 * currently: CSS only 
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>
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
// #include "util.h"

const int maxrow=10;

#define DBG 0

/** 
 * Add row=i of sparse matrix spaQ to row.
 * no range check is done 
 */ 
static inline void addto_inline(mzd_t *row, const csr_t *spaQ, const int i){
#if DBG  
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
 *  P: LDPC parity check (with row weight <= maxrow )
 */

static inline int prep_neis(const int z0, int * nei, const mzd_t * v, const mzd_t * s, const csr_t * P){ 
  int cnt=0, max=P->p[z0+1];
  for (int j=P->p[z0]; j<max; j++){
    if (!mzd_read_bit(v,0,P->i[j])) /* never flip a bit twice */
      nei[cnt++]=P->i[j]; 
  }
#if DBG
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
		     const csr_t * P, const csr_t * PT, const mzd_t * G){
  int res=0, all_zero=1;
  mzd_t *v0=NULL;
  word * rawrow = s->rows[0];  
  rci_t ns=v->ncols;
#if DBG  
  //  if(prm.debug & 512){
    printf("w=%d wgh(Pv)=%d  \nv=",w,syndrome_bit_count(v, P));
    mzd_print(v);
    printf("s=");  mzd_print(s);
    //  }
#endif
    
  for(rci_t i=0 ; i< ns ; i++){
    i=nextelement(rawrow,s->width,i);
    if (i<0)
      break; /* no more non-zero syndrome bits */
    else if(i<ns){ /* found a syndrome to fix */
      if(w>=wmax)
	return -1; /* failed to find a cluster */
      all_zero=0;
      int nei[maxrow];
      int neisize=prep_neis(i,nei,v,s,P);
      //      cout << " row="<< i<< " nei=" << nei<< endl;
      //      if(neisize==0) ERROR("neisize=0");  /* this is actually OK: just continue */
      for(int p=0; p< neisize; p++){
	int j=nei[p];
	mzd_write_bit(v,0,j,1); 
	addto_inline(s,PT,j); 
	res=start_rec(w+1,wmax,v,s,P,PT,G);
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
    //    cout<< "found v="<< v<< endl;
    v0=mzd_copy(v0,v);
    int result=do_reduce(v0,G,G->nrows);
    if (result==-1)  
      return -1; /* go back down, the found row was trivial */
#if DBG    
    if ((mzd_weight(v)!=wmax)||(0!=syndrome_bit_count(v, P)))
      ERROR(" start_rec: something is wrong %d != wmax=%d ",mzd_weight(v),wmax);
#endif     
    return 1; /* success */
  }
  return 0; /* keep going */
}

int do_dist_clus(csr_t *P, mzd_t *G, int debug, int wmax, int start){
  // input: wmax=max cluster weight,
  // start=initial position (-1 to scan all),
  // P:  LDPC parity check 
  // PT:  LDPC parity check *** TRANSPOSED ***
  // G: systematic form of the dual generator
  // these must be orthogonal (no check)

  int nc=P->cols, ns=P->rows;
  csr_t* PT=csr_transpose(NULL,P); // transpose
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
      mzd_set_ui(v,0);   mzd_set_ui(s,0); 
      mzd_write_bit(v,0,i,1);      
      //      printf("v="); mzd_print(v);
      addto_inline(s,PT,i); //  s+=col[i];
      //      printf("s=");  mzd_print(s);
      int done=start_rec(1,w,v,s,P,PT,G);
      if(done==1){
	//	cout << "distance="<< weight(v)<< endl;
	//	cout << "prod="<< P.get_H()*v<< endl;
	return w;
      }
    }
  }
  return -wmax; /* failed up to wmax */
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
  }
  if (debug & 2) 
    printf("# n=%d k=%d weimin=%d\n",n,rt-matP0->nrows,weimin);

  return weimin;
}

#ifdef DEBUG

typedef struct{
  //  int css; /* 1: css, 0: non-css -- currently not supported */
  int debug; /* debug information */ 
  int method; /* bitmap. 1: random window; 2: cluster */
  int steps; /* how many random decoding steps */
  //  char *finP, *finG; /* generators */
  int wmax; /* max cluster size to try */
  int wmin; /* min distance below which we are not interested at all */
  int seed;/* rng seed, set=0 for automatic */
  int dist; /* target distance of the code */
  int dist_max; /* distance actually checked */
  int dist_min; /* distance actually checked */
  int max_row_wgt_G; /* needed for C */
  int start;
    //  int n0;  /* code length, =n for css, (n/2) for non-css */
  //  int n; /* actual n = matrix size*/
} params_t; 


//extern params_t prm;
params_t prm={0,0,0,0,0,0,0,0};


void local_init(int argc, char **argv,int start){
  int i;
  //  memset(&prm,0,sizeof(params_t));
  int dbg=0;
  prm.debug=pio.debug;

  for(i=start; i<argc; i++){
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      if(dbg==0)
	prm.debug=0;
      else{
	prm.debug|=dbg;
	printf("# read %s, debug=%d octal=%o\n",argv[i],prm.debug,prm.debug);
      }
    }					
    else if (sscanf(argv[i],"method=%d",&dbg)==1){
      prm.method=dbg;
      if (prm.debug)
	printf("# read %s, method=%d\n",argv[i],prm.method);
      if( (prm.method<=0) || (prm.method>3))
	ERROR("Unsupported method %d",prm.method);
    }
    else if (sscanf(argv[i],"wmax=%d",&dbg)==1){
      prm.wmax=dbg;
      if (prm.debug)
	printf("# read %s, wmax=%d\n",argv[i],prm.wmax);
    }
    else if (sscanf(argv[i],"start=%d",&dbg)==1){
      prm.start=dbg;
      if (prm.debug)
	printf("# read %s, start=%d\n",argv[i],prm.start);
    }
    else if (sscanf(argv[i],"wmin=%d",&dbg)==1){
      prm.wmin=dbg;
      if (prm.debug)
	printf("# read %s, wmin=%d\n",argv[i],prm.wmin);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){
      prm.steps=dbg;
      if (prm.debug)
	printf("# read %s, steps=%d\n",argv[i],prm.steps);
    }
    else if (sscanf(argv[i],"seed=%d",&dbg)==1){
      prm.seed=dbg;
      if (prm.debug)
	printf("# read %s, seed=%d\n",argv[i],prm.seed);
    }    
    else if((strcmp(argv[i],"--help")==0)||(strcmp(argv[i],"-h")==0)){
      printf( "%s: calculate the minumum distance of a q-LDPC code\n"
	      "\tusage: %s [input parameters] [algorithm params]\n"
	      USAGE_IO
	      "Supported algorithm parameters:\n"	
	      "\tseed=0: rng seed  [0 for time(NULL)]\n"
	      "\tmethod=1: bitmap for method used. \n"
	      "\t\t1: random window (RW)\n"
	      "\t\t2: cluster algorithm (C)\n"
	      "\tsteps=1: how many RW decoding cycles\n"
	      "\twmax=5: max cluster weight (C)\n"
	      "\twmin=0: min distance of interest (RW)\n"
	      "\t-h or --help gives this help\n"
	      "\tdefault: wmax=5 seed=0\n",argv[0],argv[0]);
      exit (-1);
    }
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for options",argv[0]);
    }      
  } /* end parameter scan cycle */

  if (prm.seed==0){
    if(prm.debug)
      printf("# initializing rng from time(NULL)\n");
    srand(time(NULL));
  }
  else {
    srand(prm.seed);
    if(prm.debug)
      printf("# setting srand(%d)\n",prm.seed);
  }
  if(prm.method &1 ){ /* RW */
    if (prm.steps<=0)
      prm.steps=1; /* need to run at least once */
    if (prm.debug)
      printf("# using RW method, wmin=%d steps=%d\n",prm.wmin,prm.steps);
  }
  if(prm.method &2 ){ /* RW */
    if (prm.wmax<=0)
      prm.wmax=10; /* some reasonable lower bound */
    if (prm.debug)
      printf("# using Cluster method, wmax=%d steps=%d\n",prm.wmax,prm.steps);
  }
}

int main(int argc, char **argv)
{
  csr_t *spaP=NULL;  
  csr_t *spaG=NULL;
  int istart=local_io_init(argc,argv,&spaP,&spaG);
  local_init(argc,argv,istart); /* initialized variables */

  if (prm.method & 2)
  { /* cluster */
    prm.max_row_wgt_G=csr_max_row_wght(spaG);
    if(prm.max_row_wgt_G>maxrow)
      ERROR("main: increase maxrow=%d to %d",maxrow,prm.max_row_wgt_G);

  }
  int n=pio.n;

  
  /* convert G to standard form */
  mzp_t *piv0=mzp_init(n);  //  mzp_out(piv0);
  mzd_t *matG0=mzd_from_csr(NULL,spaG); 
  rci_t rankG=mzd_gauss_naive(matG0,piv0,1); 
  if(prm.debug & 1)
    printf("# rankG=%d\n",rankG);
  mzd_apply_p_right_trans(matG0,piv0);
  matG0->nrows=rankG;

  mzp_t *q0=perm_p_trans(NULL,piv0,1);    // permutation equiv to piv0 
  csr_t* spaP0=csr_apply_perm(NULL,spaP,q0); // permuted sparse P
  //  csr_t* spaG0=csr_apply_perm(NULL,spaG,q0); // permuted sparse G -- not needed here
  if(prm.debug & 2048)
  {
    if ((n<=150))
    {
      printf("matG0=\n");
      mzd_print(matG0);

      printf("matP0=\n");
      mzd_t *matP0=mzd_from_csr(NULL,spaP0);  
      mzd_print(matP0);
      mzd_free(matP0);
    }
    int wei=product_weight_csr_mzd(spaP0, matG0,1);    
    printf("weigt of P0*G0_T=%d\n",wei);
    if(wei>0)
      ERROR("expected zero weight product!");
  }  

  if (prm.method & 1)
  { /* RW method */
    prm.dist_max=do_dist_rnd(spaP0,matG0,prm.debug,prm.steps,prm.wmin);
    if (prm.debug) 
      printf("### RW upper bound on the distance: %d\n",prm.dist_max);
    prm.wmax=minint(prm.wmax,prm.dist_max-1);
  }
  if (prm.method & 2)
  { /* cluster method */    
    int dmin=do_dist_clus(spaP0,matG0,prm.debug,prm.wmax,prm.start);
    if (prm.debug)
      printf("### Cluster: dmin=%d\n",dmin);
    if (dmin>0)
    { 
      prm.dist_min = dmin; /* actual distance found */
      prm.dist_max = dmin;
    }
    else if (dmin<0)
    {
      if (-dmin==prm.dist_max-1)
	prm.dist_min=prm.dist_max; /* OK */
      else 
	prm.dist_min=-dmin;
    }
    //    if(prm.debug){
      if (prm.dist_min==prm.dist_max)
	printf("success d=%d\n",prm.dist_min);
      else if (prm.dist_max>prm.dist_min)
	printf("distance in the interval %d to %d\n", prm.dist_min,prm.dist_max);
      else
	printf("cluster algorithm failed up to wmax=%d\n",-dmin);
      //    }
  }
  else
  { /* just RW */
    //    if(prm.debug){
      printf("RW algorithm upper bound for the distance d=%d\n", prm.dist_max);      
      //    }
  }
  return 0;
}

#endif /* DEBUG */
