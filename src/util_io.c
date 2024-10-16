#include <unistd.h>
#include "util_io.h"

params_t prm={
  .debug=3,
  .method=0,
  .classical=0,
  .steps=1,
  .css=1,
  .wmax=5,
  .wmin=1,
  .start=-1, 
  .seed=0,
  .dist=0,
  .dist_max=0,
  .dist_min=0,
  .max_row_wgt_H =0,
  .max_col_wgt_H =0,
  .n0=0,
  .nvar=0,
  .nchk=0,
  .maxC=0,
  .finH=NULL,
  .finG=NULL,
  .finL=NULL,
  .fin="../examples/try", 
  .spaH=NULL,
  .spaG=NULL,
  .spaL=NULL
};

params_t * const p = &prm;

void var_init(int argc, char **argv, params_t * const p){
  int dbg=0;
  int swit=0;

  if(argc <= 1)
    ERROR("no command-line arguments given, " BRIEF_HELP,argv[0]);

  for (int i=1; i<argc;i++) /* scan arguments for help message */
    if((strcmp(argv[i],"--help")==0)||(strcmp(argv[i],"-h")==0)){
      printf( USAGE,argv[0],argv[0]);
      exit (-1);
    }

  for(int i=1; i<argc; i++){
    if(sscanf(argv[i],"debug=%d",& dbg)==1){/** `debug` */
      if(dbg==0)
	p->debug = 0;
      else{
        if(i==1)
          p->debug = dbg; /** just assign if in the `1st position` */
        else
          p->debug ^= dbg; /** otherwise `XOR` */
        if(p->debug&4)
	  printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);
      }
    }
    else if (sscanf(argv[i],"css=%d",&dbg)==1){
      p->css=dbg;
      if (p->debug&4)
	printf("# read %s, css=%d\n",argv[i],p->css);
    }
    else if (0==strncmp(argv[i],"finH=",5)){ /** `finH` */
      if(strlen(argv[i])>5)
        p->finH = argv[i]+5;
      else
        p->finH = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finH=%s; setting fin=\"\"\n",argv[i],p->finH);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"finL=",5)){ /** `finL` */
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finL=%s; setting fin=\"\"\n",argv[i],p->finL);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"finG=",5)){/** `finG` degeneracy generator matrix */
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finG=%s; setting fin=\"\"\n",argv[i],p->finG);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"fin=",4)){
      if(p->finH)
	ERROR("arg[%d]='%s' in conflict with finH=%s\n",i,argv[i],p->finH);
      if(p->finG)
	ERROR("arg[%d]='%s' in conflict with finG=%s\n",i,argv[i],p->finG);
      if(p->finL)
	ERROR("arg[%d]='%s' in conflict with finL=%s\n",i,argv[i],p->finL);
      if (strlen(argv[i])>4)
	p->fin = argv[i]+4;
      else{
	if (i+1 < argc)
	  p->fin = argv[i+1];
	else
	  ERROR("argv[%d]='%s', empty string for 'fin'\n",i,argv[i]);
      }
    }
    else if (sscanf(argv[i],"method=%d",&dbg)==1){
      p->method=dbg;
      if (p->debug&4)
	printf("# read %s, method=%d\n",argv[i],p->method);
      if( (p->method<=0) || (p->method>3))
	ERROR("Unsupported method %d",p->method);
    }
    else if (sscanf(argv[i],"wmax=%d",&dbg)==1){
      p->wmax=dbg;
      if (p->debug&4)
	printf("# read %s, wmax=%d\n",argv[i],p->wmax);
    }
    else if (sscanf(argv[i],"start=%d",&dbg)==1){
      p->start=dbg;
      if (p->debug&4)
	printf("# read %s, start=%d\n",argv[i],p->start);
    }
    else if (sscanf(argv[i],"wmin=%d",&dbg)==1){
      p->wmin=dbg;
      if (p->debug&4)
	printf("# read %s, wmin=%d\n",argv[i],p->wmin);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){
      p->steps=dbg;
      if (p->debug&4)
	printf("# read %s, steps=%d\n",argv[i],p->steps);
    }
    else if (sscanf(argv[i],"seed=%d",&dbg)==1){
      p->seed=dbg;
      if (p->debug&4)
	printf("# read %s, seed=%d\n",argv[i],p->seed);
      if (p->seed<=0){
	p->seed=time(NULL)+1000*p->seed+10*getpid();
	if(p->debug&4)
	  printf("# initializing rng from time(NULL), seed=%d\n",p->seed);
      }
    }    
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for options",argv[0]);
    }
  } /* end parameter scan cycle */

  if(p->method &1 ){ /* RW */
    if (p->steps<=0)
      ERROR("parameter steps=%d should be positive for RW method=%d", p->steps,p->method);
  }
  


  if((strlen(p->fin)!=0) && (!p->finH)){
    int len = strlen(p->fin);
    char *s = (char *) malloc((len+6)*sizeof(char));
    if(!s)
      ERROR("memory allocation");
    sprintf(s,"%s%s",p->fin,swit>0?"X.mtx":"Z.mtx");
    p->finG=s;
    s = (char *) malloc((len+6)*sizeof(char));
    if(!s)
      ERROR("memory allocation");
    sprintf(s,"%s%s",p->fin,swit>0?"Z.mtx":"X.mtx");
    p->finH=s;
    if (p->debug & 2)
      printf("# read 'fin=%s'; " //"since switch=%d "
	     "assigning \n# finH=%s\n# finG=%s\n",
	     p->fin,// swit,
	     p->finH,p->finG);
  }
  
  if (p->finH){
    p->spaH=csr_mm_read(p->finH,p->spaH,0);
    if(p->debug&1)
      printf("# read H <- file '%s'\n",p->finH);
    if(p->debug&32){
      if((p->spaH->cols<150)||(p->debug&2048))
	csr_print(p->spaH,"H");
    }
  }
  else
    ERROR("need to specify H=Hx input file name; use fin=[str] or finH=[str]\n");

  if(p->method &2 ){ /* CC */
    if (p->wmax<=0)
      ERROR("parameter wmax=%d should be positive for CC method=%d", p->wmax,p->method);
    if(p->wmax>=MAX_W)
      ERROR("increase MAX_W=%d defined in 'util_io.h'",MAX_W);
    for(int i=0; i<=p->wmax; i++)
      p->swei[i]=p->spaH->rows +1; 
  }
  
#if 0
  if (p->method&2){ /* cluster */
    p->max_row_wgt_H = csr_max_row_wght(p->spaH);
    if(p->max_row_wgt_H > max_row_wt)
      ERROR("increase max_row_wt=%d to %d",max_row_wt,p->max_row_wgt_H);
  }
#endif


  if((p->finG) && (p->finL))
    ERROR("either G=Hz or L=Lx matrix should be specified but not both! finG='%s' finL='%s'\n",
	  p->finG, p->finL);

  if(p->finG){
    p->classical=0;
    p->spaG=csr_mm_read(p->finG,p->spaG,0);
    if(p->debug&1)
      printf("# read G <- file '%s'\n",p->finG);
    if(csr_csr_mul_non_zero(p->spaH, p->spaG))
       ERROR("rows of H and G matrices are not orthogonal");
    if(p->debug&32){
      if((p->spaG->cols<150)||(p->debug&2048))
	csr_print(p->spaG,"G");
    }
  } 
  else if (p->finL){
    p->classical=0;
    p->spaL=csr_mm_read(p->finL,p->spaL,0);
    if(p->debug&1)
      printf("# read L <- file '%s'\n",p->finL);
    if(p->debug&32){
      if((p->spaL->cols<150)||(p->debug&2048))
	csr_print(p->spaL,"L");
    }
  } 
  else{
    p->classical=1;
    p->spaG=NULL;
  }

  srand(p->seed);

  rci_t n = (p->spaH)-> cols;
  if ((!p->classical) && ((p->spaG) && (n != (p->spaG) -> cols)))
    ERROR("Column count mismatch in H and G matrices: %d != %d",
	  (p->spaH)-> cols, (p->spaG)->cols);
  p->nvar = n; 
  p->n0 = n;
  if (p->css!=1)
    ERROR("Non-CSS codes are currently not supported, css=%d",p->css);
  
  if((p->spaG) && (p->spaL==NULL)){
    /** create `Lx` */
    /** WARNING: this does not necessarily have minimal row weights */
    p->spaL = Lx_for_CSS_code(p->spaH,p->spaG);
    p->nchk = p->spaL->rows;
  }

  if ((p->method <= 0) || (p->method > 3)){
      printf("invalid method=%d specified\n", p->method);
      ERROR(BRIEF_HELP,argv[0]);
  }
  
}

void var_kill(params_t * const p){
  if(p->spaL)
    csr_free(p->spaL);
  if(p->spaH)
    csr_free(p->spaH);
  if(p->spaG)
    csr_free(p->spaG);
#if 0
  if(strlen(p->fin) != 0){
    if(p->finH){
      printf("freeing finH=%s\n", p->finH);
      free(p->finH);
    }
    if(p->finG)
      free(p->finH);    
  }
#endif 
}
