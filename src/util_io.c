#include <unistd.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include "util_io.h"

params_t prm={
  .debug=3,
  .method=0,
  .classical=0,
  .steps=1,
  .css=1,
  .smax=20,
  .wmax=0,
  .dmax=0,
  .wmin=1,
  .noscan=0,
  .fdem=NULL,
  .pmin=0.0,
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
  .fin="", 
  .spaH=NULL,
  .spaG=NULL,
  .spaL=NULL
};

params_t * const p = &prm;

void var_init(int argc, char **argv, params_t * const p){
  int dbg=0;
  int swit=0;
  double prob=0.0;

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
	printf("# read %s, finH=%s; setting finH=\"\"\n",argv[i],p->finH);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"finL=",5)){ /** `finL` */
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finL=%s; setting finL=\"\"\n",argv[i],p->finL);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"finG=",5)){/** `finG` degeneracy generator matrix */
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finG=%s; setting finG=\"\"\n",argv[i],p->finG);
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
    else if (sscanf(argv[i],"smax=%d",&dbg)==1){
      p->smax=dbg;
      if (p->debug&4)
	printf("# read %s, smax=%d\n",argv[i],p->smax);
    }
    else if (sscanf(argv[i],"wmax=%d",&dbg)==1){
      p->wmax=dbg;
      if (p->debug&4)
	printf("# read %s, wmax=%d\n",argv[i],p->wmax);
    }
    else if (sscanf(argv[i],"dmax=%d",&dbg)==1){
      p->dmax=dbg;
      if (p->debug&4)
	printf("# read %s, dmax=%d\n",argv[i],p->dmax);
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
    }    
    else if (sscanf(argv[i],"noscan=%d",&dbg)==1){
      p->noscan=dbg;
      if (p->debug&4)
	printf("# read %s, noscan=%d\n",argv[i],p->noscan);
    }
    else if (0==strncmp(argv[i],"fdem=",5)){
      if(strlen(argv[i])>5)
        p->fdem = argv[i]+5;
      else
        p->fdem = argv[++i];
      if (p->debug&4)
	printf("# read %s, fdem=%s\n",argv[i],p->fdem);
    }
    else if (sscanf(argv[i],"pmin=%lg",&prob)==1){
      p->pmin=prob;
      if (p->debug&4)
	printf("# read %s, pmin=%g\n",argv[i],p->pmin);
    }
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for options",argv[0]);
    }
  } /* end parameter scan cycle */

  if (p->noscan && p->method != 2) {
    ERROR("noscan=1 only works with method=2");
  }

  if (!p->fdem && strlen(p->fin) == 0 && !p->finH && !p->finG && !p->finL) {
    p->fin = "../examples/try";
  }

  if (p->fdem) {
    if (p->finH || p->finG || p->finL || strlen(p->fin) > 0) {
      ERROR("Cannot specify matrix files (fin, finH, finG, finL) along with fdem");
    }
  }
  if (p->pmin != 0.0 && !p->fdem) {
    ERROR("pmin can only be used when fdem is specified");
  }

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
  
  if (p->fdem) {
    read_dem_file(p->fdem, &(p->spaH), &(p->spaL), p->pmin, p->debug);
    p->classical = 0;
    p->nvar = p->spaH->cols;
    p->n0 = p->nvar;
    p->nchk = p->spaL->rows;
  } else {
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
      p->nchk = p->spaL->rows;
    } 
    else{
      p->classical=1;
      p->spaG=NULL;
    }
  }

  if(p->method &2 ){ /* CC */
    if ((p->wmax<=0) && ((p->method & 1 )==0))
      ERROR("parameter wmax=%d should be positive for CC method=%d", p->wmax,p->method);
    if(p->wmax>=MAX_W)
      ERROR("increase MAX_W=%d defined in 'util_io.h'",MAX_W);
    for(int i=0; i<=p->smax; i++)
      p->swei[i]=p->spaH->rows +1; 
  }

  if (p->seed<=0){
    p->seed = time(NULL) - 1000 * p->seed + 10*getpid();
    if(p->debug&4)
      printf("# initializing rng from time(NULL), seed=%d\n",p->seed);
  }
  else if(p->debug&4)
    printf("# initializing rng from seed=%d\n",p->seed);

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

void read_dem_file(char *fnam, csr_t **p_spaH, csr_t **p_spaL, double pmin, int debug){
  ssize_t linelen, col=0;
  size_t lineno=0, bufsiz=0; /**< buffer size for `readline` */
  char *buf = NULL;          /** actual buffer for `readline` */
  int maxH=100, maxL=100; 
  int_pair * inH = malloc(maxH*sizeof(int_pair));
  int_pair * inL = malloc(maxL*sizeof(int_pair));
  if ((!inH)||(!inL))
    ERROR("memory allocation failed\n");

  if(debug & 1)
    printf("# opening DEM file %s\n",fnam);
  FILE *f = fopen(fnam, "r");
  if(f==NULL){
    printf("FILE I/O ERROR: %s\n", strerror(errno));     
    ERROR("can't open the (DEM) file %s for reading\n",fnam);
  }

  int r=-1, k=-1, n=0;
  int iD=0, iL=0; /** numbers of `D` and `L` entries */
  do{ /** read lines one-by-one until end of file is found *************/
    lineno++; col=0; linelen = getline(&buf, &bufsiz, f);
    if(linelen<0)
      break;
    if(debug & 32) printf("# %s",buf);
    char *c=buf;
    double prob;
    int num=0, val;
    while(isspace(*c)){ c++; col++; } /** `skip` white space */
    if((*c != '\0')&& (*c != '#') &&(col < linelen)){
      if(sscanf(c,"error( %lg ) %n",&prob,&num)){
        if((prob<=0)||(prob>=1))
          ERROR("probability should be in (0,1) exclusive p=%g\n"
                "%s:%zu:%zu: '%s'\n", prob,fnam,lineno,col+1,buf);
        c+=num; col+=num;
        
        if (prob < pmin) {
          if (debug & 32) printf("skipping error with prob %g < pmin %g\n", prob, pmin);
          continue; // Skip this error line
        }
        
        do{/** deal with the rest of the line */
          num=0;
          if(sscanf(c," D%d %n",&val, &num)){/** `D` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>=r)
              r=val+1;  /** update the number of `D` pairs */
            if(iD>=maxH){
              maxH=2*maxH;
              inH=realloc(inH,maxH*sizeof(*inH));
            }
            inH[iD].a   = val;   /** add a pair */
            inH[iD++].b = n;
            if(debug & 32) printf("n=%d iD=%d val=%d r=%d\n",n,iD,val, r);
          }
          else if(sscanf(c," L%d %n",&val, &num)){/** `L` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>=k)
              k=val+1;  /** update the number of `L` pairs */
            if(iL>=maxL){
              maxL=2*maxL;
              inL=realloc(inL,maxL*sizeof(*inL));
            }
            inL[iL].a   = val;   /** add a pair */
            inL[iL++].b = n;
            if(debug & 32) printf("n=%d iL=%d val=%d k=%d\n",n,iD,val,k);
          }
	  else if (c[0]=='^') /** entry added by `--decompose_errors` switch in `stim` */
	    c++;            /** just ignore */
          else
            ERROR("unrecognized entry %s"
		  "%s:%zu:%zu: '%s'\n",c,fnam,lineno,col+1,buf);
        }
        while((c[0]!='#')&&(c[0]!='\n')&&(c[0]!='\0')&&(col<linelen));
        n++;
      }
      else if (sscanf(c,"detector( %d %n",&val,&num)){
        /** do nothing */
      }
      else if (sscanf(c,"shift_detectors( %d %n",&val,&num)){
        /** do nothing */
      }
      else
        ERROR("unrecognized DEM entry %s"
              "%s:%zu:%zu: '%s'\n",c,fnam,lineno,col+1,buf);

    }
    /** otherwise just go to next row */
  }
  while(!feof(f));

  if(debug &1)
    printf("# read DEM %s: rows_H=%d rows_L=%d cols=%d; nz_H=%d nz_L=%d\n",fnam,r,k,n,iD,iL);
  if((r<=0)||(k<=0)||(n<=0))
    ERROR("invalid DEM file %s at line=%zu: rows_H=%d rows_L=%d cols=%d; nz_H=%d nz_L=%d\n",
	  fnam,lineno,r,k,n,iD,iL);
  
  *p_spaH = csr_from_pairs(*p_spaH, iD, inH, r, n);
  *p_spaL = csr_from_pairs(*p_spaL, iL, inL, k, n);
  
  if (buf)
    free(buf);
  free(inH);
  free(inL);
  fclose(f);
}
