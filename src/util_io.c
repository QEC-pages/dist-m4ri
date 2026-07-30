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
  .smax=5,
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
  .dW=0,
  .finC=NULL,
  .outC=NULL,
  .codewords=NULL,
  .num_cws=0,
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
  long long int dbg_ll=0;

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
    else if (0==strncmp(argv[i],"finC=",5)){
      if(strlen(argv[i])>5)
        p->finC = argv[i]+5;
      else
        p->finC = argv[++i];
      if (p->debug&4)
	printf("# read %s, finC=%s\n",argv[i],p->finC);
    }
    else if (0==strncmp(argv[i],"outC=",5)){
      if(strlen(argv[i])>5)
        p->outC = argv[i]+5;
      else
        p->outC = argv[++i];
      if (p->debug&4)
	printf("# read %s, outC=%s\n",argv[i],p->outC);
    }
    else if (sscanf(argv[i],"maxC=%lld",&dbg_ll)==1){
      p->maxC=dbg_ll;
      if (p->debug&4)
	printf("# read %s, maxC=%lld\n",argv[i],p->maxC);
    }
    else if (sscanf(argv[i],"dW=%d",&dbg)==1){
      p->dW=dbg;
      if (p->debug&4)
	printf("# read %s, dW=%d\n",argv[i],p->dW);
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
    for(int i=0; i<MAX_W; i++)
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

  cw_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    HASH_DEL(p->codewords, cw);
    free(cw);
  }
}

typedef struct {
    char **lines;
    int size;
    int capacity;
} dem_program_t;

static dem_program_t *read_dem_to_program(const char *fnam) {
  FILE *f = fopen(fnam, "r");
  if (f == NULL) {
    printf("FILE I/O ERROR: %s\n", strerror(errno));
    ERROR("can't open the (DEM) file %s for reading\n", fnam);
  }

  dem_program_t *prog = malloc(sizeof(dem_program_t));
  prog->size = 0;
  prog->capacity = 100;
  prog->lines = malloc(prog->capacity * sizeof(char *));

  char *buf = NULL;
  size_t bufsiz = 0;
  ssize_t linelen;

  while ((linelen = getline(&buf, &bufsiz, f)) >= 0) {
    if (linelen > 0 && buf[linelen - 1] == '\n') {
      buf[linelen - 1] = '\0';
    }
    if (prog->size >= prog->capacity) {
      prog->capacity *= 2;
      prog->lines = realloc(prog->lines, prog->capacity * sizeof(char *));
    }
    prog->lines[prog->size++] = strdup(buf);
  }
  if (buf) free(buf);
  fclose(f);
  return prog;
}

static void free_dem_program(dem_program_t *prog) {
  for (int i = 0; i < prog->size; i++) {
    free(prog->lines[i]);
  }
  free(prog->lines);
  free(prog);
}

static void parse_instructions(dem_program_t *prog, int *p_line_idx, 
                        int *p_iD, int_pair **p_inH, int *p_maxH, int *p_r,
                        int *p_iL, int_pair **p_inL, int *p_maxL, int *p_k,
                        int *p_n, double pmin, int *p_detector_shift, int debug) {
    while (*p_line_idx < prog->size) {
        char *line = prog->lines[*p_line_idx];
        (*p_line_idx)++;
        
        char *c = line;
        while (isspace(*c)) c++;
        
        if (*c == '\0' || *c == '#') continue;
        
        if (*c == '}') {
            return;
        }
        
        int num = 0;
        int val = 0;
        double prob = 0.0;
        
        // Parse repeat
        if (sscanf(c, "repeat %d { %n", &val, &num) == 1) {
            int start_idx = *p_line_idx;
            int temp_idx = start_idx;
            for (int r = 0; r < val; r++) {
                temp_idx = start_idx;
                parse_instructions(prog, &temp_idx, 
                                   p_iD, p_inH, p_maxH, p_r,
                                   p_iL, p_inL, p_maxL, p_k,
                                   p_n, pmin, p_detector_shift, debug);
            }
            if (val > 0) {
                *p_line_idx = temp_idx;
            } else {
                int depth = 1;
                while (*p_line_idx < prog->size && depth > 0) {
                    char *s = prog->lines[*p_line_idx];
                    (*p_line_idx)++;
                    while (isspace(*s)) s++;
                    if (strncmp(s, "repeat", 6) == 0 && strchr(s, '{')) depth++;
                    if (*s == '}') depth--;
                }
            }
            continue;
        }
        
        // Parse shift_detectors
        int shift_val = 0;
        if (sscanf(c, "shift_detectors ( %*[^)] ) %d %n", &shift_val, &num) == 1) {
            *p_detector_shift += shift_val;
            continue;
        } else if (sscanf(c, "shift_detectors %d %n", &shift_val, &num) == 1) {
            *p_detector_shift += shift_val;
            continue;
        }
        
        // Parse error
        if (sscanf(c, "error( %lg ) %n", &prob, &num) == 1) {
            if ((prob <= 0) || (prob >= 1))
                ERROR("probability should be in (0,1) exclusive p=%g\n"
                      "line %d: '%s'\n", prob, *p_line_idx, line);
            c += num;
            
            if (prob < pmin) {
                continue;
            }
            
            while (1) {
                while (isspace(c[0])) c++;
                if (c[0] == '\0' || c[0] == '#' || c[0] == '\n') break;
                
                num = 0;
                if (sscanf(c, "D%d%n", &val, &num) == 1) {
                    c += num;
                    assert(val >= 0);
                    int shifted_val = val + *p_detector_shift;
                    if (shifted_val >= *p_r)
                        *p_r = shifted_val + 1;
                    if (*p_iD >= *p_maxH) {
                        *p_maxH = 2 * (*p_maxH);
                        *p_inH = realloc(*p_inH, (*p_maxH) * sizeof(**p_inH));
                    }
                    (*p_inH)[*p_iD].a = shifted_val;
                    (*p_inH)[*p_iD].b = *p_n;
                    (*p_iD)++;
                } else if (sscanf(c, "L%d%n", &val, &num) == 1) {
                    c += num;
                    assert(val >= 0);
                    if (val >= *p_k)
                        *p_k = val + 1;
                    if (*p_iL >= *p_maxL) {
                        *p_maxL = 2 * (*p_maxL);
                        *p_inL = realloc(*p_inL, (*p_maxL) * sizeof(**p_inL));
                    }
                    (*p_inL)[*p_iL].a = val;
                    (*p_inL)[*p_iL].b = *p_n;
                    (*p_iL)++;
                } else if (c[0] == '^') {
                    c++;
                } else {
                    ERROR("unrecognized entry %s in error line %d: '%s'\n", c, *p_line_idx, line);
                }
            }
            (*p_n)++;
            continue;
        }
        
        if (strncmp(c, "detector", 8) == 0) {
            continue;
        }
        
        if (strncmp(c, "logical_observable", 18) == 0) {
            continue;
        }
        
        ERROR("unrecognized DEM entry in line %d: '%s'\n", *p_line_idx, line);
    }
}

void read_dem_file(char *fnam, csr_t **p_spaH, csr_t **p_spaL, double pmin, int debug){
  dem_program_t *prog = read_dem_to_program(fnam);
  
  int maxH=100, maxL=100; 
  int_pair * inH = malloc(maxH*sizeof(int_pair));
  int_pair * inL = malloc(maxL*sizeof(int_pair));
  if ((!inH)||(!inL))
    ERROR("memory allocation failed\n");

  int r=-1, k=-1, n=0;
  int iD=0, iL=0;
  int detector_shift = 0;
  int line_idx = 0;

  parse_instructions(prog, &line_idx, 
                     &iD, &inH, &maxH, &r,
                     &iL, &inL, &maxL, &k,
                     &n, pmin, &detector_shift, debug);

  if (line_idx < prog->size) {
      ERROR("Unmatched '}' in DEM file %s at line %d\n", fnam, line_idx);
  }

  if(debug & 1)
    printf("# read DEM %s: rows_H=%d rows_L=%d cols=%d; nz_H=%d nz_L=%d\n",fnam,r,k,n,iD,iL);
  if((r<=0)||(k<=0)||(n<=0))
    ERROR("invalid DEM file %s: rows_H=%d rows_L=%d cols=%d; nz_H=%d nz_L=%d\n",
	  fnam,r,k,n,iD,iL);
  
  *p_spaH = csr_from_pairs(*p_spaH, iD, inH, r, n);
  *p_spaL = csr_from_pairs(*p_spaL, iL, inL, k, n);
  
  free(inH);
  free(inL);
  free_dem_program(prog);
}

FILE * nzlist_w_new(const char fnam[], const char comment[]){
  FILE *f=fopen(fnam,"w");
  if(!f){
    printf("FILE I/O ERROR: %s\n", strerror(errno));
    ERROR("can't open file %s for writing",fnam);
  }
  fprintf(f,"%%%% NZLIST\n");
  if(comment)
    fprintf(f,"%% %s\n",comment);
  return f;
}

int nzlist_w_append(FILE *f, const cw_vec_t * const vec){
  assert(vec && vec->weight >0 );
  assert(f!=NULL);
  const int w=vec->weight;
  if(fprintf(f,"%d ",w)<=0)
    ERROR("can't write to `NZLIST` file");
  for(int i=0; i < w; i++)
    if(fprintf(f," %d%s", 1 + vec->arr[i], i+1 < w ? "" :"\n")<=0)
      ERROR("can't write to `NZLIST` file");
  return 0;
}

FILE * nzlist_r_open(const char fnam[], long long int *lineno){
  FILE *f=fopen(fnam,"r");
  if(!f)
    return(NULL);
  *lineno=1;
  char c=fgetc(f);
  while(c=='%'){
    do{
      c=fgetc(f);
      if(feof(f))
	return NULL;
    }
    while(c!='\n');
    (*lineno)++;
    c=fgetc(f);
  }
  ungetc(c,f); 
  return f;
}

cw_vec_t * nzlist_r_one(FILE *f, cw_vec_t * vec, const char fnam[], long long int *lineno){
  assert(f!=NULL);
  if ( ferror (f)|| feof(f) )
    return NULL;
  int w;

  char c=fgetc(f);
  while(c=='%'){
    do{
      c=fgetc(f);     
      if(feof(f))
	return NULL;
    }
    while(c!='\n');
    (*lineno)++;
    c=fgetc(f);
  }
  ungetc(c,f); 
  
  if(fscanf(f," %d",&w) != 1){
    if (feof(f)) return NULL;
    printf("%s:%lld: invalid NZLIST entry\n", fnam, *lineno);
    ERROR("expected an integer");
  }
  if ((vec!=NULL) && (vec->weight<w)){
    free(vec);
    vec=NULL;
  }
  if(vec==NULL){
    vec = calloc(sizeof(cw_vec_t)+w*sizeof(int), sizeof(char));
    if(!vec)
      ERROR("memory allocation");
  }
  vec->weight = w;
  vec->cnt = 1;
  for(int i=0; i<w; i++){
    if(fscanf(f," %d ",vec->arr + i) != 1){
      printf("%s:%lld: invalid entry of weight w=%d\n",fnam, *lineno, w);
      ERROR("expected an integer i=%d of %d",i,w);
    }    
    vec->arr[i]--;
  }

  for(int i=1; i<w; i++){
    if((vec->arr[i-1] < 0) || (vec->arr[i-1] >= vec->arr[i])){
      printf("%s:%lld: invalid entry of weight w=%d\n",fnam, *lineno, w);
      ERROR("expected strictly increasing positive entries");
    }   
  }
  (*lineno)++;
  return vec;
}

cw_vec_t * codeword_add_maybe(cw_vec_t *codewords, const int arr[], int weight, long long int *p_num_cws, long long int maxC) {
  if (maxC && *p_num_cws >= maxC) {
    return codewords;
  }
  const size_t keylen = weight * sizeof(int);
  cw_vec_t *pvec = NULL;
  HASH_FIND(hh, codewords, arr, keylen, pvec);
  if (!pvec) {
    cw_vec_t *entry = malloc(sizeof(cw_vec_t) + keylen);
    if (!entry) ERROR("memory allocation");
    entry->weight = weight;
    entry->cnt = 1;
    for (int i = 0; i < weight; i++) {
      entry->arr[i] = arr[i];
    }
    HASH_ADD(hh, codewords, arr, keylen, entry);
    (*p_num_cws)++;
  } else {
    pvec->cnt++;
  }
  return codewords;
}

long long int nzlist_read(const char fnam[], params_t *p){
  long long int count = 0, lineno;
  long long int skipped_invalid = 0;
  assert(fnam);
  FILE * f=nzlist_r_open(fnam, &lineno);
  if(!f){
    if ((p->outC ==NULL) || (strcmp(fnam,p->outC)!=0)){      
      printf("codeword input file I/O ERROR: %s, outC=%s\n", strerror(errno),p->outC);
      ERROR("can't open file %s for reading",fnam);
    }
    else
      return 0;
  }
  cw_vec_t *entry=NULL;
  while((entry=nzlist_r_one(f,NULL, fnam, &lineno))){
    if((p->maxC) && (p->num_cws >= p->maxC)) {
      free(entry);
      break;
    }
    int valid = 1;
    if (p->spaH) {
      if (sparse_syndrome_non_zero(p->spaH, entry->weight, entry->arr)) {
        valid = 0;
      }
    }
    if (valid && p->spaL) {
      if (!sparse_syndrome_non_zero(p->spaL, entry->weight, entry->arr)) {
        valid = 0;
      }
    }
    if (!valid) {
      skipped_invalid++;
      free(entry);
      continue;
    }
    if((p->wmax==0) ||((p->wmax) && (entry->weight <= p->wmax))){
      long long int old_num = p->num_cws;
      p->codewords = codeword_add_maybe(p->codewords, entry->arr, entry->weight, &(p->num_cws), p->maxC);
      if (p->num_cws > old_num) {
        count++;
      }
    }
    free(entry);
  }
  fclose(f);
  if (skipped_invalid > 0) {
    printf("# Warning: skipped %lld invalid codewords (not orthogonal to H or orthogonal to L)\n", skipped_invalid);
  }
  if(p->debug&1)
    printf("# read %lld codewords from %s, total %lld\n",count, fnam, p->num_cws);
  return count; 
}

long long int nzlist_write(const char fnam[], const char comment[], params_t *p){
  long long int count=0;
  assert(fnam);
  FILE * f = nzlist_w_new(fnam, comment);
  cw_vec_t *pvec;
  
  for(pvec = p->codewords; pvec != NULL; pvec = (cw_vec_t *)(pvec->hh.next)){
    if((p->wmax==0) ||((p->wmax) && (pvec->weight <= p->wmax))){
      count ++;
      nzlist_w_append(f,pvec);
    }
  }
  fclose(f);
  return count;
}
