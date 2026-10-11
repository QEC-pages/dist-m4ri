#include <unistd.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "util_io.h"

params_t prm={
  .debug=3,
  .method=3,
  .classical=-1,
  .steps=100000,
  .css=1,
  .smax=0,
  .wmax=0,
  .dmin=1, /* the trivial lower bound (as documented in --help); smaller values are treated as 1 */
  .dmax=0,
  .wmin=1,
  .noscan=0,
  .fdem=NULL,
  .pmin=0.0,
  .start_list=NULL,
  .start_num=0,
  .cbeg=-1,
  .cend=-1,
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
  .min_w=INT_MAX,
  .cw_max_w=0,
  .cw_cnt_w=NULL,
  .cw_hits_w=NULL,
  .finH=NULL,
  .finG=NULL,
  .finL=NULL,
  .fin="", 
  .spaH=NULL,
  .spaG=NULL,
  .spaL=NULL,
  .maskL=NULL,
  .threads=0,
  .dexp=0,
  .dstop=0,
  .timeout=60.0,
  .nothrottle=0,
  .chunk_size=0,
  .ksub=0,
  .kwin=0,
  .win_mode=0,
  .min_hits=5,
  .cov_cws=100,
  .refresh=0
};

params_t * const p = &prm;

void print_short_help(const char *prog) {
  fprintf(stderr, SHORT_HELP, prog, DIST_M4RI_VERSION, prog);
}

/** @brief Parse `start=a,b,c` (comma-separated list of CC start columns) into `p->start_list`.
 *
 * A single negative value (e.g., the legacy default `start=-1`) clears the list, i.e., all
 * columns are used.  A repeated `start=` argument replaces the previous list.  Range checks,
 * sorting, and removal of duplicates are done in `var_init()` once the matrix size is known.
 */
static void parse_start_list(const char * const str, params_t * const p){
  if (p->start_list) {
    free(p->start_list);
    p->start_list = NULL;
  }
  p->start_num = 0;
  int cnt = 1;
  for (const char *c = str; *c; c++)
    if (*c == ',')
      cnt++;
  int *list = malloc(cnt * sizeof(int));
  if (!list)
    ERROR("memory allocation");
  const char *s = str;
  for (int k = 0; k < cnt; k++) {
    char *end = NULL;
    errno = 0;
    long val = strtol(s, &end, 10);
    if ((end == s) || errno || (val > INT_MAX) || (val < INT_MIN) || ((*end != ',') && (*end != '\0')))
      ERROR("invalid start='%s': expected a comma-separated list of column indices, e.g., start=0,48", str);
    list[k] = (int) val;
    s = end + 1;
  }
  if ((cnt == 1) && (list[0] < 0)) { /* legacy default `start=-1`: all columns */
    free(list);
    return;
  }
  for (int k = 0; k < cnt; k++)
    if (list[k] < 0)
      ERROR("invalid start='%s': column indices must be non-negative", str);
  p->start_list = list;
  p->start_num = cnt;
}

static int cmp_int(const void *a, const void *b){
  const int x = *(const int *) a, y = *(const int *) b;
  return (x > y) - (x < y);
}

/** @brief Parse the value `str` of the argument `arg` = `debug=...`: a non-negative decimal or hexadecimal (0x...)
 *  integer (no octal: a leading zero is ignored) */
static int parse_debug_value(const char * const str, const char * const arg){
  const int hex = (str[0] == '0') && ((str[1] == 'x') || (str[1] == 'X'));
  const char * const digits = hex ? str + 2 : str;
  char *end = NULL;
  errno = 0;
  const long val = strtol(digits, &end, hex ? 16 : 10);
  if ((end == digits) || (*end != '\0') || errno || (val < 0) || (val > INT_MAX))
    ERROR("invalid '%s': expected a non-negative decimal or hexadecimal (0x...) integer", arg);
  return (int) val;
}

/** @brief The file name given in the argument `argv[*i]` = `name=file` (`len` is the length of `name=`), or in the
 *  space form `name= file`; then the file name is the next argument, and `*i` is advanced */
static char * file_arg_value(const int argc, char **argv, int * const i, const size_t len){
  if (strlen(argv[*i]) > len)
    return argv[*i] + len;
  if (*i + 1 >= argc)
    ERROR("argv[%d]='%s': missing file name", *i, argv[*i]);
  return argv[++(*i)];
}

/** @brief Number of nonzero entries of a CSR matrix (compressed or pair form) */
static int csr_nnz(const csr_t * const M){
  return (M->nz == -1) ? M->p[M->rows] : M->nz;
}

/** @brief Maximum column weight of a CSR matrix in compressed form (0 otherwise) */
static int csr_max_col_wght(const csr_t * const M){
  if ((M->nz != -1) || (M->cols <= 0)) return 0;
  int *cnt = calloc(M->cols, sizeof(int));
  if (!cnt)
    ERROR("memory allocation");
  int wmax = 0;
  for (int j = 0; j < M->p[M->rows]; j++)
    if (++cnt[M->i[j]] > wmax)
      wmax = cnt[M->i[j]];
  free(cnt);
  return wmax;
}

void csr_dump(FILE *stream, const csr_t * const M, const char name[], const int debug){
  if (!(debug & DBG_MATRICES) || !M || !stream) return;
  fprintf(stream, "# matrix %s: %d x %d, %d nonzeros", name, M->rows, M->cols, csr_nnz(M));
  if ((M->cols >= 150) && !(debug & DBG_LARGE)) {
    fprintf(stream, " (not shown for n >= 150 without debug bit 2048)\n");
    return;
  }
  fprintf(stream, " ('1': nonzero, '.': zero)\n");
  char *row = malloc(M->cols + 1);
  if (!row)
    ERROR("memory allocation");
  mzd_t *D = mzd_from_csr(NULL, M);
  for (int i = 0; i < M->rows; i++) {
    for (int j = 0; j < M->cols; j++)
      row[j] = mzd_read_bit(D, i, j) ? '1' : '.';
    row[M->cols] = '\0';
    fprintf(stream, "# %s\n", row);
  }
  mzd_free(D);
  free(row);
}

/** @brief One line of the input summary: size, number of nonzeros, and maximum row and column weights of a matrix */
static void print_matrix_summary(FILE *stream, const char name[], const csr_t * const M, const char note[]){
  fprintf(stream, "#   %s: %d x %d%s, %d nonzeros, max row weight %d, max column weight %d\n", name, M->rows, M->cols,
          note, csr_nnz(M), (M->nz == -1) ? csr_max_row_wght(M) : 0, csr_max_col_wght(M));
}

/** @brief DBG_SUMMARY: the code type, its length n, the input files, and the input matrices */
static void print_input_summary(FILE *stream, const params_t * const p){
  const char * const type = p->classical ? "classical code" : (p->fdem ? "detector error model" : "quantum CSS code");
  fprintf(stream, "# input: %s, n=%d (", type, p->nvar);
  if (p->fdem) {
    fprintf(stream, "fdem='%s'", p->fdem);
    if (p->pmin > 0.0)
      fprintf(stream, ", pmin=%g", p->pmin);
  } else {
    fprintf(stream, "finH='%s'", p->finH);
    if (p->finG)
      fprintf(stream, ", finG='%s'", p->finG);
    if (p->finL)
      fprintf(stream, ", finL='%s'", p->finL);
  }
  fprintf(stream, ")\n");
  print_matrix_summary(stream, "H", p->spaH, p->fdem ? " (detectors)" : "");
  if (p->spaG)
    print_matrix_summary(stream, "G", p->spaG, "");
  if (p->spaL)
    print_matrix_summary(stream, "L", p->spaL, p->fdem ? " (observables)" : (p->finL ? "" : " (from H and G)"));
}

/** @brief Rank of a binary CSR matrix (dense elimination) */
static int csr_rank(const csr_t * const M){
  mzd_t *D = mzd_from_csr(NULL, M);
  const int rank = mzd_echelonize(D, 0);
  mzd_free(D);
  return rank;
}

/** @brief DBG_PARAMS: the ranks of the input matrices and the number k of encoded (qu)bits (dense elimination) */
static void print_code_params(FILE *stream, const params_t * const p){
  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  const int n = p->nvar;
  const int rank_H = csr_rank(p->spaH);
  int k;
  char buf[128];
  if (p->classical || !p->spaL) {
    k = n - rank_H;
    snprintf(buf, sizeof(buf), "k=n-rank(H)=%d", k);
  } else if (p->spaG) { /* CSS code with H=Hx and G=Hz */
    const int rank_G = csr_rank(p->spaG);
    k = n - rank_H - rank_G;
    snprintf(buf, sizeof(buf), "rank(G)=%d, k=n-rank(H)-rank(G)=%d", rank_G, k);
  } else { /* the logical operators L independent of the rows of H */
    mzd_t *mH = mzd_from_csr(NULL, p->spaH);
    mzd_t *mL = mzd_from_csr(NULL, p->spaL);
    mzd_t *mHL = mzd_stack(NULL, mH, mL);
    const int rank_L = mzd_echelonize(mL, 0);
    k = mzd_echelonize(mHL, 0) - rank_H;
    mzd_free(mHL);
    mzd_free(mL);
    mzd_free(mH);
    snprintf(buf, sizeof(buf), "rank(L)=%d, k=rank([H;L])-rank(H)=%d", rank_L, k);
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);
  fprintf(stream, "# code parameters: n=%d, rank(H)=%d, dim ker(H)=%d, %s (%.3g s)\n", n, rank_H, n - rank_H, buf,
          (double)(t1.tv_sec - t0.tv_sec) + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec));
  if (k <= 0)
    fprintf(stream, "# Warning: k=%d, there are no non-trivial codewords\n", k);
}

void var_init(int argc, char **argv, params_t * const p){
  int dbg=0;
  int swit=0;
  double prob=0.0;
  long long int dbg_ll=0;

  for (int i = 1; i < argc; i++) /* scan arguments for version */
    if ((strcmp(argv[i], "--version") == 0) || (strcmp(argv[i], "-version") == 0)) {
      printf("dist_m4ri version %s\n", DIST_M4RI_VERSION);
      exit(0);
    }

  for (int i = 1; i < argc; i++) /* scan arguments for full help message */
    if ((strcmp(argv[i], "--morehelp") == 0) || (strcmp(argv[i], "-morehelp") == 0)
        || (strcmp(argv[i], "--more-help") == 0) || (strcmp(argv[i], "-more-help") == 0)) {
      printf(MORE_HELP, argv[0], DIST_M4RI_VERSION, argv[0]);
      exit(0);
    }

  for (int i = 1; i < argc; i++) /* scan arguments for standard help message */
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)
        || (strcmp(argv[i], "-?") == 0)) {
      printf(USAGE, argv[0], DIST_M4RI_VERSION, argv[0]);
      exit(0);
    }

  if (argc <= 1) {
    fprintf(stderr, "%s: no command-line arguments given\n\n", argv[0]);
    print_short_help(argv[0]);
    exit(-1);
  }

  int refresh_set=0;
  int steps_set=0;

  /* `debug` is parsed first, so that DBG_ARGS applies to all arguments: the first `debug=` argument replaces the
   * default, further ones are OR-combined (e.g., `debug=0` alone is silent, `debug=0 debug=4` the same as `debug=4`) */
  int debug_set=0;
  for (int i = 1; i < argc; i++)
    if (0 == strncmp(argv[i], "debug=", 6)) {
      const int val = parse_debug_value(argv[i] + 6, argv[i]);
      p->debug = debug_set ? (p->debug | val) : val;
      debug_set = 1;
    }
  if (p->debug & DBG_ARGS)
    fprintf(stderr, "# debug=%d (0x%x)%s\n", p->debug, (unsigned int) p->debug, debug_set ? "" : " (default)");

  for(int i=1; i<argc; i++){
    if (0 == strncmp(argv[i], "debug=", 6)){ /** `debug`: already parsed above */
    }
    else if (sscanf(argv[i],"css=%d",&dbg)==1){
      p->css=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, css=%d\n",argv[i],p->css);
    }
    else if (0==strncmp(argv[i],"finH=",5)){ /** `finH` */
      p->finH = file_arg_value(argc, argv, &i, 5); /**< allow space before file name */
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, finH=%s; setting fin=\"\"\n",argv[i],p->finH);
      if (strlen(p->fin) != 0)
        fprintf(stderr, "# Warning: finH='%s' given: ignoring fin='%s'\n", p->finH, p->fin);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"finL=",5)){ /** `finL` */
      p->finL = file_arg_value(argc, argv, &i, 5); /**< allow space before file name */
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, finL=%s; setting fin=\"\"\n",argv[i],p->finL);
      if (strlen(p->fin) != 0)
        fprintf(stderr, "# Warning: finL='%s' given: ignoring fin='%s'\n", p->finL, p->fin);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"finG=",5)){/** `finG` degeneracy generator matrix */
      p->finG = file_arg_value(argc, argv, &i, 5); /**< allow space before file name */
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, finG=%s; setting fin=\"\"\n",argv[i],p->finG);
      if (strlen(p->fin) != 0)
        fprintf(stderr, "# Warning: finG='%s' given: ignoring fin='%s'\n", p->finG, p->fin);
      p->fin="";
    }
    else if (0==strncmp(argv[i],"fin=",4)){
      if(p->finH)
	ERROR("arg[%d]='%s' in conflict with finH=%s\n",i,argv[i],p->finH);
      if(p->finG)
	ERROR("arg[%d]='%s' in conflict with finG=%s\n",i,argv[i],p->finG);
      if(p->finL)
	ERROR("arg[%d]='%s' in conflict with finL=%s\n",i,argv[i],p->finL);
      p->fin = file_arg_value(argc, argv, &i, 4); /**< allow space before file name */
    }
    else if (sscanf(argv[i],"method=%d",&dbg)==1){
      p->method=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, method=%d\n",argv[i],p->method);
      if( (p->method<=0) || (p->method>3)) {
        fprintf(stderr, "%s: unsupported method=%d specified\n\n", argv[0], p->method);
        print_short_help(argv[0]);
        exit(-1);
      }
    }
    else if (sscanf(argv[i],"smax=%d",&dbg)==1){
      p->smax=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, smax=%d\n",argv[i],p->smax);
    }
    else if (sscanf(argv[i],"wmax=%d",&dbg)==1){
      p->wmax=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, wmax=%d\n",argv[i],p->wmax);
    }
    else if (sscanf(argv[i],"dmin=%d",&dbg)==1){
      p->dmin=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, dmin=%d\n",argv[i],p->dmin);
    }
    else if (sscanf(argv[i],"dmax=%d",&dbg)==1){
      p->dmax=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, dmax=%d\n",argv[i],p->dmax);
    }
    else if (0==strncmp(argv[i],"start=",6)){ /** `start=a,b,c` list of CC start columns */
      parse_start_list(argv[i]+6, p);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, start_num=%d\n",argv[i],p->start_num);
    }
    else if (sscanf(argv[i],"cbeg=%d",&dbg)==1){
      p->cbeg=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, cbeg=%d\n",argv[i],p->cbeg);
    }
    else if (sscanf(argv[i],"cend=%d",&dbg)==1){
      p->cend=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, cend=%d\n",argv[i],p->cend);
    }
    else if (sscanf(argv[i],"wmin=%d",&dbg)==1){
      p->wmin=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, wmin=%d\n",argv[i],p->wmin);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){
      p->steps=dbg;
      steps_set=1;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, steps=%d\n",argv[i],p->steps);
    }
    else if (sscanf(argv[i],"seed=%d",&dbg)==1){
      p->seed=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, seed=%d\n",argv[i],p->seed);      
    }    
    else if (sscanf(argv[i],"noscan=%d",&dbg)==1){
      p->noscan=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, noscan=%d\n",argv[i],p->noscan);
    }
    else if (0==strncmp(argv[i],"fdem=",5)){
      p->fdem = file_arg_value(argc, argv, &i, 5);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, fdem=%s\n",argv[i],p->fdem);
    }
    else if (sscanf(argv[i],"pmin=%lg",&prob)==1){
      p->pmin=prob;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, pmin=%g\n",argv[i],p->pmin);
    }
    else if (0==strncmp(argv[i],"finC=",5)){
      p->finC = file_arg_value(argc, argv, &i, 5);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, finC=%s\n",argv[i],p->finC);
    }
    else if (0==strncmp(argv[i],"outC=",5)){
      p->outC = file_arg_value(argc, argv, &i, 5);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, outC=%s\n",argv[i],p->outC);
    }
    else if (sscanf(argv[i],"maxC=%lld",&dbg_ll)==1){
      p->maxC=dbg_ll;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, maxC=%lld\n",argv[i],p->maxC);
    }
    else if (sscanf(argv[i],"dW=%d",&dbg)==1){
      p->dW=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, dW=%d\n",argv[i],p->dW);
    }
    else if (sscanf(argv[i],"classical=%d",&dbg)==1){
      p->classical=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, classical=%d\n",argv[i],p->classical);
    }
    else if (sscanf(argv[i],"threads=%d",&dbg)==1){
      p->threads=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, threads=%d\n",argv[i],p->threads);
    }
    else if (sscanf(argv[i],"dexp=%d",&dbg)==1){
      p->dexp=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, dexp=%d\n",argv[i],p->dexp);
    }
    else if (sscanf(argv[i],"dest=%d",&dbg)==1){
      p->dexp=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, dest=%d (alias for dexp)\n",argv[i],p->dexp);
    }
    else if (sscanf(argv[i],"dstop=%d",&dbg)==1){ /** `dstop`: stop once CC certifies dmin >= dstop */
      if (dbg < 0)
        ERROR("invalid '%s': dstop must be non-negative", argv[i]);
      p->dstop=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, dstop=%d\n",argv[i],p->dstop);
    }
    else if (sscanf(argv[i],"timeout=%lf",&prob)==1){
      p->timeout=prob;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, timeout=%g\n",argv[i],p->timeout);
    }
    else if (sscanf(argv[i],"nothrottle=%d",&dbg)==1){
      p->nothrottle=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, nothrottle=%d\n",argv[i],p->nothrottle);
    }
    else if (strcmp(argv[i], "--no-throttle") == 0 || strcmp(argv[i], "-no-throttle") == 0
             || strcmp(argv[i], "nothrottle") == 0) {
      p->nothrottle=1;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, nothrottle=1\n",argv[i]);
    }
    else if (sscanf(argv[i],"chunk_size=%d",&dbg)==1 || sscanf(argv[i],"batch=%d",&dbg)==1){
      p->chunk_size=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, chunk_size=%d\n",argv[i],p->chunk_size);
    }
    else if (sscanf(argv[i],"ksub=%d",&dbg)==1){
      p->ksub=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, ksub=%d\n",argv[i],p->ksub);
    }
    else if (sscanf(argv[i],"kwin=%d",&dbg)==1 || sscanf(argv[i],"win=%d",&dbg)==1){
      p->kwin=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, kwin=%d\n",argv[i],p->kwin);
    }
    else if (sscanf(argv[i],"win_mode=%d",&dbg)==1){
      p->win_mode=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, win_mode=%d\n",argv[i],p->win_mode);
    }
    else if (sscanf(argv[i],"min_hits=%d",&dbg)==1 || sscanf(argv[i],"max_hits=%d",&dbg)==1){
      p->min_hits=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, min_hits=%d\n",argv[i],p->min_hits);
    }
    else if (sscanf(argv[i],"cov_cws=%d",&dbg)==1){
      p->cov_cws=dbg;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, cov_cws=%d\n",argv[i],p->cov_cws);
    }
    else if (sscanf(argv[i],"refresh=%d",&dbg)==1){
      p->refresh=dbg;
      refresh_set=1;
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read %s, refresh=%d\n",argv[i],p->refresh);
    }
    else{ /* unrecognized option */
      fprintf(stderr, "%s: unrecognized parameter \"%s\" at position %d\n\n", argv[0], argv[i], i);
      print_short_help(argv[0]);
      exit(-1);
    }
  } /* end parameter scan cycle */

  if (!refresh_set && p->ksub > 0) {
    p->refresh = 5000;
  }

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

  if (p->dmin < 0) {
    ERROR("parameter dmin=%d cannot be negative\n", p->dmin);
  }
  if (p->dmax < 0) {
    ERROR("parameter dmax=%d cannot be negative\n", p->dmax);
  }
  if (p->dmin > 0 && p->dmax > 0 && p->dmin > p->dmax) {
    ERROR("parameter dmin=%d cannot be larger than dmax=%d\n", p->dmin, p->dmax);
  }
  if (p->wmax > 0 && p->dmin > p->wmax) {
    ERROR("parameter dmin=%d cannot be larger than wmax=%d\n", p->dmin, p->wmax);
  }

  if (p->wmax > 0 && p->wmin > p->wmax) {
    ERROR("parameter wmin=%d cannot be larger than wmax=%d\n", p->wmin, p->wmax);
  }
  if (p->start_num > 0) { /* `start` list uses unlimited clusters: cannot be mixed with split runs */
    if (p->cbeg >= 0 || p->cend >= 0) {
      ERROR("Cannot specify start along with cbeg or cend\n");
    }
  }

  if (p->method == 1) {
    if (p->cbeg >= 0 || p->cend >= 0 || p->start_num > 0) {
      ERROR("Parameters start, cbeg, and cend only work with CC method (method=2 or method=3)\n");
    }
  }

  if (p->cbeg >= 0 && p->cend >= 0 && p->cbeg > p->cend) {
    ERROR("cbeg=%d cannot be larger than cend=%d\n", p->cbeg, p->cend);
  }
  if (p->noscan && p->smax > 0) {
    fprintf(stderr,
            "# WARNING: smax=%d disabled (set to 0) because noscan=1 skips small cluster weights\n",
            p->smax);
    p->smax = 0;
  } else if (p->dmin > 1 && p->smax > 0) {
    fprintf(stderr,
            "# WARNING: smax=%d disabled (set to 0) because dmin=%d skips small cluster weights\n",
            p->smax, p->dmin);
    p->smax = 0;
  }

  if (p->method == 2) {
    if (steps_set) {
      fprintf(stderr, "# WARNING: steps=%d is ignored for CC method\n", p->steps);
    }
  }

  if (p->method == 1) { /* RW */
    if (p->steps <= 0)
      ERROR("parameter steps=%d should be positive for RW method=%d", p->steps, p->method);
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
    if (p->debug & DBG_ARGS)
      fprintf(stderr, "# read 'fin=%s'; " //"since switch=%d "
	     "assigning \n# finH=%s\n# finG=%s\n",
	     p->fin,// swit,
	     p->finH,p->finG);
  }
  
  if (p->fdem) {
    read_dem_file(p->fdem, &(p->spaH), &(p->spaL), p->pmin, p->debug);
    if (p->classical == -1) p->classical = 0;
    p->nvar = p->spaH->cols;
    p->n0 = p->nvar;
    p->nchk = p->spaL->rows;
    csr_dump(stderr, p->spaH, "H", p->debug);
    csr_dump(stderr, p->spaL, "L", p->debug);
  } else {
    if (p->finH){
      p->spaH=csr_mm_read(p->finH,p->spaH,0);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read H <- file '%s'\n",p->finH);
      csr_dump(stderr, p->spaH, "H", p->debug);
    }
    else
      ERROR("need to specify H=Hx input file name; use fin=[str] or finH=[str]\n");

    if((p->finG) && (p->finL))
      ERROR("either G=Hz or L=Lx matrix should be specified but not both! finG='%s' finL='%s'\n",
	    p->finG, p->finL);

    if(p->finG){
      if (p->classical == -1) p->classical = 0;
      p->spaG=csr_mm_read(p->finG,p->spaG,0);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read G <- file '%s'\n",p->finG);
      if(csr_csr_mul_non_zero(p->spaH, p->spaG))
	 ERROR("rows of H and G matrices are not orthogonal");
      csr_dump(stderr, p->spaG, "G", p->debug);
    } 
    else if (p->finL){
      if (p->classical == -1) p->classical = 0;
      p->spaL=csr_mm_read(p->finL,p->spaL,0);
      if (p->debug & DBG_ARGS)
	fprintf(stderr, "# read L <- file '%s'\n",p->finL);
      csr_dump(stderr, p->spaL, "L", p->debug);
      p->nchk = p->spaL->rows;
    } 
    else{
      if (p->classical == -1) p->classical = 1;
      p->spaG=NULL;
    }
  }

  if(p->method & 2){ /* CC */
    if ((p->wmax<=0) && ((p->method & 1 )==0)) {
      /* the CC rounds end with the timeout, or at the latest once the bounds coincide (w = dmax - 1), or once
       * dmin >= dstop (w = dstop - 1) */
      if (p->timeout <= 0.0 && p->dmax <= 0 && p->dstop <= 0) {
        ERROR("either parameter wmax>0, dmax>0, dstop>0, or timeout>0 should be specified for CC method=%d",
              p->method);
      }
      p->wmax = MAX_W - 1;
    }
    if(p->wmax>=MAX_W)
      ERROR("increase MAX_W=%d defined in 'util_io.h'",MAX_W);
    for(int i=0; i<MAX_W; i++)
      p->swei[i]=p->spaH->rows +1; 
  }

  if (p->seed<=0){
    p->seed = time(NULL) - 1000 * p->seed + 10*getpid();
    if (p->debug & DBG_ARGS)
      fprintf(stderr, "# initializing rng from time(NULL), seed=%d\n",p->seed);
  }
  else if (p->debug & DBG_ARGS)
    fprintf(stderr, "# initializing rng from seed=%d\n",p->seed);

  srand(p->seed);

  rci_t n = (p->spaH)-> cols;
  if ((!p->classical) && ((p->spaG) && (n != (p->spaG) -> cols)))
    ERROR("Column count mismatch in H and G matrices: %d != %d",
	  (p->spaH)-> cols, (p->spaG)->cols);
  p->nvar = n; 
  p->n0 = n;
  if (p->css!=1)
    ERROR("Non-CSS codes are currently not supported, css=%d",p->css);

  if (p->cbeg >= p->nvar) {
    ERROR("cbeg=%d cannot be larger than nvar-1=%d\n", p->cbeg, p->nvar-1);
  }
  if (p->cend >= p->nvar) {
    ERROR("cend=%d cannot be larger than nvar-1=%d\n", p->cend, p->nvar-1);
  }
  if (p->start_num > 0) { /* range check, sort, and remove duplicates */
    for (int k = 0; k < p->start_num; k++)
      if (p->start_list[k] >= p->nvar)
        ERROR("start column %d cannot be larger than nvar-1=%d\n", p->start_list[k], p->nvar-1);
    qsort(p->start_list, p->start_num, sizeof(int), cmp_int);
    int num = 1;
    for (int k = 1; k < p->start_num; k++)
      if (p->start_list[k] != p->start_list[num-1])
        p->start_list[num++] = p->start_list[k];
    p->start_num = num;
  }
  
  if((p->spaG) && (p->spaL==NULL)){
    /** create `Lx` */
    /** WARNING: this does not necessarily have minimal row weights */
    p->spaL = Lx_for_CSS_code(p->spaH,p->spaG);
    p->nchk = p->spaL->rows;
    csr_dump(stderr, p->spaL, "L (from H and G)", p->debug);
  }

  if (p->classical) {
    if (p->finL != NULL || p->finG != NULL) {
      ERROR("Conflict: classical=1 specified, but finL or finG was also provided.");
    }
    if (p->spaL != NULL) {
      if (p->debug & DBG_SUMMARY) {
        fprintf(stderr, "# Warning: classical=1 specified, discarding L matrix (logical operators)\n");
      }
      p->spaL = csr_free(p->spaL);
    }
  } else {
    if (p->spaL == NULL) {
      ERROR("L matrix (logical operators) is required for quantum code (classical=0).\n"
            "Provide finL, fdem, or finG to construct it. Alternatively, set classical=1 "
            "to find the distance of the stabilizer code as a classical code.");
    }
    /* column masks of L: the check L c != 0 of a candidate codeword c takes O(|c|) operations (for up to 64 rows
     * of L) instead of O(nnz(L)) */
    p->maskL = colmask_from_csr(p->spaL);
  }

  if ((p->method <= 0) || (p->method > 3)){
      fprintf(stderr, "%s: invalid method=%d specified\n\n", argv[0], p->method);
      print_short_help(argv[0]);
      exit(-1);
  }

  if (p->debug & DBG_SUMMARY)
    print_input_summary(stderr, p);
  if (p->debug & DBG_PARAMS)
    print_code_params(stderr, p);

  if ((p->debug & DBG_SUMMARY) && !p->fdem && (p->method & 2) && p->smax == 0) {
    fprintf(stderr, "# Warning: smax=0, confinement profile is not computed\n");
  }

  /* Warnings for the expert CC options and for the experimental ksub (printed regardless of `debug`) */
  if (p->noscan) {
    const int d0 = (p->dmin > 1) ? p->dmin : 1;
    if (d0 < p->wmax)
      fprintf(stderr,
              "# WARNING: noscan=1 (expert option): CC checks only weight w=wmax=%d; weights %d..%d are not\n"
              "#   scanned, so a codeword found only gives an upper bound (dmax), and dmin=%d is not raised\n",
              p->wmax, d0, p->wmax - 1, d0);
    else
      fprintf(stderr,
              "# WARNING: noscan=1 (expert option): CC checks only weight w=wmax=%d; the result relies on\n"
              "#   the supplied dmin=%d being a certified lower bound\n", p->wmax, p->dmin);
  }
  if (p->start_num > 0) {
    fprintf(stderr, "# WARNING: start=");
    for (int k = 0; (k < p->start_num) && (k < 8); k++)
      fprintf(stderr, "%s%d", k ? "," : "", p->start_list[k]);
    fprintf(stderr,
            "%s (expert option): CC clusters are grown only from %d listed column(s)\n"
            "#   (not limited to larger column indices); dmin is valid only if a code symmetry maps every\n"
            "#   minimum-weight codeword to one containing a listed column\n",
            (p->start_num > 8) ? ",..." : "", p->start_num);
  }
  if (p->cbeg >= 0 || p->cend >= 0) {
    const int beg = (p->cbeg >= 0) ? p->cbeg : 0;
    const int end = (p->cend >= 0) ? p->cend : p->nvar - 1;
    fprintf(stderr,
            "# WARNING: cbeg/cend (expert option for split runs): CC starts only from columns [%d,%d];\n"
            "#   dmin only covers codewords whose smallest column is in this range; take the minimum of\n"
            "#   dmin (and of dmax) over runs covering all columns 0..%d\n", beg, end, p->nvar - 1);
  }
  /* The experimental RW option ksub (with or without min_hits, also if it falls back to ksub=0 for m < nu) */
  if (p->ksub > 0 && (p->method & 1)) {
    fprintf(stderr,
            "# WARNING: ksub=%d is experimental and should not be used: a RW step reduces the span of ksub rows\n"
            "#   of a fixed basis of ker(H), so that a few codewords are found much more often than others; with\n"
            "#   min_hits, RW may end early, and exp(-<n>) underestimates the probability to miss a lighter codeword\n",
            p->ksub);
  }
}

void var_kill(params_t * const p){
  if (p->start_list) {
    free(p->start_list);
    p->start_list = NULL;
    p->start_num = 0;
  }
  if(p->spaL)
    csr_free(p->spaL);
  p->maskL = colmask_free(p->maskL);
  if(p->spaH)
    csr_free(p->spaH);
  if(p->spaG)
    csr_free(p->spaG);
  if (strlen(p->fin) != 0) {
    if (p->finH) {
      free(p->finH);
      p->finH = NULL;
    }
    if (p->finG) {
      free(p->finG);
      p->finG = NULL;
    }
  }

  cw_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    HASH_DEL(p->codewords, cw);
    free(cw);
  }
  p->num_cws = 0;
  p->min_w = INT_MAX;
  p->cw_max_w = 0;
  free(p->cw_cnt_w);
  free(p->cw_hits_w);
  p->cw_cnt_w = NULL;
  p->cw_hits_w = NULL;
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

  if (debug & DBG_ARGS)
    fprintf(stderr, "# read DEM %s: rows_H=%d rows_L=%d cols=%d; nz_H=%d nz_L=%d\n",fnam,r,k,n,iD,iL);
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
    fprintf(stderr, "FILE I/O ERROR: %s\n", strerror(errno));
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
  int c=fgetc(f);
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

  int c=fgetc(f);
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
    fprintf(stderr, "%s:%lld: invalid NZLIST entry\n", fnam, *lineno);
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
      fprintf(stderr, "%s:%lld: invalid entry of weight w=%d\n",fnam, *lineno, w);
      ERROR("expected an integer i=%d of %d",i,w);
    }    
    vec->arr[i]--;
  }

  for(int i=1; i<w; i++){
    if((vec->arr[i-1] < 0) || (vec->arr[i-1] >= vec->arr[i])){
      fprintf(stderr, "%s:%lld: invalid entry of weight w=%d\n",fnam, *lineno, w);
      ERROR("expected strictly increasing positive entries");
    }   
  }
  (*lineno)++;
  return vec;
}

/* Codewords are collected (and exported with outC) with outC or maxC */
static inline int cw_collecting(const params_t * const p) {
  return (p->outC != NULL) || (p->maxC > 0);
}

/* Width of the collection window: codewords of weight up to min_w + dW are collected (with outC or maxC) */
static inline int cw_window_dw(const params_t * const p) {
  return (cw_collecting(p) && p->dW > 0) ? p->dW : 0;
}

/* Size of the representative set of lowest-weight codewords for the `min_hits` statistic (0: only the
 * minimum-weight codewords) */
static inline long long int cw_rep_size(const params_t * const p) {
  return (p->min_hits > 0 && p->cov_cws > 0) ? p->cov_cws : 0;
}

long long int codeword_window_count(const params_t * const p) {
  if (p->min_w == INT_MAX || !p->cw_cnt_w) return 0;
  const long long int top = (long long int)p->min_w + cw_window_dw(p);
  const int w_top = (top < p->cw_max_w) ? (int)top : p->cw_max_w;
  long long int num = 0;
  for (int w = p->min_w; w <= w_top; w++) num += p->cw_cnt_w[w];
  return num;
}

int codeword_maxc_reached(const params_t * const p) {
  return (p->maxC > 0) && (codeword_window_count(p) >= p->maxC);
}

int codeword_feed_limit(const params_t * const p) {
  if (p->min_hits <= 0) return 0;
  if (p->min_w == INT_MAX) return INT_MAX;
  const long long int rep = cw_rep_size(p);
  if (rep <= 0) return p->min_w + 1;
  if (p->num_cws < rep) return INT_MAX;
  return p->cw_max_w + 1;
}

/* Remove the heaviest weight classes from the hash as long as the remaining codewords keep the collection window
 * and at least `cw_rep_size()` codewords (whole weight classes are kept) */
static void cw_trim(params_t * const p) {
  if (p->min_w == INT_MAX || p->cw_max_w <= p->min_w) return;
  const long long int rep = cw_rep_size(p);
  if (rep > 0 && p->num_cws <= rep) return;
  const long long int top = (long long int)p->min_w + cw_window_dw(p);
  if (top >= p->cw_max_w) return;
  int keep = (int)top;
  if (rep > 0) { /* the smallest weight such that the codewords up to this weight form the representative set */
    long long int cum = 0;
    int w_rep = p->cw_max_w;
    for (int w = p->min_w; w < p->cw_max_w; w++) {
      cum += p->cw_cnt_w[w];
      if (cum >= rep) {
        w_rep = w;
        break;
      }
    }
    if (w_rep > keep) keep = w_rep;
  }
  if (keep >= p->cw_max_w) return;
  cw_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    if (cw->weight > keep) {
      p->cw_cnt_w[cw->weight]--;
      p->cw_hits_w[cw->weight] -= cw->cnt;
      HASH_DEL(p->codewords, cw);
      free(cw);
      p->num_cws--;
    }
  }
  p->cw_max_w = keep;
  while (p->cw_max_w > p->min_w && p->cw_cnt_w[p->cw_max_w] == 0) p->cw_max_w--;
}

cw_vec_t * codeword_add_maybe(params_t * const p, const int arr[], int weight) {
  if (weight <= 0 || weight > p->nvar) return p->codewords;
  if (!p->cw_cnt_w) { /* weight histograms, allocated on first use */
    p->cw_cnt_w = calloc((size_t)p->nvar + 2, sizeof(long long int));
    p->cw_hits_w = calloc((size_t)p->nvar + 2, sizeof(long long int));
    if (!p->cw_cnt_w || !p->cw_hits_w) ERROR("memory allocation");
  }

  const size_t keylen = weight * sizeof(int);
  cw_vec_t *pvec = NULL;
  HASH_FIND(hh, p->codewords, arr, keylen, pvec);
  if (pvec) { /* a known codeword: one more hit */
    pvec->cnt++;
    p->cw_hits_w[weight]++;
    return p->codewords;
  }

  int insert;
  if (weight < p->min_w) {
    insert = 1; /* a new minimum weight is always recorded (also with maxC) */
  } else if ((long long int)weight <= (long long int)p->min_w + cw_window_dw(p)) { /* the collection window */
    if (cw_collecting(p))
      insert = !codeword_maxc_reached(p);
    else /* weight == min_w: at most cov_cws codewords tracked for `min_hits` (all if cov_cws <= 0) */
      insert = (p->cov_cws <= 0) || (p->cw_cnt_w[weight] < p->cov_cws);
  } else { /* heavier: only to keep a representative set of cov_cws lowest-weight codewords for `min_hits` */
    const long long int rep = cw_rep_size(p);
    insert = (rep > 0) && (p->num_cws < rep || weight < p->cw_max_w);
  }
  if (!insert) return p->codewords;

  cw_vec_t *entry = malloc(sizeof(cw_vec_t) + keylen);
  if (!entry) ERROR("memory allocation");
  entry->weight = weight;
  entry->cnt = 1;
  memcpy(entry->arr, arr, keylen);
  HASH_ADD(hh, p->codewords, arr, keylen, entry);
  p->num_cws++;
  p->cw_cnt_w[weight]++;
  p->cw_hits_w[weight]++;
  if (weight > p->cw_max_w) p->cw_max_w = weight;
  if (weight < p->min_w) p->min_w = weight;

  /* gradually replace heavier codewords by lighter ones */
  cw_trim(p);
  return p->codewords;
}

void compute_min_w_hit_stats(const params_t * const p, int *min_cnt, int *max_cnt,
                             double *avg_cnt, double *stdev_cnt) {
  *min_cnt = 0;
  *max_cnt = 0;
  *avg_cnt = 0.0;
  *stdev_cnt = 0.0;
  if (!p || !p->codewords || p->min_w == INT_MAX) return;

  long long n_cws = 0;
  long long total_hits = 0;
  int c_min = INT_MAX;
  int c_max = 0;
  cw_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    if (cw->weight == p->min_w) {
      n_cws++;
      total_hits += cw->cnt;
      if (cw->cnt < c_min) c_min = cw->cnt;
      if (cw->cnt > c_max) c_max = cw->cnt;
    }
  }
  if (n_cws <= 0) return;

  double avg = (double)total_hits / (double)n_cws;
  double sum_sq = 0.0;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    if (cw->weight == p->min_w) {
      double diff = (double)cw->cnt - avg;
      sum_sq += diff * diff;
    }
  }
  *min_cnt = c_min;
  *max_cnt = c_max;
  *avg_cnt = avg;
  *stdev_cnt = sqrt(sum_sq / (double)n_cws);
}

/* Hit statistics of the codewords of one weight in hash */
typedef struct {
  long long int n;    /* number of codewords */
  long long int hits; /* total hits */
  double sumsq;       /* sum of squared hit counts */
  int c_min, c_max;   /* smallest and largest hit count */
} cw_class_stats_t;

/* Hit statistics of the codewords in hash of each weight w_lo..w_hi (one pass over the hash; free() the result) */
static cw_class_stats_t *cw_class_stats(const params_t * const p, const int w_lo, const int w_hi) {
  const int nw = w_hi - w_lo + 1;
  cw_class_stats_t *cs = calloc(nw > 0 ? nw : 1, sizeof(cw_class_stats_t));
  if (!cs) ERROR("memory allocation");
  for (int i = 0; i < nw; i++) cs[i].c_min = INT_MAX;
  cw_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    if (cw->weight < w_lo || cw->weight > w_hi) continue;
    cw_class_stats_t * const c = &cs[cw->weight - w_lo];
    c->n++;
    c->hits += cw->cnt;
    c->sumsq += (double)cw->cnt * (double)cw->cnt;
    if (cw->cnt < c->c_min) c->c_min = cw->cnt;
    if (cw->cnt > c->c_max) c->c_max = cw->cnt;
  }
  return cs;
}

void codeword_hit_stats(const params_t * const p, cw_hit_stats_t * const st) {
  memset(st, 0, sizeof(*st));
  if (p->min_w == INT_MAX || !p->cw_cnt_w || p->cw_cnt_w[p->min_w] <= 0) return;
  const long long int rep = cw_rep_size(p);
  st->w_lo = st->w_hi = p->min_w;
  st->min_cws = p->cw_cnt_w[p->min_w];
  st->min_hits = p->cw_hits_w[p->min_w];
  for (int w = p->min_w; w <= p->cw_max_w; w++) { /* the lowest weight classes with at least `rep` codewords */
    if (p->cw_cnt_w[w] <= 0) continue;
    st->set_cws += p->cw_cnt_w[w];
    st->set_hits += p->cw_hits_w[w];
    st->w_hi = w;
    if (st->set_cws >= rep) break;
  }
  const double avg_min = (double)st->min_hits / (double)st->min_cws;
  const double avg_set = (double)st->set_hits / (double)st->set_cws;
  st->avg = (avg_set < avg_min) ? avg_set : avg_min;
}

void print_codeword_stats(FILE *stream, const params_t * const p) {
  if (!stream || !p) return;
  if (p->num_cws <= 0 || !p->codewords || p->min_w == INT_MAX) {
    fprintf(stream, "# codewords accumulated: total=0\n");
    return;
  }

  const int w_lo = p->min_w;
  const int w_hi = (p->cw_max_w > w_lo) ? p->cw_max_w : w_lo;
  cw_class_stats_t *cs = cw_class_stats(p, w_lo, w_hi);
  /* without DBG_LARGE, at most 4 heavier weight classes are shown individually, the others in one line */
  int shown = 0, rest_classes = 0, rest_lo = 0, rest_hi = 0, rest_min = INT_MAX, rest_max = 0;
  long long rest_cws = 0, rest_hits = 0;
  for (int w = w_lo; w <= w_hi; w++) {
    const cw_class_stats_t * const c = &cs[w - w_lo];
    if (c->n <= 0) continue;
    const double avg = (double)c->hits / (double)c->n;
    const double var = c->sumsq / (double)c->n - avg * avg;
    const double stdev = (var > 0.0) ? sqrt(var) : 0.0;
    if (w == w_lo) {
      fprintf(stream,
              "# codewords accumulated: total=%lld, min_w=%d: cws=%lld, total_hits=%lld, "
              "hits min=%d, max=%d, avg=%.2f, stdev=%.2f",
              p->num_cws, w, c->n, c->hits, c->c_min, c->c_max, avg, stdev);
      if (p->min_hits > 0) {
        cw_hit_stats_t st;
        codeword_hit_stats(p, &st);
        fprintf(stream, ", <n>=%.2f over %lld cws of w=%d..%d (min_hits=%d, cov_cws=%d)",
                st.avg, st.set_cws, st.w_lo, st.w_hi, p->min_hits, p->cov_cws);
      }
      fprintf(stream, "\n");
    } else if ((p->debug & DBG_LARGE) || (shown < 4)) {
      shown++;
      fprintf(stream,
              "# codewords w=%d: cws=%lld, total_hits=%lld, "
              "hits min=%d, max=%d, avg=%.2f, stdev=%.2f\n",
              w, c->n, c->hits, c->c_min, c->c_max, avg, stdev);
    } else {
      if (rest_classes++ == 0) rest_lo = w;
      rest_hi = w;
      rest_cws += c->n;
      rest_hits += c->hits;
      if (c->c_min < rest_min) rest_min = c->c_min;
      if (c->c_max > rest_max) rest_max = c->c_max;
    }
  }
  if (rest_classes > 0)
    fprintf(stream, "# codewords w=%d..%d: cws=%lld in %d weight classes, total_hits=%lld, hits min=%d, max=%d "
            "(each class with debug bit 2048)\n", rest_lo, rest_hi, rest_cws, rest_classes, rest_hits, rest_min,
            rest_max);
  free(cs);
}

int check_min_hits_convergence(const params_t * const p) {
  if (p->min_hits <= 0) return 0;
  cw_hit_stats_t st;
  codeword_hit_stats(p, &st);
  return (st.set_cws > 0) && (st.avg >= (double)p->min_hits);
}

/* Upper tail probability P(X >= x) of the chi-square distribution with `df` degrees of freedom (exact for df <= 2,
 * Wilson-Hilferty approximation otherwise) */
static double chi2_upper_tail(const double x, const long long int df) {
  if (!(x > 0.0)) return 1.0;
  if (df == 1) return erfc(sqrt(0.5 * x));
  if (df == 2) return exp(-0.5 * x);
  const double s = 2.0 / (9.0 * (double)df);
  const double z = (cbrt(x / (double)df) - (1.0 - s)) / sqrt(s);
  return 0.5 * erfc(z / sqrt(2.0));
}

/* Mean mu of the Poisson distribution whose zero-truncated version has the mean m > 1, i.e., m = mu / (1 - e^-mu):
 * Newton iteration for the convex function g(mu) = mu - m (1 - e^-mu), starting to the right of the root */
static double ztp_mu(const double m) {
  double mu = m;
  for (int it = 0; it < 200; it++) {
    const double e = exp(-mu);
    const double dg = 1.0 - m * e;
    if (!(dg > 0.0)) break;
    const double step = (mu - m * (1.0 - e)) / dg;
    mu -= step;
    if (fabs(step) <= 1e-12 * (1.0 + mu)) break;
  }
  return mu;
}

int codeword_hit_check(FILE *stream, const params_t * const p) {
  if (!p || p->min_w == INT_MAX || !p->cw_cnt_w) return 0;
  cw_hit_stats_t st;
  codeword_hit_stats(p, &st);
  if (st.set_cws <= 0) return 0;
  cw_class_stats_t *cs = cw_class_stats(p, st.w_lo, st.w_hi);
  int num_bad = 0;
  for (int w = st.w_lo; w <= st.w_hi; w++) {
    const cw_class_stats_t * const c = &cs[w - st.w_lo];
    if (c->n < 2) continue;
    const double n = (double)c->n;
    const double m = (double)c->hits / n;
    if (m <= 1.0 + 1e-9) continue; /* every codeword found once: no information */
    double s2 = (c->sumsq - n * m * m) / (n - 1.0); /* sample variance of the hit counts */
    if (s2 < 0.0) s2 = 0.0;
    /* with equal hit rates, a hit count is Poisson(mu) conditioned on >= 1 (zero-truncated), with mean m and
     * variance v0; codewords tracked only after their first hit have a smaller variance (the test is conservative) */
    const double mu = ztp_mu(m);
    const double v0 = m * (1.0 + mu - m);
    if (!(mu > 0.0) || !(v0 > 0.0)) continue;
    const double pval = chi2_upper_tail((n - 1.0) * s2 / v0, c->n - 1); /* dispersion test */
    const double cv2 = (s2 - v0) / (mu * mu); /* squared relative spread of the hit rates */
    if (!(pval < 1e-3) || !(cv2 >= 1.0)) continue; /* significant and strong: spread of hit rates >= their mean */
    /* a lighter codeword with a gamma-distributed hit rate of the same relative spread is missed with the
     * probability E[exp(-lambda)] = (1 + cv2 mu)^(-1/cv2), instead of exp(-mu) with equal hit rates: warn if this
     * is at least 3 times larger */
    const double p_unif = exp(-mu);
    const double p_spread = pow(1.0 + cv2 * mu, -1.0 / cv2);
    if (p_spread >= 3.0 * p_unif) {
      num_bad++;
      if (stream) {
        fprintf(stream,
                "# Warning: non-uniform RW hit counts of the %lld codewords of weight %d: hits min=%d, max=%d, "
                "avg=%.2f, stdev=%.2f (expected %.2f; p=%.1e)\n"
                "#   some codewords are found much more often than others of the same weight (relative spread of "
                "hit rates %.2f):\n"
                "#   a lighter codeword may be missed with probability ~%.2g rather than exp(-%.2f)=%.2g; "
                "consider a larger min_hits\n",
                c->n, w, c->c_min, c->c_max, m, sqrt(s2), sqrt(v0), pval, sqrt(cv2), p_spread, mu, p_unif);
      }
    }
  }
  free(cs);
  return num_bad;
}

/* Logarithm of the binomial coefficient C(a, b), 0 <= b <= a */
static double log_binom(const double a, const double b) {
  return lgamma(a + 1.0) - lgamma(b + 1.0) - lgamma(a - b + 1.0);
}

double rw_infoset_find_prob(const int n, const int rank, const int w) {
  const int k = n - rank; /* number of non-pivot columns */
  if (n <= 0 || rank < 0 || k <= 0 || w < 1 || w > n || w - 1 > rank) return 0.0;
  /* exactly one of the w positions among the k non-pivot columns */
  const double p1 = exp(log((double)w) + log_binom(n - w, k - 1) - log_binom(n, k));
  return (p1 < 1.0) ? p1 : 1.0;
}

/* Format a duration in seconds (s, min, h, or days) */
static void format_duration(char * const buf, const size_t size, const double t) {
  if (t < 120.0) snprintf(buf, size, "%.3g s", t);
  else if (t < 7200.0) snprintf(buf, size, "%.3g min", t / 60.0);
  else if (t < 172800.0) snprintf(buf, size, "%.3g h", t / 3600.0);
  else snprintf(buf, size, "%.3g days", t / 86400.0);
}

void print_rw_infoset_estimate(FILE *stream, const int n, const int rank, const int w_lo, const int w_hi,
                               const long steps, const long steps_total, const double t_step) {
  if (!stream || w_lo < 1 || w_hi < w_lo) return;
  int w_hard = w_hi; /* the weight hardest to find (smallest P1) */
  double p1 = rw_infoset_find_prob(n, rank, w_hi);
  for (int w = w_lo; w < w_hi; w++) {
    const double pw = rw_infoset_find_prob(n, rank, w);
    if (pw < p1) {
      p1 = pw;
      w_hard = w;
    }
  }
  if (!(p1 > 0.0)) {
    fprintf(stream, "#   information-set estimate: a codeword of weight %d cannot be found by RW (n=%d, rank(H)=%d)\n",
            w_hard, n, rank);
    return;
  }
  const double lq = (p1 < 1.0) ? log1p(-p1) : -INFINITY; /* log of the miss probability per step */
  const double p_miss = (steps > 0) ? exp((double)steps * lq) : 1.0;
  /* RW steps for a 1% miss probability: uniform steps, scaled to all steps with the fraction of uniform steps */
  const double unif_steps_1pc = (p1 < 1.0) ? ceil(log(0.01) / lq) : 1.0;
  const double steps_1pc = (steps > 0 && steps_total > steps) ? ceil(unif_steps_1pc * steps_total / steps)
                                                                : unif_steps_1pc;
  char steps_str[32], time_str[64] = "";
  snprintf(steps_str, sizeof(steps_str), (steps_1pc < 1e9) ? "%.0f" : "%.3g", steps_1pc);
  if (t_step > 0.0) {
    char dur[32];
    format_duration(dur, sizeof(dur), steps_1pc * t_step);
    snprintf(time_str, sizeof(time_str), ", ~%s", dur);
  }
  fprintf(stream,
          "#   information-set estimate (uniform random information sets, n=%d, rank(H)=%d): a codeword of weight %d\n"
          "#   is found with probability %.3g per RW step, missed in %ld uniform steps with probability %.2g "
          "(1%% after %s RW steps in total%s)\n",
          n, rank, w_hard, p1, steps, p_miss, steps_str, time_str);
}

long long int nzlist_read(const char fnam[], params_t *p){
  long long int count = 0, lineno;
  long long int skipped_invalid = 0;
  assert(fnam);
  FILE * f=nzlist_r_open(fnam, &lineno);
  if(!f){
    if ((p->outC ==NULL) || (strcmp(fnam,p->outC)!=0)){      
      fprintf(stderr, "codeword input file I/O ERROR: %s, outC=%s\n", strerror(errno),p->outC);
      ERROR("can't open file %s for reading",fnam);
    }
    else
      return 0;
  }
  cw_vec_t *entry=NULL;
  while((entry=nzlist_r_one(f,NULL, fnam, &lineno))){
    if (codeword_maxc_reached(p)) {
      free(entry);
      break;
    }
    int valid = (entry->weight > 0) && (entry->arr[entry->weight - 1] < p->nvar); /* column indices in range */
    if (valid && p->spaH) {
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
      const size_t keylen = entry->weight * sizeof(int);
      cw_vec_t *known = NULL, *added = NULL;
      HASH_FIND(hh, p->codewords, entry->arr, keylen, known);
      p->codewords = codeword_add_maybe(p, entry->arr, entry->weight);
      if (!known) {
        HASH_FIND(hh, p->codewords, entry->arr, keylen, added);
        if (added) count++;
      }
    }
    free(entry);
  }
  fclose(f);
  if (skipped_invalid > 0) {
    fprintf(stderr,
            "# Warning: skipped %lld invalid codewords (column out of range, not orthogonal to H, or orthogonal "
            "to L)\n",
            skipped_invalid);
  }
  if (p->debug & DBG_SUMMARY) {
    fprintf(stderr, "# read %lld codewords from %s, total %lld", count, fnam, p->num_cws);
    if (p->min_w != INT_MAX)
      fprintf(stderr, ", min weight %d", p->min_w);
    fprintf(stderr, "\n");
  }
  if ((p->debug & DBG_CODEWORDS) && (p->min_w != INT_MAX)) { /* a lightest codeword (the hash holds only finC ones) */
    cw_vec_t *cw, *tmp;
    HASH_ITER(hh, p->codewords, cw, tmp) {
      if (cw->weight == p->min_w) {
        print_codeword_support(stderr, "# finC: lightest codeword", cw->arr, cw->weight, p->debug);
        break;
      }
    }
  }
  return count; 
}

void print_codeword_support(FILE *stream, const char *prefix, const int arr[], const int weight, const int debug){
  if (!stream) return;
  const int max = ((debug & DBG_LARGE) || (weight <= DBG_CW_MAX)) ? weight : DBG_CW_MAX;
  fprintf(stream, "%s of weight %d (1-based columns):", prefix, weight);
  for (int i = 0; i < max; i++)
    fprintf(stream, " %d", arr[i] + 1);
  if (max < weight)
    fprintf(stream, " ... (%d more)", weight - max);
  fprintf(stream, "\n");
}

long long int nzlist_write(const char fnam[], const char comment[], params_t *p){
  long long int count=0;
  assert(fnam);
  FILE * f = nzlist_w_new(fnam, comment);
  /* only the collection window (weight up to min_w + dW): the hash may also hold heavier codewords, which are kept
   * for the `min_hits` statistic only */
  const long long int w_top = (p->min_w == INT_MAX) ? 0 : (long long int)p->min_w + cw_window_dw(p);
  cw_vec_t *pvec;
  
  for(pvec = p->codewords; pvec != NULL; pvec = (cw_vec_t *)(pvec->hh.next)){
    if (pvec->weight > w_top) continue;
    count ++;
    nzlist_w_append(f,pvec);
  }
  fclose(f);
  return count;
}
