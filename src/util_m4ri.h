#ifndef UTIL_M4RI_H
#define UTIL_M4RI_H

/************************************************************************ 
 * helper functions for use with m4ri library, including binary sparse
 * matrices and conversion utilities 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 * some code borrowed from various sources 
 ************************************************************************/

#define SWAPINT(a,b) do{ int t=a; a=b; b=t; } while(0)

#define ERROR(fmt,...)                                                 \
  do{                                                                  \
    fprintf (stderr, "%s:%d: *** ERROR in function '%s()' ***\n", __FILE__, __LINE__, __FUNCTION__); \
    printf("[31;1m " fmt "[0m\n",##__VA_ARGS__); \
    exit(-1);                                                          \
  }                                                                    \
  while(0)



/**
 * macros from nauty.h
 * SETWD(pos) gives the setword in which pos is located
 * SETBT(pos) gives the location of bit pos in a setword
 */
#define SETWD(pos) ((pos)>>6)
#define SETBT(pos) ((pos)&0x3F)
#define TIMESWORDSIZE(w) ((w)<<6)    /* w*WORDSIZE */

#define FIRSTBIT(x) __builtin_ctzll(x) // number of trailing zeros 

// #ifdef __POPCNT__
#if 1

/* 
 * optimized code copied verbatim from https://danluu.com/assembly-intrinsics/ 
 */
/* uint32_t builtin_popcnt_unrolled_errata_manual(const uint64_t* buf, int len) { */
/*   assert(len % 4 == 0); */
/*   uint64_t cnt[4]; */
/*   for (int i = 0; i < 4; ++i) { */
/*     cnt[i] = 0; */
/*   } */

/*   for (int i = 0; i < len; i+=4) { */
/*     __asm__( */
/*         "popcnt %4, %4  \n\t" */
/*         "add %4, %0     \n\t" */
/*         "popcnt %5, %5  \n\t" */
/*         "add %5, %1     \n\t" */
/*         "popcnt %6, %6  \n\t" */
/*         "add %6, %2     \n\t" */
/*         "popcnt %7, %7  \n\t" */
/*         "add %7, %3     \n\t" // +r means input/output, r means intput */
/*         : "+r" (cnt[0]), "+r" (cnt[1]), "+r" (cnt[2]), "+r" (cnt[3]) */
/*         : "r"  (buf[i]), "r"  (buf[i+1]), "r"  (buf[i+2]), "r"  (buf[i+3])); */
/*   } */
/*   return cnt[0] + cnt[1] + cnt[2] + cnt[3]; */
/* } */

static inline int m4ri_bitcount(word w){
  return __builtin_popcountll(w);  
}

#else /* no __POPCNT__ */

#define MASK(c)    (((uint64_t)(-1)) / (__M4RI_TWOPOW(__M4RI_TWOPOW(c)) + 1))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (__M4RI_TWOPOW(c))) & MASK(c))

static inline int m4ri_bitcount(word w)  {
   uint64_t n = __M4RI_CONVERT_TO_UINT64_T(w);
   n = COUNT(n, 0);
   n = COUNT(n, 1);
   n = COUNT(n, 2);
   n = COUNT(n, 3);
   n = COUNT(n, 4);
   n = COUNT(n, 5);
   return (int)n;
}

static inline int std_bitcount ( uint64_t x) {
  uint64_t c1 = UINT64_C (0x5555555555555555 );
  uint64_t c2 = UINT64_C (0x3333333333333333 );
  uint64_t c4 = UINT64_C (0x0F0F0F0F0F0F0F0F );
  x -= (x >> 1) & c1;
  x = (( x >> 2) & c2) + (x & c2);
  x = ( x + (x >> 4) ) & c4;
  x *= UINT64_C (0x0101010101010101 );
  return (int) (x >> 56);
}

#endif /* __POPCNT__ */


/**
 * sparse binary matrix in compressed-row form (CSR, nz=-1) or 
 * List-Of-Pairs (nz pairs).
 * use mzp_compress() to convert from LOP to CSR. 
 */
typedef struct{    /*  */
  int rows ;	    /* number of rows */
  int cols ;	    /* number of columns */
  int nz ;	    /* # of entries in triplet matrix */
  int nzmax ;	    /* # allocated size */
  int *p ;	    /* row pointers (size rows+1) OR row indices */
  int *i ;	    /* col indices, size nzmax */
} csr_t ;


typedef struct { int a; int b; } int_pair;


#if defined(__cplusplus) && !defined (_MSC_VER)
extern "C" {
#endif

  /** 
   * number of set bits in a matrix.  
   * TODO: use the built-in version for the entire matrix
   */
  size_t mzd_weight(const mzd_t *A);
  size_t mzd_weight_naive(const mzd_t *A);
  /**
   * number of set bits in the row i of a matrix 
   */  
  size_t mzd_weight_row(const mzd_t *A, rci_t i);

  /**
   * nextelement(set1,m,pos) = the position of the first element in set set1   
   * which occupies a position greater or equal than pos.  If no such element exists,   
   * the value is -1.  pos can have any value less than n, including negative  
   * values.                                                                   
   *  
   * near verbatim copy from naututil.c (Nauty library by Brendan McKay)
   */
  //  int nextelement(word *set1, int m, int pos);

  /** 
   * return first non-zero bit in raw set-word vector set1 
   * of length m, starting with position pos.
   * with all zero bits, return -1 or number outside the range
   */
  static inline int nextelement(const word * const set1, const int m, const int pos){
    word setwd;
    int w;
#if 1
    if (pos < 0){
      w = 0;
      setwd = set1[0];
    }
    else
#endif 
      //    {
      w = SETWD(pos);
    setwd = set1[w] & (m4ri_ffff<< SETBT(pos));
    //  }

    for (;;){
      if (setwd != 0) return  TIMESWORDSIZE(w) + FIRSTBIT(setwd);
      if (++w == m) return -1;
      setwd = set1[w];
    }
  }


  /**
   * Copy of mzd_gauss_delayed from mzd.c (m4ri package) except additionally 
   * returns the list of pivot columns in second argument 
   */
  rci_t mzd_gauss_naive(mzd_t *M, mzp_t *q, int full);

  /** 
   * return max row weight of CSR matrix p
   * TODO: add code for List of Pairs 
   */
  int csr_max_row_wght(const csr_t * const p);
  
  /** 
   * transpose compressed CSR matrix, 
   * (re) allocate if needed 
   * return resulting matrix
   * TODO: add code for List of Pairs 
   */
  csr_t * csr_transpose(csr_t *dst, const csr_t * const p);

  
  /**
   * Convert CSR sparse binary matrix to MZD
   * allocate dst if needed (must be correct size or NULL)
   */
  mzd_t *mzd_from_csr(mzd_t *dst, const csr_t *p);

  /**
   * Convert a sparse binary matrix CSR into a standard form [ I C ],
   * with some col permutations if needed, create the dense generator
   * [ CT I ], and permute the cols back.  
   * (re)allocate G if needed.
   */
  mzd_t *mzd_generator_from_csr(mzd_t *G, const csr_t * const H);

  /**
   * sparse-S by dense B multiplication
   * C=C+S*B; allocate C if needed.
   */
  mzd_t * csr_mzd_mul(mzd_t *C, const csr_t *S, const mzd_t *B, int clear);

  /** 
   * Calculate the syndrome vector change: syndrome=syndrome +row.spaQ
   * optionally clear the destination 
   */ 
  mzd_t * syndrome_vector(mzd_t *syndrome, mzd_t *row, csr_t *spaQ, int clear);


  /**
   * helper function to compute the weight of the product 
   * A*B (transpose == 0) or A*B^T (transpose == 1)
   * with A sparse, B dense binary matrices
   */
  size_t product_weight_csr_mzd(const csr_t *A, const mzd_t *B, int transpose);

  /**
   * return uniformly distributed random number in the range [0,...,max-1] 
   * uses RAND internally 
   * \todo Replace by a better generator 
   */
  int rand_uniform(const int max);

  /**
   * replace pivot q with a random pivot smaller or equal length 
   * (second parameter); remaining positions in place, 
   * *** note: LAPACK style pivot permutations! ***
   * return pointer to q.
   * input: perm -- existing permutation
   */ 
  mzp_t * mzp_rand_len(mzp_t *q, rci_t length);

  /* same as above but of equal length */
  static inline mzp_t * mzp_rand(mzp_t *q){
    if (q==NULL){
      printf("mzp_rand: permutation must be initialized!");
      exit(-1);
    }
    return mzp_rand_len(q,q->length);
  }


  /**
   * print out the permutation (only needed under windows)
   */
  void mzp_out(mzp_t const *p);

  /**
   * apply pivot p to permutation q in place from start; 
   * initialize q to identity permutation if NULL
   * return q 
   */
  mzp_t *perm_p(mzp_t *q, const mzp_t *p,rci_t start);

  /**
   * apply pivot p (transposed) to permutation q in place from start; 
   * initialize q to identity permutation if NULL
   * return q 
   */
  mzp_t *perm_p_trans(mzp_t *q, const mzp_t *p,const rci_t start);


  /**
   * kill a CSR matrix 
   */
  csr_t *csr_free(csr_t *p);

  /**
   * initialize a CSR matrix 
   * check existing size and (re)allocate if  needded 
   */
  csr_t *csr_init(csr_t *mat, int rows, int cols, int nzmax);

  /**
   *  compress a CSR matrix  
   */ 
  void csr_compress(csr_t *mat);

  /**
   *  output a CSR matrix  
   */ 
  void csr_out(const csr_t *mat);
  void csr_print(const csr_t * const smat, const char str[]);
  /**
   * read sparse matrix into a (binary) CSR (all entries default to 1)
   * (re)allocate mat if needed
   * use transpose=1 to transpose.
   */
  csr_t *csr_mm_read(char *fin, csr_t *mat, int transpose);

  /** 
   * Permute columns of a CSR matrix with permutation perm.
   */
  csr_t *csr_apply_perm(csr_t *dst, const csr_t * const src, const mzp_t * const perm);


  /**
   * \brief Flip the bit at position M[row,col].
   *
   * \param M Matrix
   * \param row Row index
   * \param col Column index
   *
   * \note No bounds checks whatsoever are performed.
   *
   */

  static inline void mzd_flip_bit(mzd_t const *M, rci_t const row, rci_t const col ) {
    __M4RI_FLIP_BIT(M->rows[row][col/m4ri_radix], col%m4ri_radix);
  }


/**
 * @brief one step of gauss on column `idx` of matrix `M`
 * @param M the matrix
 * @param idx index of the column of `M` to deal with
 * @param begrow row to start with
 * @return number of pivot points found, `0` or `1` only
 */
static inline int gauss_one(mzd_t *M, const int idx, const int begrow){
  /** note: force-inlining actually slows it down (`???`) */
  rci_t startrow = begrow;
  rci_t pivots = 0;
  const rci_t i = idx;
  //  for (rci_t i = startcol; i < endcol ; ++i) {
  for(rci_t j = startrow ; j < M->nrows; ++j) {
    if (mzd_read_bit(M, j, i)) {
      mzd_row_swap(M, startrow, j);
      ++pivots;
      for(rci_t ii = 0 ;  ii < M->nrows; ++ii) {
        if (ii != startrow) {
          if (mzd_read_bit(M, ii, i)) {
            mzd_row_add_offset(M, ii, startrow,0);
          }
        }
      }
      startrow = startrow + 1;
      break;
    }
  }
  //  }
  return pivots; /** 0 or 1 only */
  // if one, need to update the current pivot list
}

/**
 * @brief return 1 if syndrome is non-zero
 * 
 * @param ee vector with `cnt` ordered set bit coordinates 
 */

static inline int sparse_syndrome_non_zero(const csr_t * const H, const int cnt, const int ee[]){
  for(int ir=0; ir < H->rows; ir++){
    int nz=0;
    for(int iL = H->p[ir], iE = 0; iL < H->p[ir+1]; iL++){
      int ic = H->i[iL];
      while((iE < cnt) && (ee[iE] < ic))
	iE++;
      if(iE >= cnt)
	break;
      if(ee[iE]==ic)
	nz ^= 1;      
    }
    if (nz)
      return 1;
  }
  return 0;
}

  /** @brief return 1 if matrix product A*B^T is non-zero 
   * @param A first matrix
   * @param B second matrix 
   * */  
  int csr_csr_mul_non_zero(const csr_t * const A, const csr_t * const B);
  
  /** 
   * Check whether syndrome is zero or not 
   */ 
  int syndrome_bit_count(const mzd_t * const row, const csr_t * const spaQ);

  /** 
   * Check if row is linearly dependent with the rows of matP0
   * which is assumed to be in standard form.
   * rankP0 is the number of non-zero rows in rankP0.
   */

  int do_reduce(mzd_t *row, const mzd_t *matP0, const rci_t rankP0);

  /**
   * generate binary error vector with error probability p 
   */
  void make_err(mzd_t *row, double p);

  int do_dist_clus(const csr_t * const P, const mzd_t * const G, int debug, int wmax, int start, const int rankG);

  csr_t * Lx_for_CSS_code(const csr_t * const Hx, const csr_t *const Hz);
  
#if defined(__cplusplus) && !defined (_MSC_VER)
}
#endif
  
#endif /* UTIL_M4RI_H */
