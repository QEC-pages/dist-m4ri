#ifndef UTIL_M4RI_H
#define UTIL_M4RI_H

/************************************************************************ 
 * helper functions for use with m4ri library, including binary sparse
 * matrices and conversion utilities 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 * some code borrowed from various sources 
 ************************************************************************/

#define SWAPINT(a,b) do{ int t=a; a=b; b=t; } while(0)

#define ERROR(fmt,...) \
  do{ \
    fprintf (stderr, "%s:%d: *** ERROR in function '%s()' ***\n", __FILE__, __LINE__, __FUNCTION__); \
    fprintf(stderr, " [31;1m " fmt " [0m\n",##__VA_ARGS__); \
    exit(-1); \
  } \
  while(0)

/** `kludge` to work around the differences between old and new m4ri libraries */
static inline word const * mzd_row_cons(const mzd_t * mat, const int row){
  //  return mat->rows[row] ;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdiscarded-qualifiers"
  return mzd_row(mat,row);
#pragma GCC diagnostic pop
}


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
   * @brief Compute the number of set bits (Hamming weight) in a dense matrix A.
   * @param A Pointer to the dense matrix.
   * @return Hamming weight of the matrix.
   */
  size_t mzd_weight(const mzd_t *A);

  /** 
   * @brief Naive implementation to compute the Hamming weight of a dense matrix A.
   * @param A Pointer to the dense matrix.
   * @return Hamming weight of the matrix.
   */
  size_t mzd_weight_naive(const mzd_t *A);

  /**
   * @brief Compute the Hamming weight of a specific row in a dense matrix.
   * @param A Pointer to the dense matrix.
   * @param i Index of the row.
   * @return Hamming weight of row i.
   */  
  size_t mzd_weight_row(const mzd_t *A, rci_t i);

  /** 
   * @brief Find the position of the next non-zero bit in a raw word vector.
   * 
   * Searches the word array set1 of length m for the first set bit at or after pos.
   *
   * @param set1 Pointer to the array of words representing the bit vector.
   * @param m Length of the set1 array in words.
   * @param pos Bit position to start searching from.
   * @return Position of the next set bit, or -1 if none found.
   */
  static inline int nextelement(const word * const set1, const int m, const int pos){
    word setwd;
    int w;
    if (pos < 0){
      w = 0;
      setwd = set1[0];
    }
    else{
	w = SETWD(pos);
	if (w >= m) return -1;
	setwd = set1[w] & (m4ri_ffff<< SETBT(pos));
      }

    for (;;){
      if (setwd != 0)
	return TIMESWORDSIZE(w) + FIRSTBIT(setwd);
      if (++w == m)
	return -1;
      setwd = set1[w];
    }
  }


  /**
   * @brief Perform Gaussian elimination (naive) and return pivot column list.
   * 
   * Performs row reduction on M. If full is 1, performs full row reduction (RREF),
   * otherwise only upper triangular. Stores the pivot column indices in q.
   *
   * @param M Matrix to reduce (modified in place).
   * @param q Permutation structure to store pivot columns.
   * @param full Flag for full reduction (1) or triangular (0).
   * @return Rank of the matrix.
   */
  rci_t mzd_gauss_naive(mzd_t *M, mzp_t *q, int full);

  /** 
   * @brief Get the maximum row weight of a CSR sparse matrix.
   * @param p Pointer to the CSR sparse matrix.
   * @return Maximum row weight.
   */
  int csr_max_row_wght(const csr_t * const p);
  
  /** 
   * @brief Transpose a compressed CSR sparse matrix.
   * 
   * Transposes matrix p and stores the result in dst. (Re)allocates dst if needed.
   *
   * @param dst Destination sparse matrix (can be NULL, in which case it is allocated).
   * @param p Source sparse matrix.
   * @return Pointer to the transposed sparse matrix (dst or newly allocated).
   */
  csr_t * csr_transpose(csr_t *dst, const csr_t * const p);

  
  /**
   * @brief Convert a CSR sparse matrix to an MZD dense matrix.
   * 
   * @param dst Destination dense matrix (can be NULL or must have correct dimensions).
   * @param p Source sparse matrix.
   * @return Pointer to the dense matrix.
   */
  mzd_t *mzd_from_csr(mzd_t *dst, const csr_t *p);

  /**
   * @brief Construct the generator matrix from a parity check matrix in CSR form.
   * 
   * Permutes columns of H to bring it to standard form [ I C ], computes
   * the generator matrix G = [ C^T I ], and permutes columns of G back to align
   * with the original H.
   *
   * @param G Destination dense matrix for generator (can be NULL).
   * @param H Parity check matrix in CSR format.
   * @return Pointer to the generator matrix G.
   */
  mzd_t *mzd_generator_from_csr(mzd_t *G, const csr_t * const H);

  /**
   * @brief Multiply a sparse matrix by a dense matrix.
   * 
   * Computes C = S * B (if clear is 1) or C = C + S * B (if clear is 0).
   *
   * @param C Destination dense matrix (can be NULL).
   * @param S Source sparse matrix.
   * @param B Source dense matrix.
   * @param clear Flag to clear C before adding the product.
   * @return Pointer to the resulting dense matrix C.
   */
  mzd_t * csr_mzd_mul(mzd_t *C, const csr_t *S, const mzd_t *B, int clear);

  /** 
   * @brief Calculate syndrome vector change.
   * 
   * Computes syndrome = syndrome + row * spaQ (or sets it if clear is 1).
   *
   * @param syndrome Destination syndrome vector.
   * @param row Error vector (dense, 1 row).
   * @param spaQ Sparse check matrix.
   * @param clear Flag to clear syndrome before accumulation.
   * @return Pointer to the updated syndrome vector.
   */ 
  mzd_t * syndrome_vector(mzd_t *syndrome, mzd_t *row, csr_t *spaQ, int clear);


  /**
   * @brief Compute the weight of the product of a sparse matrix and a dense matrix.
   * 
   * Computes Hamming weight of A * B (if transpose is 0) or A * B^T (if transpose is 1).
   *
   * @param A Source sparse matrix.
   * @param B Source dense matrix.
   * @param transpose Flag to transpose B.
   * @return Hamming weight of the product.
   */
  size_t product_weight_csr_mzd(const csr_t *A, const mzd_t *B, int transpose);

  /**
   * @brief Return a uniformly distributed random integer in the range [0, max-1].
   * @param max Upper bound (exclusive).
   * @return Random integer.
   */
  int rand_uniform(const int max);

  /**
   * @brief Generate a random permutation of a given length in-place.
   * 
   * Modifies the permutation q by randomizing the first 'length' positions.
   * Uses LAPACK style pivot permutations.
   *
   * @param q Permutation structure to modify.
   * @param length Number of positions to permute.
   * @return Pointer to the modified permutation q.
   */ 
  mzp_t * mzp_rand_len(mzp_t *q, rci_t length);

  /**
   * @brief Generate a random permutation of the full length of q in-place.
   * @param q Permutation structure to modify.
   * @return Pointer to the modified permutation q.
   */
  static inline mzp_t * mzp_rand(mzp_t *q){
    if (q==NULL){
      printf("mzp_rand: permutation must be initialized!");
      exit(-1);
    }
    return mzp_rand_len(q,q->length);
  }

  /**
   * @brief Print the permutation to stdout.
   * @param p Pointer to the permutation structure.
   */
  void mzp_out(mzp_t const *p);

  /**
   * @brief Apply pivot permutation p to permutation q in-place starting from index.
   * 
   * @param q Destination permutation (initialized to identity if NULL).
   * @param p Source pivot permutation.
   * @param start Index to start applying the permutation from.
   * @return Pointer to the updated permutation q.
   */
  mzp_t *perm_p(mzp_t *q, const mzp_t *p,rci_t start);

  /**
   * @brief Apply transposed pivot permutation p to permutation q in-place.
   * 
   * @param q Destination permutation (initialized to identity if NULL).
   * @param p Source pivot permutation.
   * @param start Index to start applying the permutation from.
   * @return Pointer to the updated permutation q.
   */
  mzp_t *perm_p_trans(mzp_t *q, const mzp_t *p,const rci_t start);

  /**
   * @brief Free memory allocated for a CSR sparse matrix.
   * @param p Pointer to the CSR sparse matrix to free.
   * @return Always returns NULL.
   */
  csr_t *csr_free(csr_t *p);

  /**
   * @brief Initialize or reallocate a CSR sparse matrix to target size.
   * 
   * Checks dimensions and nzmax of mat, reallocating arrays if they are too small.
   *
   * @param mat Existing matrix structure (can be NULL, in which case it is allocated).
   * @param rows Target number of rows.
   * @param cols Target number of columns.
   * @param nzmax Target capacity (number of non-zero entries).
   * @return Pointer to the initialized/reallocated CSR matrix.
   */
  csr_t *csr_init(csr_t *mat, int rows, int cols, int nzmax);

  /**
   * @brief Compress CSR matrix (sort indices and remove duplicates/zeros).
   * 
   * Compresses the representation of the sparse matrix and sorts column indices
   * within each row to ensure consistency (critical for binary search / merging).
   *
   * @param mat Pointer to the CSR matrix to compress in-place.
   */ 
  void csr_compress(csr_t *mat);

  /**
   * @brief Construct a CSR sparse matrix from an array of coordinate pairs.
   * 
   * @param mat Destination CSR matrix (can be NULL).
   * @param nz Number of coordinate pairs.
   * @param prs Array of row-column index pairs.
   * @param nrows Number of rows in the matrix.
   * @param ncols Number of columns in the matrix.
   * @return Pointer to the constructed CSR matrix.
   */
  csr_t * csr_from_pairs(csr_t *mat, const int nz, int_pair * const prs, const int nrows, const int ncols);

  /**
   * @brief Print raw contents of a CSR matrix to stdout.
   * @param mat Pointer to the CSR matrix.
   */ 
  void csr_out(const csr_t *mat);

  /**
   * @brief Print a CSR matrix with a label string (formatted output).
   * @param smat Pointer to the CSR matrix.
   * @param str Label/name of the matrix.
   */
  void csr_print(const csr_t * const smat, const char str[]);

  /**
   * @brief Read a sparse matrix from a Matrix Market (.mtx) file into CSR format.
   * 
   * @param fin Path to the Matrix Market file.
   * @param mat Destination CSR matrix (can be NULL).
   * @param transpose Set to 1 to transpose the matrix while reading.
   * @return Pointer to the loaded CSR matrix.
   */
  csr_t *csr_mm_read(char *fin, csr_t *mat, int transpose);

  /** 
   * @brief Permute columns of a CSR sparse matrix.
   * 
   * @param dst Destination sparse matrix (can be NULL).
   * @param src Source sparse matrix.
   * @param perm Column permutation to apply.
   * @return Pointer to the permuted sparse matrix.
   */
  csr_t *csr_apply_perm(csr_t *dst, const csr_t * const src, const mzp_t * const perm);


  /**
   * @brief Flip the bit at position M[row,col] in a dense matrix.
   * @param M Dense matrix.
   * @param row Row index.
   * @param col Column index.
   * @note No bounds checks are performed.
   */

static inline void mzd_flip_bit(mzd_t * const M, rci_t const row, rci_t const col ) {
  word * const rawrow = mzd_row(M,row);
  __M4RI_FLIP_BIT(rawrow[col/m4ri_radix], col%m4ri_radix);
}

/**
 * @brief Perform one step of Gaussian elimination on column idx of M.
 * @param M Dense matrix.
 * @param idx Column index to eliminate.
 * @param begrow Starting row index.
 * @return Number of pivots found (0 or 1).
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
 * @brief Compute the syndrome of a sparse error vector against H and check if non-zero.
 * 
 * Relies on sorted column indices in H (must be compressed).
 *
 * @param H Parity check matrix (CSR format).
 * @param cnt Weight of the error vector.
 * @param ee Array of sorted indices representing the error support.
 * @return 1 if the syndrome is non-zero, 0 if it is zero.
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

  /** 
   * @brief Check if the product of two sparse matrices A * B^T is non-zero.
   * @param A First sparse matrix.
   * @param B Second sparse matrix.
   * @return 1 if the product is non-zero (matrices are not orthogonal), 0 otherwise.
   */  
  int csr_csr_mul_non_zero(const csr_t * const A, const csr_t * const B);
  
  /** 
   * @brief Compute the Hamming weight of the syndrome of a dense row vector against spaQ.
   * @param row Dense row vector (1 row).
   * @param spaQ Sparse check matrix.
   * @return Hamming weight of the syndrome.
   */ 
  int syndrome_bit_count(const mzd_t * const row, const csr_t * const spaQ);

  /** 
   * @brief Reduce a row vector against a matrix in standard form.
   * 
   * Checks if row is linearly dependent on the rows of matP0 (using rankP0 rows).
   * Modifies row in place.
   *
   * @param row Dense row vector to reduce (modified in place).
   * @param matP0 Dense matrix in standard form.
   * @param rankP0 Rank of matP0 (number of active rows).
   * @return 1 if row was reduced to zero (linearly dependent), 0 otherwise.
   */

  int do_reduce(mzd_t *row, const mzd_t *matP0, const rci_t rankP0);

  /**
   * @brief Generate a random binary error vector with error probability p.
   * @param row Dense row vector to store the error (modified in place).
   * @param p Bit-flip probability.
   */
  void make_err(mzd_t *row, double p);

  /**
   * @brief (Internal) Run cluster distance algorithm on matrices.
   */
  int do_dist_clus(const csr_t * const P, const mzd_t * const G, int debug, int wmax, int start, const int rankG);

  /**
   * @brief Construct the logical matrix Lx for a quantum CSS code given Hx and Hz.
   * @param Hx Sparse X-check matrix.
   * @param Hz Sparse Z-check matrix (passed as Hz).
   * @return Pointer to the constructed sparse logical matrix Lx.
   */
  csr_t * Lx_for_CSS_code(const csr_t * const Hx, const csr_t *const Hz);
  
#if defined(__cplusplus) && !defined (_MSC_VER)
}
#endif
  
#endif /* UTIL_M4RI_H */
