/* 
*   Matrix Market I/O library for ANSI C
*
*   See http://math.nist.gov/MatrixMarket for details.
*
*
*/

#ifndef MM_IO_H
#define MM_IO_H

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Convert a Matrix Market typecode to a human-readable string.
 * @param matcode The Matrix Market typecode.
 * @return Dynamically allocated string describing the typecode.
 */
char *mm_typecode_to_str(MM_typecode matcode);

/**
 * @brief Read the Matrix Market banner from a file.
 * @param f File pointer.
 * @param matcode Pointer to store the read typecode.
 * @return 0 on success, or an error code on failure.
 */
int mm_read_banner(FILE *f, MM_typecode *matcode);

/**
 * @brief Read dimensions and number of non-zero entries for a coordinate matrix.
 * @param f File pointer.
 * @param M Pointer to store the number of rows.
 * @param N Pointer to store the number of columns.
 * @param nz Pointer to store the number of non-zero entries.
 * @return 0 on success, or an error code on failure.
 */
int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);

/**
 * @brief Read dimensions for an array (dense) matrix.
 * @param f File pointer.
 * @param M Pointer to store the number of rows.
 * @param N Pointer to store the number of columns.
 * @return 0 on success, or an error code on failure.
 */
int mm_read_mtx_array_size(FILE *f, int *M, int *N);

/**
 * @brief Write the Matrix Market banner to a file.
 * @param f File pointer.
 * @param matcode The typecode to write.
 * @return 0 on success, or an error code on failure.
 */
int mm_write_banner(FILE *f, MM_typecode matcode);

/**
 * @brief Write dimensions and number of non-zero entries for a coordinate matrix.
 * @param f File pointer.
 * @param M Number of rows.
 * @param N Number of columns.
 * @param nz Number of non-zero entries.
 * @return 0 on success, or an error code on failure.
 */
int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);

/**
 * @brief Write dimensions for an array (dense) matrix.
 * @param f File pointer.
 * @param M Number of rows.
 * @param N Number of columns.
 * @return 0 on success, or an error code on failure.
 */
int mm_write_mtx_array_size(FILE *f, int M, int N);


#ifdef __cplusplus
     }
#endif

/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode)	((typecode)[0]=='M')

#define mm_is_sparse(typecode)	((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)	((typecode)[1]=='A')
#define mm_is_array(typecode)	((typecode)[1]=='A')

#define mm_is_complex(typecode)	((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)	((typecode)[3]=='G')
#define mm_is_skew(typecode)	((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

int mm_is_valid(MM_typecode matcode);		/* too complex for a macro */


/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode)	((*typecode)[0]='M')
#define mm_set_coordinate(typecode)	((*typecode)[1]='C')
#define mm_set_array(typecode)	((*typecode)[1]='A')
#define mm_set_dense(typecode)	mm_set_array(typecode)
#define mm_set_sparse(typecode)	mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)	((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')


#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)	((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
									(*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17


/******************** Matrix Market internal definitions ********************

   MM_matrix_typecode: 4-character sequence

				    ojbect 		sparse/   	data        storage 
						  		dense     	type        scheme

   string position:	 [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
						        A(array)	C(omplex)   H(ermitian)
											P(attern)   S(ymmetric)
								    		I(nteger)	K(kew)

 ***********************************************************************/

#define MM_MTX_STR		"matrix"
#define MM_ARRAY_STR	"array"
#define MM_DENSE_STR	"array"
#define MM_COORDINATE_STR "coordinate" 
#define MM_SPARSE_STR	"coordinate"
#define MM_COMPLEX_STR	"complex"
#define MM_REAL_STR		"real"
#define MM_INT_STR		"integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR		"symmetric"
#define MM_HERM_STR		"hermitian"
#define MM_SKEW_STR		"skew-symmetric"
#define MM_PATTERN_STR  "pattern"


/*  high level routines */
#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Write a coordinate matrix to a file.
 * @param fname File name.
 * @param M Number of rows.
 * @param N Number of columns.
 * @param nz Number of non-zero entries.
 * @param I Array of row indices.
 * @param J Array of column indices.
 * @param val Array of values (can be NULL for pattern matrices).
 * @param matcode Matrix Market typecode.
 * @return 0 on success, or an error code on failure.
 */
int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
		 double val[], MM_typecode matcode);

/**
 * @brief Read coordinate matrix data from a file into arrays.
 * @param f File pointer.
 * @param M Number of rows.
 * @param N Number of columns.
 * @param nz Number of non-zero entries.
 * @param I Array to store row indices.
 * @param J Array to store column indices.
 * @param val Array to store values.
 * @param matcode Matrix Market typecode.
 * @return 0 on success, or an error code on failure.
 */
int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
		double val[], MM_typecode matcode);

/**
 * @brief Read a single coordinate entry from a file.
 * @param f File pointer.
 * @param I Pointer to store row index.
 * @param J Pointer to store column index.
 * @param real Pointer to store real part of the value.
 * @param img Pointer to store imaginary part of the value.
 * @param matcode Matrix Market typecode.
 * @return 0 on success, or an error code on failure.
 */
int mm_read_mtx_crd_entry(FILE *f, int *I, int *J, double *real, double *img,
			MM_typecode matcode);

/**
 * @brief High-level routine to read an unsymmetric sparse matrix in coordinate format.
 * @param fname File name.
 * @param M_ Pointer to store number of rows.
 * @param N_ Pointer to store number of columns.
 * @param nz_ Pointer to store number of non-zero entries.
 * @param val_ Pointer to store allocated array of values.
 * @param I_ Pointer to store allocated array of row indices.
 * @param J_ Pointer to store allocated array of column indices.
 * @return 0 on success, or an error code on failure.
 */
int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
                double **val_, int **I_, int **J_);

#ifdef __cplusplus
     }
#endif


#endif
