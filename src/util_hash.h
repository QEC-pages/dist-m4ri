#ifndef UTIL_HASH_H
#define UTIL_HASH_H
/**
 * @file util_hash.h
 *
 * @brief Utility functions for use with `uthash.h`
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2022 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 */
#ifdef __cplusplus
extern "C"{
#endif
#include "util_io.h"
#include "util_m4ri.h"

  
  /** hash storage helper functions *** use `uthash.h` *************************************/

  /**< @brief structure to hold two sparse vectors (syndrome,error) in a hash */
  typedef struct TWO_VEC_T {
    UT_hash_handle hh;
    int *err; /** error vector, pointer to arr+w_s */
    //    int cnt;  /** how many times this vector was encountered */
    int w_e; /**< weight of error (non-negative) */
    int w_s; /**< weight of syndrome (non-negative) */
    //  int w_tot; /** `= w_e+w_s` */
    //  size_t len; /** `weight*sizeof(int)` (is this really needed?) */
    int syn[0]; /** array of `w_e+w_s` integers, the actual key  */
  } two_vec_t;

  typedef struct ONE_VEC_T{
    int wei; /** current weight */
    //    int max; /** allocated */
    int vec[0];
  } one_vec_t;


  /** @brief print entire `one_vec_t` structure by pointer */
  void one_vec_print(const one_vec_t * const pvec);
  
  static inline two_vec_t * two_vec_init(const one_vec_t * const syn, const one_vec_t * const err){
    two_vec_t * ans = malloc(sizeof(two_vec_t)+sizeof(int)*(syn->wei + err->wei));
    if(!ans)
      ERROR("memory allocation fail");
    for(int i=0; i < syn->wei; i++)
      ans->syn[i] = syn->vec[i];
    ans->w_e = err->wei;
    ans->w_s = syn->wei;
    ans->err = &( ans->syn[syn->wei] );
    for(int j=0 ; j < err->wei; j++)
      ans->err[j] = err->vec[j];
    
    return ans;
  } 
  
  static inline void two_vec_print(const two_vec_t * const it){
    if(!it){
      printf("two_vec_print(): null structure!\n");
      return;
    }
    printf("# two_vec: w_e=%d e=[",it->w_e);
    for(int i=0; i < it->w_e; i++)
      printf(" %d",it->err[i]);
    printf(" ]%s w_s=%d s=[",it->w_s>20?"\n#":" ",it->w_s);
    for(int i=0; i < it->w_s; i++)
      printf(" %d",it->syn[i]);
    printf(" ]\n");
  }

  /** @brief compare two `two_vec_t` structures by syndrome */
static inline int by_syndrome(void *a, void *b){
  const two_vec_t * const pa = (two_vec_t *) a;
  const two_vec_t * const pb = (two_vec_t *) b;
  if (pa->w_s < pb->w_s)
    return -1;
  else if (pa->w_s > pb->w_s)
    return +1;
  else{ /** Wa == Wb */
    for(int i=0; i < pa->w_s; i++){
      if (pa->syn[i] < pb->syn[i])
	return -1;
      else if (pa->syn[i] > pb->syn[i])
	return +1;
    }
  }
  return 0;
}

  /** @brief compare two `two_vec_t` structures by error vectors */
  static inline int by_error(void *a, void *b){
    const two_vec_t * const pa = (two_vec_t *) a;
    const two_vec_t * const pb = (two_vec_t *) b;
    if (pa->w_e < pb->w_e)
      return -1;
    else if (pa->w_e > pb->w_e)
      return +1;
    else{ /** Wa == Wb */
      for(int i=0; i < pa->w_e; i++){
	if (pa->err[i] < pb->err[i])
	  return -1;
	else if (pa->err[i] > pb->err[i])
	  return +1;
      }
    }
    return 0;
  }

  /** @brief search hash `errors` for `syn` and, if not found, insert `syn`,`err` 
   * @param syn input syndrome, must be ordered
   * @param err input error vector, must be ordered
   * @param errors input/output hash
   * @param p_swei[] array to store error weights 
   * @return errors (updated hash)
   */
  static inline two_vec_t * hash_add_maybe(one_vec_t * syn, const one_vec_t * const err,
			     two_vec_t * errors, int p_swei[],
					   __attribute__ ((unused)) const int debug){
    two_vec_t *pvec, *entry;
    const size_t keylen = syn->wei * sizeof(int);
    HASH_FIND(hh, errors, syn->vec, keylen, pvec);
    if(!pvec){ /** syndrome not found, inserting */
      entry = two_vec_init(syn, err);		
      HASH_ADD(hh, errors, syn, keylen, entry); /** add to `hash` */
      if(p_swei[err->wei] > syn->wei){
	//#ifndef NDEBUG
	if(debug&64){
	  printf("# swei[%d]=%d -> %d change\n# err: ",
		 err->wei,p_swei[err->wei],syn->wei);
	  one_vec_print(err);
	  printf("# syn: ");
	  one_vec_print(syn);
	}
	//#endif
	p_swei[err->wei]=syn->wei;      
      }
    }
    else{  /** we construct small-weight vectors first, thus should not
	       worry about replacing error vectors already in the hash */
    //#ifndef NDEBUG  
      if(debug&128){
	printf("err: ");
	one_vec_print(err);
	printf("syn: ");
	one_vec_print(syn);
	printf("already in the hash:\n");
	two_vec_print(pvec);
      }
    //#endif   
    }
    return errors;
  }
  

#ifdef __cplusplus
}
#endif

#endif /* UTIL_HASH_H */
  
