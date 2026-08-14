#ifndef DIST_CC_H
#define DIST_CC_H

/************************************************************************ 
 * @file dist_cc.h
 * @brief Connected Cluster (CC) algorithm routines and helpers
 * 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>, Weilei Zeng
 ************************************************************************/

#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <m4ri/m4ri.h>

#include "mmio.h"
#include "uthash.h"
#include "util_hash.h"
#include "util_m4ri.h"
#include "util_io.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief v1[:] = v0[:] + mat[row,:] */
static inline int one_csr_row_combine(one_vec_t * const v1, const one_vec_t * const v0,
				  const csr_t * const mat, const int row){
#ifndef NDEBUG
  if ((!v1) || (!v0) || (!mat))
    ERROR("all arguments must be allocated: v1=%p v0=%p mat=%p\n",v1,v0,mat);
  if(v1 == v0)
    ERROR("the two vectors should not be the same !");
  if((row<0) || (row >= mat->rows))
    ERROR("this should not happen\n");
#endif
  int iM, i0=0, i1=0; /** iterators */
  for (iM = mat->p[row]; iM < mat->p[row+1]; iM++){
    const int ic = mat->i[iM];
    while((i0 < v0->wei) && (v0->vec[i0] < ic))
      v1->vec[i1++] = v0->vec[i0++];
    if(i0 >= v0->wei)
      break;
    if(v0->vec[i0]==ic)
      i0++; /** `1+1=0` just skip this position */
    else
      v1->vec[i1++] = ic;
  }
  if(i0 >= v0->wei) /** remaining `mat[row]` entries */
    for (                      ; iM < mat->p[row+1]; iM++){
      const int ic = mat->i[iM];
      v1->vec[i1++] = ic;
    }
  else /** remaining `v0` entries */
    while(i0 < v0->wei)
      v1->vec[i1++] = v0->vec[i0++];
  
  v1->wei = i1;
  return i1; /** weight of the out vector */
}

/** @brief insert `j` (originally absent) into ordered array, return position */
static inline int one_ordered_ins(one_vec_t * const err, const int j){
  int pos=err->wei-1;
  while(j < err->vec[pos]){
    err->vec[pos+1] = err->vec[pos];
    pos--;
  }
#ifndef NDEBUG  
  if (j == err->vec[pos]) 
    ERROR("Unexpected! vec[%d]=%d is already present!",pos,j);
#endif   
  err->vec[pos+1]=j;
#ifndef NDEBUG
  for(int i=0; i < err->wei; i++)
    if(err->vec[i] >= err->vec[i+1]){
      printf("check ordering at i=%d! ",i);
      one_vec_print(err);
      ERROR("unexpected");
    }
#endif   
  err->wei ++;
  return pos+1;
}

/** @brief find `val` in ordered array 
 *  @return position of `val` if found, -1 otherwise 
 */
static inline int one_ordered_search(one_vec_t * const err, const int val){
  /** binary search for pos of `val` */
  int bot=0, top=err->wei , mid=0;
#ifndef NDEBUG  
  if (!top)
    return -1;
#endif   
  while(top - bot > 1){
    mid = (top+bot) >> 1;
#ifndef NDEBUG  
  if (mid>=err->wei)
    ERROR("this should not happen");
#endif   
    if (err->vec[mid] <= val)
      bot = mid;
    else
      top = mid;
  }
  if ( err->vec[bot] == val)
    return bot;
  else
    return -1;
} 

/** @brief delete `val` (if originally present) from ordered array 
 *  @return 1 if `val` was found, 0 otherwise 
 */
static inline int one_ordered_find_del(one_vec_t * const err, const int val){
  /** binary search for pos of `j` */
  int bot=0, top=err->wei, mid=0;
#ifndef NDEBUG  
  if (!top)
    return 0;
#endif   
  while(top - bot > 1){
    mid = (top+bot) >> 1;
    if (err->vec[mid] <= val)
      bot = mid;
    else
      top = mid;
  }
  if ( err->vec[mid] != val)
    return 0;
  
  for(int i=mid; i < err->wei; i++)
      err->vec[i] = err->vec[i+1];
  err->wei --;
  return 1;
}

/** @brief delete `val` in known position `pos` from ordered array */
static inline void one_ordered_pos_del(one_vec_t * const err, _maybe_unused const int val, const int pos){
#ifndef NDEBUG
  if ((pos<0) || (pos >= err->wei) || (err->wei == 0) || (err->vec[pos] != val))
    ERROR("this should not happen!");
#endif
  err->wei --; 
  for(int i=pos; i < err->wei; i++)
      err->vec[i] = err->vec[i+1];
}

int start_CC_recurs(one_vec_t *err, one_vec_t *urr, one_vec_t * const syn[],
		    const int w_limit, const int max_col_wt, 
		    const csr_t * const mH, const csr_t * const mHT,
                    params_t * const p);

int do_CC_dist(params_t * const p);

#ifdef __cplusplus
}
#endif

#endif /* DIST_CC_H */
