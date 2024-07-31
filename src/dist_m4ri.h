#ifndef DIST_M4RI
/************************************************************************ 
 * qLDPC code input utility routines for distance/decoder package               
 *                                                                      
 * currently: CSS only 
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 ************************************************************************/
#define DIST_M4RI



static const int maxrow = 10;

  /** @brief helper function to sort `int`
   *  use `qsort(array, len, sizeof(rci_t), cmp_rci_t);`
   */
  static inline int cmp_rci_t(const void *a, const void *b){
    const int va= *((int *) a);
    const int vb= *((int *) b);
    return va-vb;
  }


#endif /* DIST_M4RI */
