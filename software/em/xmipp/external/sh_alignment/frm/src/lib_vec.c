/*********************************************************************
*                           L I B _ V E C                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Pablo Chacon and Willy Wriggers, 2001-2003                     *
**********************************************************************
*                                                                    *
* Creating, resetting, copying arrays and other data structures.     * 
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "situs.h"

/* function list, functions are more ore less self-explanatory  */

void zero_vect(double *, unsigned long);
void do_vect(double **, unsigned long);
void zero_mat(double **, unsigned long, unsigned long);
void do_mat(double ***, unsigned long, unsigned long);
void cp_vect(double **, double **, unsigned long);
void cp_vect_destroy(double **, double **, unsigned long);
void add_scaled_vect(double *, double *, double, unsigned long);


/*====================================================================*/
void zero_vect(double *vect, unsigned long len) {
  unsigned long i;
  for(i=0;i<len;++i) vect[i] = 0.0;
}

/*====================================================================*/
void do_vect(double **vect, unsigned long len) {
  *vect = (double *) malloc(len*sizeof(double));
  if (*vect == NULL) {
    fprintf(stderr, "lib_vec> Error: Unable to satisfy memory allocation request [e.c. 18010]\n"); 
    exit(18010);
  }
  zero_vect(*vect,len);
}

/*====================================================================*/
void zero_mat(double **mat,unsigned long len_i,unsigned long len_j) {
  unsigned long i;
  for(i=0;i<len_i;++i) zero_vect(mat[i],len_j);
}

/*====================================================================*/
void do_mat(double ***pmat,unsigned long len_i,unsigned long len_j) {
  unsigned long i;
  *pmat = (double **) malloc(len_i*sizeof(double *));
  if (*pmat == NULL) {
    fprintf(stderr, "lib_vec> Error: Unable to satisfy memory allocation request [e.c. 18020]\n"); 
    exit(18020);
  }
  for(i=0;i<len_i;i++) do_vect(&((*pmat)[i]),len_j);  
}

/*====================================================================*/
void cp_vect(double **vect1,double **vect2,unsigned long len) {
  memcpy(*vect1,*vect2, len*sizeof(double));
}

/*====================================================================*/
void cp_vect_destroy(double **vect1,double **vect2,unsigned long len) {
/* destroys memory allocated to vect2 after copying */  
  free(*vect1); 
  do_vect(vect1,len); 
  cp_vect(vect1,vect2,len); 
  free(*vect2);
}

/*====================================================================*/
void add_scaled_vect(double *to_vect, double *from_vect, double scalar, unsigned long len) { 
  unsigned long i;
  for(i=0;i<len;++i) to_vect[i] += scalar * from_vect[i];
}



