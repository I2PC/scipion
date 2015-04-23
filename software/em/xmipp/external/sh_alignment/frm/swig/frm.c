/***************************************************************************
  **************************************************************************

                     Python interface of Fast Rotational Matching

   AUTHOR:
   Yuxiang Chen 
   Contact: chenyxkk@hotmail.com 

   Copyright 2011 Yuxiang Chen

   HOW TO USE:
   1. run mkswig.sh, which should generate shared library called _swig_frm.so
   2. start python, import swig_frm and enjoy!

  ************************************************************************
  ************************************************************************/

#include "situs.h"
#include "fftw3.h"
#define SQT2 sqrt(2.0)

#include "cospmls.h"
#include "FST_semi_memo.h"

/* external general library funtions */
extern void do_vect (double **, unsigned long);

/* functions defined in this file */
static void wigner(int, double, double *);
static int prime(int, int);
static double distance_eulers(unsigned long, unsigned long, double *[]);

static void fourier_corr(int, double *, double *, double *, double *, double *, double *);

/* Given two spherical function f and g, return the best correlation values and Euler angles. */
int frm(double *f, int dim1, double *g, int dim2, double *res, int dim3)
{
  if(f==NULL || g==NULL || res==NULL)
  {
	 fprintf(stderr, "Either no input data given or forget to allocate the memory for output!");
    return -1;
  }
  
  int bw = sqrt(dim1)/2;	/* bandwidth for sampling on the sphere */
  if(dim1!=dim2 || (dim3%4)!=0)
  {
	 fprintf(stderr,"Inputs' dimensions are wrong!");
	 return -1;
  }

  long i, j, k, l, m,  ind1, ind2, ind3, size2;
  long ncorr, npeaks;
  

  int log2s;       /* log2 of size */
  int size;        /* no. of points in each coordinate = 2*bw */
  int cut_naive;   /* switch from semi-naive to naive at order bw - cut_naive */

  int ifit;

  fftw_plan plan_fftw;     /* plan for the FFTW routine */
  
  double temp, tempi, tempj, tempk;
  double f000, f100, f200, f010, f020, f001, f002;

  double phi, th, psi;        /* Euler angles */
  double phe, phc, ohm;       /* angular parameters */
  
  double dist_cut = 3.0;            /* cutoff distance for eliminating nearby peaks */

  double *seminaive_naive_tablespace, *workspace;
  double **seminaive_naive_table;

  double *idata;      /* imaginary data (0) */

  double *coeffr_hi, *coeffi_hi;  /* spherical harmonics coefficients for high res. map */
  double *coeffr_lo, *coeffi_lo;  /* spherical harmonics coefficients for low res. map */

  double *ddd;        /* piramid array that will contain the d-functions */

  double *coef_corr;  /* Fourier coefficients of correlation function */

  double *corr[4];    /* correlation peaks and corresponding Euler angles */



 /* ==================== FRM PRECOMPUTATIONS =======================*/  
//  fprintf(stderr, "frmr> Doing precomputations for FST and d-matrices...\n");

  size=2*bw;
  size2=size*size;
  cut_naive=bw;

  do_vect(&seminaive_naive_tablespace, Reduced_Naive_TableSize(bw,cut_naive) +
               Reduced_SpharmonicTableSize(bw,cut_naive));

  do_vect(&workspace,(8*bw*bw+29*bw)); // different with the remmendation 29, the larger the better?

  seminaive_naive_table = (double **)SemiNaive_Naive_Pml_Table(bw, cut_naive,
                            seminaive_naive_tablespace,
                            workspace);
  
  do_vect(&ddd, bw*(4*bw*bw-1)/3);
  wigner(bw, 0.5*M_PI, ddd);    

 /* ==================== FRM COMPUTATIONS =======================*/  

  do_vect(&idata,size*size); /* set imagnary part to 0*/

  
  /* compute coeff for hi */
//  fprintf(stderr, "frmr> Computing coefficients for high resolution map...\n");
  do_vect(&coeffr_hi,bw*bw);
  do_vect(&coeffi_hi,bw*bw);
  FST_semi_memo(f, idata,
          coeffr_hi, coeffi_hi,
          size,
          seminaive_naive_table,
          workspace,
          1, /* 1 - real sample; 0 - complex sample */
          cut_naive);
//  for(l=0;l<bw;l++)              /* correction factor for m=0 coefficients. Dont need anymore! */
//    coeffr_hi[l] *= SQT2;

  
  /* compute coeff for lo*/
  do_vect(&coeffr_lo,bw*bw);
  do_vect(&coeffi_lo,bw*bw); 
  FST_semi_memo(g, idata,
          coeffr_lo, coeffi_lo,
          size,
          seminaive_naive_table,
          workspace,
          1,
          cut_naive);
//  for(l=0;l<bw;l++)              /* correction factor for m=0 coefficients. Dont need anymore! */
//    coeffr_lo[l] *= SQT2;


  /* compute FT of correlation function */
//  fprintf(stderr, "frmr> Computing Fourier coefficients of the correlation function...\n");
  do_vect(&coef_corr,2*size*size2);
  fourier_corr(bw, coeffr_hi, coeffi_hi, coeffr_lo, coeffi_lo, ddd, coef_corr);  /* compute T(p,q,r) */

  /* free space */
  free(coeffr_hi); free(coeffi_hi);
  free(coeffr_lo); free(coeffi_lo);
  free(workspace); free(ddd);
  free(seminaive_naive_tablespace);
  free(seminaive_naive_table);
  free(idata);

  
  /* compute the correlation function using FFTW */
//  fprintf(stderr, "frmr> Computing inverse FFT...\n");
  /* FFTW initialization */
  plan_fftw=fftw_plan_dft_3d(size, size, size, (fftw_complex *)coef_corr, (fftw_complex *)coef_corr, FFTW_BACKWARD,
                               FFTW_ESTIMATE);
  fftw_execute(plan_fftw);
  fftw_destroy_plan(plan_fftw);


 /*=========================== DETECT AND SORT PEAKS  ============================*/ 

 /* this block searches the relative maxima of the correlation function */
 /* and stores them in the array corr */
 /* j ranges over half the values due to periodicity */
 /* the correlation value is actually interpolated! */
  
 for(i=0;i<4;i++){
    corr[i]  = (double *) malloc(1 *sizeof(double));
    if(corr[i] == NULL){
      perror("Error in allocating memory -1\n");
      return -1;
    }
  }

// fprintf(stderr, "frmr> Searching peaks...\n");

 npeaks=0;                      
 for(i=0;i<size;i++){           
   ind1=i*size2;
   for(j=0;j<=bw;j++){         
     ind2=ind1+j*size;
     for(k=0;k<size;k++){
       ind3=ind2+k;
        
       if ((f020=coef_corr[2*(ind3+size)])<=(f000=coef_corr[2*ind3])) {       
     if (i==0) f100=coef_corr[2*(ind3+(size-1)*size2)];
     else f100=coef_corr[2*(ind3-size2)];
         if (f100<=f000) {         
       if (i==size-1) f200=coef_corr[2*(j*size+k)];
       else f200=coef_corr[2*(ind3+size2)];
       if (f200<=f000) {         
         if (j==0) f010=coef_corr[2*(ind3+size2-size)];
         else f010=coef_corr[2*(ind3-size)];
         if (f010<=f000) {
           if(k==0) f001=coef_corr[2*(ind3+size-1)];
           else f001=coef_corr[2*(ind3-1)];
           if (f001<=f000) {
         if(k==size-1) f002=coef_corr[2*ind2];
         else f002=coef_corr[2*(ind3+1)];
         if (f002<=f000) {
     
           tempi=0.0; tempj=0.0; tempk=0.0;
      
           if((temp=f100+f200-2*f000)!=0)
             tempi=(f100-f200)/(2*temp);         
           phe=(i+tempi)*M_PI/bw;

           if((temp=f010+f020-2*f000)!=0)
             tempj=(f010-f020)/(2*temp);
           phc=(j+tempj)*M_PI/bw;

           if((temp=f001+f002-2*f000)!=0)
             tempk=(f001-f002)/(2*temp);
           ohm=(k+tempk)*M_PI/bw;
           
          
           /* reallocate memory */   
                   
           for(l=0;l<4;l++)
             corr[l]=(double *) realloc(corr[l],(npeaks+1)*sizeof(double));
           phi=M_PI-ohm; phi -= 2.0*M_PI*floor(phi/2.0/M_PI);
           th=M_PI-phc;
           psi=-phe;     psi -= 2.0*M_PI*floor(psi/2.0/M_PI);
           
               corr[0][npeaks]=(0.5*(f100+f200)-f000)*tempi*tempi+(0.5*(f010+f020)-f000)*tempj*tempj+
             (0.5*(f001+f002)-f000)*tempk*tempk+0.5*(f200-f100)*tempi+0.5*(f020-f010)*tempj+
             0.5*(f002-f001)*tempk+f000;;
           corr[1][npeaks]=psi;
           corr[2][npeaks]=th;
           corr[3][npeaks]=phi;
           npeaks++;
         }
           }
         }
       }
     }
       }
     } 
   }
 }
 free(coef_corr);



 /* sort correlation values in decreasing order */
 for(i=0;i<npeaks;i++)      
   for(j=i+1;j<npeaks;j++)
     if(corr[0][i]<corr[0][j])
       for(k=0;k<4;k++)  SWAPPING( corr[k][i], corr[k][j], double );


 /* eliminate close-by solutions */
 for(i=0;i<npeaks;i++)      
   for(j=i+1;j<npeaks;j++)
     if(distance_eulers(i, j, corr)<dist_cut*M_PI/bw) { 
       npeaks--;
       for(k=0;k<4;k++) corr[k][j]=corr[k][npeaks];       
       j--; 
     }


  /* sort again correlation values in decreasing order */
  for(i=0;i<npeaks;i++)      
    for(j=i+1;j<npeaks;j++)
      if(corr[0][i]<corr[0][j])
        for(k=0;k<4;k++)  SWAPPING(corr[k][i], corr[k][j],double );
 

  /* The result here is already zxz convention. Transfer to degrees! */
  int max_peaks = dim3/4; /* maximal return peaks */
  for(i=0; i<npeaks && i<max_peaks; i++)
  {
    res[i*4] = corr[0][i]; /* correlation value*/
    res[i*4+1] = corr[1][i]*180/M_PI; /* psi in deg */
    res[i*4+2] = corr[2][i]*180/M_PI; /* the in deg */
    res[i*4+3] = corr[3][i]*180/M_PI; /* phi in deg */
  }


  for(i=0;i<4;i++)
    free(corr[i]);  
/*
  free(corr);
*/

  return 0;
}


/* Direct return correlation volume */
int frm_corr(double *f, int dim1, double *g, int dim2, double *coef_corr, int dim3)
{
  if(f==NULL || g==NULL || coef_corr==NULL)
  {
	fprintf(stderr, "Either no input data given or forget to allocate the memory for output!");
    return -1;
  }
  
  int bw = sqrt(dim1)/2;	/* bandwidth for sampling on the sphere */
  if(dim1!=dim2 || dim1!=4*bw*bw || dim3!=16*bw*bw*bw || !(bw && !(bw & (bw - 1))))
  {
	fprintf(stderr,"Inputs' dimensions are wrong!");
	return -1;
  }

  long i, j, k, l, m,  ind1, ind2, ind3, size2;
  long ncorr, npeaks;
  

  int log2s;       /* log2 of size */
  int size;        /* no. of points in each coordinate = 2*bw */
  int cut_naive;   /* switch from semi-naive to naive at order bw - cut_naive */

  int ifit;

  fftw_plan plan_fftw;     /* plan for the FFTW routine */
  
  double temp, tempi, tempj, tempk;
  double f000, f100, f200, f010, f020, f001, f002;

  double phi, th, psi;        /* Euler angles */
  double phe, phc, ohm;       /* angular parameters */
  
  // double dist_cut = 3.0;            /* cutoff distance for eliminating nearby peaks */

  double *seminaive_naive_tablespace, *workspace;
  double **seminaive_naive_table;

  double *idata;      /* imaginary data (0) */

  double *coeffr_hi, *coeffi_hi;  /* spherical harmonics coefficients for high res. map */
  double *coeffr_lo, *coeffi_lo;  /* spherical harmonics coefficients for low res. map */

  double *ddd;        /* piramid array that will contain the d-functions */



 /* ==================== FRM PRECOMPUTATIONS =======================*/  
//  fprintf(stderr, "frmr> Doing precomputations for FST and d-matrices...\n");

  size=2*bw;
  size2=size*size;
  cut_naive=bw;

  do_vect(&seminaive_naive_tablespace, Reduced_Naive_TableSize(bw,cut_naive) +
               Reduced_SpharmonicTableSize(bw,cut_naive));

  do_vect(&workspace,(8*bw*bw+29*bw));

  seminaive_naive_table = (double **)SemiNaive_Naive_Pml_Table(bw, cut_naive,
                            seminaive_naive_tablespace,
                            workspace);
  
  do_vect(&ddd, bw*(4*bw*bw-1)/3);
  wigner(bw, 0.5*M_PI, ddd);    

 /* ==================== FRM COMPUTATIONS =======================*/  

  do_vect(&idata,size*size); /* set imagnary part to 0*/

  
  /* compute coeff for hi */
//  fprintf(stderr, "frmr> Computing coefficients for high resolution map...\n");
  do_vect(&coeffr_hi,bw*bw);
  do_vect(&coeffi_hi,bw*bw);
  FST_semi_memo(f, idata,
          coeffr_hi, coeffi_hi,
          size,
          seminaive_naive_table,
          workspace,
          1,
          cut_naive);
  // for(l=0;l<bw;l++)              /* correction factor for m=0 coefficients */
  //   coeffr_hi[l] *= SQT2;

  
  /* compute coeff for lo*/
  do_vect(&coeffr_lo,bw*bw);
  do_vect(&coeffi_lo,bw*bw); 
  FST_semi_memo(g, idata,
          coeffr_lo, coeffi_lo,
          size,
          seminaive_naive_table,
          workspace,
          1,
          cut_naive);
  // for(l=0;l<bw;l++)              /* correction factor for m=0 coefficients */
  //   coeffr_lo[l] *= SQT2;


  /* compute FT of correlation function */
//  fprintf(stderr, "frmr> Computing Fourier coefficients of the correlation function...\n");
//  do_vect(&coef_corr,2*size*size2);
  fourier_corr(bw, coeffr_hi, coeffi_hi, coeffr_lo, coeffi_lo, ddd, coef_corr);  /* compute T(p,q,r) */

  /* free space */
  free(coeffr_hi); free(coeffi_hi);
  free(coeffr_lo); free(coeffi_lo);
  free(workspace); free(ddd);
  free(seminaive_naive_tablespace);
  free(seminaive_naive_table);
  free(idata);

  
  /* compute the correlation function using FFTW */
//  fprintf(stderr, "frmr> Computing inverse FFT...\n");
  plan_fftw=fftw_plan_dft_3d(size, size, size, (fftw_complex *)coef_corr, (fftw_complex *)coef_corr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan_fftw);
  fftw_destroy_plan(plan_fftw);

  /* copy the real part for return */
  return 0;
}


/* Calculate the FRM in Fourier space */
int frm_fourier_corr(double *fr, int dim1, double *fi, int dim2, double *gr, int dim3, double *gi, int dim4, double *coef_corr, int dim5)
{
  if(fr==NULL || fi==NULL || gr==NULL || gi==NULL || coef_corr==NULL)
  {
	fprintf(stderr, "Either no input data given or forget to allocate the memory for output!");
    return -1;
  }
  
  int bw = sqrt(dim1)/2;	/* bandwidth for sampling on the sphere */
  if(dim1!=dim2 || dim3!=dim4 || dim1!=dim3 || dim5!=16*bw*bw*bw)
  {
	fprintf(stderr,"Inputs' dimensions are wrong!");
	return -1;
  }

  long i, j, k, l, m,  ind1, ind2, ind3, size2;
  long ncorr, npeaks;
  

  int log2s;       /* log2 of size */
  int size;        /* no. of points in each coordinate = 2*bw */
  int cut_naive;   /* switch from semi-naive to naive at order bw - cut_naive */

  int ifit;

  fftw_plan plan_fftw;     /* plan for the FFTW routine */
  
  double temp, tempi, tempj, tempk;
  double f000, f100, f200, f010, f020, f001, f002;

  double phi, th, psi;        /* Euler angles */
  double phe, phc, ohm;       /* angular parameters */
  
  // double dist_cut = 3.0;            /* cutoff distance for eliminating nearby peaks */

  double *seminaive_naive_tablespace, *workspace;
  double **seminaive_naive_table;

//  double *idata;      /* imaginary data (0) */

  double *coeffr_hi, *coeffi_hi;  /* spherical harmonics coefficients for high res. map */
  double *coeffr_lo, *coeffi_lo;  /* spherical harmonics coefficients for low res. map */

  double *ddd;        /* piramid array that will contain the d-functions */



 /* ==================== FRM PRECOMPUTATIONS =======================*/  
//  fprintf(stderr, "frmr> Doing precomputations for FST and d-matrices...\n");

  size=2*bw;
  size2=size*size;
  cut_naive=bw;

  do_vect(&seminaive_naive_tablespace, Reduced_Naive_TableSize(bw,cut_naive) +
               Reduced_SpharmonicTableSize(bw,cut_naive));

  do_vect(&workspace,(8*bw*bw+29*bw));

  seminaive_naive_table = (double **)SemiNaive_Naive_Pml_Table(bw, cut_naive,
                            seminaive_naive_tablespace,
                            workspace);
  
  do_vect(&ddd, bw*(4*bw*bw-1)/3);
  wigner(bw, 0.5*M_PI, ddd);    

 /* ==================== FRM COMPUTATIONS =======================*/  

//  do_vect(&idata,size*size); /* set imagnary part to 0*/

  
  /* compute coeff for hi */
//  fprintf(stderr, "frmr> Computing coefficients for high resolution map...\n");
  do_vect(&coeffr_hi,bw*bw);
  do_vect(&coeffi_hi,bw*bw);
  FST_semi_memo(fr, fi,
          coeffr_hi, coeffi_hi,
          size,
          seminaive_naive_table,
          workspace,
          0,  // complex numbers!!!
          cut_naive);
  // for(l=0;l<bw;l++)              /* correction factor for m=0 coefficients */
  //   coeffr_hi[l] *= SQT2;
  
  /* compute coeff for lo*/
  do_vect(&coeffr_lo,bw*bw);
  do_vect(&coeffi_lo,bw*bw); 
  FST_semi_memo(gr, gi,
          coeffr_lo, coeffi_lo,
          size,
          seminaive_naive_table,
          workspace,
          0,  // complex numbers!!!
          cut_naive);
  // for(l=0;l<bw;l++)              /* correction factor for m=0 coefficients */
  //   coeffr_lo[l] *= SQT2;


  /* compute FT of correlation function */
//  fprintf(stderr, "frmr> Computing Fourier coefficients of the correlation function...\n");
//  do_vect(&coef_corr,2*size*size2);
  fourier_corr(bw, coeffr_hi, coeffi_hi, coeffr_lo, coeffi_lo, ddd, coef_corr);  /* compute T(p,q,r) */

  /* free space */
  free(coeffr_hi); free(coeffi_hi);
  free(coeffr_lo); free(coeffi_lo);
  free(workspace); free(ddd);
  free(seminaive_naive_tablespace);
  free(seminaive_naive_table);
//  free(idata);

  
  /* compute the correlation function using FFTW */
//  fprintf(stderr, "frmr> Computing inverse FFT...\n");
  /* FFTW initialization */
  plan_fftw=fftw_plan_dft_3d(size, size, size, (fftw_complex *)coef_corr, (fftw_complex *)coef_corr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan_fftw);
  fftw_destroy_plan(plan_fftw);

  return 0;
}

/* Do the element-wise multiplication of FST coefficients from two function. */
static void fourier_corr(int bw, double *coeffr_hi, double *coeffi_hi,
                  double *coeffr_lo, double *coeffi_lo, double *ddd, double *coef_corr){
/* computes the Fourier coefficients of the correlation function */
/* transposing the blocks of data */
  
  int l, p, q, r, tr, zi, zii, size, tsize, isigp, isigr;
  double temp_r, temp_i, tempd;
  unsigned long index, indexn, init, ind1, ind3, indexa, indexan, indexb, indexbn,
                indexmod, indexnmod, size2, tsize2;

  size=2*bw; tsize=2*size;
  size2=size*size;
  tsize2=2*size2;

  for(l=0;l<bw;l++){                /* compute T(p,q,r) for p>=0 and q>=0 */
    init=l*(4*l*l-1)/3;
    for(p=0;p<=l;p++){
      ind1=seanindex(p,l,bw);
      indexa=p*size2;
      for(r=-l;r<0;r++){            /* r<0 */
        ind3=seanindex(r,l,bw);
        indexb=indexa+r+size;
        temp_r = coeffr_lo[ind1]*coeffr_hi[ind3]+coeffi_lo[ind1]*coeffi_hi[ind3];
        temp_i =-coeffr_lo[ind1]*coeffi_hi[ind3]+coeffi_lo[ind1]*coeffr_hi[ind3];
        for(q=0;q<=l;q++){
          tempd = ddd[init+(p+l)*(2*l+1)+(q+l)] * ddd[init+(q+l)*(2*l+1)+(r+l)];
          index=2*(indexb+q*size);
          *(coef_corr + index)     += temp_r * tempd;
          *(coef_corr + index + 1) += temp_i * tempd;
        }
      }
      for(r=0;r<=l;r++){            /* r>=0 */
        ind3=seanindex(r,l,bw);
        indexb=indexa+r;
        temp_r = coeffr_lo[ind1]*coeffr_hi[ind3]+coeffi_lo[ind1]*coeffi_hi[ind3];
        temp_i =-coeffr_lo[ind1]*coeffi_hi[ind3]+coeffi_lo[ind1]*coeffr_hi[ind3];
        for(q=0;q<=l;q++){
          tempd = ddd[init+(p+l)*(2*l+1)+(q+l)] * ddd[init+(q+l)*(2*l+1)+(r+l)];
          index=2*(indexb+q*size);
          *(coef_corr + index)     += temp_r * tempd;
          *(coef_corr + index + 1) += temp_i * tempd;
        }
      }
    }
  }

  isigp=-1;
  for(p=0,indexa=0; p<bw; p++,indexa+=size2){   /* fill in the region p>=0,q<0 of T(p,q,r) */
    isigp = -isigp;                     /* parity of p */
    zi=(indexa+bw+size2)<<1;
    zii=zi-tsize2;
    for(q=-bw+1;q<0;q++){
      isigr=1;
      index=indexmod=zi+q*tsize;
      indexn=indexnmod=zii-q*tsize;
      for(tr=2;tr<size;tr+=2){          /* for r=-(bw-1),...,-1  ("tr=-2r") */
        isigr=-isigr;                   /* parity of r */
        index+=2;
        indexn+=2;
        /* this is coef_corr[indexn]*(-1)^(p+r) */
        *(coef_corr + index)     = isigp * isigr * *(coef_corr + indexn);
        *(coef_corr + index + 1) = isigp * isigr * *(coef_corr + indexn + 1);
      }
      index=indexmod-size;
      indexn=indexnmod-size;
      for(r=0;r<bw;r++){
        isigr=-isigr;                   /* parity of r */
        /* this is coef_corr[indexn]*(-1)^(p+r) */
        *(coef_corr + index)     = isigp * isigr * *(coef_corr + indexn);
        *(coef_corr + index + 1) = isigp * isigr * *(coef_corr + indexn + 1);
        index+=2;
        indexn+=2;
      }
    }
  }

  for(p=-bw+1;p<0;p++){                 /* fill in the region p<0 of T(p,q,r) */
    indexa=(p+size)*size2;
    indexan=-p*size2;
    for(q=-bw+1;q<bw;q++){
      indexb=indexa+prime(q, size)*size;
      indexbn=indexan+prime(-q, size)*size;
      for(r=-bw+1;r<bw;r++){
        index=(indexb+prime(r, size))<<1;
        indexn=(indexbn+prime(-r, size))<<1;
        *(coef_corr + index)     =  *(coef_corr + indexn);
        *(coef_corr + index + 1) = -*(coef_corr + indexn + 1);
      }
    }
  }
}





/*====================================================================*/
static void wigner(int bw, double theta, double *ddd){
  /* computes the d-functions for argument theta, for degrees 0 through bw-1 */
  /* the matrices are returned in the (3D) piramid ddd */
  /* Reference: T. Risbo, Journal of Geodesy (1996) 70:383-396 */

  double p, q, pc, qc, temp, *d, *dd;
  int size, index, i, j2, k, l, hdeg;
  unsigned long init;
  double fact1, fact2, fact3, fact4;

  size=2*bw;
  
  d  = (double *) malloc(size*size*sizeof(double));
  dd = (double *) malloc(size*size*sizeof(double));

  if( (d == NULL) || (dd == NULL) ){
      perror("Error in allocating memory");
      exit(1);
  }

  p=sin(theta/2); q=cos(theta/2);   /* Cayley-Klein parameters */
  pc=p; qc=q;

  *(ddd+0)=1.0;                     /* d-matrix of degree 0 */
  
  /* d-matrix of degree 1/2 */
  *(d+0)=q;
  *(d+1)=p;
  *(d+1*size+0)=-pc;
  *(d+1*size+1)=qc;
  
  for(l=1;l<bw;l++){       /* l is the degree of the next d-matrix to be saved in ddd */
    j2=(l<<1)-1;           /* j2 will be twice the degree of the d-matrix about to be computed */
    for(hdeg=0;hdeg<2;hdeg++){
      fact1=q/ ++j2;
      fact2=pc/j2;
      fact3=p/j2;
      fact4=qc/j2;
      for(i=0;i<=j2+1;i++)
        for(k=0;k<=j2+1;k++)
          *(dd+i*size+k)=0.0;
      for(i=0;i<j2;i++)
        for(k=0;k<j2;k++){
          index=i*size+k;
          temp=*(d+index);
          *(dd+index)        += sqrt((j2-i)*(j2-k))*temp*fact1;
      *(dd+index+size)   += -sqrt((i+1)*(j2-k))*temp*fact2;
      *(dd+index+1)      += sqrt((j2-i)*(k+1)) *temp*fact3;
      *(dd+index+size+1) += sqrt((i+1) *(k+1)) *temp*fact4;
        }
      for(i=0;i<=j2;i++)        /* move dd to d */
        for(k=0;k<=j2;k++)
          *(d+i*size+k)=*(dd+i*size+k);
      if(hdeg==0){              /* if degree is integer, copy d to the proper location in ddd */
        init=l*(4*l*l-1)/3;
        for(i=0;i<=j2;i++)
          for(k=0;k<=j2;k++)
            *(ddd+init+i*(j2+1)+k) = *(d+i*size+k);
      }
      if(l==bw-1) break;     /* for l=bw-1 do hdeg=0 only */
    }
  }
  free(d); free(dd);
}


/*====================================================================*/
static int prime(int n, int s){
  /* returns n if it is >=0, or n+s otherwise */

  if(n>=0) return n;
    return n+s;
}

/*====================================================================*/
static double distance_eulers(unsigned long i, unsigned long j, double *corr[]){
  /* returns the distance between the rotations given by the Euler angles in
     corr[1,2,3][i] and corr[1,2,3][j] */

  extern void get_rot_matrix ();
  double matrix1[3][3],matrix2[3][3]; 
  int k, l;
  double temp, trace;
  
  get_rot_matrix(matrix1,corr[1][i],corr[2][i],corr[3][i]);
  get_rot_matrix(matrix2,corr[1][j],corr[2][j],corr[3][j]);
  

  trace=0.0;             /* trace(matrix1*transpose(matrix2)) */
  for(k=0;k<3;k++)
    for(l=0;l<3;l++)
      trace += matrix1[k][l]*matrix2[k][l];

  temp=0.5*(trace-1.0);
  if(temp>=1.0)
    return 0.0;
  if(temp<=-1.0)
    return M_PI;
  return acos(temp);
}

double angle_distance(double a1, double a2, double a3, double b1, double b2, double b3)
/* psi, the, phi order */
{
  extern void get_rot_matrix ();
  double matrix1[3][3],matrix2[3][3]; 
  int k, l;
  double temp, trace;
  
  get_rot_matrix(matrix1,a1,a2,a3);
  get_rot_matrix(matrix2,b1,b2,b3);
  

  trace=0.0;             /* trace(matrix1*transpose(matrix2)) */
  for(k=0;k<3;k++)
    for(l=0;l<3;l++)
      trace += matrix1[k][l]*matrix2[k][l];

  temp=0.5*(trace-1.0);
  if(temp>=1.0)
    return 0.0;
  if(temp<=-1.0)
    return M_PI;
  return acos(temp);
}

/*====================================================================*/
int find_topn_angles(double *coef_corr, int dim1, int bw, double *res, int dim2, double dist_cut)
{
    /*=========================== DETECT AND SORT PEAKS  ============================*/ 

    /* this block searches the relative maxima of the correlation function */
    /* and stores them in the array corr */
    /* j ranges over half the values due to periodicity */
    /* the correlation value is actually interpolated! */

    if(8*bw*bw*bw != dim1)
    {
      fprintf(stderr,"Inputs' dimensions are wrong!");
      return -1;
    }

    double *corr[4];
    long i, j, k, l, m,  ind1, ind2, ind3;
    long npeaks;
    int size=2*bw;        /* no. of points in each coordinate = 2*bw */
    long size2=size*size;
    
    double temp, tempi, tempj, tempk;
    double f000, f100, f200, f010, f020, f001, f002;

    double phi, th, psi;        /* Euler angles */
    double phe, phc, ohm;       /* angular parameters */

    for(i=0;i<4;i++){
      corr[i]  = (double *) malloc(1 *sizeof(double));
      if(corr[i] == NULL){
        perror("Error in allocating memory -1\n");
        return -1;
      }
    }

   npeaks=0;                      
   for(i=0;i<size;i++){           
     ind1=i*size2;
     for(j=0;j<=bw;j++){         
       ind2=ind1+j*size;
       for(k=0;k<size;k++){
         ind3=ind2+k;
          
         if ((f020=coef_corr[ind3+size])<=(f000=coef_corr[ind3])) {       
           if (i==0) f100=coef_corr[ind3+(size-1)*size2];
           else f100=coef_corr[ind3-size2];
           if (f100<=f000) {         
             if (i==size-1) f200=coef_corr[j*size+k];
             else f200=coef_corr[ind3+size2];
             if (f200<=f000) {         
               if (j==0) f010=coef_corr[ind3+size2-size];
               else f010=coef_corr[ind3-size];
               if (f010<=f000) {
                 if(k==0) f001=coef_corr[ind3+size-1];
                 else f001=coef_corr[ind3-1];
                 if (f001<=f000) {
                   if(k==size-1) f002=coef_corr[ind2];
                   else f002=coef_corr[ind3+1];
                   if (f002<=f000) {
       
             tempi=0.0; tempj=0.0; tempk=0.0;
        
             if((temp=f100+f200-2*f000)!=0)
               tempi=(f100-f200)/(2*temp);         
             phe=(i+tempi)*M_PI/bw;

             if((temp=f010+f020-2*f000)!=0)
               tempj=(f010-f020)/(2*temp);
             phc=(j+tempj)*M_PI/bw;

             if((temp=f001+f002-2*f000)!=0)
               tempk=(f001-f002)/(2*temp);
             ohm=(k+tempk)*M_PI/bw;
            
             /* reallocate memory */   
             for(l=0;l<4;l++)
               corr[l]=(double *) realloc(corr[l],(npeaks+1)*sizeof(double));

             phi=M_PI-ohm; phi -= 2.0*M_PI*floor(phi/2.0/M_PI);
             th=M_PI-phc;
             psi=-phe;     psi -= 2.0*M_PI*floor(psi/2.0/M_PI);
             
             corr[0][npeaks]=(0.5*(f100+f200)-f000)*tempi*tempi+(0.5*(f010+f020)-f000)*tempj*tempj+
               (0.5*(f001+f002)-f000)*tempk*tempk+0.5*(f200-f100)*tempi+0.5*(f020-f010)*tempj+
               0.5*(f002-f001)*tempk+f000;
             corr[1][npeaks]=psi;
             corr[2][npeaks]=th;
             corr[3][npeaks]=phi;
             npeaks++;
           }
             }
           }
         }
       }
         }
       } 
     }
   }

   /* sort correlation values in decreasing order */
   for(i=0;i<npeaks;i++)      
     for(j=i+1;j<npeaks;j++)
       if(corr[0][i]<corr[0][j])
         for(k=0;k<4;k++)  SWAPPING( corr[k][i], corr[k][j], double );

   /* find the top N and eliminate close-by solutions*/
   int max_peaks = dim2/4;
   res[0] = corr[0][0]; /* correlation value*/
   res[1] = corr[1][0]; /* psi in deg */
   res[2] = corr[2][0]; /* the in deg */
   res[3] = corr[3][0]; /* phi in deg */
   int current_npeaks = 1;
   int bSmall;
   for(i=1; i<npeaks; i++)
   {
     bSmall = 0;
     for(j=0; j<current_npeaks; j++)
     {
       if(angle_distance(corr[1][i],corr[2][i],corr[3][i],res[j*4+1],res[j*4+2],res[j*4+3])<dist_cut*M_PI/bw)
       {
         bSmall = 1;
         break;
       }
     }
     if(bSmall == 0)
     {
       res[current_npeaks*4] = corr[0][i];
       res[current_npeaks*4+1] = corr[1][i];
       res[current_npeaks*4+2] = corr[2][i];
       res[current_npeaks*4+3] = corr[3][i];
       current_npeaks++;
     }
     if(current_npeaks >= max_peaks)
       break;
   }
   for(i=0;i<current_npeaks;i++)
   {
     res[i*4+1] = res[i*4+1]*180/M_PI;
     res[i*4+2] = res[i*4+2]*180/M_PI;
     res[i*4+3] = res[i*4+3]*180/M_PI;
   }

    for(i=0;i<4;i++)
      free(corr[i]);

    return 0;
}


/*====================================================================*/
int get_3d_index(int nx, int ny, int nz, int x, int y, int z)
{
  return x*ny*nz+y*nz+z; // C order
}

double get_3d_entry(double *data, int nx, int ny, int nz, int x, int y, int z)
{
  if(x>=nx) // force the index to be inside the range
    x = nx-1;
  if(y>=ny)
    y = ny-1;
  if(z>=nz)
    z = nz-1;

  return data[get_3d_index(nx, ny, nz, x, y, z)];
}

/**
  Enlarge the input data into a twice big output data array.
  */
int enlarge2(double *in, int dim1, int nx, int ny, int nz, double *out, int dim2)
{
  if(dim1*8 != dim2)
  {
    fprintf(stderr, "Input and output dimensions are wrong!");
    return -1;
  }

  int x,y,z;
  for(x=0; x<nx; x++)
  {
    for(y=0; y<ny; y++)
    {
      for(z=0; z<nz; z++)
      {
        out[get_3d_index(2*nx,2*ny,2*nz,2*x,2*y,2*z)] = get_3d_entry(in, nx,ny,nz,x,y,z);
        out[get_3d_index(2*nx,2*ny,2*nz,2*x,2*y,2*z+1)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                           get_3d_entry(in, nx,ny,nz,x,y,z+1))/2;
        out[get_3d_index(2*nx,2*ny,2*nz,2*x,2*y+1,2*z)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                           get_3d_entry(in, nx,ny,nz,x,y+1,z))/2;
        out[get_3d_index(2*nx,2*ny,2*nz,2*x,2*y+1,2*z+1)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                             get_3d_entry(in, nx,ny,nz,x,y,z+1)+
                                                             get_3d_entry(in, nx,ny,nz,x,y+1,z)+
                                                             get_3d_entry(in, nx,ny,nz,x,y+1,z+1))/4;
        out[get_3d_index(2*nx,2*ny,2*nz,2*x+1,2*y,2*z)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                           get_3d_entry(in, nx,ny,nz,x+1,y,z))/2;
        out[get_3d_index(2*nx,2*ny,2*nz,2*x+1,2*y,2*z+1)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                             get_3d_entry(in, nx,ny,nz,x,y,z+1)+
                                                             get_3d_entry(in, nx,ny,nz,x+1,y,z)+
                                                             get_3d_entry(in, nx,ny,nz,x+1,y,z+1))/4;
        out[get_3d_index(2*nx,2*ny,2*nz,2*x+1,2*y+1,2*z)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                             get_3d_entry(in, nx,ny,nz,x+1,y,z)+
                                                             get_3d_entry(in, nx,ny,nz,x,y+1,z)+
                                                             get_3d_entry(in, nx,ny,nz,x+1,y+1,z))/4;
        out[get_3d_index(2*nx,2*ny,2*nz,2*x+1,2*y+1,2*z+1)] = (get_3d_entry(in, nx,ny,nz,x,y,z)+
                                                               get_3d_entry(in, nx,ny,nz,x,y,z+1)+
                                                               get_3d_entry(in, nx,ny,nz,x,y+1,z)+
                                                               get_3d_entry(in, nx,ny,nz,x,y+1,z+1)+
                                                               get_3d_entry(in, nx,ny,nz,x+1,y,z)+
                                                               get_3d_entry(in, nx,ny,nz,x+1,y,z+1)+
                                                               get_3d_entry(in, nx,ny,nz,x+1,y+1,z)+
                                                               get_3d_entry(in, nx,ny,nz,x+1,y+1,z+1))/8;
      }
    }
  }

  return 0;
}


// typedef struct
// {
//   double real;
//   double imag;
// } cmplx;

// cmplx c_mul(cmplx a, cmplx b)
// {
//   cmplx res = {a.real*b.real-a.imag*b.imag, a.real*b.imag+a.imag+b.real};
//   return res;
// }

// int fourier_shift_sf(double *sf, int dim1, double r, int sizex, int sizey, int sizez, double shiftX,double shiftY,double shiftZ, double *res, int dim2)
// {
//   // assume the zero frequency of the input volume has already been shifted in the center
//   // validate the input vol
//   int b = sqrt(dim1/8);
//   if(dim2 != 8*b*b)
//     return -1;
  
//   double the = 0.0;
//   double phi = 0.0;
//   double x = 0.0;
//   double y = 0.0;
//   double z = 0.0;

//   int jj = 0, kk = 0;
//   for(; jj < b; jj++) // only loop half sphere
//   {
//     for(; kk < 2*b; kk++)
//     {
//       the = M_PI * (2*jj+1) / (4*b);
//       phi = M_PI * kk / b;
//       x = r*cos(phi)*sin(the);
//       y = r*sin(phi)*sin(the);
//       z = r*cos(the);

//       cmplx shiftx = {cos(2*M_PI/sizex*x*shiftX), -sin(2*M_PI/sizex*x*shiftX)};
//       cmplx shifty = {cos(2*M_PI/sizey*y*shiftY), -sin(2*M_PI/sizey*y*shiftY)};
//       cmplx shiftz = {cos(2*M_PI/sizez*z*shiftZ), -sin(2*M_PI/sizez*z*shiftZ)};
//       cmplx c = {sf[2*(jj*2*b+kk)], sf[2*(jj*2*b+kk)+1]};
//       c = c_mul(c_mul(c_mul(shiftx, shifty), shiftz), c);

//       res[2*(jj*2*b+kk)] = c.real;
//       res[2*(jj*2*b+kk)+1] = c.imag;
//       int idx = (2*b-jj-1)*2*b+(kk<b?kk+b:kk-b); // set the other half
//       res[2*idx] = c.real;
//       res[2*idx+1] = -c.imag;
//     }
//   }

//   return 0;
// }


/* *******************************
 * For constrained angular search
 * *******************************
 */

void idx2angle(int bw, int i, int j, int k, double *angle)
{
  double phe, phc, ohm, phi, the, psi;

  phe=i*M_PI/bw;
  phc=j*M_PI/bw;
  ohm=k*M_PI/bw;

  phi = M_PI-ohm; phi -= 2.0*M_PI*floor(phi/2.0/M_PI);
  the=M_PI-phc;
  psi=-phe; psi -= 2.0*M_PI*floor(psi/2.0/M_PI);

  angle[0] = psi;
  angle[1] = the;
  angle[2] = phi;

  return;
}

int get_constraint_vol(double *cv, int dim1, int bw, double phi, double psi, double the, double nearby)
{
  if(dim1!=8*bw*bw*bw)
  {
    fprintf(stderr,"Input's dimension not consistent!");
    return -1;
  }

  int i,j,k;
  double angle[3];
  for(i=0; i<2*bw; i++)
  {
    for(j=0; j<2*bw; j++)
    {
      for(k=0; k<2*bw; k++)
      {
        idx2angle(bw, i,j,k, angle);
        if(angle_distance(angle[0], angle[1], angle[2], psi, the, phi) < nearby)
          cv[i*4*bw*bw+j*2*bw+k] = 1;
        else
          cv[i*4*bw*bw+j*2*bw+k] = 0;
      }
    }
  }

  return 1;
}

