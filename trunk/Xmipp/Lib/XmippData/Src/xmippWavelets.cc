/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Antonio Jose Rodriguez Sanchez (ajr@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
/* ------------------------------------------------------------------------- */
/* WAVELETS                                                                  */
/* ------------------------------------------------------------------------- */

#include "../xmippWavelets.hh"
#include "../xmippArgs.hh"
#include "../xmippHistograms.hh"
#include "../xmippMasks.hh"

/* Wavelet ----------------------------------------------------------------- */

// Set the DWT type --------------------------------------------------------
void set_DWT_type(int DWT_type) {
   pwtset(DWT_type);
}

// DWT ---------------------------------------------------------------------

template <class T>
   void DWT(const matrix1D<T> &v,
   matrix1D<double> &result, int isign)
   {
   unsigned long int nn[1];
   unsigned long int *ptr_nn=nn-1;

   type_cast(v, result);
   nn[0] = XSIZE(result);
   double *ptr_result=MULTIDIM_ARRAY(result)-1;
   wtn(ptr_result, ptr_nn, 1, isign, pwt);
   }

template <class T>
   void DWT(const matrix2D<T> &v,
   matrix2D<double> &result, int isign)
   {
   unsigned long int nn[2];
   unsigned long int *ptr_nn=nn-1;

   type_cast(v, result);
   nn[1] = YSIZE(result);
   nn[0] = XSIZE(result);
   double *ptr_result=MULTIDIM_ARRAY(result)-1;
   wtn(ptr_result, ptr_nn, 2, isign, pwt);
   }

template <class T>
   void DWT(const matrix3D<T> &v,
   matrix3D<double> &result, int isign)
   {
   unsigned long int nn[2];
   unsigned long int *ptr_nn=nn-1;

   type_cast(v, result);
   nn[2] = ZSIZE(result);
   nn[1] = YSIZE(result);
   nn[0] = XSIZE(result);
   double *ptr_result=MULTIDIM_ARRAY(result)-1;
   wtn(ptr_result, ptr_nn, 3, isign, pwt);
   }

// IDWT --------------------------------------------------------------------
void IDWT(const matrix1D<double> &v, matrix1D<double> &result) {
   DWT(v,result,-1);
}

void IDWT(const matrix2D<double> &v, matrix2D<double> &result) {
   DWT(v,result,-1);
}

void IDWT(const matrix3D<double> &v, matrix3D<double> &result) {
   DWT(v,result,-1);
}

// Instantiation -----------------------------------------------------------
template <class T>
   void instantiate_XmippWavelets1(T &t) {
      matrix1D<T> v1;
      matrix2D<T> v2;
      matrix3D<T> v3;
      matrix1D<double> v1d;
      matrix2D<double> v2d;
      matrix3D<double> v3d;
      DWT(v1,v1d);
      DWT(v2,v2d);
      DWT(v3,v3d);
      int x;
      SelectDWTBlock(0,v1,"0",x,x);
      SelectDWTBlock(0,v2,"00",x,x,x,x);
      SelectDWTBlock(0,v3,"000",x,x,x,x,x,x);
   }

void instantiate_XmippWavelets() {
   double d; instantiate_XmippWavelets1(d);
   int    i; instantiate_XmippWavelets1(i);
}

// Lowpass DWT -------------------------------------------------------------
void DWT_lowpass(const matrix2D<double> &v, matrix2D<double> &result) {
   matrix2D<double> dwt, aux;
   result.init_zeros(YSIZE(v),XSIZE(v)/2);
   DWT(v,dwt);
   int Nx=Get_Max_Scale(XSIZE(v));
   for (int s=0; s<Nx; s++) {
      // Perform the inverse DWT transform of the low pass
      dwt.resize(XSIZE(dwt)/2,YSIZE(dwt)/2);
      IDWT(dwt,aux);
      // Copy the result to the 01 quadrant of the result
      int x1,y1,x2,y2,x,y,i,j;
      SelectDWTBlock(s,v,"01",x1,x2,y1,y2);
      for (y=y1, i=0; y<=y2; y++, i++)
         for (x=x1, j=0; x<=x2; x++, j++)
	    result(y,x)=aux(i,j);
   }
}

// Select block ------------------------------------------------------------
#define DWT_Imin(s,smax,l) (int)((l=='0')?0:pow(2.0,smax-s-1))
#define DWT_Imax(s,smax,l) (int)((l=='0')?pow(2.0,smax-s-1)-1:pow(2.0,smax-s)-1)
template <class T>
void SelectDWTBlock(int scale, const matrix1D<T> &I,
const string &quadrant,
int &x1, int &x2) {   
   double Nx = Get_Max_Scale(XSIZE(I));
   I.physical2logical(DWT_Imin(scale,Nx,quadrant[0]),x1);
   I.physical2logical(DWT_Imax(scale,Nx,quadrant[0]),x2);
}

template <class T>
void SelectDWTBlock(int scale, const matrix2D<T> &I,
const string &quadrant,
int &x1, int &x2, int &y1, int &y2) {   
   double Nx = Get_Max_Scale(XSIZE(I));
   double Ny = Get_Max_Scale(YSIZE(I));
   x1=DWT_Imin(scale,Nx,quadrant[0]);
   y1=DWT_Imin(scale,Ny,quadrant[1]);
   x2=DWT_Imax(scale,Nx,quadrant[0]);
   y2=DWT_Imax(scale,Ny,quadrant[1]);
   I.physical2logical(y1,x1,y1,x1);
   I.physical2logical(y2,x2,y2,x2);
}

template <class T>
void SelectDWTBlock(int scale, const matrix3D<T> &I,
const string &quadrant,
int &x1, int &x2, int &y1, int &y2, int &z1, int &z2) {   
   double Nx = Get_Max_Scale(XSIZE(I));
   double Ny = Get_Max_Scale(YSIZE(I));
   double Nz = Get_Max_Scale(ZSIZE(I));
   x1=DWT_Imin(scale,Nx,quadrant[0]);
   y1=DWT_Imin(scale,Ny,quadrant[1]);
   z1=DWT_Imin(scale,Nz,quadrant[2]);
   x2=DWT_Imax(scale,Nx,quadrant[0]);
   y2=DWT_Imax(scale,Ny,quadrant[1]);
   z2=DWT_Imax(scale,Nz,quadrant[2]);
   I.physical2logical(z1,y1,x1,z1,y1,x1);
   I.physical2logical(z2,y2,x2,z2,y2,x2);
}

// Quadrant .---------------------------------------------------------------
string Quadrant2D(int q) {
   switch (q) {
      case 0: return "00"; break;
      case 1: return "01"; break;
      case 2: return "10"; break;
      case 3: return "11"; break;
   }
}

string Quadrant3D(int q) {
   switch (q) {
      case 0: return "000"; break;
      case 1: return "001"; break;
      case 2: return "010"; break;
      case 3: return "011"; break;
      case 4: return "100"; break;
      case 5: return "101"; break;
      case 6: return "110"; break;
      case 7: return "111"; break;
   }
}

// Provide block -----------------------------------------------------------
#define DWT_Scale(i,smax) ((int)((i==0)?smax-1:(ABS((CEIL(log10(i+1)/log10(2.0))-smax)))))
#define DWT_Quadrant1D(i,s,smax) ((s!=smax-1)?'1':((i==0)?'0':'1'))
#define DWT_QuadrantnD(i,s,sp,smax) \
  ((s!=sp)?'0':DWT_Quadrant1D(i,s,smax))
   
void Get_Scale_Quadrant(int size_x, int x,
   int &scale, string &quadrant) {
   double Nx = Get_Max_Scale(size_x);
   quadrant="x";
   scale=DWT_Scale(x,Nx);
   quadrant[0]=DWT_Quadrant1D(x,scale,Nx);
}

void Get_Scale_Quadrant(int size_x, int size_y, int x, int y,
   int &scale, string &quadrant) {
   double Nx = Get_Max_Scale(size_x);
   double Ny = Get_Max_Scale(size_y);
   quadrant="xy";
   double scalex=DWT_Scale(x,Nx);
   double scaley=DWT_Scale(y,Ny);
   scale = (int)(MIN(scalex,scaley));
   quadrant[1]=DWT_QuadrantnD(y,scaley,scale,Ny);
   quadrant[0]=DWT_QuadrantnD(x,scalex,scale,Nx);
}

void Get_Scale_Quadrant(int size_x, int size_y, int size_z,
   int x, int y, int z,
   int &scale, string &quadrant) {
   double Nx = Get_Max_Scale(size_x);
   double Ny = Get_Max_Scale(size_y);
   double Nz = Get_Max_Scale(size_z);
   quadrant="xyz";
   double scalex=DWT_Scale(x,Nx);
   double scaley=DWT_Scale(y,Ny);
   double scalez=DWT_Scale(z,Nz);
   scale = (int)(MIN(scalez,MIN(scalex,scaley)));
   quadrant[2]=DWT_QuadrantnD(z,scalez,scale,Nz);
   quadrant[1]=DWT_QuadrantnD(y,scaley,scale,Ny);
   quadrant[0]=DWT_QuadrantnD(x,scalex,scale,Nx);
}

// Clean quadrant ----------------------------------------------------------
void clean_quadrant(matrix2D<double> &I, int scale, const string &quadrant) {
   int x1, y1, x2, y2;
   matrix1D<int> corner1(2), corner2(2);
   matrix1D<double> r(2);
   SelectDWTBlock(scale, I, quadrant, XX(corner1), XX(corner2),
      YY(corner1), YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) I(r)=0;
}

void clean_quadrant(matrix3D<double> &I, int scale, const string &quadrant) {
   int x1, y1, z1, x2, y2, z2;
   SelectDWTBlock(scale, I, quadrant, x1, x2, y1, y2, z1, z2);
   matrix1D<int> corner1(3), corner2(3);
   matrix1D<double> r(3);
   SelectDWTBlock(scale, I, quadrant, XX(corner1), XX(corner2),
      YY(corner1), YY(corner2), ZZ(corner1), ZZ(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) I(r)=0;
}

// Soft thresholding -------------------------------------------------------
void soft_thresholding(matrix2D<double> &I, double th) {
   FOR_ALL_ELEMENTS_IN_MATRIX2D(I) 
      if (ABS(I(i,j))>th) 
         if (I(i,j)>0) I(i,j)-=th; else I(i,j)+=th;
      else I(i,j)=0;
}

void soft_thresholding(matrix3D<double> &I, double th) {
   FOR_ALL_ELEMENTS_IN_MATRIX3D(I) 
      if (ABS(I(k,i,j))>th) 
         if (I(k,i,j)>0) I(k,i,j)-=th; else I(k,i,j)+=th;
      else I(k,i,j)=0;
}

// Adaptive soft thresholding ----------------------------------------------
void adaptive_soft_thresholding_block(matrix2D<double> &I, int scale,
   const string &quadrant, double sigma) {
   // Compute block variance
   matrix1D<int> corner1(2), corner2(2);
   matrix1D<double> r(2);
   SelectDWTBlock(scale, I, quadrant,
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   double dummy, avg, stddev;
   I.compute_stats(avg,stddev, dummy, dummy,corner1,corner2);

   // Now denoise
   double th=sigma*sigma/stddev;
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      if (ABS(I(r))>th) 
         if (I(r)>0) I(r)-=th; else I(r)+=th;
      else I(r)=0;
   }
}

double compute_noise_power(matrix2D<double> &I) {
   // Compute histogram of the absolute values of the DWT coefficients
   // at scale=0
   histogram1D hist;
   double avg, stddev, min_val, max_val;
   I.compute_stats(avg,stddev,min_val,max_val);
   hist.init(0,MAX(ABS(min_val),ABS(max_val)),100);   

   matrix1D<int> corner1(2), corner2(2);
   matrix1D<double> r(2);
   SelectDWTBlock(0, I, "01",
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2)
      hist.insert_value(ABS(I(r)));

   SelectDWTBlock(0, I, "10",
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2)
      hist.insert_value(ABS(I(r)));

   SelectDWTBlock(0, I, "11",
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2)
      hist.insert_value(ABS(I(r)));

   return hist.percentil(50)/0.6745;
}

void adaptive_soft_thresholding(matrix2D<double> &I, int scale) {
   double sigma=compute_noise_power(I);
   for (int s=0; s<=scale; s++) {
      adaptive_soft_thresholding_block(I,s,"01",sigma);
      adaptive_soft_thresholding_block(I,s,"10",sigma);
      adaptive_soft_thresholding_block(I,s,"11",sigma);
   }
}

// Keep central part -------------------------------------------------------
void DWT_keep_central_part(matrix2D<double> &I, double R) {
   Mask_Params mask(INT_MASK);
   mask.type=BINARY_DWT_CIRCULAR_MASK;
   if (R==-1) mask.R1=(double)XSIZE(I)/2+1;
   else       mask.R1=R;
   mask.smin=0; mask.smax=Get_Max_Scale(XSIZE(I));
   mask.quadrant="xx";
   mask.resize(I);
   mask.generate_2Dmask();
   mask.apply_mask(I,I);
}


#ifdef NEVER_DEFINED
// CO: Doesn't work very fine
// Bayesian Wiener filtering -----------------------------------------------
#define DEBUG
void estimate_gaussian_for_scale_with_limits(const matrix2D<double> &I,
   int scale, double min_val, double max_val,
   double &alpha, double &sigma) {
   matrix1D<int> corner1(2), corner2(2);
   matrix1D<double> r(2);
   #ifdef DEBUG
      histogram1D hist;
      hist.init(min_val,max_val,20);
   #endif
   
   double sum=0, sum2=0;
   int    N_accounted=0, N=0;

   // "01"
   SelectDWTBlock(scale, I, "01",
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      if (I(r)>min_val && I(r)<max_val) {
         N_accounted++;
	 sum+=I(r);
	 sum2+=I(r)*I(r);
	 #ifdef DEBUG
	    hist.insert_value(I(r));
	 #endif
      }
      N++;
   }

   // "10"
   SelectDWTBlock(scale, I, "10",
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      if (I(r)>min_val && I(r)<max_val) {
         N_accounted++;
	 sum+=I(r);
	 sum2+=I(r)*I(r);
	 #ifdef DEBUG
	    hist.insert_value(I(r));
	 #endif
      }
      N++;
   }

   // "11"
   SelectDWTBlock(scale, I, "11",
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      if (I(r)>min_val && I(r)<max_val) {
         N_accounted++;
	 sum+=I(r);
	 sum2+=I(r)*I(r);
	 #ifdef DEBUG
	    hist.insert_value(I(r));
	 #endif
      }
      N++;
   }

   if (N_accounted!=0) {
      sum2/=N_accounted;
      sum/=N_accounted;
      sigma=sqrt(sum2-sum*sum);
      alpha=(double)N_accounted/(double)N;
   } else {sigma=alpha=0;}
   #ifdef DEBUG
      hist/=hist.sampleNo();
      hist.write((string)"Hist_"+ItoA(scale));;
   #endif
}
#undef DEBUG

void estimate_gaussian_for_scale(const matrix2D<double> &I,
   int scale, double &alpha, double &sigma) {
   double min_val, max_val;
   I.compute_double_minmax(min_val,max_val);
   estimate_gaussian_for_scale_with_limits(I,scale,min_val,max_val,alpha,sigma);
   double alpha_ant, sigma_ant;
   int N=1;
   do {
      alpha_ant=alpha;
      sigma_ant=sigma;
      max_val=3*sigma;
      min_val=-max_val;
      estimate_gaussian_for_scale_with_limits(I,scale,min_val,max_val,alpha,sigma);
      N++;
   } while (ABS(sigma-sigma_ant)/sigma>0.02 && N<10);
}

// See Bijaoui Signal Processing 2002, 82:709-712
// for the variable notation
void bayesian_wiener_filtering_block(const matrix2D<double> &I,
   int scale, const string &quadrant, double N) {
   // Compute block average
   double avg=0, Nb=0;
   matrix1D<int> corner1(2), corner2(2);
   matrix1D<double> r(2);
   SelectDWTBlock(scale, I, quadrant,
      XX(corner1),XX(corner2),YY(corner1),YY(corner2));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      avg+=I(r); Nb++;
   }
   avg/=Nb;
   
   // Compute second and fourth moments
   double M2=0, M4=0;
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      double d=(I(r)-avg)*(I(r)-avg);
      M2+=d; M4+=d*d;
   }
   M2/=Nb;
   M4/=Nb;
   
   // Compute a and S
   double a, S;
   if (M2-N<0 || M4/3-N*N<0) {a=0; S=M2;}
   else {
      a=MIN(1,(M2-N)*(M2-N)/(M4/3-N*N));
      S=(M2-N)/a;
   }

   cout << scale << " " << quadrant << " M2=" << M2 << " M4=" << M4
        << " a=" << a << " S=" << S << " N=" << N << endl;

   // Now denoise
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1,corner2) {
      double G_y_S_N=gaussian1D(I(r),S+N);
      double G_y_N  =gaussian1D(I(r),N);
      I(r)*=(a*S/(S-N)*G_y_S_N)/((1-a)*G_y_N+a*G_y_S_N);
   }
}

//#define DEBUG
void bayesian_wiener_filtering(matrix2D<double> &I, int scale) {
   // Compute the amount of noise
   double alpha,N;
   estimate_gaussian_for_scale(I,0,alpha,N);
   N*=N; // compute the power of the noise
   
   // Denoise all scales up to the given one
   for (int s=0; s<=scale; s++) {
      bayesian_wiener_filtering_block(I,s,"01",N);
      bayesian_wiener_filtering_block(I,s,"10",N);
      bayesian_wiener_filtering_block(I,s,"11",N);
   }
}
#undef DEBUG
#endif
