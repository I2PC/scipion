/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Antonio Jose Rodriguez Sanchez (ajr@cnb.uam.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
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

#include "../xmippImages.hh"

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
#define DWT_Scale(i,smax) ((int)((i==0)?smax-1:(ABS((CEIL(log10((double)(i+1))/log10(2.0))-smax)))))
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

// Bayesian Wiener filtering -----------------------------------------------
//DWT_Bijaoui_denoise_LL -- Bijaoui denoising at a perticular scale.
void DWT_Bijaoui_denoise_LL(matrix2D<double> &WI,int scale,
      	             	   const string &orientation,
			   double mu,double S,double N){
   matrix1D<int> x0(2), xF(2), r(2);
   SelectDWTBlock(scale, WI, orientation, XX(x0), XX(xF), YY(x0), YY(xF));
  
   double SN=S+N;
   double S_N=S/SN;
   if (S<1e-6 && N<1e-6) {
      FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(x0,xF) WI(r)=0;
   } else {
      FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(x0,xF){
	 double y=WI(r);
	 double ymu=y-mu;
	 double ymu2=-0.5*ymu*ymu;
	 double expymu2SN=exp(ymu2/SN);
	 double den=exp(ymu2/N)+expymu2SN;
	 if (den>1e-10) WI(r)=S_N*expymu2SN/den*y;
      }
   }   
}

void DWT_Bijaoui_denoise_LL(matrix3D<double> &WI,int scale,
      	             	   const string &orientation,
			   double mu,double S,double N){
   matrix1D<int> x0(3), xF(3), r(3);
   SelectDWTBlock(scale, WI, orientation, XX(x0), XX(xF), YY(x0), YY(xF),
      ZZ(x0), ZZ(xF));
  
   double SN=S+N;
   double S_N=S/SN;
   if (S<1e-6 && N<1e-6) {
      FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(x0,xF) WI(r)=0;
   } else {
      FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(x0,xF){
	 double y=WI(r);
	 double ymu=y-mu;
	 double ymu2=-0.5*ymu*ymu;
	 double expymu2SN=exp(ymu2/SN);
	 double den=exp(ymu2/N)+expymu2SN;
	 if (den>1e-10) WI(r)=S_N*expymu2SN/den*y;
      }
   }   
}

//#define DEBUG
void bayesian_solve_eq_system(
   const matrix1D<double> &power,
   const matrix1D<double> &average,
   const matrix1D<double> &Ncoefs,
   double SNR0,
   double SNRF,
   double powerI,
   double power_rest,
   bool white_noise,
   int tell,
   matrix1D<double> &estimatedS) {

   int scale_dim=XSIZE(power);

   matrix2D<double> A;
   int extra_constraints=0;
   if (white_noise) extra_constraints+=(scale_dim-1);
   A.init_zeros(2*(scale_dim-1)+2*scale_dim+2+extra_constraints,2*scale_dim);
   for(int i=1;i<scale_dim;i++){
      A(i-1,i-1)=1;
      A(i-1,i)=-1;
      A(i-1+scale_dim-1,i-1+scale_dim)=1;
      A(i-1+scale_dim-1,i+scale_dim)=-1;
   }
   for(int i=0;i<2*scale_dim;i++)
      A(i+2*(scale_dim-1),i)=-1;
       
   // Constraints on the SNR
   matrix1D<double> aux0coefs(scale_dim);
   for(int j=0;j<scale_dim;j++)
      aux0coefs(j)=Ncoefs(j)*SNR0;
   matrix1D<double> auxFcoefs(scale_dim);
   for(int j=0;j<scale_dim;j++)
      auxFcoefs(j)=Ncoefs(j)*SNRF;
   
   //initializing the second last row of A
   for(int j=0;j<scale_dim;j++)
      A(2*(scale_dim-1)+2*scale_dim,j)=(-1)*auxFcoefs(j);
   for(int j=scale_dim;j<2*scale_dim;j++)
      A(2*(scale_dim-1)+2*scale_dim,j)=Ncoefs(j-scale_dim);
   
   //initializing the last row of A
   for(int j=0;j<scale_dim;j++)
      A(2*(scale_dim-1)+2*scale_dim+1,j)=aux0coefs(j);
   for(int j=scale_dim;j<2*scale_dim;j++)
      A(2*(scale_dim-1)+2*scale_dim+1,j)=(-1)*Ncoefs(j-scale_dim);
      
   // White noise constraints
   if (white_noise)
      for (int i=0; i<scale_dim-1; i++) {
          A(YSIZE(A)-(scale_dim-1)+i,i)  =-1.01;
          A(YSIZE(A)-(scale_dim-1)+i,i+1)= 1;
      }
   
   //initialize the matrix b
   matrix1D<double> b(A.RowNo());
   
   // Initialize Aeq matrix
   matrix2D<double> Aeq;
   Aeq.init_zeros(1,2*scale_dim);
   for(int j=0;j<scale_dim;j++){
      Aeq(0,j)=Ncoefs(j);
      Aeq(0,j+scale_dim)=Ncoefs(j);
   }
   
   //initialize beq matrix
   matrix1D<double> beq;
   beq.init_zeros(1);
   beq(0)=powerI-power_rest;
   
   //initialization of Matrix C (cost matrix)
   matrix2D<double> C;
   C.init_zeros(scale_dim,2*scale_dim);
   for(int j=0;j<scale_dim;j++){
      C(j,j)=1;
      C(j,j+scale_dim)=1;
   }
   
   // initialise the estimatedS which will contain the solution vector
   estimatedS.init_zeros(2*scale_dim);
   #ifdef DEBUG
      //Writing the matrices to ASCII files for comparing with matlab
      C.write("./matrices/C.txt");
      power.write("./matrices/power.txt");
      A.write("./matrices/A.txt");
      b.write("./matrices/b.txt");
      Aeq.write("./matrices/Aeq.txt");
      beq.write("./matrices/beq.txt");
      
      cout << "Equation system Cx=d\n"
           << "C=\n" << C << endl
	   << "d=" << (power/Ncoefs).transpose() << endl
	   << "Constraints\n"
	   << "Ax<=b\n"
	   << "A=\n" << A << endl
	   << "b=" << b.transpose() << endl
	   << "Aeq x=beq\n"
	   << "Aeq=\n" << Aeq << endl
	   << "beq=" << beq.transpose() << endl;
   #endif
      
   // Solve the system
   matrix1D<double> bl,bu;
   lsqlin(C,power/Ncoefs,A,b,Aeq,beq,bl,bu,estimatedS);
   // COSS
   estimatedS/=2;
   
   #ifdef DEBUG
      cout<<"scale_dim :: "<<scale_dim<<endl;
      cout<<"--------estimatedS -------- \n";
      cout<<estimatedS;
      cout<<"--------------------------- \n";
      cout<<"Inequality constraints agreement"<<endl
          <<(A*estimatedS).transpose() << endl;
      cout<<"Equality constraints agreement"<<endl
          <<(Aeq*estimatedS).transpose()<<endl;
      cout<<"Goal function value: " <<(C*estimatedS).transpose() <<endl;
   #endif
}
#undef DEBUG

//#define DEBUG
matrix1D<double> bayesian_wiener_filtering(matrix2D<double> &WI, int allowed_scale,
   double SNR0, double SNRF, bool white_noise, int tell, bool denoise) {
   /*Calculate the power of the wavelet transformed image */
   double powerI=WI.sum2();
         
   /*Number of pixels and some constraints on SNR*/
   int Ydim=0,Xdim=0;
   WI.get_dim(Ydim,Xdim);
   int Ncoef_total=Ydim*Xdim;
   int max_scale=ROUND(log(double(Xdim))/log(2.0));  

   #ifdef DEBUG
      cout<<"powerI= "<<powerI<<endl;
      double powerWI=WI.sum2();
      cout<<"powerWI= "<<powerWI<<endl;
      cout<<"Ydim = "<<Ydim<<"  Xdim = "<<Xdim<<"\n";
      cout<<"Ncoef_total= "<<Ncoef_total<<endl;
      cout<<"max_scale = "<<max_scale<<"\n";
   #endif
   
   /*Calculate the power at each band*/
   //init the scale vector
   matrix1D<int> scale(MIN(allowed_scale+1,max_scale-1));
   FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i)=i;
   int scale_dim=scale.get_dim();
   
     //define some vectors
   matrix1D<double> power(scale_dim), average(scale_dim), Ncoefs(scale_dim);
   matrix1D<int> x0(2), xF(2), r(2);
   vector<string> orientation;
   orientation.push_back("01");
   orientation.push_back("10");
   orientation.push_back("11");
   for (int j=0;j<scale.get_dim();j++) {
      for (int k=0; k<orientation.size(); k++) {
	 SelectDWTBlock(scale(j), WI, orientation[k],
	    XX(x0), XX(xF), YY(x0), YY(xF));
	 FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(x0,xF){
      	    power(j)+=WI(r)*WI(r);
	    average(j)+=WI(r);
	 }
      }
      Ncoefs(j)=(int)pow(2.0,2*(max_scale-scale(j)-1))*orientation.size();
      average(j)=average(j)/Ncoefs(j);
   }
   
   /*Evaluate the power of the unconsidered part of the image */
   double power_rest=0.0;
   int Ncoefs_rest=0;
   SelectDWTBlock(scale(scale_dim-1), WI, "00", XX(x0), XX(xF), YY(x0), YY(xF));
   FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(x0,xF)
     power_rest+=WI(r)*WI(r);
   Ncoefs_rest=(int)pow(2.0,2*(max_scale-1-scale(scale_dim-1)));
   
   if (tell) {
      cout<<"scale = "<<endl<<scale<<endl;
      cout<<"power= "<<endl<<power<<"\n";
      cout<<"average= "<<endl<<average<<"\n";
      cout<<"Ncoefs= "<<endl<<Ncoefs<<"\n";
      cout<<"power_rest= "<<power_rest<<"\n";
      cout<<"Ncoefs_rest= "<<Ncoefs_rest<<"\n";
      cout<<"powerI= "<<powerI<<endl;
      cout<<"Total sum of powers = "<<power.sum()+power_rest<<endl;
      cout<<"Difference powers = "<<powerI-power.sum()-power_rest<<endl;
   }
   
   /*Solve the Equation System*/
   matrix1D<double> estimatedS;
   bayesian_solve_eq_system(power, average, Ncoefs,
      SNR0, SNRF, powerI, power_rest, white_noise, tell, estimatedS);
   
   if (tell) {
      cout<<"estimatedS =\n"<< estimatedS <<endl;
      double S=0,N=0;
      for (int i=0; i<scale_dim; i++) {
         N+=Ncoefs(i)*estimatedS(i);
         S+=Ncoefs(i)*estimatedS(scale_dim+i);
      }
      cout << "SNR value=" << S/N << endl << endl;
   }

   /* Apply the Bijaoui denoising to all scales >= allowed_scale */
   if (denoise)
      bayesian_wiener_filtering(WI,allowed_scale,estimatedS);

   return estimatedS;
}
#undef DEBUG

void bayesian_wiener_filtering(matrix2D<double> &WI,
   int allowed_scale, matrix1D<double> &estimatedS) {
   vector<string> orientation;
   orientation.push_back("01");
   orientation.push_back("10");
   orientation.push_back("11");

   int max_scale=ROUND(log(double(XSIZE(WI)))/log(2.0));  
   matrix1D<int> scale(MIN(allowed_scale+1,max_scale-1));
   FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i)=i;

   for(int i=0;i<XSIZE(scale);i++){
      double N=estimatedS(i);
      double S=estimatedS(i+XSIZE(scale));
      for (int k=0; k<orientation.size(); k++)
          DWT_Bijaoui_denoise_LL(WI,scale(i),orientation[k],0,S,N);
   }
}

//#define DEBUG
matrix1D<double> bayesian_wiener_filtering(matrix3D<double> &WI, int allowed_scale,
   double SNR0, double SNRF, bool white_noise, int tell, bool denoise) {
   /*Calculate the power of the wavelet transformed image */
   double powerI=WI.sum2();
         
   /*Number of pixels and some constraints on SNR*/
   int Zdim,Ydim,Xdim;
   WI.get_dim(Zdim, Ydim,Xdim);
   int Ncoef_total=Zdim*Ydim*Xdim;
   int max_scale=ROUND(log(double(Xdim))/log(2.0));  

   #ifdef DEBUG
      cout<<"powerI= "<<powerI<<endl;
      double powerWI=WI.sum2();
      cout<<"powerWI= "<<powerWI<<endl;
      cout<<"Zdim= " << Zdim << " Ydim = "<<Ydim<<"  Xdim = "<<Xdim<<"\n";
      cout<<"Ncoef_total= "<<Ncoef_total<<endl;
      cout<<"max_scale = "<<max_scale<<"\n";
   #endif
   
   /*Calculate the power at each band*/
   //init the scale vector
   matrix1D<int> scale(allowed_scale+1);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i)=i;
   int scale_dim=scale.get_dim();
   
     //define some vectors
   matrix1D<double> power(scale_dim), average(scale_dim), Ncoefs(scale_dim);
   matrix1D<int> x0(3), xF(3), r(3);
   vector<string> orientation;
   orientation.push_back("001");
   orientation.push_back("010");
   orientation.push_back("011");
   orientation.push_back("100");
   orientation.push_back("101");
   orientation.push_back("110");
   orientation.push_back("111");
   for (int j=0;j<scale.get_dim();j++) {
      for (int k=0; k<orientation.size(); k++) {
	 SelectDWTBlock(scale(j), WI, orientation[k],
	    XX(x0), XX(xF), YY(x0), YY(xF), ZZ(x0), ZZ(xF));
	 FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(x0,xF){
      	    power(j)+=WI(r)*WI(r);
	    average(j)+=WI(r);
	 }
      }
      Ncoefs(j)=(int)pow(2.0,3*(max_scale-scale(j)-1))*orientation.size();
      average(j)=average(j)/Ncoefs(j);
   }
   
   /*Evaluate the power of the unconsidered part of the image */
   double power_rest=0.0;
   int Ncoefs_rest=0;
   SelectDWTBlock(scale(scale_dim-1), WI, "000", XX(x0), XX(xF), YY(x0), YY(xF),
      ZZ(x0), ZZ(xF));
   FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(x0,xF)
     power_rest+=WI(r)*WI(r);
   Ncoefs_rest=(int)pow(2.0,3*(max_scale-1-scale(scale_dim-1)));
   
   if (tell) {
      cout<<"scale = "<<endl<<scale<<endl;
      cout<<"power= "<<endl<<power<<"\n";
      cout<<"average= "<<endl<<average<<"\n";
      cout<<"Ncoefs= "<<endl<<Ncoefs<<"\n";
      cout<<"power_rest= "<<power_rest<<"\n";
      cout<<"Ncoefs_rest= "<<Ncoefs_rest<<"\n";
      cout<<"powerI= "<<powerI<<endl;
      cout<<"Total sum of powers = "<<power.sum()+power_rest<<endl;
      cout<<"Difference powers = "<<powerI-power.sum()-power_rest<<endl;
   }
   
   /*Solve the Equation System*/
   matrix1D<double> estimatedS;
   bayesian_solve_eq_system(power, average, Ncoefs,
      SNR0, SNRF, powerI, power_rest, white_noise, tell, estimatedS);
   if (tell) {
      cout<<"estimatedS =\n"<< estimatedS <<endl;
      double S=0,N=0;
      for (int i=0; i<scale_dim; i++) {
         N+=Ncoefs(i)*estimatedS(i);
         S+=Ncoefs(i)*estimatedS(scale_dim+i);
      }
      cout << "SNR value=" << S/N << endl << endl;
   }
   
   /* Apply the Bijaoui denoising to all scales >= allowed_scale */
   if (denoise)
      bayesian_wiener_filtering(WI,allowed_scale,estimatedS);
   return estimatedS;
}
#undef DEBUG

void bayesian_wiener_filtering(matrix3D<double> &WI,
   int allowed_scale, matrix1D<double> &estimatedS) {
   vector<string> orientation;
   orientation.push_back("001");
   orientation.push_back("010");
   orientation.push_back("011");
   orientation.push_back("100");
   orientation.push_back("101");
   orientation.push_back("110");
   orientation.push_back("111");

   int max_scale=ROUND(log(double(XSIZE(WI)))/log(2.0));  
   matrix1D<int> scale(MIN(allowed_scale+1,max_scale-1));
   FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i)=i;

   for(int i=0;i<XSIZE(scale);i++){
      double N=estimatedS(i);
      double S=estimatedS(i+XSIZE(scale));
      for (int k=0; k<orientation.size(); k++)
          DWT_Bijaoui_denoise_LL(WI,scale(i),orientation[k],0,S,N);
   }
}
