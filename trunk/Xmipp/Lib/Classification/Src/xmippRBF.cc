/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#include "../xmippRBF.hh"

#define Vij(v,i,j) v.itemAt(i)[j]

// Design Matrix -----------------------------------------------------------
void RBF_design_matrix(xmippCTVectors &C, matrix1D<double> &r, 
   xmippCTVectors &X, matrix2D<double> &H) _THROW {
   int p=X.size(); // Number of points to evaluate
   int m=C.size(); // Number of centers
   H.init_zeros(p,m);

   if (p==m && m==0) return;
   if (X.dimension()!=C.dimension()) 
      REPORT_ERROR(1,"RBF_design_matrix: Input vectors not of the "
         "same size");
   if (X.dimension()!=XSIZE(r)) 
      REPORT_ERROR(1,"RBF_design_matrix: Input radius not of the "
         "same size than the input vector");
   int vsize=X.dimension();
   for (int i=0; i<p; i++)
       for (int j=0; j<m; j++) {
          // Compute distance between the given point and the center
	  double z=0;
	  for (int l=0; l<vsize; l++) {
	     //cout << "X(" << i << "," << l << ")=" << Vij(X,i,l) << " " 
	     //     << "C(" << j << "," << l << ")=" << Vij(C,j,l) << " " 
	     //	    << "r=" << r(l) << endl;
	     z+=(Vij(X,i,l)-Vij(C,j,l))*(Vij(X,i,l)-Vij(C,j,l))/(r(l)*r(l));
	  }
      	  H(i,j)=exp(-z);
	  // cout << "H(" << i << "," << j << ")=" << H(i,j) << endl;
       }
}

// Training Best scale -----------------------------------------------------
void RBF_train_best_scale(xmippCTVectors &candidate_C,  xmippCTVectors &X,
   vector<double> &y, double minscale,
   double maxscale, double scalestep, xmippRBF &RBF,
   double &error, vector<double> &y_predicted) _THROW {
   int p=X.size(); // Number of points to evaluate
   int m=candidate_C.size(); // Number of centers
   if (p==m && m==0) return;
   if (X.dimension()!=candidate_C.dimension()) 
      REPORT_ERROR(1,"RBF_train_best_scale: Input vectors not of the "
         "same size");

   // Select the default radius 
   int vsize=X.dimension();
   RBF.r.resize(vsize);
   matrix1D<double> maxX(vsize), minX(vsize);
   maxX.init_zeros();
   minX.init_zeros();
   for (int i=0; i<p; i++)
       for (int l=0; l<vsize; l++) {
           maxX(l)=MAX(maxX(l),Vij(X,i,l));
           minX(l)=MIN(minX(l),Vij(X,i,l));
       }
   RBF.r=maxX-minX;
   
   // Compute model for each scale
   double  aux_error, best_error=1e38;
   vector<int> aux_idx, best_idx;
   matrix1D<double> aux_r, best_r;
   matrix1D<double> aux_w, best_w;
   
   for (double scale=minscale; scale<=maxscale; scale+=scalestep) {
      aux_r=RBF.r; aux_r*=scale;
      RBF_train(candidate_C, X, y, RBF.r, scale, aux_idx, aux_r, aux_w,
         aux_error);
      if (aux_error<best_error) {
         best_idx   = aux_idx;
	 best_r     = aux_r;
	 best_w     = aux_w;
	 best_error = aux_error;
      }
   }
   
   RBF.r=best_r;
   RBF.w=best_w;

   // Get selected centers
   RBF.C.clear();
   int imax=best_idx.size();
   RBF.C.theItems.resize(imax);
   RBF.C.theTargets.resize(imax);
   for (int i=0; i<imax; i++)
      RBF.C.itemAt(i)=candidate_C.itemAt(best_idx[i]);

   // Predict the values at X
   // The error is the standard deviation of the prediction
   // with respect to the tru value
   RBF_predict(RBF,X,y_predicted);
   double sum=0, sum2=0;
   for (int i=0; i<p; i++) {
      sum +=(y_predicted[i]-y[i]);
      sum2+=(y_predicted[i]-y[i])*(y_predicted[i]-y[i]);
   }
   error=sqrt(sum2/p-(sum/p)*(sum/p));
}

// Train a single scale ----------------------------------------------------
//#define DEBUG
void RBF_train(xmippCTVectors &C, xmippCTVectors &X,
   vector<double> &y, matrix1D<double> &r, double scale,
   vector<int> &idx_out, matrix1D<double> &r_out,
   matrix1D<double> &w_out, double &error) {
   
   int d=X.dimension();   // Dimension of vectors
   int p=X.size();        // Number of points
   int M=C.size();        // Number of centers
   r_out=r; r_out*=scale; // Multiply input radii

   // Compute design matrix
   matrix2D<double> F;
   RBF_design_matrix(C, r_out, X, F);

   // Some initialisation
   int m=0;          	   // Number of regressors
   bool finished=false;
   double yy=0;
   for (int i=0; i<p; i++) yy+=y[i]*y[i];
   double sse=yy;
   idx_out.clear();
   double lam=1e-9;   
   
   matrix2D<double> Hm, Hn, Fm, U;
   matrix1D<double> Hy, FmFm, Fty, err, hh, aux1D;
   double fy, ff, msc, min_msc;
   int age;

   // Resize all matrices to the maximum possible size
   Hm.init_zeros(YSIZE(F),M);
   Hn.init_zeros(YSIZE(F),M);
   Hy.init_zeros(M);
   hh.resize(M); 
   U.init_zeros(M,M);

   while (!finished) {
      m++; // Increment the number of regressors
      
      if (m==1) {
         // It is the first regressor ......................................
	 // Compute the change in the cost function due to the
	 // first regressor
	 matrix1D<double> FF,Fty,err; // Sum of the squares of F by columns
	 FF.init_zeros(XSIZE(F));
	 Fty.init_zeros(XSIZE(F));
	 // MATLAB: FF=sum(F.*F,1)';
	 // MATLAB: err=((F'*y).^2.*(2*lam+FF))./(lam+FF).^2
	 FOR_ALL_ELEMENTS_IN_MATRIX2D(F) {
	    FF(j)+=F(i,j)*F(i,j);
	    Fty(j)+=F(i,j)*y[i];
	 }
	 err.init_zeros(XSIZE(F));
      	 FOR_ALL_ELEMENTS_IN_MATRIX1D(err) 
	    err(i)=Fty(i)*Fty(i)*(2*lam+FF(i))/((lam+FF(i))*(lam+FF(i)));

      	 // Select the maximum
	 int jmax; err.max_index(jmax);
	 idx_out.push_back(jmax);
	 
	 // Initialize Hm, the current orthogonalised design matrix
	 // and Hn, the same but with normalised columns
	 // MATLAB: f=F(:,j)
	 // MATLAB: ff=f'*f;
	 // MATLAB: fy=fy'*y;
	 // MATLAB: Hm=f;
	 // MATLAB: Hn=f/ff;
	 // MATLAB: Hy=fy;
	 // MATLAB: hh=ff;
	 fy=ff=0;
	 for (int i=0; i<YSIZE(F); i++) {
	    ff+=F(i,jmax)*F(i,jmax);
	    fy+=F(i,jmax)*y[i];
	    Hm(i,m-1)=Hn(i,m-1)=F(i,jmax);
	 }
	 for (int i=0; i<YSIZE(F); i++) Hn(i,m-1)/=ff;
	 Hy(m-1)=fy;
	 hh(m-1)=ff;
	 
      	 // Initialize Fm ready for the second iteration
	 // MATLAB: Fm=F-Hn*(Hm'*F);
	 // MATLAB: FmFm=sum(Fm.*Fm,1);
      	 Fm=F;
	 aux1D.init_zeros(XSIZE(F));
	 // Compute Hm'*F
	 for (int i=0; i<XSIZE(F); i++)
	     for (int j=0; j<YSIZE(Hm); j++)
	        aux1D(i)+=Hm(j,m-1)*F(j,i);
	 // Compute F-Hn*(Hm'*F)
	 FOR_ALL_ELEMENTS_IN_MATRIX2D(Fm)
	    Fm(i,j)-=Hn(i,m-1)*aux1D(j);
	 //Fm=F-Hn*(Hm.transpose()*F);
      	 FmFm.init_zeros(XSIZE(F)); // Squared sum of Fm by rows
	 FOR_ALL_ELEMENTS_IN_MATRIX2D(Fm)
	    FmFm(j)+=Fm(i,j)*Fm(i,j);
	 
	 // Initialize Upper triangular matrix U
	 // MATLAB: U=1
	 U(0,0)=1;
      } else {
         // It is not the first regressor ..................................
	 // Compute marginal error
	 // MATLAB: numerator=(Fm'*y).^2.*(2*lam+FmFm)
	 // MATLAB: denominator=(lam+FmFm).^^2;
	 // MATLAB: denominator(subset)=ones(m-1,1)
	 // MATLAB: err=numerator./denominator;
	 // MATLAB: err(subset)=zeros(m-1,1);
	 Fty.init_zeros(XSIZE(Fm));
	 FOR_ALL_ELEMENTS_IN_MATRIX2D(Fm)
	    Fty(j)+=Fm(i,j)*y[i];
      	 int imax=idx_out.size();
	 err.resize(XSIZE(Fm)); err.init_constant(1);
	 for (int i=0; i<imax; i++)
	    err(idx_out[i])=0;
      	 FOR_ALL_ELEMENTS_IN_MATRIX1D(err) 
	    if (err(i))
	       err(i)=Fty(i)*Fty(i)*(2*lam+FmFm(i))/((lam+FmFm(i))*(lam+FmFm(i)));

      	 // Select the maximum
	 int jmax; err.max_index(jmax);
	 idx_out.push_back(jmax);
	 #ifdef DEBUG
            cout << "j=" << jmax+1 << endl;
	 #endif
	 
	 // Collect next columns of Hm and Hn
	 // MATLAB: f=F(:,j)
	 // MATLAB: ff=f'*f;
	 // MATLAB: fy=fy'*y;
	 // MATLAB: Hm=f;
	 // MATLAB: Hn=f/ff;
	 // MATLAB: Hy=fy;
	 // MATLAB: hh=ff;
	 fy=ff=0;
	 for (int i=0; i<YSIZE(Fm); i++) {
	    ff+=Fm(i,jmax)*Fm(i,jmax);
	    fy+=Fm(i,jmax)*y[i];
	    Hm(i,m-1)=Hn(i,m-1)=Fm(i,jmax);
	 }
	 for (int i=0; i<YSIZE(Fm); i++) {Hn(i,m-1)/=ff;}
	 Hy(m-1)=fy;
	 hh(m-1)=ff;

      	 // Recompute Fm ready for the next computation
      	 // MATLAB: Fm=Fm-Hn(:,m)*(Hm(:,m)'*Fm)
	 // Fm-=Hncol*(Hmcol_t*Fm);
	 aux1D.init_zeros(XSIZE(Fm));
	 // Compute Hm(:,m)'*Fm
	 for (int i=0; i<XSIZE(Fm); i++)
	     for (int j=0; j<YSIZE(Hm); j++)
	        aux1D(i)+=Hm(j,m-1)*Fm(j,i);
	 // Compute Fm-Hn(:,m)*(Hm(:,m)'*Fm)
	 FOR_ALL_ELEMENTS_IN_MATRIX2D(Fm)
	    Fm(i,j)-=Hn(i,m-1)*aux1D(j);
      	 FmFm.init_zeros(XSIZE(Fm)); // Squared sum of Fm by rows
	 FOR_ALL_ELEMENTS_IN_MATRIX2D(Fm)
	    FmFm(j)+=Fm(i,j)*Fm(i,j);

      	 // Update U
	 U(m-1,m-1)=1;
	 for (int i=0; i<m-1; i++)
	    for (int j=0; j<YSIZE(F); j++)
	       U(i,m-1)+=Hn(j,i)*F(j,jmax);
      }
      #ifdef DEBUG
	 cout << "U="  << U.compute_avg()  << endl
              << "Hm=" << Hm.compute_avg() << endl
              << "Hn=" << Hn.compute_avg() << endl
              << "Hy=" << Hy.compute_avg() << endl
              << "hh=" << hh.compute_avg() << endl;
      #endif

      // Update the sum of squared errors
      sse-=fy*fy*(2*lam+ff)/((lam+ff)*(lam+ff));
      #ifdef DEBUG
         cout << "SSE="<< sse << endl;
      #endif

      // Compute the current MSC
      double gam=0;
      FOR_ALL_ELEMENTS_IN_MATRIX1D(hh) gam+=hh(i)/(lam+hh(i));
      msc=p*sse/((p-gam)*(p-gam));
      #ifdef DEBUG
	 cout << "MSC=" << msc << endl
              << "Gam=" << gam << endl;
      #endif

      // Are we ready to terminate
      #define conf_thresh 10000
      #define conf_wait       2
      if (m==M) { // Run out of candidates
	 finished=true;
	 if (msc<min_msc) age=0;
	 else             age++;
      } else if (sse==0) { // Reduced error to 0
         finished=true;
	 age=0;
      } else if (sse<0) { // Numerical error
         finished=true;
	 age=1;
      } else if (FmFm.compute_max()<2.3e-16) { // No more significant regressors
      	 finished=true;
	 age=0;
      } else {
         if (m==1) {min_msc=msc; age=0;}
	 else if (msc<min_msc) {
	    if (msc<min_msc*(conf_thresh-1)/conf_thresh) {
	       // Yes, it is small enough to be a new minimum
	       min_msc=msc; age=0;
	    } else {
	       // Not small enough
	       age++;
	    }
	 } else {
	    // Age old minimum
	    age++;
	 }
	 if (age>=conf_wait) finished=true;
      }
      #ifdef DEBUG
	 cout << "Finished=" << finished << endl
              << "Age=" << age << endl;

      #endif
   } // while (!finished)

   // Don't include last regressors which aged the minimum
   m-=age;
   #ifdef DEBUG
      cout << "m=" << m << endl;
   #endif
   idx_out.resize(m);

   // Also truncate recursive OLS productions
   U.resize(m,m);
   hh.resize(m);
   Hy.resize(m);
   
   // Weight vector
   FOR_ALL_ELEMENTS_IN_MATRIX1D(Hy) Hy(i)/=lam+hh(i);
   w_out=U.inv()*Hy;

   // Error measure
   error=msc;
   #ifdef DEBUG
      cout << "w=" << w_out.compute_avg() << endl;
   #endif
}
#undef DEBUG

/* Predict ----------------------------------------------------------------- */
void RBF_predict(xmippRBF &RBF,  xmippCTVectors &X,
   vector<double> &y_predicted) {
   matrix2D<double> H;
   RBF_design_matrix(RBF.C, RBF.r, X, H);
   int p=X.size();
   y_predicted.resize(p);
   for (int i=0; i<p; i++) y_predicted[i]=0;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(H) 
      y_predicted[i]+=H(i,j)*RBF.w(j);
}

/* Show and read ----------------------------------------------------------- */
ostream & operator << (ostream &out, xmippRBF &rbf) {
   out << "Model dimension: " << XSIZE(rbf.w)      << endl;
   out << "Radii:           " << XSIZE(rbf.r) << " " << rbf.r.transpose() << endl;
   out << "Weights:         " << rbf.w.transpose() << endl;
   out << "Centers:\n"        << rbf.C             << endl;
   return out;
}

istream & operator >> (istream &in, xmippRBF &rbf) {
   int p, dim;
   in.scan("Model dimension: %d",&p);
   rbf.w.resize(p);
   in.scan("Radii: %d",&dim);
   rbf.r.resize(dim); in >> rbf.r;
   in.scan("Weights:"); in >> rbf.w; rbf.w.self_transpose();
   in.scan("Centers:"); in >> rbf.C;
   return in;
}
