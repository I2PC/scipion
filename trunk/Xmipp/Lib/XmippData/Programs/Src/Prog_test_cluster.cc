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

#include "../Prog_test_cluster.hh"
#include <XmippData/Src/NumericalRecipes.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippHistograms.hh>
#include <XmippData/xmippMasks.hh>

// Empty constructor -------------------------------------------------------
Test_cluster_parameters::Test_cluster_parameters() {
   fn_selfile="";
   fn_out="";
   distance=MAHALANOBIS;
}

// Read from command line --------------------------------------------------
void Test_cluster_parameters::read(int argc, char **argv) _THROW {
   fn_selfile = get_param(argc,argv,"-i");
   fn_out     = get_param(argc,argv,"-o","");
   mask.read(argc,argv);
   if (check_param(argc,argv,"-dist")) {
      string strdist=get_param(argc,argv,"-dist");
      if      (strdist=="euclidean")   distance=EUCLIDEAN;
      else if (strdist=="correlation") distance=CORRELATION;
      else if (strdist=="mahalanobis") distance=MAHALANOBIS;
   }
   produce_side_info();
}

// Produce side info -------------------------------------------------------
void Test_cluster_parameters::produce_side_info() _THROW {
   SF_in.read(fn_selfile);
   int Xdim, Ydim; SF_in.ImgSize(Ydim,Xdim);
   mask.resize(Ydim,Xdim);
   mask.generate_2Dmask();
   
   if (distance==MAHALANOBIS) build_covariance_matrix();
}

// Show --------------------------------------------------------------------
void Test_cluster_parameters::show() {
   cout << "Selfile: " << fn_selfile << endl
        << "Output:  " << fn_out     << endl;
   cout << "Distance:";
   switch (distance) {
      case EUCLIDEAN:   cout << "Euclidean\n"; break;
      case CORRELATION: cout << "Correlation\n"; break;
      case MAHALANOBIS: cout << "Mahalanobis\n"; break;
   }
   mask.show();
}


// Usage -------------------------------------------------------------------
void Test_cluster_parameters::usage() {
    cerr << "Test_cluster\n"
         << "   -i <selfile>      : Selfile with the cluster images\n"
	 << "  [-o <txt file>]    : By default, the histogram is printed on the screen\n"
	 << "  [-dist <distance=mahalanobis>]:"
	 << "                     : Valid distances are euclidean, correlation, mahalanobis\n"
    ;
    mask.usage();	 
}

// Build covariance matrix -------------------------------------------------
void Test_cluster_parameters::build_covariance_matrix() {
   SF_in.go_first_ACTIVE();
   bool first=true;
   int N=0;
   // Compute the mean of the population
   while (!SF_in.eof()) {
      matrix1D<double> v;
      FileName fn_img=SF_in.NextImg();
      ImageXmipp I; I.read(fn_img);
      mask.produce_vector(I(),v);
      
      if (first) {
         covariance.init_zeros(XSIZE(v),XSIZE(v));
	 mean.init_zeros(XSIZE(v));
	 first=false;
      }
      
      mean+=v;
      N++;
   }
   if (N==0) return;
   mean/=N;
   
   // Compute the covariance
   SF_in.go_first_ACTIVE();
   while (!SF_in.eof()) {
      matrix1D<double> v;
      FileName fn_img=SF_in.NextImg();
      ImageXmipp I; I.read(fn_img);
      mask.produce_vector(I(),v);
      v=-mean;
      
      // v*v' (compute only the upper triangle)
      for (int i=0; i<XSIZE(v); i++) 
         for (int j=i; j<XSIZE(v); j++) {
	    covariance(i,j)+=v(i)*v(j);
	 }
   }
   covariance/=N;
   for (int i=0; i<YSIZE(covariance); i++) 
      for (int j=i+1; j<XSIZE(covariance); j++)
	 covariance(j,i)=covariance(i,j);
}

// Do the work -------------------------------------------------------------
void Test_cluster_parameters::test_cluster() {
   // Compute all distances between image pairs
   int imax=SF_in.ImgNo();
   matrix1D<double> dij(imax*(imax-1)/2);
   int p=0;
   SF_in.go_first_ACTIVE();
   FileName fni, fnj;
   cerr << "Computing all vs. all distances ...\n";
   init_progress_bar(XSIZE(dij));
   int significative=0, positive=0;
   for (int i=0; i<imax; i++) {
       fni=SF_in.get_file_number(i);
       ImageXmipp Ii; Ii.read(fni);
       matrix1D<double> vi; mask.produce_vector(Ii(),vi);
       for (int j=i+1; j<imax; j++) {
          fnj=SF_in.get_file_number(j);
          ImageXmipp Ij; Ij.read(fnj);
          matrix1D<double> vj; mask.produce_vector(Ij(),vj);
	  
          // Prepare the vectors for the distance computation
	  double avgi, avgj, stddevi, stddevj, dummy;
	  switch (distance) {
      	     case EUCLIDEAN:
	        vi-=vj;
		vj=vi;
		break;
	     case CORRELATION:
	        vi.compute_stats(avgi,stddevi,dummy,dummy);
	        vj.compute_stats(avgj,stddevj,dummy,dummy);
		vi-=avgi; vi/=stddevi;
		vj-=avgj; vj/=stddevj;
		break;
	     case MAHALANOBIS:
	        vi-=vj;
		vj=vi;
		vj=covariance*vj;
	        break;
	  }
	  
          double dist=0;
          for (int k=0; k<XSIZE(vi); k++)
	     dist+=vi(k)*vj(k);
	  dist/=XSIZE(vi);
	  dij(p++)=dist;
	  
          // Check if the correlation is significative
	  if (distance==CORRELATION) {
             // D. Sheskin
	     // Handbook of parametric and nonparametric statistical procedures
	     // Page 953
	     int n=XSIZE(vi);
	     int df=n-2;
	     float t=(float)(dist*sqrt(df/(1-dist*dist)));
	     float p_val=student_up_to_t0(t,df);
	     if (p_val>0.95) significative++;
	     
	     if (dist>=0) positive++;
	  }
	  
	  if (p%50==0) progress_bar(p);
       }
   }
   progress_bar(XSIZE(dij));

   // Compute the histogram   
   histogram1D hist;
   compute_hist(dij,hist,100);
   if (fn_out=="") cout << hist;
   else            hist.write(fn_out);
   
   // Show the number of significative correlations
   if (distance==CORRELATION) {
      cout << "There are " << significative << " correlations out of "
           << XSIZE(dij) << " that are significantly greater than 0\n";
	   
      int expected=XSIZE(dij)/2;
      double chi2=(positive-expected)*(positive-expected)/expected+
                  (XSIZE(dij)-positive-expected)*(XSIZE(dij)-positive-expected)/expected;
      double p_val=chi2_from_t0(chi2,1);
      cout << chi2 << " " << p_val << endl;
      if (p_val<0.05) {
         // It is significant
	 if (positive>expected) 
	    cout << "There is a significantly positive correlation: "
	         << positive << " out of " << XSIZE(dij) << endl;
	 else
	    cout << "There is a significantly negative correlation: "
	         << positive << " out of " << XSIZE(dij) << endl;
      } else {
         // It is not significant
         cout << "The hypothesis that positive and negative correlations are equally distributed\n"
	      << "cannot be rejected with a confidence of 95%\n";
      }
   }
}

