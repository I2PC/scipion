/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Javier Ángel Velázquez Muriel (javi@cnb.uam.es)
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

#include <XmippData/xmippArgs.hh>
#include "../Prog_SpARMA.hh"
#include <XmippData/xmippFFT.hh> // auto_correlation_matrix()


double    global_tolerance;

/**************************************************************************

   NAME:          read
   
   DESCRIPTION:   Read parameters from command line

    PARAMETERS:   arc and **argv
				
/**************************************************************************/
void ARMA_parameters::read(int argc, char **argv)
{
      fn_in     = get_param(argc,argv,"-i");
      fn_filter = get_param(argc,argv,"-o");                 
      N_AR         = AtoI(get_param(argc,argv,"-N_AR"));
      M_AR         = AtoI(get_param(argc,argv,"-M_AR","0"));
         if (M_AR==0) M_AR=N_AR;
      N_MA         = AtoI(get_param(argc,argv,"-N_MA"));
      M_MA         = AtoI(get_param(argc,argv,"-M_MA","0"));
         if (M_MA==0) M_MA=N_MA;

   global_tolerance         = AtoF(get_param(argc,argv,"-tol","1e-20"));
}


/**************************************************************************

   NAME:          read
   
   DESCRIPTION:   Read parameters from file

    PARAMETERS:   arc and **argv
				
/**************************************************************************/
void ARMA_parameters::read(const FileName &InputFile) _THROW
{
   // Read the parameters file to get every one
   FILE *file;
   if ((file = fopen(InputFile.c_str(), "r")) == NULL)
   	 REPORT_ERROR(1,(string)"ARMA_parameters::read: There is a problem "
            "opening the file "+InputFile);
			  	         
      fn_in        = get_param(file,"image", 0,"");
      fn_filter    = get_param(file,"ARMAfile",0,"");                 
      N_AR         = AtoI(get_param(file,"N_AR",0));
      M_AR         = AtoI(get_param(file,"M_AR",0,"0"));
         if (M_AR==0) M_AR=N_AR;
      N_MA         = AtoI(get_param(file,"N_MA",0));
      M_MA         = AtoI(get_param(file,"M_MA",0,"0"));
         if (M_MA==0) M_MA=N_MA;
   fclose(file);
}

// Write to a file =========================================================
void ARMA_parameters::write(const FileName &fn_prm, bool rewrite) _THROW {
   ofstream fh_param;
   if (!rewrite) fh_param.open(fn_prm.c_str(),ios::app);
   else          fh_param.open(fn_prm.c_str(),ios::out);
   if (!fh_param)
      REPORT_ERROR(1,(string)"ARMA_parameters::write: There is a problem "
            "opening the file "+fn_prm);
   fh_param << "# ARMA parameters\n";
   if (fn_in!="")
      fh_param << "image=" << fn_in     << endl;
   if (fn_filter!="")
      fh_param << "ARMAfile=" << fn_filter << endl;
   fh_param << "N_AR="     << N_AR      << endl
            << "M_AR="     << M_AR      << endl
            << "N_MA="     << N_MA      << endl
            << "M_MA="     << M_MA      << endl
   ;
   fh_param << endl;
   fh_param.close();
}
      	 	  	  
/**************************************************************************

   NAME:          ComputeDenominator
   
   DESCRIPTION:   This function computes the denominator of the ARMA filter.
                 (see ARMAFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
		  ARParameters - Parameters of the AR model to perform the 
		                 calculus.
				
   OUTPUT: 	The value of the denominator, a complex number
   
   DATE:        26-1-2001
   
/**************************************************************************/
complex<double> ComputeDenominator(matrix1D<double> &dDigitalFreq,
                          matrix2D<double> &ARParameters)
{
    complex<double> A(0,0);

   long NumberOfARParameters=ARParameters.RowNo();
   
   for (long n=0 ; n<NumberOfARParameters; n++)
   {
      int p=(int)ARParameters(n,0);	  
      int q=(int)ARParameters(n,1);
	  complex<double> c(0,(-2)*PI*(p*YY(dDigitalFreq)+q*XX(dDigitalFreq)));
      complex<double> d(MAT_ELEM(ARParameters,n,2),0);
      A+=d*(exp(c)+exp(-c));
	  
   }
 
   return (1.0-A);
}

/**************************************************************************

   NAME:          ComputeDenominator_fast
   
   DESCRIPTION:   This function computes the denominator of the ARMA filter.
                 (see ARMAFilter function), but using real arithmetics

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
		  ARParameters - Parameters of the AR model to perform the 
		                 calculus.
				
   OUTPUT: 	The value of the denominator, a complex number
   
   DATE:        26-1-2001
   
/**************************************************************************/
complex<double> ComputeDenominator_fast(matrix1D<double> &dDigitalFreq,
                          matrix2D<double> &ARParameters)
{
    double A=0;

   long NumberOfARParameters=ARParameters.RowNo();
   
   for (long n=0 ; n<NumberOfARParameters; n++)
   {
      int p=(int)ARParameters(n,0);	  
      int q=(int)ARParameters(n,1);
      A+=MAT_ELEM(ARParameters,n,2)*2*
	     cos((-2)*PI*(p*YY(dDigitalFreq)+q*XX(dDigitalFreq)));
   }
 
   complex<double> An(A,0);
   return (1.0-An);
}


/**************************************************************************

   NAME:          ComputeNumerator
   
   DESCRIPTION:   This function computes the numerator of the ARMA filter.
                 (see ARMAFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
		  ARParameters - Parameters of the AR part of the model to perform the 
		                 calculus.
				
   OUTPUT: 	The value of the numerator, a complex number
   
   DATE:        26-1-2001
   
/**************************************************************************/
complex<double> ComputeNumerator(matrix1D<double> &dDigitalFreq, 
                                 matrix2D<double> &MAParameters)
{
    complex<double> B(0,0);
         
   long NumberOfMAParameters=MAParameters.RowNo();
   
   for (long n=0 ; n<NumberOfMAParameters; n++)
   {
      int p=(int)MAParameters(n,0);	  
      int q=(int)MAParameters(n,1);
	  complex<double> c(0,(-2)*PI*(p*YY(dDigitalFreq) + q*XX(dDigitalFreq)));
      complex<double> d(MAT_ELEM(MAParameters,n,2),0);
      B+=d*(exp(c)+exp(-c));
   }
 
   return (1.0+B);   
}

/**************************************************************************

   NAME:          ComputeNumerator_fast
   
   DESCRIPTION:   This function computes the numerator of the ARMA filter.
                 (see ARMAFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
		  ARParameters - Parameters of the AR part of the model to perform the 
		                 calculus.
				
   OUTPUT: 	The value of the numerator, a complex number
   
   DATE:        26-1-2001
   
/**************************************************************************/
complex<double> ComputeNumerator_fast(matrix1D<double> &dDigitalFreq, 
                                 matrix2D<double> &MAParameters)
{
        
   double B=0;
   long NumberOfMAParameters=MAParameters.RowNo();
   
   for (long n=0 ; n<NumberOfMAParameters; n++)
   {
      int p=(int)MAParameters(n,0);	  
      int q=(int)MAParameters(n,1);
      B+=MAT_ELEM(MAParameters,n,2)*2*
	              cos((-2)*PI*(p*YY(dDigitalFreq) + q*XX(dDigitalFreq)));
   }
 
   complex<double> Bn(B,0);
  
   return (1.0+Bn);   
}


void First_Quadrant_Neighbors(int N, int M,matrix2D<double> &Neighbors)
{
   long NumberOfPoints=(N+1)*M;
   int n;  
   Neighbors.resize(NumberOfPoints,3);
   
   n=0; // Number of neighbors found so far. 

   for (int p=N; p>=0; p--)
   {
       for (int q=M; q>0; q--)
       {
			Neighbors(n,0)=(double)p;
			Neighbors(n,1)=(double)q;
			n++;
       }
   }
}

void Second_Quadrant_Neighbors(int N, int M,matrix2D<double> &Neighbors)
{
   long NumberOfPoints=N*(M+1);
   int n;  
   Neighbors.resize(NumberOfPoints,3);
   
   n=0; // Number of neighbors found so far. 

   for (int p=N; p>=1; p--)
   {
       for (int q=0; q>=-M; q--)
       {
			Neighbors(n,0)=(double)p;
			Neighbors(n,1)=(double)q;
			n++;
       }
   }
}

#undef DEBUG



/**************************************************************************

   NAME:          advert_of_ill_condition
   
   DESCRIPTION:   This function is called when an ARMA model is suspicious
                  of being ill-conditioned

    PARAMETERS:   
				
   OUTPUT: 	
   
   DATE:        26-1-2001
   
/**************************************************************************/
void advert_of_ill_condition(void)
{
   cerr << "Warning: ARMA filter for this image is ill-conditioned.\n" 
        << "         Some values of the filter couldn't not be correct.\n";
}

// ****************************************************************************************
// See the .hh file for a description of this function
//#define DEBUG
double CausalARMA(matrix2D<double> &Img, int N_AR, int M_AR,
                  int N_MA, int M_MA, matrix2D<double> &ARParameters,
				  matrix2D<double> &MAParameters)
{
   double dSigma; // To store de sigma coeficient of the model
    
   // Calculate the autocorrelation matrix
   matrix2D<double> R;
   auto_correlation_matrix(Img,R);
   R.set_Xmipp_origin();

   /**********************************************************************/	
   // Set equation system for AR part of the model
   /**********************************************************************/
      cerr << "Generating AR part of the model ...\n";

   matrix2D<double> Coeficients;
   matrix1D<double> Indep_terms, ARcoeficients;
   matrix2D<double> N3;

   // Assign the support region for the AR part of the model (N1)
   First_Quadrant_Neighbors(N_AR,M_AR,ARParameters); 
   // Assign the support region for the MA part of the model (N2)
   Second_Quadrant_Neighbors(N_MA,M_MA,MAParameters); 
   // Assign the support region for the AR equations (N3)
   // Here is the same of N1, but it hasn´t to be
   First_Quadrant_Neighbors(N_AR,M_AR,N3); 


   long NumberOfARParameters=ARParameters.RowNo();
   long NumberOfMAParameters=MAParameters.RowNo();
   

   Coeficients.resize(NumberOfARParameters,NumberOfARParameters);
   STARTINGX(Coeficients)=0;
   STARTINGY(Coeficients)=0;

   Indep_terms.resize(NumberOfARParameters);
   STARTINGX(Indep_terms)=0;

   ARcoeficients.resize(NumberOfARParameters);
   STARTINGX(ARcoeficients)=0;   
   
   // Generate matrix (eq stands for equation number and co for coeficents)
   for (long eq=0 ; eq<NumberOfARParameters; eq++)
   {
      // take the independet term from the correlation matrix (or calculate it
	  // if it was not calculated before).
      int l=(int)N3(eq,0);	  
      int m=(int)N3(eq,1);
	  Indep_terms(eq)=MAT_ELEM(R,l,m);
   
      // take the coeficients
      for (long co=0 ; co<NumberOfARParameters; co++)
      {
		 // Take the pertinent coeficient form the correlation matrix (or calculate it)
		 int alpha1= (int)(N3(eq,0) - ARParameters(co,0)); 
		 int alpha2= (int)(N3(eq,1) - ARParameters(co,1));
		 int beta1 = (int)(N3(eq,0) + ARParameters(co,0)); 
		 int beta2 = (int)(N3(eq,1) + ARParameters(co,1));
         MAT_ELEM(Coeficients,eq,co)= R(alpha1,alpha2) + R(beta1,beta2);					 
	  }  
   }
   
   /**********************************************************************/	
   // Solve the equation system to determine the AR model coeficients and sigma.
   /**********************************************************************/
	
   
   // Check for posible ill-definition of the  linear system to solve the 
   // AR coeficients
   if(Coeficients.det()==0) 
   {
      advert_of_ill_condition();
   }
   
   cerr << "Solving AR part of the model ...\n";   
   solve_by_svd(Coeficients, Indep_terms, ARcoeficients,global_tolerance);

/*****************************************/
   #ifdef DEBUG
     ofstream fichero("coeficientes.txt");
     double suma=0;
	 FOR_ALL_ELEMENTS_IN_MATRIX1D(ARcoeficients)
	    suma+=ARcoeficients(i);
	 fichero.precision(15); 
	 fichero << Indep_terms << endl
	         << "Determinante: " << Coeficients.det() << endl
             << "Coeficientes AR por descomposicion LU: " << endl << ARcoeficients 
			 << "Suma del coeficientes AR " << suma << endl;
			 
   // Get info about the coeficients matrix with svd
   matrix2D<double> u,v;
   matrix1D<double> w;

   svdcmp(Coeficients,u,w,v);      
   double min=1e10;
	 FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
	    if(w(i)<min) min=w(i);

   solve_by_svd(Coeficients, Indep_terms, ARcoeficients,global_tolerance);
   suma=0;
	 FOR_ALL_ELEMENTS_IN_MATRIX1D(ARcoeficients)
	    suma+=ARcoeficients(i);
    fichero << "Autovalores de la matriz : " << endl << w << endl
           << "Coeficientes AR por descomposicion SV: " << endl << ARcoeficients
		   << "Suma del coeficientes AR " << suma << endl;
		   
   fichero.close();   
		   
   cout << N_AR << " " << N_MA << " " << min << " " << suma << endl;   
   exit(0);
        
   #endif
/******************************************/


   // Assing the results to the matrix given as parameter
   for (long n=0 ; n<NumberOfARParameters; n++)
   {
     MAT_ELEM(ARParameters,n,2)=VEC_ELEM(ARcoeficients,n);
   }

   /**********************************************************************/
   // Determine the sigma coeficient from the equation for (p,q)=(0,0)
   /**********************************************************************/
   double dSum=0;
   for (long n=0 ; n<NumberOfARParameters; n++)
   {
      int p=(int)ARParameters(n,0);	  
      int q=(int)ARParameters(n,1);
	  dSum+=MAT_ELEM(ARParameters,n,2) * MAT_ELEM(R,p,q);
   }

   	// And calculate sigma
	dSigma=(MAT_ELEM(R,0,0)-2*dSum); 
	double idSigma=1.0/dSigma;
	
    /**********************************************************************/
	// Determine the MA parameters of the model using the AR parameters and
	// sigma
    /**********************************************************************/

   cerr << "Solving MA part of the model ...\n";

   for (long n=0 ; n<NumberOfMAParameters; n++)
   {
   
       dSum=0;
       double MAn0=MAT_ELEM(MAParameters,n,0);
       double MAn1=MAT_ELEM(MAParameters,n,1);
       for (long m=0 ; m<NumberOfARParameters; m++)
	   {
	   double ARm0=MAT_ELEM(ARParameters,m,0);
	   double ARm1=MAT_ELEM(ARParameters,m,1);
    	   int alpha1= (int)(MAn0 - ARm0); 
	   int alpha2= (int)(MAn1 - ARm1);
           int beta1 = (int)(MAn0 + ARm0); 
	   int beta2 = (int)(MAn1 + ARm1);
	   dSum += MAT_ELEM(ARParameters,m,2) * (
	       MAT_ELEM(R,alpha1,alpha2) + MAT_ELEM(R,beta1,beta2));
	   }
 
       int p=(int)MAn0;	  
       int q=(int)MAn1;
       MAT_ELEM(MAParameters,n,2)=(MAT_ELEM(R,p,q)-dSum)*idSigma;
   }	      

   // return the sigma coeficient
   return dSigma;
}
#undef DEBUG


//****************************************************************************************
// See the .hh file for a description of this function
//#define DEBUG
void ARMAFilter(matrix2D<double> &Img, matrix2D< complex<double> > &Filter, 
                matrix2D<double> &ARParameters, matrix2D<double> &MAParameters,
				double dSigma)
{
#define MINIMUM_PERMITED_VALUE_OF_AR_DENOMINATOR    0

   complex<double> A,B;              /* Two Terms involved in calculation 
                                     of the filter defined by the ARMA model */
									// A is the denominator
                                    // B is the numerator

   matrix1D<double> dDigitalFreq(2); // digital frequency corresponding to and
   				                     // index

   // Resize de Filter to the image dimensions
   Filter.resize(Img.RowNo(),Img.ColNo());
   Filter.init_zeros();

/*****************************************/
   #ifdef DEBUG
     ofstream fichero("coeficientes_filtro.txt");
   #endif
/******************************************/

   // Compute the filter (only half the values are computed)
   // The other half is computed based in symmetry. 
    int sizeX=Img.ColNo();
   int sizeY=Img.RowNo();   
   for(int i=0;i<sizeY;i++)
       for(int j=0;j<(sizeX/2);j++)
	   {
		   // Compute dDigitalFreq
		   XX(dDigitalFreq)=j/(double)sizeX;
		   YY(dDigitalFreq)=i/(double)sizeY;
		   // Compute terms A and B for Filter
		   B=ComputeNumerator_fast(dDigitalFreq,MAParameters);	   
		   A=ComputeDenominator_fast(dDigitalFreq,ARParameters);      

   // Check for posible problems calculating the value of A that could lead
   // to ARMA filters with erroneous values due to a value of A near 0.
   if(A.real()<MINIMUM_PERMITED_VALUE_OF_AR_DENOMINATOR)
   {
      advert_of_ill_condition();
   }
    
/*****************************/
   #ifdef DEBUG
     fichero << "A " << A << " B " << B << " ARMA " << sqrt(dSigma*B/A) << endl ;
   #endif
/*****************************/

		  // This is the filter proposed by original paper
		  // As the filter values are symmetrical respect the center (frequency zero),
		  // take advantage of this fact.
		  Filter(sizeY-1-i,sizeX-1-j)=Filter(i,j)=sqrt(dSigma*B/A);
	   }
	
	
   
/*****************************************/
   #ifdef DEBUG
     fichero.close();
   #endif
/******************************************/

}

#undef DEBUG


#ifdef NEVER_DEFINED
/**************************************************************************

   NAME:          Causal_Nearest_Neighbors
   
   DESCRIPTION:   This function determines the nearest points (i,j) to (0,0) that 
                  meet this conditions
				  - Both i and j are integers
				  - (i,j) doesn´t belong to the NSHP defined by N_MA, M_MA 

				  As many points as those in the NSHP defined by N_AR, M_AR are
				  dtermined.
				  
    PARAMETERS:  N_AR,M_AR - Orders for rows and columns of the AR part of the model
	             N_MA,M_MA - Orders for rows and columns of the MA part of the model
				
   OUTPUT: 	     A matrix with the neighbors
   
   DATE:        26-3-2001
   
/**************************************************************************/
//#define DEBUG
void Causal_Nearest_Neighbors(int N_AR, int M_AR,int N_MA,int M_MA, matrix2D<double> &Neighbors)
{
   long NumberOfPoints=N_AR*(2*M_AR+1)+M_AR;
   int r, i,Top_i,j,n ;
   
   Neighbors.resize(NumberOfPoints,3);
   
   n=0; // Number of neighbors found so far. 
   
   // Go through all posible distances
   for(r=1;;r++)
   {
      // Go through every posible j (column) coordinate
	  for(j=-r;j<=r;j++)
	  {
	     // Determine top of posibilities for i (row)
		 Top_i=(int)floor(sqrt(r*r-j*j));

		 // Take every point into cosideration descending from i
		 // if a point is at distance r-1, it should have been taken before
		 // since every r is checked.
		 for(i=Top_i;((i*i+j*j)>((r-1)*(r-1))) && i>=0;i--)
		 {

            // If the point considered is not in the NSPH of the MA part
			if(( i>N_MA || j<-M_MA || j>M_MA ) || (i==0 && j<0))
			{
			    Neighbors(n,0)=(double)i;
			    Neighbors(n,1)=(double)j;
			    n++;
			    // Check if all required points have been determined
				if(n==NumberOfPoints) 
                {
                    return;
                }
			}
			
		 }
	  }
   }

}

void Causal_Neighbors(int N, int M,matrix2D<double> &Neighbors)
{
   long NumberOfPoints=N*(2*M+1)+M;
   int n;  
   Neighbors.resize(NumberOfPoints,3);
   
   n=0; // Number of neighbors found so far. 

   for (int p=N; p>=0; p--)
   {
       for (int q=M; q>=-M; q--)
       {
		   // This line is to avoid points out the NSHP.
		   if(p==0 && q<=0)
		   {
	    	  break;
		   }
		   else
		   {
			    Neighbors(n,0)=(double)p;
			    Neighbors(n,1)=(double)q;
				n++;
		   }
       }
   }
}

void AntiCausal_Neighbors(int N, int M,matrix2D<double> &Neighbors)
{
   long NumberOfPoints=N*(2*M+1)+M;
   int n;  
   Neighbors.resize(NumberOfPoints,3);
   
   n=0; // Number of neighbors found so far. 

   for (int p=0; p>=-N; p--)
   {
       for (int q=M; q>=-M; q--)
       {
		   // This line is to avoid points out the lower NSHP.
		   if(p==0 && q>=0)
		   {
	    	  continue;
		   }
		   else
		   {
			    Neighbors(n,0)=(double)p;
			    Neighbors(n,1)=(double)q;
				n++;
		   }
       }
   }
}

/**************************************************************************

   NAME:          NonCausal_Nearest_Neighbors
   
   DESCRIPTION:   This function determines the nearest points (i,j) to (0,0) that 
                  meet this conditions
				  - Both i and j are integers
				  - (i,j) doesn´t belong to the plane defined by N_MA, M_MA 

				  As many points as those in the plane defined by N_AR, M_AR are
				  determined.
				  
    PARAMETERS:  N_AR,M_AR - Orders for rows and columns of the AR part of the model
	             N_MA,M_MA - Orders for rows and columns of the MA part of the model
				
   OUTPUT: 	     A matrix with the neighbors, Neignbors
   
   DATE:        26-3-2001
   
/**************************************************************************/
//#define DEBUG
void NonCausal_Nearest_Neighbors(int N_AR, int M_AR,int N_MA,int M_MA, matrix2D<int> &Neighbors)
{
   long NumberOfPoints=(2*N_AR+1)*(2*M_AR+1)-1;
   int r, i,Top_i,j,n ;
   
   Neighbors.resize(NumberOfPoints,3);
   
   n=0; // Number of neighbors found so far. 
   
   // Go through all posible distances
   for(r=1;;r++)
   {
      // Go through every posible j (column) coordinate
	  for(j=-r;j<=r;j++)
	  {
	     // Determine top of posibilities for i (row)
		 Top_i=(int)floor(sqrt(r*r-j*j));
		 // Take every point into cosideration descending from i
		 // if a point is at distance r-1, it should have been taken before
		 // since every r is checked.
		 for(i=Top_i;((i*i+j*j)>((r-1)*(r-1))) && i>=0;i--)
		 {
		    // If the point considered is not in the plane of the MA part
			if( i>N_MA || j<-M_MA || j>M_MA )
			{
			    Neighbors(n,0)=i;
			    Neighbors(n,1)=j;
			    n++;
			    // Check if all required points have been determined
				if(n==NumberOfPoints) 
                {
                    return;
                }
			    
                // Add the symmetric point too, if i is not zero
                if(i!=0)
                {
                	Neighbors(n,0)=(-i);
			    	Neighbors(n,1)=j;
			    	n++;
			    	// Check if all required points have been determined
					if(n==NumberOfPoints)
                    {
                        return;
                    }
                }
			}
		 }
	  }
   }
}

void NonCausal_Neighbors(int N_AR, int M_AR,matrix2D<int> &Neighbors)
{
   long NumberOfPoints=(2*N_AR+1)*(2*M_AR+1)-1;
   int n;  
   Neighbors.resize(NumberOfPoints,2);
   
   n=0; // Number of neighbors found so far. 

   for (int p=N_AR; p>=-N_AR; p--)
   {
       for (int q=M_AR; q>=-M_AR; q--)
       {
		   // This line is to avoid points out the NSHP.
		   if(p==0 && q==0)
		   {
	    	  continue;
		   }
		   else
		   {
			    Neighbors(n,0)=p;
			    Neighbors(n,1)=q;
				n++;
		   }
       }
   }
}
#undef DEBUG
#endif
