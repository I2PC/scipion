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

#ifdef _HAVE_VTK

#include <XmippData/xmippFilters.hh>
#include "../Prog_SpAR.hh"

/**************************************************************************

   NAME:          ComputeTermA
   
   DESCRIPTION:   This function computes term called A. (see ARFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
		  ARParameters - Parameters of the AR model to perform the 
		                 calculus.
				
   OUTPUT: 	The value of A
   
   DATE:        26-1-2001
   
/**************************************************************************/
double ComputeTermA(matrix1D<double> &dDigitalFreq, matrix2D<double> &ARParameters)
{
    double A=0;
    
    
    for(int p=STARTINGY(ARParameters);p<=FINISHINGY(ARParameters);p++)
    {
        for(int q=STARTINGX(ARParameters);q<=FINISHINGX(ARParameters);q++)
        {
	   // The term for (p,q)=(0,0) is not part of the AR model. It
	   // contains sigma.
	   if(!(p==0 && q==0)) 
           {
	      A+=MAT_ELEM(ARParameters,p,q) *
	            cos((-2)*PI*(p*YY(dDigitalFreq)+q*XX(dDigitalFreq)));
           }		  
        }
    }
    return A;
}


/**************************************************************************

   NAME:          ComputeTermB
   
   DESCRIPTION:   This function computes term called B. (see ARFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
		  ARParameters - Parameters of the AR model to perform the 
		                 calculus.
				
   OUTPUT: 	The value of B
   
   DATE:        26-1-2001
   
/**************************************************************************/
double ComputeTermB(matrix1D<double> &dDigitalFreq, matrix2D<double> &ARParameters)
{
    double B=0;
        
    for(int p=STARTINGY(ARParameters);p<=FINISHINGY(ARParameters);p++)
    {
        for(int q=STARTINGX(ARParameters);q<=FINISHINGX(ARParameters);q++)
        {
	   // The term for (p,q)=(0,0) is not part of the AR model. It
	   // contains sigma.
	   if(!(p==0 && q==0)) 
           {
	      B+=MAT_ELEM(ARParameters,p,q) *
	             sin((-2)*PI*(p*YY(dDigitalFreq)+q*XX(dDigitalFreq)));
           }		  
        }
    }
    return B;
}



/**************************************************************************

   NAME:          CausalAR
   
   DESCRIPTION:   This function determines the coeficients of an 2D - AR model
                                   pf  qf
   		  Img(y,x)=(-1.0)*sum(sum( AR(p,q)*Img(y-p,x-q)) + sigma^2 * h(y,x)
		                  p=0 q=q0
		  (except for the case where p=0 and q=0)		  		  
		  		  	  
    		  that adjust to the matrix provided as argument. To do 
		  that work, it solves the Yule-Walker equations for espectral
		  estimation, involving the calculus of autocorrelation factors.
		  It returns the AR model coeficients in a matrix ARParameters(p,q)
		  ARParameters(0,0) contains the sigma^2 value, and it's also 
		  returned by the function.
		  
		  In this program the region of support considered can be the upper NSHP 
		  (NonSymmetric Half Plane) or lower NSHP.
		  {p= 1, ...,pf; q=q0, ..., 0, ...,qF} U {p=0; q=0, ...,qF}
		  (See Ultramicroscopy 68 (1997), pp. 276)
		  
		  For more details: 
		  Dudgeon "Multidimensional DSP", 
		  Prentice Hall, signal proc. series, pp. 325 

    PARAMETERS:   Img - The matrix - Here it's supposed that it comes from
    		        an image
		  ordeR, orderC - The order in Rows and Columns directions of the AR
		                  model.
		  ARParameters - The matrix to store the resulting parameteres
		  		for the AR model.
				
   OUTPUT: 	The function stores the AR parameters into ARParameters
   		value for (0,0) is sigma form the model.
		Sigma is also returned by the function.
   
   DATE:        19-1-2001
  
/**************************************************************************/
// #define DEBUG
double CausalAR(matrix2D<double> &Img,
   int orderR,int orderC, matrix2D<double> &ARParameters)
{

   int l0;  // initial Rows coeficient for correlations
   int lF;  // final Rows coeficient for correlations
   int m0;  // initial Columns coeficient for correlations
   int mF;  // final Columns coeficient for correlations      
   int p0;  // initial Rows coeficient for AR parameters
   int pF;  // final Rows coeficient for AR parameters 
   int q0;  // initial Columns coeficient for AR parameters
   int qF;  // final Columns coeficient for AR parameters

   // Choose region of the image that will affect the pixel considered (quadrants)
   // The region considered is the upper NSHP when orderR is greater than zero      
   if(orderR>=0) 
   {
       l0=0;        lF=orderR; 
       m0=-orderC;  mF=orderC;
       p0=0;        pF=orderR; 
       q0=-orderC;  qF=orderC;   
   }
   else
   // The region considered is the lower NSHP when orderR is greater than zero      
   {
       l0=orderR;   lF=0; 
       m0=-orderC;  mF=orderC;
       p0=orderR;   pF=0; 
       q0=-orderC;  qF=orderC;   
   }

   int eq,co; // auxiliary indexes for equation and coeficient

   // Compute correlation matrix (we'll name it R) 
   matrix2D<double> R((lF-p0)-(l0-pF)+1,(mF-q0)-(m0-qF)+1);
   R.init_zeros();
   STARTINGY(R)=l0-pF;
   STARTINGX(R)=m0-qF;
   cerr << "Generating correlation coefficients ...\n";
   FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
      MAT_ELEM(R,i,j)=correlation(Img,Img,NULL,i,j);

   // Set equation system for AR model
   matrix2D<double> Coeficients;
   matrix1D<double> Indep_terms, ARcoeficients;

   Coeficients.resize((lF-l0+1)*(mF-m0+1)-mF,(lF-l0+1)*(mF-m0+1)-mF);

   STARTINGX(Coeficients)=0;
   STARTINGY(Coeficients)=0;

   Indep_terms.resize((lF-l0+1)*(mF-m0+1)-mF);
   STARTINGX(Indep_terms)=0;

   ARcoeficients.resize((lF-l0+1)*(mF-m0+1)-mF);
   STARTINGX(ARcoeficients)=0;

   // Generate matrix
   eq=0; // equation number      
   for (int l=lF; l>=l0; l--)
   {   
     for (int m=mF;m>=m0; m--)
	 {   
	   // This line is included to avoid points not in the NSHP.
       if(l==0 && m!=0 && SGN(m)!=SGN(orderR)) 
	   {
	      continue; 
	   } 
	   else
	   {

          // take the independet terms from the correlation matrix
          Indep_terms(eq)=(-1.0)*R(l,m);

	      co=0; // coeficient number
          // take the coeficients
          for (int p=pF; p>=p0; p--)
          {
              for (int q=qF; q>=q0; q--) 
         	  {	      
		         // This line is included to avoid points not in the NSHP.
                 if(p==0 && q!=0 && SGN(q)!=SGN(orderR))
			     {
			         continue;
			     }
			     else
			     {
			         // in the site for a(0,0) coeficient we put the sigma coeficient  
                        if(p==0 && q==0)
                        {
                           // The coeficient for sigma, de std. dev. of the random process
                           // asociated with the AR model is determined here.
                           // It's supposed that the filter asociated to the random process
                           // has a response h(0,0)=1 and h(i,j)=0 elsewhere.
                           if(l==0 && m==0) 
                             MAT_ELEM(Coeficients,eq,co)=-1;
                           else
                             MAT_ELEM(Coeficients,eq,co)=0;
                        }		        
                        else
                        {
                           MAT_ELEM(Coeficients,eq,co)=R(l-p,m-q);
                        }
		          // increment the coeficient counter into an equation		     
		          co++;
	              }
		      } 
           }

	       // take the next equation
	       eq++;
	    }
     }   
  }

   // Solve the equation system to determine the AR model coeficients and sigma.
   cerr << "Solving AR model ...\n";
/*****************************************/
   #ifdef DEBUG
     ofstream fichero("coeficientes.txt");
     fichero << Coeficients << endl << Indep_terms << endl ;
     fichero.close();
   #endif
/******************************************/
   
   solve(Coeficients,Indep_terms,ARcoeficients);

   // Put the ARcoeficients into the matrix given as parameter
   ARParameters.resize(pF-p0+1,qF-q0+1);
   STARTINGY(ARParameters)=p0;
   STARTINGX(ARParameters)=q0;

   co=0;
   for (int p=pF; p>=p0; p--)
       for (int q=qF; q>=q0; q--)
       {
	   // This line is to avoid points out the NSHP.
	   if(p==0 && q!=0 && SGN(q)!=SGN(orderR))
	   {
	      continue;
	   }
	   else
	   {
	      MAT_ELEM(ARParameters,p,q)=VEC_ELEM(ARcoeficients,co);
	      co++;
	   }
       }

   // return the sigma coeficient
   return MAT_ELEM(ARParameters,0,0);
}
#undef DEBUG

/* Non causal AR ----------------------------------------------------------- */
double NonCausalAR(matrix2D<double> &Img,
   int orderR,int orderC, matrix2D<double> &ARParameters)
{

   int l0;  // initial Rows coeficient for correlations
   int lF;  // final Rows coeficient for correlations
   int m0;  // initial Columns coeficient for correlations
   int mF;  // final Columns coeficient for correlations      
   int p0;  // initial Rows coeficient for AR parameters
   int pF;  // final Rows coeficient for AR parameters 
   int q0;  // initial Columns coeficient for AR parameters
   int qF;  // final Columns coeficient for AR parameters

       l0=-orderR;  lF=orderR; 
       m0=-orderC;  mF=orderC;
       p0=-orderR;  pF=orderR; 
       q0=-orderC;  qF=orderC;   

   int eq,co; // auxiliary indexes for equation and coeficient

   // Compute correlation matrix (we'll name it R) 
   matrix2D<double> R((lF-p0)-(l0-pF)+1,(mF-q0)-(m0-qF)+1);
   R.init_zeros();
   STARTINGY(R)=l0-pF;
   STARTINGX(R)=m0-qF;
   cerr << "Generating correlation coefficients ...\n";
   FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
      MAT_ELEM(R,i,j)=correlation(Img,Img,NULL,i,j);

   // Set equation system for AR model
   matrix2D<double> Coeficients;
   matrix1D<double> Indep_terms, ARcoeficients;

   Coeficients.resize((lF-l0+1)*(mF-m0+1),(lF-l0+1)*(mF-m0+1));

   STARTINGX(Coeficients)=0;
   STARTINGY(Coeficients)=0;

   Indep_terms.resize((lF-l0+1)*(mF-m0+1));
   STARTINGX(Indep_terms)=0;

   ARcoeficients.resize((lF-l0+1)*(mF-m0+1));
   STARTINGX(ARcoeficients)=0;

   // Generate matrix
   eq=0; // equation number      
   for (int l=lF; l>=l0; l--)
   {   
       for (int m=mF;m>=m0; m--)
	   {   	      
          // take the independet terms from the correlation matrix
          Indep_terms(eq)=(-1.0)*R(l,m);

	      co=0; // coeficient number
          // take the coeficients
          for (int p=pF; p>=p0; p--)
          {
              for (int q=qF; q>=q0; q--) 
         	  {	      
      		     // in the site for a(0,0) coeficient we put the sigma coeficient  
                 if(p==0 && q==0)
                 {
                     // The coeficient for sigma, de std. dev. of the random process
                     // It's supposed that the filter asociated to the random process
                     // has a response h(0,0)=1 and h(i,j)=0 elsewhere.
                     if(l==0 && m==0) 
                      {
                          MAT_ELEM(Coeficients,eq,co)=-1;
                      }
                      else
                      {
                          MAT_ELEM(Coeficients,eq,co)=0;
                       }                       
                  }		        
                   else
                  {
                      MAT_ELEM(Coeficients,eq,co)=R(l-p,m-q);
                  }
		       // increment the coeficient counter into an equation		     
		       co++;
	           }
		   } 

	       // take the next equation
	       eq++;

		}
	}

   // Solve the equation system to determine the AR model coeficients and sigma.
   cerr << "Solving AR model ...\n";
   solve(Coeficients,Indep_terms,ARcoeficients);

   // Put the ARcoeficients into the matrix given as parameter
   ARParameters.resize(pF-p0+1,qF-q0+1);
   STARTINGY(ARParameters)=p0;
   STARTINGX(ARParameters)=q0;

   co=0;
   for (int p=pF; p>=p0; p--)
       for (int q=qF; q>=q0; q--)
	   {
	      MAT_ELEM(ARParameters,p,q)=VEC_ELEM(ARcoeficients,co);
	      co++;
	   }

   // return the sigma coeficient
   return MAT_ELEM(ARParameters,0,0);
}

/* AR Filter --------------------------------------------------------------- */
#define DEBUG
void ARFilter(matrix2D<double> &Img, vtkImageData *&vtkFilter, 
   matrix2D<double> &ARParameters) {

   int iDimensions[2];              // dimensions of an image
   double A,B;                      /* Two Terms involved in calculation 
                                        of the filter defined by the AR model */
					
   // The next lines get a image prepared to receive the FFT 
   // components of the filter.	
   SPEED_UP_vtk;   
   vtkFilter=NULL;
   vtkFilter=vtkImageData::New();   
   vtkFilter->SetDimensions(XSIZE(Img),YSIZE(Img),1);
   // vtkFilter->SetUpdateExtent(0,Img.ColNo()-1,0,Img.RowNo()-1,0,0);
   vtkFilter->SetNumberOfScalarComponents(2);   
   vtkFilter->SetScalarType(VTK_FLOAT);
   vtkFilter->AllocateScalars();

   /**** Then, the Fourier Transform of the filter defined by the AR model is done ****/
   
   // Get image dimensions
   vtkFilter->GetDimensions(iDimensions); 
 
   // Compute the filter  
   float *vtk=(float *)vtkFilter->GetScalarPointer();	
   matrix1D<int>    iIndex(2);       // index in a VTK image
   matrix1D<double> dDigitalFreq(2); // digital frequency corresponding to and
   				     // index

/*****************************************/
   #ifdef DEBUG
     ofstream fichero("coeficientes.txt");
   #endif
/******************************************/

   FOR_ALL_ELEMENTS_IN_VTK(vtkFilter)
   {
     // Compute dDigitalFreq
     XX(dDigitalFreq)=j/(double)iDimensions[0];
     YY(dDigitalFreq)=i/(double)iDimensions[1];
     // Compute terms A and B for Filter
     A=ComputeTermA(dDigitalFreq,ARParameters);	   
     B=ComputeTermB(dDigitalFreq,ARParameters);      

   double sigma=sqrt(MAT_ELEM(ARParameters,0,0));
    *vtk=sigma*(1+A)/((1+A)*(1+A)+B*B);    // Real part of the filter coeficient
    *(vtk+1)=sigma*(-B)/((1+A)*(1+A)+B*B); // Imag part of the filter coeficient

/*****************************************/
   #ifdef DEBUG
     fichero << "A " << A << " B " << B << " po " << (*vtk)*(*vtk)+(*(vtk+1))*(*(vtk+1)) << endl ;
   #endif
/******************************************/
    // Increment the pointer
    vtk+=2;
   }
/*****************************************/
   #ifdef DEBUG
     fichero.close();
   #endif
/******************************************/
  
}

#undef DEBUG
/* Combine AR filters ------------------------------------------------------ */
void combineARFilters(vtkImageData *&vtkFilter1,
		      vtkImageData *&vtkFilter2,
		      vtkImageData *&vtkFilter,
	              string method)
{
   // The next lines get a new image prepared to receive data
   SPEED_UP_vtk;   
   vtkFilter=vtkImageData::New();
   vtkFilter->CopyStructure(vtkFilter1);
   vtkFilter->SetScalarType(VTK_FLOAT);
   vtkFilter->SetNumberOfScalarComponents(2);
   vtkFilter->AllocateScalars();

   float *vtk1=(float *)vtkFilter1->GetScalarPointer();
   float *vtk2=(float *)vtkFilter2->GetScalarPointer();
   float *vtk=(float *)vtkFilter->GetScalarPointer();
   
   FOR_ALL_ELEMENTS_IN_VTK(vtkFilter)
   {
	   // Get the coeficients from the matrices
	   complex<double> c1(*vtk1,*(vtk1+1));
	   complex<double> c2(*vtk2,*(vtk2+1));
	   	    
	   if(method=="arithmetic_mean")
	   {
	   // Take the arithmetic mean of the two filters
	   complex<double> c=(c1+c2)/2.0;	  
	   *vtk=c.real();
	   *(vtk+1)=c.imag();
	   }
	   else if(method=="armonic_mean")
	   {
	   // Take the armonic mean of the two filters
	   complex<double> c=2.0/((1.0/c1)+(1.0/c2));
	   *vtk=c.real();
	   *(vtk+1)=c.imag();
	   }
	   else if(method=="geometric_mean")
	   {
	   // Take the geometric mean of the two filters
	   complex<double> c=sqrt(c1*c2);
	   *vtk=c.real();
	   *(vtk+1)=c.imag();
	   }
	   else if(method=="armonic_spectrum")
	   {
	   // Take the armonic mean of the spectrum. 
	   *vtk=sqrt(2/(1/(abs(c1)*abs(c1))+1/(abs(c2)*abs(c2))));
	   *(vtk+1)=0;
	   }
	   else
	   {
	   // by default take the arithmetic mean of the two filters
	   complex<double> c=(c1+c2)/2.0;
	   *vtk=c.real();
	   *(vtk+1)=c.imag();
	   }
	   
	   // Increment the pointers 
	   vtk1+=2; vtk2+=2; vtk+=2;	  	   	   
   }

}

#endif
