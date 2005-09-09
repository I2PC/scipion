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

#ifndef _PROG_SPARMA_HH
   #define _PROG_SPARMA_HH
   #define AR          0   // To distinguish between an AR parameter or a MA   
   #define MA	       1   // parameter, or sigma in the output ARMAParameters matrix
   #define SIGMA       2   // returned by functions.
						   
#include <XmippData/xmippMatrices2D.hh>

/**@name Spectrum modelling by ARMA filters */
//@{
/** CausalARMA.
   
   This function determines the coeficients of an 2D - ARMA model
   		    
   Img(y,x)=sum( AR(p,q)*Img(y+p,x+q)) + sqrt(sigma) * e(y,x)

    
   Where:
        (p,q) is a vector belonging to a support region N1
		AR(p,q) is the matrix with the AR part coeficients of the model
		sigma is the variance of the random correlated input e(x,y)
		e(x,y) is random correlated input with zero mean and cross-spectral
		density:
		
		E[e(x,y)*Img(x+a,y+b)]= 
		                           sqrt(sigma)            if (p,q)=(0,0)
								   sqrt(sigma)*MA(p,q)    if (p,q) belongs to N2
								   0                      otherwise 
		                          
		
		
		N1 - The support region for the AR part is the first quadrant region
		     defined by N_AR and M_AR
		N2 - The support region the the MA part  is the second quadrant region
		     defined by N_MA and M_MA 
			
   This model is determined following the method depicted in:
   
   R. L. Kashyap, "Characterization and Estimation of Two-Dimensional
   ARMA Models, " IEEE Transactions on Information Theory, Vol. IT-30, No.
   5, pp. 736-745, September 1984.


   				  

    PARAMETERS:   
	      Img - The matrix - Here it's supposed that it comes from an image
		  N_AR, M_AR - The order in Rows and Columns directions of the AR part
		               of the model.
		  N_AR, M_AR - The order in Rows and Columns directions of the MA part
		               of the model.
		  ARParameters - The matrix to return the resulting parameteres
		  		           for the AR part of the model.
		  MAParameters - The matrix to return the resulting parameteres
		  		           for the MA part of the model.
				
   OUTPUT: 	The function stores the ARMA parameters into ARParameters and
            MAParameters
            Every row of this output matrices has 3 values:
			1nd and 2nd- indicate the support point (p,q) for the coeficient
			3th column - indicates the value of the coeficient
 
   DATE:        26-3-2001
*/
double CausalARMA(matrix2D<double> &Img, int N_AR, int M_AR,
                  int N_MA, int M_MA, matrix2D<double> &ARParameters,
				  matrix2D<double> &MAParameters);


/** ARMAFilter.
   This function returns the ARMA Filter associated to an ARMA model.
   This filter should be applied to a random matrix of independent identically 
   distributed variables. In Xmipp the program to do that is Fourierfilter

    PARAMETERS:   Img - The matrix - Here it's supposed that it comes from
                        an input image
                  Filter - The matrix that will contain the filter.
		  ARParameters - The matrix with the AR model coeficients, as
			   is returned by CausalARMA or NonCausalARMA
		  MAParameters - The matrix with the MA model coeficients, as
                           is returned by CausalARMA or NonCausalARMA
		  dSigma        - The Sigma Coeficient for the ARMA model
				
   OUTPUT: 	The function stores the output in Filter.
   
   DATE:        26-3-2001
*/
void ARMAFilter(matrix2D<double> &Img, matrix2D< double > &Filter, 
                matrix2D<double> &ARParameters, matrix2D<double> &MAParameters,
				double dSigma);
//@}

class ARMA_parameters 
{
public:
   FileName 	fn_in;		    // Name of input image
   FileName     fn_filter;          // Name of filter image     
   int      	N_AR;               // order in the Rows direction of the AR part of the model
   int      	M_AR;               // order in the Cols direction of the AR part of the model
   int      	N_MA;               // order in the Rows direction of the MA part of the model
   int      	M_MA;               // order in the Cols direction of the MA part of the model

public:
   /// Read parameters from command line
   void read(int argc, char **argv);
   /// Read parameters from command line
   void read(const FileName &file) _THROW;
   /// Write to a file
   void write(const FileName &fn, bool rewrite=TRUE) _THROW;
};

#endif
