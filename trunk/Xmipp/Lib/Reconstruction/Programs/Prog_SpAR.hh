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

#ifndef _PROG_SPAR_HH
   #define _PROG_SPAR_HH

#include <XmippData/xmippMatrices2D.hh>

/**@name Spectrum modelling by AR filters */
//@{
/** CausalAR.
   
   This function determines the coeficients of an 2D - AR model
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
*/
double CausalAR(matrix2D<double> &Img,
   int orderR,int orderC, matrix2D<double> &ARParameters);

/** NonCausalAR.
   
   Does the same work of CausalAR but in this function the
   region of support considered is all the pixels that surround
	 	   the actual one considered. 

   	   For more details: 
   	   Huang, T.  "Advances in computer vision and image processing", 
   	   Vol.2 (1986), pp. 328 (Noncausal prediction) 

    PARAMETERS:   Img - The matrix: Here it's supposed that it comes from
    		            an image
		     
		          ordeR, orderC - The order in Rows and Columns directions of the AR
		                          model.
		          ARParameters - The matrix to store the resulting parameteres
		  		                 for the AR model.
				
   OUTPUT: 	      The function stores the AR parameters into ARParameters
   		          value for (0,0) is sigma form the model.
		          Sigma is also returned by the function.
   
   DATE:          21-2-2001
*/
double NonCausalAR(matrix2D<double> &Img,
   int orderR,int orderC, matrix2D<double> &ARParameters);

/** ARFilter.
   This function returns the AR Filter associated to an AR model. 

    PARAMETERS:   Img - The matrix - Here it's supposed that it comes from
    		        an input image
		  Filter - The matrix if Fourier space
		           that will contain the filter. It has double
		           size in columns to store complex numbers. 		
		  ARParameters - The matrix with the AR model coeficients.
				
   OUTPUT: 	The function stores the output in Filter.
   
   DATE:        6-1-2001
*/
void ARFilter(matrix2D<double> &Img, matrix2D< complex<double> > &Filter, 
   matrix2D<double> &ARParameters);

/** Combine AR Filters.
   
   This function combines two Filters in the Fourier plane
   doing their arithmetic, geometric, or armonic mean, or
   doing the armonic mean of their magnitudes. 

    PARAMETERS:   
		  Filter1 - The matrix with coeficients of the 1st filter,
		            in Fourier space
		  Filter2 - The matrix with coeficients of the 2nd filter,
		            in Fourier space
		  Filter  - The matrix with coeficients of the resulting filter,
		            in Fourier space	
 	          method _ The method employed to combine the filters
		  
   OUTPUT: 	The function stores the output in Filter.
   
   DATE:        6-1-2001
*/
void combineARFilters( const matrix2D< complex<double> > &Filter1,
		       const matrix2D< complex<double> > &Filter2,
		       matrix2D< complex<double> > &Filter,
	               const string &method);

//@}
#endif
