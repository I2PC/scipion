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
#include <XmippData/xmippImages.hh>
#include <Reconstruction/Programs/Prog_SpAR.hh>

/**************************************************************************

	Prototypes
  
/**************************************************************************/
void Usage(void);

/**************************************************************************

   NAME:          main
   
   DESCRIPTION:   This function takes the parameters needed from the 
   		  command line, call the function that does the actual work
		  and writes the output volume.
   
   DATE:          18-1-2001
  
/**************************************************************************/
int main(int argc, char **argv) 
{
   FileName 	fn_in;		        // Name of input image
   FileName     fn_filter;              // Name of filter image     
   int      	N;                      // order in the Rows direction of the AR models
   int      	M;                      // order in the Cols direction of the AR models
   ImageXmipp 	Img;                    // image where to perform the AR model
   FourierImageXmipp Filter_Out;        // image to store the filter
   string       method;                 // method of combine AR fliters
   double	dSigma;                 // std. deviation of the AR model.
      
   // Obtain command line parameters form the program 
   try
   {
      fn_in     = get_param(argc,argv,"-i");
      fn_filter = get_param(argc,argv,"-o");                 
      N         = AtoI(get_param(argc,argv,"-N"));
      M         = AtoI(get_param(argc,argv,"-M"));
      method     = get_param(argc,argv,"-method");
	      
   } catch (Xmipp_error XE) {cout << XE; Usage(); exit(1);}
   
   
   try 
   {
      // read the input image
      Img.read(fn_in);
      
      matrix2D<double> ARParameters;
     // If a causal filter is required by the user
      if(method=="causal")
         // Determine the AR model for the upper NSHP
	     dSigma=CausalAR(Img(),N,M, ARParameters);
      else if(method=="anticausal")
    	 // Determine the AR model for the lower NSHP
	     dSigma=CausalAR(Img(),-N,M, ARParameters);	 
      else if(method=="noncausal")
    	 // Determine the AR model for the whole plane
	     dSigma=NonCausalAR(Img(),N,M, ARParameters);	 
      else
         EXIT_ERROR(1,"No valid method specified.\n");
      cout << "Noise power = " << dSigma << endl;

      cerr << "Generating Fourier mask ...\n";
      // Get the AR filter coeficients in the fourier plane
      ARFilter(Img(),Filter_Out(),ARParameters);	
      	            
      // Write the output filter
      Filter_Out.write(fn_filter);
   } catch (Xmipp_error XE) {cout << XE; Usage(); exit(1);}
 
 
}       

/**************************************************************************

   NAME:          usage
   
   DESCRIPTION:   This function displays how to use the program                  
   
   DATE:          19-1-2001
  
/**************************************************************************/
void Usage() {
     cout << "This program takes a Xmipp image obtained from a electron microscope,\n "
             "computes an spectral estimation of its noise by means of an 2-D AR\n"
	     "model and produces an output (complex) image with the filter of the AR model.\n."
	     "The AR filter can then be applied to a random image to get simulated noise.\n\n"
	  << "  -i input file		   : name of the image file to determine the AR model\n"
	  << "  -o output file  	   : Output Xmipp Fourier Image with AR filter\n"
	  << "  -N			   : Order of the AR model in Rows direction\n"
	  << "  -M			   : Order of the AR model in Columns direction\n"
	  << "  -method causal  	   : If you want a causal filter (support is upper NSHP)\n"
	  << "          anticausal         : If you want an anticausal filter (support is lower NSHP)\n"
          << "          noncausal          : If you want a non causal filter  (support is all plane)\n"
     ;
}
