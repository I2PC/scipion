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
 
#include <XmippData/xmippImages.hh>
#include <Reconstruction/Programs/Prog_SpARMA.hh>

//#define DEBUG  // when DEBUG is defined, the depuration code is compilated

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
   ARMA_parameters    prm;                    // Parameters for the program
   double	          dSigma;                 // std. deviation in the ARMA model.    
   ImageXmipp 	      Img;                    // image where to perform the ARMA model
   FourierImageXmipp Filter_Out;              // image to store the filter
      
   // Obtain command line parameters form the program 
   try
   {
      prm.read(argc,argv); 	      
   }
   catch (Xmipp_error XE) {cout << XE; Usage(); exit(1);}
   
   
   try 
   {
      // read the input image
      Img.read(prm.fn_in);
     matrix2D<double> ARParameters,MAParameters;
     // Determine the ARMA model
	 dSigma=CausalARMA(Img(),prm.N_AR,prm.M_AR,prm.N_MA,prm.M_MA,
	                                           ARParameters,MAParameters);
	 cerr << "Generating Fourier mask ...\n";
      // Get the AR filter coeficients in the fourier plane
      ARMAFilter(Img(),Filter_Out(),ARParameters,MAParameters,dSigma);	
      	            
     
     cout << "Potencia de ruido " << dSigma << endl;
/************************************************/
#ifdef DEBUG

	  ofstream fichero("prueba.txt");
	  
      fichero << "AR Parameters" << endl << ARParameters << endl;
	  fichero << "dSigma" << endl << dSigma << endl;
	  fichero << "MA Parameters" << endl << MAParameters << endl;
 	  fichero << "Image" << endl << Img() << endl << endl;
      fichero << "Filter" << endl <<Filter_Out() << endl << endl;
	  fichero.close();
	  
#endif
/************************************************/

	  // Write the output filter
      Filter_Out.write(prm.fn_filter);
   } catch (Xmipp_error XE) {cout << XE; exit(1);}
 

#undef DEBUG

}       

/**************************************************************************

   NAME:          usage
   
   DESCRIPTION:   This function displays how to use the program                  
   
   DATE:          19-1-2001
  
/**************************************************************************/
void Usage() {
     cout << "This program takes a Xmipp image obtained from a electron microscope,\n "
             "computes an spectral estimation of its noise by means of an 2-D ARMA\n"
	     "model and produces an output (complex) image with the filter of the ARMA model.\n."
	     "The ARMA filter can then be applied to a random image to get simulated noise.\n\n"
	  << "  -i <input file>            : name of the image file to determine the ARMA model\n"
	  << "  -o <output file>           : Output Xmipp Fourier Image with ARMA filter\n"
	  << "  -N_AR <N_AR>               : Order of the AR part of the model in Rows direction\n"
	  << " [-M_AR <M_AR=N_AR>]         : Order of the AR part of the model in Columns direction\n"
	  << "  -N_MA <N_MA>               : Order of the MA part of the model in Rows direction\n"
	  << " [-M_MA <M_MA=N_MA>]         : Order of the MA part of the model in Columns direction\n"
          << " [-tol  <tol=1e-20>]         : SVD tolerance\n"
     ;
}



/* Colimate menu =========================================================== */
/*Colimate:

   PROGRAM SpARMA 
   {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/SpARMA/Help/spARMA.html";
      help="Spectral Estimation by means of ARMA models. Microscope noise generation";
      OPEN MENU SpARMA;
      COMMAND LINES
	  {
         + usual: xmipp_spARMA -i $FILEIN -o $FILEOUT
                  -N_AR $N_AR
				  -M_AR $M_AR
                  -N_MA $N_MA
				  -M_MA $M_MA
				  
      }
      PARAMETER DEFINITIONS
	  {
         $FILEIN
		 {
            label="Input micrograph";
            help="file in Spider format";
            type=file existing;
         }
         $FILEOUT 
		 {
            label="Ouput ARMA filter";
            help="File containig the filter generated by the ARMA model";
            type=file;
         }
		 
		 $N_AR
		 {
		     label="N_AR";
			 help="Order of the AR part of the model in Rows direction";
			 type=NATURAL;
		 }

		 $M_AR
		 {
		     label="M_AR";
			 help="Order of the AR part of the model in Columns direction";
			 type=NATURAL;
		 }

		 $N_MA
		 {
		     label="N_MA";
			 help="Order of the MA part of the model in Rows direction";
			 type=NATURAL;
		 }

		 $M_MA
		 {
		     label="M_MA";
			 help="Order of the MA part of the model in Columns direction";
			 type=NATURAL;
		 }
      }
   }
   MENU SpARMA {
      "I/O parameters"
      $FILEIN
      $FILEOUT 
      "ARMA model parameters"
       $N_AR
	   $M_AR
       $N_MA
	   $M_MA
   }
*/

