/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
#include "../Prog_microscope.hh"
#include <XmippData/xmippArgs.hh>

/* Read parameters --------------------------------------------------------- */
void Prog_Microscope_Parameters::read(int argc, char **argv) _THROW {
   Prog_parameters::read(argc,argv);
   fn_ctf=get_param(argc,argv,"-ctf","");
   sigma=AtoF(get_param(argc,argv,"-noise","0"));
   low_pass_before_CTF=AtoF(get_param(argc,argv,"-low_pass","0"));
   fn_after_ctf=get_param(argc,argv,"-after_ctf","");
   
   produce_side_info();
}

/* Usage ------------------------------------------------------------------- */
void Prog_Microscope_Parameters::usage() {
   Prog_parameters::usage();
   cerr << "  [-ctf <CTF file>]         : a Xmipp Fourier Image\n"
	<< "  [-low_pass <w=0>]         : low pass filter for noise before CTF\n"
        << "  [-noise <stddev=0>]       : noise to be added\n"
	<< "  [-after_ctf <spectrum>]   : a Xmipp Fourier Image with\n"
	<< "                              the root squared spectrum of noise\n"
	<< "                              after the CTF\n";
}

/* Show -------------------------------------------------------------------- */
void Prog_Microscope_Parameters::show() {
   Prog_parameters::show();
   cout << "CTF file: " << fn_ctf << endl
        << "Noise: " << sigma << endl
        << "Noise before: " << sigma_before_CTF << endl
	<< "Noise after: " << sigma_after_CTF << endl
	<< "Low pass freq: " << low_pass_before_CTF << endl
	<< "After CTF noise spectrum: " << fn_after_ctf << endl;
}


/* Produce side information ------------------------------------------------ */
//#define DEBUG
void Prog_Microscope_Parameters::produce_side_info() {
   int Zdim, Ydim, Xdim;
   get_input_size(Zdim,Ydim,Xdim);

   double before_power=0, after_power=0;
   
   if (fn_ctf!="") {
      ctf.FilterBand=ctf.FilterShape=FROM_FILE;
      ctf.fn_mask=fn_ctf;
      ctf.generate_mask(NULL);
      #ifdef DEBUG
	 ctf.write_amplitude("PPP.xmp");
      #endif
      ctf.resize_mask(Ydim,Xdim);
      #ifdef DEBUG
	 ctf.write_amplitude("QQQ.xmp");
      #endif
      before_power=ctf.mask_power();
   }

   if (low_pass_before_CTF!=0) {
      lowpass.FilterBand=LOWPASS;
      lowpass.FilterShape=RAISED_COSINE;
      lowpass.w1=low_pass_before_CTF;
   }   

   if (fn_after_ctf!="") {
      after_ctf.FilterBand=after_ctf.FilterShape=FROM_FILE;
      after_ctf.fn_mask=fn_after_ctf;
      after_ctf.generate_mask(NULL);
      after_ctf.resize_mask(Ydim,Xdim);
      after_power=after_ctf.mask_power();
   }

   // Compute noise balance
   if (after_power!=0 || before_power!=0) {
      double p=after_power/(after_power+before_power);
      sigma_after_CTF=sqrt(p)*sigma;
      sigma_before_CTF=sqrt(1-p)*sigma;
   }
}

/* Apply ------------------------------------------------------------------- */
void Prog_Microscope_Parameters::apply(matrix2D<double> &I) {
   // Add noise before CTF
   matrix2D<double> noisy;
   noisy.resize(I);
   noisy.init_random(0,sigma_before_CTF,"gaussian");
   if (low_pass_before_CTF!=0) lowpass.apply_mask(noisy);
   I += noisy;

   // Apply CTF
   if (fn_ctf!="") ctf.apply_mask(I);

   // Add noise after CTF
   noisy.init_random(0,sigma_after_CTF,"gaussian");
   if (fn_after_ctf!="") after_ctf.apply_mask(noisy);
   I += noisy;
}
#endif
