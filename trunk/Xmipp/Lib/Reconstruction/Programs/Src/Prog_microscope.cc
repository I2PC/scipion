/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2001 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

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
   defocus_change=AtoF(get_param(argc,argv,"-defocus_change","0"));
   
   produce_side_info();
}

/* Usage ------------------------------------------------------------------- */
void Prog_Microscope_Parameters::usage() {
   Prog_parameters::usage();
   cerr << "  [-ctf <CTF file>]         : a Xmipp Fourier Image or a CTF description\n"
	<< "  [-low_pass <w=0>]         : low pass filter for noise before CTF\n"
        << "  [-noise <stddev=0>]       : noise to be added\n"
	<< "  [-after_ctf <spectrum>]   : a Xmipp Fourier Image or a CTF description with\n"
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
	<< "After CTF noise spectrum: " << fn_after_ctf << endl
        << "Defocus change: " << defocus_change << endl
   ;
}


/* Produce side information ------------------------------------------------ */
//#define DEBUG
void Prog_Microscope_Parameters::produce_side_info() {
   int Zdim;
   get_input_size(Zdim,Ydim,Xdim);
   matrix2D<double> aux;
   vtkImageData *vtkaux=NULL;

   double before_power=0, after_power=0;
   
   if (fn_ctf!="") {
      if (Is_FourierImageXmipp(fn_ctf)) {
         ctf.FilterBand=ctf.FilterShape=FROM_FILE;
         ctf.fn_mask=fn_ctf;
         ctf.generate_mask(NULL);
         ctf.resize_mask(Ydim,Xdim);
      } else {
         ctf.FilterBand=CTF;
         ctf.ctf.read(fn_ctf);
         ctf.ctf.enable_CTFnoise=FALSE;
         ctf.ctf.Produce_Side_Info();
         aux.init_zeros(Ydim,Xdim);
         xmippArray2VTK(aux, vtkaux,2);
         ctf.generate_mask(vtkaux);
      }

      #ifdef DEBUG
	 ctf.write_amplitude("PPP.xmp");
      #endif
      before_power=ctf.mask_power();
   }

   if (low_pass_before_CTF!=0) {
      lowpass.FilterBand=LOWPASS;
      lowpass.FilterShape=RAISED_COSINE;
      lowpass.w1=low_pass_before_CTF;
   }   

   if (fn_after_ctf!="") {
      if (Is_FourierImageXmipp(fn_ctf)) {
         after_ctf.FilterBand=after_ctf.FilterShape=FROM_FILE;
         after_ctf.fn_mask=fn_after_ctf;
         after_ctf.generate_mask(NULL);
         after_ctf.resize_mask(Ydim,Xdim);
      } else {
         after_ctf.FilterBand=CTF;
         after_ctf.ctf.read(fn_after_ctf);
         after_ctf.ctf.enable_CTF=FALSE;
         after_ctf.ctf.Produce_Side_Info();
         after_ctf.generate_mask(vtkaux);
      }
      #ifdef DEBUG
	 after_ctf.write_amplitude("PPPafter.xmp");
      #endif
      after_power=after_ctf.mask_power();
   }
   if (vtkaux!=NULL) vtkaux->Delete();

   // Compute noise balance
   if (after_power!=0 || before_power!=0) {
      double p=after_power/(after_power+before_power);
      sigma_after_CTF=sqrt(p)*sigma;
      sigma_before_CTF=sqrt(1-p)*sigma;
   }
}
#undef DEBUG

/* Apply ------------------------------------------------------------------- */
//#define DEBUG
void Prog_Microscope_Parameters::apply(matrix2D<double> &I) {
   // Add noise before CTF
   matrix2D<double> noisy;
   noisy.resize(I);
   noisy.init_random(0,sigma_before_CTF,"gaussian");
   if (low_pass_before_CTF!=0) lowpass.apply_mask(noisy);
   I += noisy;

   // Check if the mask is a defocus changing CTF
   // In that case generate a new mask with a random defocus
   if (defocus_change!=0) {
      double old_DefocusU=ctf.ctf.DeltafU;
      double old_DefocusV=ctf.ctf.DeltafV;
      matrix2D<double> aux;
      vtkImageData *vtkaux=NULL;
      ctf.ctf.DeltafU*=rnd_unif(1-defocus_change/100,1+defocus_change/100);
      ctf.ctf.DeltafV*=rnd_unif(1-defocus_change/100,1+defocus_change/100);
      aux.init_zeros(Ydim,Xdim);
      xmippArray2VTK(aux, vtkaux,2);
      ctf.generate_mask(vtkaux);
      vtkaux->Delete();
      ctf.ctf.DeltafU=ctf.ctf.DeltafU;
      ctf.ctf.DeltafV=ctf.ctf.DeltafV;
      #ifdef DEBUG
         ctf.write_amplitude("PPP_particular.xmp");
         char c; cout << "Press any key\n"; cin >> c;
      #endif
   }

   // Apply CTF
   if (fn_ctf!="") ctf.apply_mask(I);

   // Add noise after CTF
   noisy.init_random(0,sigma_after_CTF,"gaussian");
   if (fn_after_ctf!="") after_ctf.apply_mask(noisy);
   I += noisy;
}
#undef DEBUG
#endif
