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
   #include "../Prog_CorrectPhase.hh"
   #include <XmippData/xmippArgs.hh>

/* Read parameters from command line. -------------------------------------- */
void CorrectPhase_Params::read(int argc, char **argv) {
   if (check_param(argc,argv,"-ctf_descr")) {
      fn_ctf=get_param(argc,argv,"-ctf_descr");
      CTF_description_file=TRUE;
   } else {
      fn_ctf=get_param(argc,argv,"-ctf");
      CTF_description_file=FALSE;
   }
   
   epsilon=AtoF(get_param(argc,argv,"-small","0"));
   string aux; aux=get_param(argc,argv,"-method","");
   if      (aux=="remove")           method=CORRECT_SETTING_SMALL_TO_ZERO;
   else if (aux=="leave" || aux=="") method=CORRECT_LEAVING_SMALL;
   else if (aux=="divide")           method=CORRECT_AMPLIFYING_NOT_SMALL;
   multiple_CTFs=!Is_FourierImageXmipp(fn_ctf);
}

/* Show -------------------------------------------------------------------- */
void CorrectPhase_Params::show() {
   cout << "CTF: " << fn_ctf << endl
        << "Small is under " << epsilon << endl;
   cout << "Correcting method: ";
   switch (method) {
      case CORRECT_SETTING_SMALL_TO_ZERO:
         cout << "Set small values to 0\n"; break;
      case CORRECT_LEAVING_SMALL:
         cout << "Leave small values as they are\n"; break;
      case CORRECT_AMPLIFYING_NOT_SMALL:
         cout << "Correct amplitude except for the small values\n";
         break;
   }
}

/* Usage ------------------------------------------------------------------- */
void CorrectPhase_Params::usage() {
   cerr << "   -ctf <Xmipp Fourier file or selfile>: It must not be centered\n"
        << "  [-small <epsilon=0>]      : Values under epsilon are small\n"
	<< "  [-method <mth=leave>]     : Valid methods are: remove, leave\n"
	<< "                              divide\n";
   ;
}

/* Produce Side information ------------------------------------------------ */
void CorrectPhase_Params::produce_side_info() {
   if (!multiple_CTFs && !CTF_description_file) {
      ctf.FilterShape=ctf.FilterBand=FROM_FILE;
      ctf.fn_mask=fn_ctf;
      ctf.generate_mask(NULL);
   } else if (multiple_CTFs) {
      ctf.FilterShape=ctf.FilterBand=FROM_FILE;
      SF_CTF.read(fn_ctf);
   } else {
      ctf.FilterBand=CTF;
      ctf.ctf.read(fn_ctf);
      ctf.ctf.Produce_Side_Info();
   }
}

/* Look for a CTF file ----------------------------------------------------- */
FileName CorrectPhase_Params::CTF_filename(const FileName &fn) {
   SF_CTF.go_first_ACTIVE();
   FileName searched_fn_root=fn.remove_all_extensions();
   searched_fn_root=searched_fn_root.remove_directories();
   while (!SF_CTF.eof()) {
      FileName fn_CTF=SF_CTF.NextImg();
      FileName fn_root=fn_CTF.remove_until_prefix("ctf-");
      fn_root=fn_root.remove_all_extensions();
      if (fn_root==searched_fn_root) return fn_CTF;
   }
   return "";
}

/* Correct a single image -------------------------------------------------- */
//#define DEBUG
void CorrectPhase_Params::correct(vtkImageData *v) _THROW {
   if (ctf.mask==NULL && CTF_description_file)
      ctf.generate_mask(v); // This is the first time
                            // that the CTF is applied
   if (!same_shape(v,ctf.mask))
      REPORT_ERROR(1,"CorrectPhase::correct: ctf and input FT do not have"
         " the same shape");
   if (v->GetScalarType()!=VTK_FLOAT || v->GetNumberOfScalarComponents()!=2)
      REPORT_ERROR(1,"CorrectPhase::correct: input array do not seem to be"
         " a Fourier Transform");
   #ifdef DEBUG
      cout << "New image ----------------------------\n";
   #endif
   
   int dim[3]; v->GetDimensions(dim);
   float *mi=(float *) ctf.mask->GetScalarPointer();
   float *vi=(float *)        v->GetScalarPointer();
   for (int k=0; k<dim[2]; k++)
       for (int i=0; i<dim[1]; i++)
      	   for (int j=0; j<dim[0]; j++) {
	      if (*(mi+1)!=0) 
	         REPORT_ERROR(1,"CorrectPhase::correct: CTF is not real\n");
	      #ifdef DEBUG
	         cout << "CTF at (" << j << "," << i << "," << k << ")="
		      << *mi << " Value there " << *vi << "+i" << *(vi+1);
	      #endif
	      switch (method) {
		 case CORRECT_SETTING_SMALL_TO_ZERO:
		    if (*mi<0)
		       if (*mi<-epsilon) {*vi=-(*vi); *(vi+1)=-(*(vi+1));}
		       else              {*vi=*(vi+1)=0;}
		    break;
		 case CORRECT_LEAVING_SMALL:
		    if (*mi<-epsilon) {*vi=-(*vi); *(vi+1)=-(*(vi+1));}
		    break;
		 case CORRECT_AMPLIFYING_NOT_SMALL:
		    if (ABS(*mi)>epsilon) {(*vi)/=(*mi); (*(vi+1))/=(*mi);}
		    break;
	      }
	      #ifdef DEBUG
	         cout << " Final value " << *vi << "+i" << *(vi+1) << endl;
	      #endif
	      vi+=2; mi+=2;
	   }
}
#undef DEBUG

/* Correct a set of images ------------------------------------------------- */
void CorrectPhase_Params::correct(SelFile &SF) {
   vtkImageData * fft=NULL;
   SF.go_first_ACTIVE();
   cerr << "Correcting CTF phase ...\n";
   int istep=CEIL((double)SF.ImgNo()/60.0);
   init_progress_bar(SF.ImgNo());
   int i=0;
   while (!SF.eof()) {
      ImageXmipp I;
      I.read(SF.NextImg());
      FFT_VTK(I(),fft,TRUE);
      correct(fft);
      IFFT_VTK(fft,I(),TRUE);
      I.write();
      if (i++%istep==0) progress_bar(i);
   }
   progress_bar(SF.ImgNo());
   if (fft!=NULL) fft->Delete();
}
#endif
