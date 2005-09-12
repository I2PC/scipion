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

#include "../Prog_CorrectPhase.hh"
#include <XmippData/xmippArgs.hh>

/* Read parameters from command line. -------------------------------------- */
void CorrectPhase_Params::read(int argc, char **argv) {
   fn_ctf=get_param(argc,argv,"-ctf");
   
   epsilon=AtoF(get_param(argc,argv,"-small","0"));
   string aux; aux=get_param(argc,argv,"-method","");
   if      (aux=="remove")           method=CORRECT_SETTING_SMALL_TO_ZERO;
   else if (aux=="leave" || aux=="") method=CORRECT_LEAVING_SMALL;
   else if (aux=="divide")           method=CORRECT_AMPLIFYING_NOT_SMALL;
   multiple_CTFs=fn_ctf.get_extension()=="sel";
   cout << fn_ctf << " " << fn_ctf.get_extension() << endl;
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
   cerr << "   -ctf <CTF descr file or selfile> : It must not be centered\n"
        << "  [-small <epsilon=0>]              : Values under epsilon are small\n"
	<< "  [-method <mth=leave>]             : Valid methods are: remove, leave\n"
	<< "                                      divide\n";
   ;
}

/* Produce Side information ------------------------------------------------ */
void CorrectPhase_Params::produce_side_info() {
   ctf.FilterBand=CTF;
   ctf.ctf.enable_CTFnoise=FALSE;
   if (multiple_CTFs) {
      SF_CTF.read(fn_ctf);
      SF_CTF.go_first_ACTIVE();
   } else {
      ctf.ctf.read(fn_ctf);
      ctf.ctf.Produce_Side_Info();
   }
}

/* Correct a single image -------------------------------------------------- */
//#define DEBUG
void CorrectPhase_Params::correct(matrix2D< complex<double> > &v) {
   if (XSIZE(ctf.mask2D)==0 || multiple_CTFs)
      ctf.generate_mask(v);
   #ifdef DEBUG
      cout << "New image ----------------------------\n";
   #endif

   FOR_ALL_ELEMENTS_IN_MATRIX2D(v) {
     complex<double> m=ctf.mask2D(i,j);
     if (m.imag()!=0) 
	REPORT_ERROR(1,"CorrectPhase::correct: CTF is not real\n");
     #ifdef DEBUG
	cout << "CTF at (" << j << "," << i << ")="
	     << m << " Value there " << v(i,j);
     #endif
     switch (method) {
	case CORRECT_SETTING_SMALL_TO_ZERO:
	   if (m.real()<0)
	      if (v(i,j).real()<-epsilon) v(i,j)*=-1;
	      else                        v(i,j)=0;
	   break;
	case CORRECT_LEAVING_SMALL:
	   if (m.real()<-epsilon) v(i,j)*=-1;
	   break;
	case CORRECT_AMPLIFYING_NOT_SMALL:
	   if (ABS(m.real())>epsilon) v(i,j)/=m.real();
	   break;
     }
     #ifdef DEBUG
	cout << " Final value " << v(i,j) << endl;
     #endif
  }
}
#undef DEBUG

/* Correct a set of images ------------------------------------------------- */
void CorrectPhase_Params::correct(SelFile &SF) {
   matrix2D< complex<double> > fft;
   SF.go_first_ACTIVE();
   cerr << "Correcting CTF phase ...\n";
   int istep=CEIL((double)SF.ImgNo()/60.0);
   init_progress_bar(SF.ImgNo());
   int i=0;
   while (!SF.eof()) {
      ImageXmipp I;
      I.read(SF.NextImg());
      FourierTransform(I(),fft);
      correct(fft);
      InverseFourierTransform(fft,I());
      I.write();
      if (i++%istep==0) progress_bar(i);
   }
   progress_bar(SF.ImgNo());
}
