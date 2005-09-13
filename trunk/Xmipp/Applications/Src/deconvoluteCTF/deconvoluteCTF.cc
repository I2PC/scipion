/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)
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
                                                                                               
/* INCLUDES ---------------------------------------------------------------- */
#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippFFT.hh>
#include <Reconstruction/CTF.hh>
#include <Reconstruction/Programs/Prog_FourierFilter.hh>


/* PROTOTYPES -------------------------------------------------------------- */
void Usage();
                                                                                              
class deconvoluteCTF_parameters: public Prog_parameters {
public:
  XmippCTF    ctf;
  FileName    fn_ctf;
  bool        first;
  int         iter;
  matrix2D< complex<double> > ctfmask;

  void read(int argc, char **argv) {
    Prog_parameters::read(argc,argv);
    fn_ctf=get_param(argc,argv,"-ctf");
    iter=AtoI(get_param(argc,argv,"-iter","1000"));
    ctf.read(fn_ctf);
    ctf.Produce_Side_Info();    
    first=TRUE;
  }
  void show() {
    Prog_parameters::show();
    cout << "CTF parameters= " << fn_ctf << endl;
    cout << "Number of iterations= " << iter << endl;
  }
  void usage() {
    Prog_parameters::usage();
    cerr << "   -ctf <filename>          : CTF-parameter file\n";
    cerr << "  [-iter <int>  ]           : number of iterations\n";
  }

};
 
bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   deconvoluteCTF_parameters  *eprm=(deconvoluteCTF_parameters *) prm;
   matrix2D< complex<double> > fft,Fnew,CCFnew;
   FourierTransform(img(), fft);
   fft.set_Xmipp_origin();
   if (eprm->first) {
     eprm->ctf.Generate_CTF(img().RowNo(),img().ColNo(),eprm->ctfmask); 
     eprm->ctfmask.set_Xmipp_origin();
     eprm->first=FALSE;
   }
   fft*=eprm->ctfmask;
   Fnew=fft;
   for (int i=0; i < eprm->iter; i++) {
     CCFnew=Fnew;
     CCFnew*=eprm->ctfmask;
     CCFnew*=eprm->ctfmask;
     Fnew+=fft-CCFnew;
   }
   InverseFourierTransform(Fnew,img());
   return TRUE;
}


void process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
}                                                                                                    

int main (int argc, char **argv) {
   deconvoluteCTF_parameters prm;
   prm.each_image_produces_an_output=TRUE;
   // Set default action for application of header transformation
   prm.apply_geo=FALSE;   
   SF_main(argc, argv, &prm, (void*)&process_img, NULL);
}
                                                                                                                    
