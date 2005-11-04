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

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippInterface/xmippSpider.hh>
#include <Reconstruction/radon.hh>

class Radon_parameters: public Prog_parameters {
public:
   double Delta_rot;
   double Delta_tilt;
   int    output_size;
   bool   Fourier;
   double low_pass;
   double temperature;
   bool   use_Xmipp;

   void read(int argc, char **argv) {
      Prog_parameters::read(argc,argv);
      Delta_rot   = AtoF(get_param(argc,argv,"-delta_rot","2"));
      Delta_tilt  = AtoF(get_param(argc,argv,"-delta_tilt","2"));
      output_size = AtoI(get_param(argc,argv,"-output_size","-1"));
      Fourier     =    check_param(argc,argv,"-fourier");
      low_pass    = AtoF(get_param(argc,argv,"-low_pass","0"));
      temperature = AtoF(get_param(argc,argv,"-temperature","0.2"));
      use_Xmipp   = check_param(argc,argv,"-use_Xmipp");

      if (Fourier && low_pass==0)
         REPORT_ERROR(1,"Radon_transform: Some low_pass parameter needed");
   }
   
   void show() {
      Prog_parameters::show();
      cout << "Delta_rot      = " << Delta_rot   << endl
           << "Delta_tilt     = " << Delta_tilt  << endl
           << "Output size    = " << output_size << endl
           << "Fourier output = " << Fourier     << endl
           << "Low pass       = " << low_pass    << endl
           << "Temperature    = " << temperature << endl
	   << "Use Xmipp      = " << use_Xmipp   << endl
      ;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "  [-delta_rot  <D=2>]       : Sampling rate for rotational angle\n"
           << "  [-delta_tilt <D=2>]       : Sampling rate for tilting angle\n"
           << "                              Only valid for 3D\n"
           << "  [-output_size <s=-1>]     : Radon output size\n"
           << "                              By default, 1.5*SIZE(input)\n"
           << "  [-fourier]                : Produce Fourier outputs\n"
           << "  [-low_pass <w=0>]         : Normalized to 0.5\n"
           << "  [-temperature <t=0.02>]   : Fermi filter paramter\n"
	   << "  [-use_Xmipp]              : Only for images, without FFT\n"
      ;
   }
};

bool process_img(ImageXmipp &img, const FileName &fn_out,
   const Prog_parameters *prm) {
   Radon_parameters *eprm=(Radon_parameters *) prm;
   if (!eprm->use_Xmipp) {
      // Do it via Spider
      radon_transform(img, fn_out, eprm->Delta_rot, eprm->output_size);
      if (eprm->Fourier)
	 Fourier_transform_of_Radon_transform(fn_out,fn_out,
            eprm->low_pass, eprm->temperature);
   } else {
      // Do it via Xmipp
      matrix2D<double> RT;
      Radon_Transform(img(),eprm->Delta_rot,RT);
      img()=RT;
      img.write(fn_out);
   }
   return true;
}

bool process_vol(VolumeXmipp &vol, const FileName &fn_out,
   const Prog_parameters *prm) {
   Radon_parameters *eprm=(Radon_parameters *) prm;
   radon_transform(vol, fn_out, eprm->Delta_rot, eprm->Delta_tilt,
      eprm->output_size);
   if (eprm->Fourier)
      Fourier_transform_of_Radon_transform(fn_out,fn_out,
         eprm->low_pass, eprm->temperature);
   return true;
}

int main (int argc, char **argv) {
   Radon_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol, IMAGE2FILE);
}
