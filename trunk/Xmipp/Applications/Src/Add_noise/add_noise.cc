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

class Add_noise_parameters: public Prog_parameters {
public:
   double noise_min, noise_max;
   double noise_avg, noise_stddev;
   bool   gaussian,  uniform;

   void read(int argc, char **argv) {
      Prog_parameters::read(argc,argv);
      gaussian=uniform=false;
      if (check_param(argc,argv,"-gaussian")) {
         gaussian=true;
         int i=position_param(argc,argv,"-gaussian");
         if (i+1>=argc) REPORT_ERROR(1,"Not enough parameters after -gaussian");
         noise_stddev=AtoF(argv[i+1]);
         if (i+2<argc) {noise_avg=AtoF(argv[i+2]);} else noise_avg=0;
      } else if (check_param(argc,argv,"-uniform")) {
         uniform=true;
         int i=position_param(argc,argv,"-uniform");
         if (i+2>=argc) REPORT_ERROR(1,"Not enough parameters after -uniform");
         noise_min=AtoF(argv[i+1]);
         noise_max=AtoF(argv[i+2]);
      } else
         REPORT_ERROR(1,"Unknown noise type");
   }

   void show() {
      Prog_parameters::show();
      if (gaussian) 
         cout << "Noise avg=" << noise_avg << endl
              << "Noise stddev=" << noise_stddev << endl;
      else if (uniform)
         cout << "Noise min=" << noise_min << endl
              << "Noise max=" << noise_max << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "  [-gaussian <stddev> [<avg>=0]] : Gaussian noise parameters\n"
           << "  [-uniform  <min> <max>]   : Uniform noise parameters\n";
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
  Add_noise_parameters *eprm=(Add_noise_parameters *) prm;
  if (eprm->gaussian)
     img().add_noise(eprm->noise_avg, eprm->noise_stddev,"gaussian");
  else if (eprm->uniform)
     img().add_noise(eprm->noise_min, eprm->noise_max,"uniform");
   return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
  Add_noise_parameters *eprm=(Add_noise_parameters *) prm;
  if (eprm->gaussian)
     vol().add_noise(eprm->noise_avg, eprm->noise_stddev,"gaussian");
  else if (eprm->uniform)
     vol().add_noise(eprm->noise_min, eprm->noise_max,"uniform");
   return true;
}

int main (int argc, char **argv) {
   Add_noise_parameters prm;
   randomize_random_generator();
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Add_noise {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Add_noise/Help/add_noise.html";
      help="Add noise to volumes and images";
      OPEN MENU menu_add_noise;
      COMMAND LINES {
	+ usual: xmipp_add_noise
               #include "prog_line.mnu"
               [-gaussian $STDDEV [$AVG]]
               [-uniform  $MIN $MAX]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        OPT(-gaussian) {label="Add gaussian noise";}
           $STDDEV {label="Standard deviation"; type=float;}
           $AVG    {label="Average";            type=float;}
        OPT(-uniform) {label="Add uniform noise";}
           $MIN    {label="Minimum value";      type=float;}
           $MAX    {label="Maximum value";      type=float;}
      }
   }

   MENU menu_add_noise {
      #include "prog_menu.mnu"
      "Noise parameters"
      OPT(-gaussian)
      OPT(-uniform)
   }
*/
