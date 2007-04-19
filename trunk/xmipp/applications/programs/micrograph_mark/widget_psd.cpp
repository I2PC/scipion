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

#include <graphics/show_2d.h>
#include <graphics/show_selfile.h>

#include "widget_psd.h"

// Set the assign CTF file -------------------------------------------------
void QtWidgetPSD::set_assign_CTF_file(Micrograph &m,
   const FileName &_fn_assign_CTF) {

   // Read the input assign CTF file
   fn_assign_CTF=_fn_assign_CTF;
   assign_ctf_prm.read(fn_assign_CTF);
   FileName fn_root=assign_ctf_prm.image_fn.remove_all_extensions();

   // Check if the CTF is computed at each particle
   if (assign_ctf_prm.compute_at_particle) {
      cerr << "QtWidgetPSD::set_assign_CTF_file: The PSDs and CTFs cannot be shown "
              "for individual particles. See file " << fn_assign_CTF << endl;
      return;
   }

   // Generate a random selfile
   FileName fn_random; fn_random.init_random(15);
   fn_random=(string)"PPP"+fn_random+".sel";
   files_to_remove.push_back(fn_random);

   // To shape the selfile
   int div_NumberX, div_NumberY;

   // Check if it is a micrograph average
   if (assign_ctf_prm.micrograph_averaging) {
      FileName fn_avg;
      if (!ctf_mode) {
         if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
              fn_avg=fn_root+"_ARMAavg.psd";
         else fn_avg=fn_root+"_Periodogramavg.psd";
      } else {
         if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
              fn_avg=fn_root+"_ARMAavg.ctfmodel";
         else fn_avg=fn_root+"_Periodogramavg.ctfmodel";
      }
      system(((string)"xmipp_do_selfile "+fn_avg+" > "+fn_random).c_str());
      div_NumberX=div_NumberY=1;
   } else {
      // The CTF must have been computed by pieces
      FileName PSDfn_root;
      if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
           PSDfn_root=fn_root+"_ARMA";
      else PSDfn_root=fn_root+"_Periodogram";
      string command;
      if (!ctf_mode)
         command=(string)"xmipp_do_selfile \""+PSDfn_root+"?????.psd\" > "+fn_random;
      else
         command=(string)"xmipp_do_selfile \""+PSDfn_root+"?????.ctfmodel\" > "+fn_random;
      system(command.c_str());

      int Ydim, Xdim; m.size(Xdim, Ydim);
      div_NumberX=CEIL((double)Xdim/assign_ctf_prm.N_horizontal);
      div_NumberY=CEIL((double)Ydim/assign_ctf_prm.N_vertical  );
   }

   // Show the selfile
   ShowSel *showsel = new ShowSel;
   showsel->apply_geo=false;
   showsel->showonlyactive = false;
   if (!ctf_mode)
      showsel->initWithFile(div_NumberY, div_NumberX, fn_random.c_str(), ShowSel::PSD_mode);
   else {
      showsel->initWithFile(div_NumberY, div_NumberX, fn_random.c_str(), ShowSel::CTF_mode);
      showsel->setAssignCTFfile(fn_assign_CTF);
   }
   showsel->show();
}

QtWidgetPSD::~QtWidgetPSD() {
   for (int i=0; i<files_to_remove.size(); i++)
      system(((string)"rm "+files_to_remove[i]).c_str());
}
