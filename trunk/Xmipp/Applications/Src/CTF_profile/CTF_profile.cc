/***************************************************************************
 *
 * Authors:     Carlos Oscar Sánchez Sorzano (coss@cnb.uam.es)
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
#include <Reconstruction/CTF.hh>
#include <XmippData/xmippArgs.hh>

void Usage();

/* Main -------------------------------------------------------------------- */
int main(int argc, char **argv) {
   FileName         fn_ctf;
   matrix1D<double> w_dir;
   double           w_step;

   // Read parameters
   try {
      fn_ctf=get_param(argc,argv,"-i");
      if (check_param(argc,argv,"-dir")) 
         w_dir=get_vector_param(argc,argv,"-w_dim",2);
      else w_dir=vector_R2(1,0);
      w_step=AtoF(get_param(argc,argv,"-w_step","0.001"));
   } catch (Xmipp_error XE) {cout << XE; Usage(); exit(-1);}

   // Show parameters
   cout << "# File: " << fn_ctf << endl
        << "# Direction: " << w_dir.transpose() << endl;

   // Show profile
   try {
      // Read CTF
      XmippCTF CTF;
      CTF.enable_CTF=CTF.enable_CTFnoise=true;
      CTF.read(fn_ctf);
      CTF.Produce_Side_Info();
      
      // Compute CTF at each frequency
      cout << "# Dig.Freq Cont.Freq Cont.Freq(A) CTF_pure CTF_noise CTF_total\n";
      w_dir/=w_dir.module();
      for (double w=0; w<=0.5; w+=w_step) {
         matrix1D<double> current_w=w/CTF.Tm*w_dir;
         double CTF_pure=CTF.CTFpure_at(XX(current_w),YY(current_w));
         double CTF_noise=CTF.CTFnoise_at(XX(current_w),YY(current_w));
         double cont_freq=current_w.module();
         double cont_freq_A=(cont_freq==0) ? 0:1/cont_freq;
         cout << w << " " << cont_freq << " " << cont_freq_A << " "
              << CTF_pure << " " << CTF_noise << " " << CTF_pure*CTF_pure+CTF_noise
              << endl;
      }
   } catch (Xmipp_error XE) {cout << XE;}
}

// Usage -------------------------------------------------------------------
void Usage() {
   cerr << "Usage: CTF_profile [options]\n"
        << "   -i <CTF description file>      : CTF to plot\n"
        << "  [-w_dir \"[X=1,Y=0]\"             : Frequency direction\n"
        << "  [-w_step <step=0.001>]          : Normalized to 0.5\n";
}
