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

#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippHistograms.hh>
#include <XmippData/xmippGeometry.hh>
#include <XmippInterface/xmippSpider.hh>

void Usage();

int main (int argc, char **argv) {
   FileName fn_vol, fn_fr, fn_img, fn_original;
   bool apply, apply_to_original, omit_check;
   double rot0,  rotF,  rot_step;
   double tilt0, tiltF, tilt_step;
   double psi0,  psiF,  psi_step;
   double max_shift, shift_step;
   double first_ring, last_ring;
   char   wrapper_mode; // There are two wrapper modes
                        // R: Radon
                        // M: projection matching

   // Get input parameters -------------------------------------------------
   try {
      if (check_param(argc,argv,"-apmq")) wrapper_mode='M';
      else wrapper_mode='R';
      fn_vol=get_param(argc,argv,"-i");
      fn_img=get_param(argc,argv,"-img");
      apply=check_param(argc,argv,"-apply");
      apply_to_original=check_param(argc,argv,"-apply_to_original");
      omit_check=check_param(argc,argv,"-omit_check");
      if (apply_to_original)
         fn_original=get_param(argc,argv,"-apply_to_original");
      max_shift=AtoF(get_param(argc,argv,"-max_shift","2"));
      
      if (check_param(argc,argv,"-rot")) {
         int i=position_param(argc,argv,"-rot");
         if (i+3>=argc)
            REPORT_ERROR(1,
               "Angular_refinement: Not enough parameters after -rot");
         rot0=AtoF(argv[i+1]);
         rotF=AtoF(argv[i+2]);
         rot_step=AtoF(argv[i+3]);
      } else {rot0=0; rotF=360; rot_step=10;}
      if (check_param(argc,argv,"-tilt")) {
         int i=position_param(argc,argv,"-tilt");
         if (i+3>=argc)
            REPORT_ERROR(1,
               "Angular_refinement: Not enough parameters after -tilt");
         tilt0=AtoF(argv[i+1]);
         tiltF=AtoF(argv[i+2]);
         tilt_step=AtoF(argv[i+3]);
      } else {tilt0=0; tiltF=180; tilt_step=10;}
      if (check_param(argc,argv,"-psi")) {
         int i=position_param(argc,argv,"-psi");
         if (i+3>=argc)
            REPORT_ERROR(1,
               "Angular_refinement: Not enough parameters after -psi");
         psi0=AtoF(argv[i+1]);
         psiF=AtoF(argv[i+2]);
         psi_step=AtoF(argv[i+3]);
      } else {psi0=0; psiF=360; psi_step=10;}

      if (wrapper_mode=='R') {
         fn_fr =get_param(argc,argv,"-fr");
      } else {
         shift_step=AtoF(get_param(argc,argv,"-shift_step","1"));
         first_ring=AtoI(get_param(argc,argv,"-first_ring","0"));
         last_ring= AtoI(get_param(argc,argv,"-last_ring","-1"));
         tilt_step= AtoF(get_param(argc,argv,"-tilt_step","10"));
      }
   } catch (Xmipp_error XE) {cout << XE; Usage(); exit(0);}

   // Perform refinement ---------------------------------------------------
   try {
      FileName fn_root=fn_img.without_extension();
      FileName fn_iter;
      int i=fn_root.find("_iter_");
      int refinement_step=1;
      if (i!=-1) {
         refinement_step=AtoI(fn_root.substr(i+6,2))+1;
         fn_root=fn_root.substr(0,i);
      }
      fn_iter="_iter_"+ItoA(refinement_step,2);

      // Get original angles ...............................................
      DocFile angles_in;
      SelFile SF_img, SF_original;
      SF_img.read(fn_img);
      if (apply_to_original) SF_original.read(fn_original);
      extract_angles(SF_img,angles_in);
      angles_in.go_beginning();
      angles_in.remove_current(); // Remove angles comment
      SF_img.go_first_ACTIVE();

      // Perform refinement ................................................
      if (wrapper_mode=='R') {
         Angular_refinement_Radon(fn_vol, fn_fr, fn_root+fn_iter+"_orfsfs",
           rot0, rotF, rot_step, tilt0, tiltF, tilt_step,
           psi0, psiF, psi_step, max_shift);
         system(((string)"mv "+fn_root+fn_iter+"_orfsfs."+fn_vol.get_extension()+" "+
            fn_root+fn_iter+"_orfsfs.txt").c_str());
      } else {
         Angular_refinement_Matching(fn_vol, fn_img, fn_root+fn_iter+"_orfsfs",
           tilt_step, max_shift, shift_step, first_ring, last_ring);
      }

      // Get report ........................................................
      DocFile report;
      report.read(fn_root+fn_iter+"_orfsfs.txt",TRUE);
      
      // Process report ....................................................
      matrix1D<double> rot_diff, tilt_diff, psi_diff, shift_diff;
      rot_diff  .resize(report.dataLineNo());
      tilt_diff .resize(rot_diff);
      psi_diff  .resize(rot_diff);
      shift_diff.resize(rot_diff);
      
      report.go_first_data_line();
      angles_in.go_first_data_line();
      angles_in.insert_comment("Refinement summary");
      angles_in.insert_comment(
         "old_rot old_tilt old_psi corr new_rot new_tilt new_psi new_X new_Y "
         "rot_diff tilt_diff psi_diff shift_diff");

      SelFile SF_out;
      FOR_ALL_ELEMENTS_IN_MATRIX1D(rot_diff) {
         double corr=report(0);
         double rot =report(1);
         double tilt=report(2);
         double psi =report(3);
         double X   =report(4);
         double Y   =report(5);
         rot_diff(i)  =realWRAP(angles_in(0)-rot,-180,180);
         tilt_diff(i) =realWRAP(angles_in(1)-tilt,-180,180);
         psi_diff(i)  =realWRAP(angles_in(2)-psi,-180,180);
         shift_diff(i)=sqrt(X*X+Y*Y);

         // Check that the new angles really meets the 
         // subsearch constraints
         if (!omit_check &&
             (0.9*rot_diff(i) <rot0  || 0.9*rot_diff(i) >rotF  ||
              0.9*tilt_diff(i)<tilt0 || 0.9*tilt_diff(i)>tiltF ||
              0.9*psi_diff(i) <psi0  || 0.9*psi_diff(i) >psiF)) {
            double other_rot, other_tilt, other_psi;
            Euler_another_set(rot,tilt,psi,
               other_rot, other_tilt, other_psi);
            rot_diff(i)  =realWRAP(angles_in(0)-other_rot,-180,180);
            tilt_diff(i) =realWRAP(angles_in(1)-other_tilt,-180,180);
            psi_diff(i)  =realWRAP(angles_in(2)-other_psi,-180,180);
            if (0.9*rot_diff(i) <rot0  || 0.9*rot_diff(i) >rotF  ||
                0.9*tilt_diff(i)<tilt0 || 0.9*tilt_diff(i)>tiltF ||
                0.9*psi_diff(i) <psi0  || 0.9*psi_diff(i) >psiF) {
               cout << "Image " << SF_img.get_current_file() << " moves from "
                    << "(" << angles_in(0) << "," << angles_in(1) << "," << angles_in(2) << ") to "
                    << "(" << rot << "," << tilt << "," << psi << ") "
                    << "that is too much for the subsearch constraints\n"
                    << "Angles are kept as the originals\n"
                    << "The Euler matrix for the original set is\n";
               matrix2D<double> E;
               Euler_angles2matrix(angles_in(0),angles_in(1),angles_in(2),E);
               cout << E << endl << "While for the new one\n";
               Euler_angles2matrix(rot,tilt,psi,E);
               cout << E << endl;

               rot_diff(i)  =0; rot=angles_in(0);
               tilt_diff(i) =0; tilt=angles_in(1);
               psi_diff(i)  =0; psi=angles_in(2);
               shift_diff(i)=0; X=Y=0;
               corr=0;
            } else {
               rot=other_rot;
               tilt=other_tilt;
               psi=other_psi;
            }
         }

         EULER_CLIPPING(rot,tilt,psi);
         if (apply) {
            ImageXmipp I;
            I.read(SF_img.get_current_file());
            I.set_eulerAngles(rot,tilt,psi);
            I().translate(vector_R2(X,Y));

            FileName fn_out=I.name().get_root();
            int n=fn_out.find("_iter_");
            if (n!=-1) fn_out=fn_out.substr(0,n);

            FileName fn_output;
            fn_output.compose(fn_out+fn_iter+"_",I.name().get_number(),
               I.name().get_extension());
            I.write(fn_output);
            SF_out.insert(fn_output);
         }

         if (apply_to_original) {
            ImageXmipp I;
            I.read(SF_original.get_current_file());
            I.set_eulerAngles(rot,tilt,psi);
            I().translate(vector_R2(X,Y));
            I.write();
            SF_original.NextImg();
         }

         angles_in.set( 3,corr);
         angles_in.set( 4,rot);
         angles_in.set( 5,tilt);
         angles_in.set( 6,psi);
         angles_in.set( 7,X);
         angles_in.set( 8,Y);
         angles_in.set( 9,rot_diff(i));
         angles_in.set(10,tilt_diff(i));
         angles_in.set(11,psi_diff(i));
         angles_in.set(12,shift_diff(i));

         report.next_data_line();
         angles_in.next_data_line();
         SF_img.NextImg();
      }
      if (apply) SF_out.write(fn_root+fn_iter+".sel");

      // Effectively report ................................................
      angles_in.write(fn_root+fn_iter+"_summary.txt");
      double avg,stddev,min_val,max_val;
      rot_diff.compute_stats(avg,stddev,min_val,max_val);
      cout << "Rotational changes: " << avg << "+-" << stddev
           << " [" << min_val << "," << max_val << "]" << endl;
      tilt_diff.compute_stats(avg,stddev,min_val,max_val);
      cout << "Tilting changes:    " << avg << "+-" << stddev
           << " [" << min_val << "," << max_val << "]" << endl;
      psi_diff.compute_stats(avg,stddev,min_val,max_val);
      cout << "Psi changes:        " << avg << "+-" << stddev
           << " [" << min_val << "," << max_val << "]" << endl;
      shift_diff.compute_stats(avg,stddev,min_val,max_val);
      cout << "Shift changes:      " << avg << "+-" << stddev
           << " [" << min_val << "," << max_val << "]" << endl;

      // Compute histograms ................................................
      histogram1D hist;

      compute_hist(rot_diff,hist,50);
      hist.write(fn_root+fn_iter+"_rot_histogram.txt");
      
      compute_hist(tilt_diff,hist,50);
      hist.write(fn_root+fn_iter+"_tilt_histogram.txt");
      
      compute_hist(psi_diff,hist,50);
      hist.write(fn_root+fn_iter+"_psi_histogram.txt");
      
      compute_hist(shift_diff,hist,50);
      hist.write(fn_root+fn_iter+"_shift_histogram.txt");
   } catch (Xmipp_error XE) {cout << XE;}
}

// Usage -------------------------------------------------------------------
void Usage() {
   cerr << "Usage: angular_refinement\n"
        << "If Radon wrapper --------------------------------------------------------------------\n"
        << "   -i <FR Volume>                        : Fourier Radon Volume\n"
        << "   -fr <FR selfile>                      : Fourier Radon images\n"
        << "   -img <selfile>                        : Real space images\n"
        << "  [-rot  <rot0=0>  <rotF=360> <step=10>] : Rotational angle range\n"
        << "  [-tilt <tilt0=0> <tiltF=180> <step=10>]: Tilting angle range\n"
        << "  [-psi  <psi0=0>  <psiF=360> <step=10>] : In-plane angle range\n"
        << "  [-max_shift <s=2>]                     : in pixels\n"
        << "If Projection Matching wrapper ------------------------------------------------------\n"
        << "   -apmq                                 : This flag sets the Projection Matching mode\n"
        << "   -i <Xmipp Volume>                     : referenece volume\n"
        << "   -img <selfile>                        : Real space images\n"
        << "  [-tilt_step <step=10>]                 : Tilting angle range\n"
        << "  [-max_shift <s=2>]                     : in pixels\n"
        << "  [-shift_step <s=1>]                    : max_shift must be a multiple\n"
        << "                                           of shift_step\n"
        << "  [-first ring <r=0>]                    : First ring to evaluate\n"
        << "  [-last_ring <r=-1>]                    : Last ring to evaluate\n"
        << "Both --------------------------------------------------------------------------------\n"
        << "  [-omit_check]                          : Check that the angles are within\n"
        << "                                           the provided limits\n"
        << "  [-apply]                               : Apply shift and angles to input images\n"
        << "  [-apply_to_original <selfile>]         : Apply shift and angles to original images.\n"
        << "                                           They are rewritten\n";
   ;
}
