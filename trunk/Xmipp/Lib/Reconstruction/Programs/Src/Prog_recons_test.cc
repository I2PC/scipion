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

#include "../Prog_recons_test.hh"

#include "../Prog_random_phantom.hh"
#include "../Prog_project.hh"
#include "../Prog_art.hh"
#include "../Prog_surface.hh"
#include "../Prog_symmetrize.hh"
#include "../Prog_normalize.hh"
#include "../Prog_range_adjust.hh"
#include "../Prog_FourierFilter.hh"
#include "../../recons_spider.hh"
#include "../Prog_CorrectPhase.hh"
#include "../Prog_IDR_art.hh"
#include "../Prog_microscope.hh"
#include <XmippData/Programs/Prog_segment.hh>
#include <XmippData/xmippMorphology.hh>

/* Read Reconstruction test parameters from file =========================== */
void Recons_test_Parameters::read(const FileName &fn_test_params) _THROW {
   FILE    *fh_param;
   string  str;
   char *  auxstr;

   // Clean old lambda and no_it
   lambda0.clear();
   lambdaF.clear();
   no_it0.clear();
   no_itF.clear();
   only_structural=FALSE;

   // Open file
   if ((fh_param = fopen(fn_test_params.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,
         (string)"Recons_test_Parameters::read: There is a problem "
         "opening the file "+fn_test_params);

   // Read parameters
   try {
      // Reconstruction method
      str=get_param(fh_param,"reconstruction method",0,NULL,
         3007,"Recons_test_Parameters::read: Reconstruction method not found");
      if      (str=="ART")         recons_method=use_ART;
      else if (str=="SIRT")        recons_method=use_SIRT;
      else if (str=="WBP")         recons_method=use_WBP;
      else if (str=="SIRT_Spider") recons_method=use_SIRT_Spider;
      else
         REPORT_ERROR(3007,(string)"Recons_test_Parameters::read: "
            "reconstruction mode " + str + " not supported");
      random_sort=check_param(fh_param,"random sort");
      sort_last_N=AtoI(get_param(fh_param,"sort last",0,"2"));
         
      // Several filenames and parameters
      fn_random_phantom=get_param(fh_param,"phantom family",0,"",
         3007,"Recons_test_Parameters::read: Random Phantom filename not found");
      fn_proj_params=get_param(fh_param,"projection parameters",0,NULL,
         3007,"Recons_test_Parameters::read: Projection parameters "
            "filename not found");
      fn_voxel_phantom=get_param(fh_param,"voxel phantom",0,"");
      fn_crystal=get_param(fh_param,"crystal parameters",0,"");
      fn_sym    =get_param(fh_param,"symmetry file",0,"");
      force_sym =AtoI(get_param(fh_param,"force symmetry",0,"0"));
      do_not_use_symproj =check_param(fh_param,"no projsym");
      fn_final_sym=get_param(fh_param,"final symmetry file",0,"");
      fn_CTF=get_param(fh_param,"CTF",0,"");
      defocus_change=AtoF(get_param(fh_param,"defocus change",0,"0"));
      sigma=AtoF(get_param(fh_param,"noise stddev",0,"0"));
      low_pass_before_CTF=AtoF(get_param(fh_param,"noise lowpass before CTF",0,"0"));
      fn_after_CTF=get_param(fh_param,"noise spectrum after CTF",0,"");

      MeasNo=AtoI(get_param(fh_param,"measurement number",0,"-1"));
      accuracy=AtoF(get_param(fh_param,"accuracy",0,"-1"));
      unluckiness=AtoF(get_param(fh_param,"unluckiness",0,"0.01"));
      global_radius=AtoF(get_param(fh_param,"global radius",0,"-1"));
      max_resolution=AtoF(get_param(fh_param,"max resolution",0,"-1"));

      // Surface mask
      probe_radius=AtoF(get_param(fh_param,"surface top",0,"0.5"));
      str=get_param(fh_param,"surface top",0,"");
      if (str!="") {
         enable_top_surface=TRUE;
         top0=AtoF(first_word(str),3007,
            "Recons_test_Parameters::read: top0 is not a true number");
         auxstr=next_token();
         // Is it a range?
         if (auxstr==NULL) topF=top0;
         else 
            topF=AtoF(auxstr,3007,
               "Recons_test_Parameters::read: topF is not a true number");
      } else enable_top_surface=FALSE;
      
      str=get_param(fh_param,"surface bottom",0,"");
      if (str!="") {
         enable_bottom_surface=TRUE;
         bottom0=AtoF(first_word(str),3007,
            "Recons_test_Parameters::read: bottom0 is not a true number");
         auxstr=next_token();
         // Is it a range?
         if (auxstr==NULL) bottomF=bottom0;
         else 
            bottomF=AtoF(auxstr,3007,
               "Recons_test_Parameters::read: bottomF is not a true number");
      } else enable_bottom_surface=FALSE;

      run_also_without_constraints=check_param(fh_param,"run also without constraints");

      // Normalization
      enable_normalization=check_param(fh_param,"enable normalization");
      if (enable_normalization) {
      	 a_avg=AtoF(get_param(fh_param,"a avg"));
      	 a_stddev=AtoF(get_param(fh_param,"a stddev"));
      	 b_avg=AtoF(get_param(fh_param,"b avg"));
      	 b_stddev=AtoF(get_param(fh_param,"b stddev"));
	 str=get_param(fh_param,"normalizing method");
	 if      (str=="OldXmipp")      normalizing_method=OLDXMIPP;
	 else if (str=="Near_OldXmipp") normalizing_method=NEAR_OLDXMIPP;
	 else if (str=="NewXmipp")      normalizing_method=NEWXMIPP;
	 else if (str=="NewXmipp2")     normalizing_method=NEWXMIPP2;
	 else if (str=="Michael")       normalizing_method=MICHAEL;
	 else if (str=="None")          normalizing_method=NONE;
	 else REPORT_ERROR(1,"Normalize: Unknown normalizing method");
	 bg_radius=AtoI(get_param(fh_param,"background radius",0,"0"));
      }

      // CTF correction
      correct_phase=check_param(fh_param,"CTF phase method");
      str=get_param(fh_param,"-method",0,"");
      if      (str=="remove")           phase_correction_method=CORRECT_SETTING_SMALL_TO_ZERO;
      else if (str=="leave" || str=="") phase_correction_method=CORRECT_LEAVING_SMALL;
      else if (str=="divide")           phase_correction_method=CORRECT_AMPLIFYING_NOT_SMALL;
      phase_correction_param=AtoF(get_param(fh_param,"CTF phase small",0,"0"));
      correct_amplitude=check_param(fh_param,"correct CTF amplitude");
      idr_iterations=AtoI(get_param(fh_param,"idr_iterations",0,"1"));
      mu0_list=get_vector_param(fh_param,"mu",-1);
      if (XSIZE(mu0_list)==0)
	 {mu0_list.resize(1); mu0_list.init_constant(1);}
      muF_list=get_vector_param(fh_param,"muF",-1);
      unmatched=check_param(fh_param,"unmatched");

      // Only valid for ART and SIRT
      str=get_param(fh_param,"blob type",0,"big");
      if      (str=="big")   blob_type=BIG_BLOB;
      else if (str=="small") blob_type=SMALL_BLOB;
      else    REPORT_ERROR(3007,
         "Recons_test_Parameters::read: unknown blob type, valid types big or small");
      stop_at=AtoI(get_param(fh_param,"stop at",0,"0"));
      succesive_params=check_param(fh_param,"succesive parameters");
      POCS_positivity=check_param(fh_param,"POCS positivity");
      reconstruction_radius=AtoF(get_param(fh_param,"reconstruction radius",0,"-1"));

      // Segmented surface
      enable_segmented_surface=check_param(fh_param,"segmented surface");
      if (enable_segmented_surface)
         threshold_surface_segment=AtoF(get_param(fh_param,"segmented surface",0));

      // Starting volume
      start_from_phantom=check_param(fh_param,"start from phantom");
      if (start_from_phantom) {
         starting_low_pass=AtoF(get_param(fh_param,"starting lowpass",0,"",
            3007,"Recons_test_Parameters::read: starting lowpass is missing"));
      }

      segmented_dilation=AtoI(get_param(fh_param,"segmented dilation",0,"0"));
      mass=AtoF(get_param(fh_param,"mass",0,"-1"));

      // If ART ..., read Iterative parameters
      if (recons_method==use_ART || recons_method==use_SIRT ||
	  recons_method==use_SIRT_Spider ) {
         int skip=0;
         do {
            str=get_param(fh_param,"iterative parameters",skip,"");
            if (str!="") {
               lambda0.push_back(AtoF(first_word(str),3007,
                  "Recons_test_Parameters::read: lambda0 is not a true number"));
               auxstr=next_word(3007,"Recons_test_Parameters::read: "
                  "number of iterations not found");
               no_it0.push_back(AtoI(auxstr,3007,
                  "Recons_test_Parameters::read: no_it0 is not a true number"));
               auxstr=next_token();
               // Is it a range?
               if (auxstr==NULL) {
                  lambdaF.push_back(lambda0.back());
                  no_itF. push_back(no_it0. back());
               } else {
                  lambdaF.push_back(AtoF(auxstr,3007,
                     "Recons_test_Parameters::read: lambdaF is not a true number"));
                  auxstr=next_word(3007,"Recons_test_Parameters::read: "
                     "number of iterations not found");
                  no_itF.push_back(AtoI(auxstr,3007,
                     "Recons_test_Parameters::read: no_itF is not a true number"));
               }
               skip++;
            } else if (skip==0)
               REPORT_ERROR(3007,"Recons_test_Parameters::read: There are"
                  " no iterative parameters");
            else break;
         } while (TRUE);

      // If WBP, read list of thresholds
      } else if (recons_method==use_WBP) {
         int skip=0;
         do {
            str=get_param(fh_param,"threshold",skip,"");
            if (str!="") {
               WBP_threshold.push_back(AtoF(first_word(str),3007,
                  "Recons_test_Parameters::read: WBP threshold is not a true number"));
               skip++;
            } else if (skip==0)
               REPORT_ERROR(3007,"Recons_test_Parameters::read: There are"
                  " no threshold parameters");
            else break;
         } while (TRUE);
      }
   } catch (Xmipp_error XE) {
      cout << XE << endl;
      REPORT_ERROR(3007,(string)"There is an error reading "+fn_test_params);
   }
   fclose(fh_param);
}

/* Show parameters ========================================================= */
ostream & operator << (ostream &out, const Recons_test_Parameters &prm) {
   out << "Reconstruction Test Parameters ===================\n";
   out << "   Reconstruction method=";
   switch (prm.recons_method) {
      case use_ART:         out << "ART\n"; break;
      case use_SIRT:        out << "SIRT\n"; break;
      case use_WBP:         out << "WBP\n"; break;
      case use_SIRT_Spider: out << "SIRT Spider\n"; break;
   }
   out << "   Phantom family: " << prm.fn_random_phantom << endl;
   out << "   Voxel Phantom: "  << prm.fn_voxel_phantom  << endl;
   out << "   Projection parameters: " << prm.fn_proj_params << endl;
   out << "   Crystal parameters: " << prm.fn_crystal << endl;
   out << "   Random Sort: "; print(out,prm.random_sort); out << endl;
   out << "   Sort with last: " << prm.sort_last_N << endl;
   out << "   Measurements: " << prm.MeasNo << endl;
   out << "   Accuracy: " << prm.accuracy << endl;
   out << "   Unluckiness: " << prm.unluckiness << endl;
   out << "   Only structural: " << prm.only_structural << endl;
   out << "   Global radius: " << prm.global_radius << endl;
   out << "   Max resolution: " << prm.max_resolution << endl;
   out << "   Symmetry file: " << prm.fn_sym << endl;
   out << "   Final Symmetry file: " << prm.fn_final_sym << endl;
   out << "   CTF file: " << prm.fn_CTF << endl;
   out << "   Defocus change: " << prm.defocus_change << endl;
   out << "   Noise stddev: " << prm.sigma << endl;
   out << "   Noise lowpass before CTF: " << prm.low_pass_before_CTF << endl;
   out << "   Noise spectrum after CTF: " << prm.fn_after_CTF << endl;
   out << "   Probe radius: " << prm.probe_radius << endl;
   out << "   Top surface: ";
   if (prm.enable_top_surface)
      out << "TRUE z0=" << prm.top0 << " zF=" << prm.topF << endl;
   else out << "FALSE\n";
   out << "   Bottom surface: ";
   if (prm.enable_bottom_surface)
      out << "TRUE z0=" << prm.bottom0 << " zF=" << prm.bottomF << endl;
   else out << "FALSE\n";
   out << "   Start from phantom: "; print(out,prm.start_from_phantom);
      out << endl;
   out << "   Start from lowpass filter: " << prm.starting_low_pass << endl;
   out << "   Stop at: " << prm.stop_at << endl;
   out << "   Reconstruction radius: " << prm.reconstruction_radius << endl;
   out << "   Run also without constraints: ";
      print(out,prm.run_also_without_constraints); out << endl;
   if (prm.enable_normalization) {
      out << "   Y=AX+B: A follows N(" << prm.a_avg << "," << prm.a_stddev
          << ") and B follows N(" << prm.b_avg << ","
	  << prm.b_stddev << ")\n"
	  << "   Normalizing method: ";
      switch (prm.normalizing_method) {
	 case OLDXMIPP:      out << "OldXmipp\n"; break;
	 case NEAR_OLDXMIPP: out << "Near_OldXmipp\n"; break;
	 case NEWXMIPP:      out << "NewXmipp\n"; break;
	 case NEWXMIPP2:     out << "NewXmipp2\n"; break;
	 case MICHAEL:       out << "Michael\n"; break;
	 case NONE:          out << "None\n"; break;
      }
      if (prm.normalizing_method==NEWXMIPP || 
	  prm.normalizing_method==NEWXMIPP2 ||
	  prm.normalizing_method==NEAR_OLDXMIPP ||
	  prm.normalizing_method==MICHAEL)
	  cout << "   Background mode: Circle, radius " << prm.bg_radius << endl;
   }
   if (prm.correct_phase) {
      out << "   Correcting CTF phase\n"
          << "   Small is under " << prm.phase_correction_param << endl
          << "   Correcting method: ";
      switch (prm.phase_correction_method) {
	 case CORRECT_SETTING_SMALL_TO_ZERO:
            out << "Set small values to 0\n"; break;
	 case CORRECT_LEAVING_SMALL:
            out << "Leave small values as they are\n"; break;
	 case CORRECT_AMPLIFYING_NOT_SMALL:
            out << "Correct amplitude except for the small values\n";
            break;
      }
   }
   
   if (prm.correct_amplitude)
      out << "   Correcting CTF amplitude\n"
          << "   No. IDR iterations: " << prm.idr_iterations << endl
	  << "   IDR relaxation parameters FROM: " << prm.mu0_list.transpose() << endl
	  << "                               TO: " << prm.muF_list.transpose() << endl;
   if (prm.unmatched)
      out << "   Unmatched CTF correction\n";
   if (prm.recons_method==use_ART || prm.recons_method==use_SIRT) {
      if (prm.blob_type==BIG_BLOB) out << "   Blob type: big\n";
      else                         out << "   Blob type: small\n";
      out << "   Succesive parameters: "; print(out,prm.succesive_params);
         out << endl;
      if (prm.POCS_positivity) out << "   Positivity constraint allowed\n";
      if (prm.force_sym!=-1)   out << "   Symmetry forced=" << prm.force_sym << "\n";
      if (prm.do_not_use_symproj) out << "   Do not use symmetrized projections\n";
      if (prm.enable_segmented_surface)
         out << "   Segmented surface threshold: " << prm.threshold_surface_segment << endl;
      out << "   Dilation for segmented volumes/surfaces: " << prm.segmented_dilation << endl;
      out << "   Constrained mass: " << prm.mass << endl;
   }
   if (prm.recons_method==use_WBP)
      for (int i=0; i<prm.WBP_threshold.size(); i++)
          out << "   Threshold=" << prm.WBP_threshold[i] << endl;
   else
      for (int i=0; i<prm.lambda0.size(); i++)
          out << "   Lambda0=" << prm.lambda0[i]
              << " LambdaF=" << prm.lambdaF[i]
              << " No It0 =" << prm.no_it0[i] 
              << " No ItF =" << prm.no_itF[i] 
              << endl;
   return out;
}

/* Make a single measure on scL2 =========================================== */
void single_measure_on_FOM(Recons_test_Parameters &prm,
   int i, int &nvol,
   double &training_avg, double &training_stddev, double &training_N,
   EVALUATE_results &results, const string &training_FOM) {
   double min, max;
   bool accuracy_mode=prm.MeasNo==-1;
   int sample_size=(accuracy_mode)? 3:prm.MeasNo;
   matrix1D<double> training_FOMs(sample_size);
   prm.only_structural=TRUE;
   for (int k=0; k<sample_size; k++) {
      cout << "Making measure number: " << k+1 << endl;
      single_recons_test(prm, i, nvol, results);
      if        (training_FOM=="scL20") training_FOMs(k)=results.scL2_FOMs(0);
      else if   (training_FOM=="scL2")  training_FOMs(k)=results.scL2_FOM;
      else if   (training_FOM=="scL2w") training_FOMs(k)=results.scL2w_FOM;
      else if   (training_FOM=="scL10") training_FOMs(k)=results.scL1_FOMs(0);
      else if   (training_FOM=="scL1")  training_FOMs(k)=results.scL1_FOM;
      else if   (training_FOM=="scL1w") training_FOMs(k)=results.scL1w_FOM;
      else if   (training_FOM=="scL21") {
         matrix1D<double> aux=results.scL2_FOMs; aux.window(1,XSIZE(aux)-1);
         training_FOMs(k)=aux.compute_avg();
      } else if (training_FOM=="scL11") {
         matrix1D<double> aux=results.scL1_FOMs; aux.window(1,XSIZE(aux)-1);
         training_FOMs(k)=aux.compute_avg();
      }
      if (nvol!=-1) nvol++;
      if (accuracy_mode && k>0) {
         matrix1D<double> aux=training_FOMs; aux.window(0,k);
         aux.compute_stats(training_avg, training_stddev, min, max);
         double t=student_outside_probb(prm.unluckiness,k+1);
         double estimated_sample_size=
            t*training_stddev/(prm.accuracy/100*training_avg);
         cout << "tFOM values=" << aux.transpose() << endl
              << estimated_sample_size << " samples will be needed\n";
         if (sample_size<estimated_sample_size && k==sample_size-1)
            {sample_size++; training_FOMs.resize(sample_size);}
      }
   }
   training_FOMs.compute_stats(training_avg, training_stddev, min, max);
   training_N=sample_size;
}

/* Make a single measure on all FOMs ======================================= */
void single_measure_on_all_FOMs(Recons_test_Parameters &prm, int i,
   int &nvol, FOMs &foms_mean, FOMs &foms_stddev, EVALUATE_results &results) {
   prm.only_structural=FALSE;

   FOMs foms(prm.MeasNo);
   for (int k=0; k<XSIZE(foms.scL2); k++) {
      cout << "Making measure number: " << k << endl;
      single_recons_test(prm, i, nvol, results);
      if (nvol!=-1) nvol++;
      foms.set_FOMs(k,results);
   }

   compute_FOMs_stats(foms,i,foms_mean,foms_stddev);
   cout << foms;
}

/* Make a single test ====================================================== */
void single_recons_test(const Recons_test_Parameters &prm,
   int i, int nvol, EVALUATE_results &results) {
// Get Filename root -------------------------------------------------------
   Prog_Project_Parameters        Prog_proj_prm;
   Projection_Parameters          proj_prm;
   Crystal_Projection_Parameters  crystal_proj_prm;

   
   Prog_proj_prm.fn_proj_param=prm.fn_proj_params;
   proj_prm.from_prog_params(Prog_proj_prm);
   if (prm.fn_crystal!="") crystal_proj_prm.read(prm.fn_crystal);
   FileName fn_root, fn_recons_root;
   fn_root=proj_prm.fn_projection_seed;
   if (nvol!=-1) fn_recons_root=fn_root+"exp"+ItoA(nvol,2);
   else          fn_recons_root=fn_root;
   FileName fn_ext=proj_prm.fn_projection_extension;


// Generate random phantom -------------------------------------------------
   Prog_Random_Phantom_Parameters rp_prm;
   Phantom realization;
   FileName fn_phantom;
   
   if (prm.fn_random_phantom!="") {
      rp_prm.fn_random=prm.fn_random_phantom;
      fn_phantom=rp_prm.fn_output=fn_recons_root+".descr";
      rp_prm.min_vol=0;

      ROUT_random_phantom(rp_prm,realization);
   } else 
      fn_phantom=prm.fn_voxel_phantom;

// Read phantom in memory? -------------------------------------------------
    VolumeXmipp vol_phantom;
    if (prm.enable_segmented_surface || prm.start_from_phantom)
       if (prm.fn_random_phantom!="") realization.draw_in(&vol_phantom);
       else                           vol_phantom.read(fn_phantom);
       vol_phantom().set_Xmipp_origin();

// Generate projections ----------------------------------------------------
   Projection Proj;
   SelFile SF;

   if (prm.recons_method!=use_WBP) Prog_proj_prm.fn_sel_file=fn_root+".sel";
   else Prog_proj_prm.fn_sel_file=fn_root.substr(0,2)+"t.sel";

   // Read projection parameters and produce side information
   proj_prm.from_prog_params(Prog_proj_prm);
   if (prm.fn_random_phantom!="") proj_prm.fn_phantom=fn_phantom;
   proj_prm.fn_projection_seed=fn_root;

   PROJECT_Side_Info side;
   side.produce_Side_Info(proj_prm);

   PROJECT_Effectively_project(proj_prm, side, crystal_proj_prm, Proj, SF);
   SF.write(Prog_proj_prm.fn_sel_file);
   SF.go_first_ACTIVE();

// Adding microscope effect ------------------------------------------------
   if (prm.sigma!=0 ||
       prm.low_pass_before_CTF!=0 ||
       prm.fn_CTF!="" || prm.fn_after_CTF!="") {
         Prog_Microscope_Parameters prm_micro;
         prm_micro.fn_in=Prog_proj_prm.fn_sel_file;
	 prm_micro.fn_ctf=prm.fn_CTF;
         prm_micro.defocus_change=prm.defocus_change;
	 prm_micro.sigma=prm.sigma;
	 prm_micro.low_pass_before_CTF=prm.low_pass_before_CTF;
	 prm_micro.fn_after_ctf=prm.fn_after_CTF;
	 prm_micro.produce_side_info();
	 
      	 cerr << "Applying microscope simulation ...\n";
	 init_progress_bar(SF.ImgNo());
	 int i=0;
	 while (!SF.eof()) {
      	    FileName fn_proj=SF.NextImg();
	    ImageXmipp I; I.read(fn_proj); I().set_Xmipp_origin();

      	    prm_micro.apply(I());

	    I.write(); i++;
	    if (i%20==0) progress_bar(i);
	 }
	 progress_bar(SF.ImgNo());
	 SF.go_first_ACTIVE();
   }

// Normalize ---------------------------------------------------------------
   randomize_random_generator();
   if (prm.enable_normalization) {
      Normalize_parameters norm_prm;
      norm_prm.fn_in=SF.name();
      norm_prm.normalizing_method=prm.normalizing_method;
      norm_prm.background_mode=CIRCLE;
      norm_prm.r=prm.bg_radius;
      norm_prm.produce_side_info();
      cerr << "Applying linear transformation and normalizing ...\n";
      init_progress_bar(SF.ImgNo());
      int n=0;
      while (!SF.eof()) {
      	 FileName fn_proj=SF.NextImg();
	 ImageXmipp I; I.read(fn_proj); I().set_Xmipp_origin();

	 double a=rnd_gaus(prm.a_avg, prm.a_stddev);
	 double b=rnd_gaus(prm.b_avg, prm.b_stddev);
         FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
	    IMGPIXEL(I,i,j)=a*IMGPIXEL(I,i,j)+b;

	 norm_prm.apply(&I);
	 I.write();n++;
	 if (n%20==0) progress_bar(n);
      }
      progress_bar(SF.ImgNo());
      SF.go_first_ACTIVE();
   }

// Correct phase -----------------------------------------------------------
   FileName fn_applied_CTF;
   FourierMask ctf;
   if (prm.correct_phase) {
      CorrectPhase_Params correct;
      correct.fn_ctf=prm.fn_CTF;
      correct.multiple_CTFs=FALSE;
      correct.CTF_description_file=!Is_FourierImageXmipp(prm.fn_CTF);
      correct.method=prm.phase_correction_method;
      correct.epsilon=prm.phase_correction_param;
      correct.produce_side_info();
      ctf=correct.ctf;
      // Correct images
      correct.correct(SF);
      // Correct the CTF itself
      correct.correct(ctf.mask2D);
      // Save corrected CTF
      fn_applied_CTF=prm.fn_CTF.insert_before_extension("_phase_corrected");
      ctf.write_mask(fn_applied_CTF,2);
   } else
      fn_applied_CTF=prm.fn_CTF;

// Generate surface --------------------------------------------------------
   Prog_Surface_Parameters prm_surface;
   FileName fn_mask;
   // In AFM mode
   if ((prm.enable_top_surface || prm.enable_bottom_surface) &&
       (prm.recons_method==use_ART || prm.recons_method==use_SIRT)) {
      if (prm.fn_random_phantom=="")
         REPORT_ERROR(1,
           "Recons_test: Cannot use surface option without a mathematical phantom");
      cerr << "Generating surface ...\n";
      prm_surface.probe_radius=prm.probe_radius;
      prm_surface.fn_phantom=fn_phantom;
      prm_surface.phantom=realization;
      prm_surface.zdim=realization.zdim;
      if (prm.enable_top_surface) {
         prm_surface.enable_ztop=TRUE;
         prm_surface.ztop=rnd_unif(prm.top0,prm.topF);
         prm_surface.fn_top=fn_recons_root+"_top_surface";
            prm_surface.fn_top.add_extension(fn_ext);
      }
      if (prm.enable_bottom_surface) {
         prm_surface.enable_zbottom=TRUE;
         prm_surface.zbottom=rnd_unif(prm.bottom0,prm.bottomF);
         prm_surface.fn_bottom=fn_recons_root+"_bottom_surface";
            prm_surface.fn_bottom.add_extension(fn_ext);
      }
      prm_surface.fn_mask=fn_mask=fn_recons_root+"_mask.vol";
      ROUT_surface(prm_surface);
   // In segmented mode
   } else if (prm.enable_segmented_surface &&
       (prm.recons_method==use_ART || prm.recons_method==use_SIRT)) {
       VolumeXmipp aux;
       vol_phantom().threshold("below",prm.threshold_surface_segment,
          prm.threshold_surface_segment);
       vol_phantom().binarize(prm.threshold_surface_segment);
       aux().resize(vol_phantom());
       if (prm.segmented_dilation!=0)
          dilate3D(vol_phantom(),aux(),18,0,prm.segmented_dilation);
       aux()*=-1; // Invert mask
       aux()+=1;
       fn_mask=fn_recons_root+"_mask.vol";
       aux.write(fn_mask);
   }

// Reconstruct -------------------------------------------------------------
   VolumeXmipp vol_recons;
   FourierMask Filter;
   if (prm.recons_method==use_ART || prm.recons_method==use_SIRT) {
      Basic_ART_Parameters art_prm;
      Plain_ART_Parameters plain_art_prm;
      GridVolume           vol_blobs;
      GridVolume           *vol_blobs_var=NULL;

      art_prm.default_values();
      // art_prm.tell |= TELL_SHOW_ERROR;
      // art_prm.tell |= TELL_SAVE_AT_EACH_STEP;
      if (prm.blob_type==BIG_BLOB) {
         art_prm.blob.alpha=3.6;
         art_prm.grid_relative_size=2.26;
      } else {
         art_prm.blob.alpha=10.4;
         art_prm.grid_relative_size=1.41;
      }
      art_prm.fn_surface_mask="";
      art_prm.fn_sym="";
      art_prm.stop_at=prm.stop_at;
      art_prm.R=prm.reconstruction_radius;
      art_prm.fn_sel=Prog_proj_prm.fn_sel_file;
      art_prm.proj_ext=0;
      art_prm.max_tilt=1e7;
      if (!prm.succesive_params) {
         art_prm.lambda_list.resize(1);
	 art_prm.lambda_list(0)=rnd_log(prm.lambda0[i],prm.lambdaF[i]);
      	 art_prm.no_it=(int)rnd_log(prm.no_it0[i],prm.no_itF[i]);
      } else {
         art_prm.no_it=prm.lambda0.size();
         art_prm.lambda_list.resize(art_prm.no_it);
         for (int j=0; j<prm.lambda0.size(); j++)
      	     art_prm.lambda_list(j)=rnd_log(prm.lambda0[j],prm.lambdaF[j]);
      }
      art_prm.random_sort=prm.random_sort;
      art_prm.sort_last_N=prm.sort_last_N;
      if (prm.recons_method==use_SIRT)
         art_prm.parallel_mode=Basic_ART_Parameters::SIRT;
      if (prm.POCS_positivity) art_prm.positivity=TRUE;
      if (prm.unmatched) {
         art_prm.unmatched=true;
         art_prm.fn_ctf=fn_applied_CTF;
      }

      cout << "Selected: Lambda= " << art_prm.lambda_list.transpose() << endl
           << " No_it= " << art_prm.no_it << endl;

      if (prm.run_also_without_constraints) {
         art_prm.fn_root=fn_recons_root+"_wos";
         Basic_ROUT_Art(art_prm,plain_art_prm,vol_recons,vol_blobs,
            vol_blobs_var);
      }

      // Extra conditions
      art_prm.fn_root=fn_recons_root;
      art_prm.fn_sym=prm.fn_sym;
      art_prm.force_sym=prm.force_sym;
      art_prm.do_not_use_symproj=prm.do_not_use_symproj;
      art_prm.known_volume=prm.mass;
      if (prm.enable_top_surface || prm.enable_bottom_surface ||
          prm.enable_segmented_surface)
         art_prm.fn_surface_mask=fn_mask;

      if (prm.start_from_phantom) {
         cerr << "Filtering phantom ...\n";
         VolumeXmipp starting_vol;
         starting_vol()=vol_phantom();
         Filter.FilterShape=RAISED_COSINE;
         Filter.FilterBand=LOWPASS;
         Filter.w1=prm.starting_low_pass;
         Filter.raised_w=0.02;
         Filter.show();
         Filter.apply_mask_Space(starting_vol());
         starting_vol.write(fn_recons_root+"_starting.vol");

         cerr << "Converting phantom to blobs ...\n";
         voxels2blobs(&starting_vol, art_prm.blob, vol_blobs, BCC,
            art_prm.grid_relative_size, 0.25, NULL, NULL, 0.01, 0,
            CEIL(XSIZE(starting_vol())/2));
         art_prm.fn_start=fn_recons_root+"_starting.blob";
         vol_blobs.write(art_prm.fn_start);
         art_prm.fn_root=fn_recons_root;
      }

      if (!prm.correct_amplitude)
	 // Do not correct 
	 Basic_ROUT_Art(art_prm,plain_art_prm,vol_recons, vol_blobs,
            vol_blobs_var);
      else {
	 Prog_IDR_ART_Parameters idr_prm;
	 idr_prm.art_prm=&art_prm;
      	 idr_prm.idr_iterations=prm.idr_iterations;
	 idr_prm.it=0;
      	 idr_prm.dont_rewrite=FALSE;
      	 idr_prm.mu0_list=prm.mu0_list;
      	 idr_prm.muF_list=prm.muF_list;
      	 idr_prm.fn_ctf=fn_applied_CTF;
         idr_prm.max_resolution=prm.max_resolution;
         idr_prm.fn_final_sym=prm.fn_final_sym;
	 idr_prm.produce_side_info();
	 Basic_ROUT_IDR_Art(idr_prm, vol_recons);
	 fn_recons_root=vol_recons.name().without_extension();
      }
   } else if (prm.recons_method==use_WBP) {
      // Estaría bien que WBP se tradujese y no hubiese que llamarlo
      // desde el sistema

      // Change to the old recvol file image format
      cerr << "Adapting angles to old recvol format ...\n";
      while (!SF.eof()) {
         ImageXmipp P;
         P.read(SF.NextImg());
         P.old_rot()=-P.rot();
         P.rot()=0;
         P.write();
      }
      SF.go_first_ACTIVE();

      string command_line=(string)"recvol "+fn_root[0]+fn_root[1]+
         " "+fn_recons_root+".vol "+ItoA((int)(proj_prm.proj_Xdim/2))+
         " -t"+FtoA(prm.WBP_threshold[i],0);
      cerr << "Reconstructing with WBP ...\n";
      system(command_line.c_str());
   } else if (prm.recons_method==use_SIRT_Spider) {
      double lambda=rnd_log(prm.lambda0[i],prm.lambdaF[i]);
      double no_it=(int)rnd_log(prm.no_it0[i],prm.no_itF[i]);
      int radius=(int)(proj_prm.proj_Xdim/2)-2;
      cout << "Selected: Lambda= " << lambda
           << " No_it= " << no_it << endl;
      SIRT_Spider(SF, lambda, no_it, radius, fn_root, fn_ext, fn_recons_root,
         "b73");
   }

// Filter result -----------------------------------------------------------
   if (prm.max_resolution!=-1) {
      cerr << "Filtering result ...\n";
      vol_recons.read(fn_recons_root+".vol");
      Filter.FilterShape=RAISED_COSINE;
      Filter.FilterBand=LOWPASS;
      Filter.w1=prm.max_resolution;
      Filter.raised_w=0.02;
      Filter.apply_mask_Space(vol_recons());
      vol_recons.write();
   }

// Symmetrize --------------------------------------------------------------
   if (prm.fn_final_sym!="") {
      Symmetrize_Parameters sym_prm;
      sym_prm.fn_in=fn_recons_root+".vol";
      sym_prm.fn_out="";
      sym_prm.fn_sym=prm.fn_final_sym;
      sym_prm.wrap=TRUE;
      ROUT_symmetrize(sym_prm);
   }

// Evaluate ----------------------------------------------------------------
   Prog_Evaluate_Parameters eval_prm;
   eval_prm.default_values();
   eval_prm.fit_gray_scales=TRUE;
   if (prm.only_structural) eval_prm.tell = ONLY_STRUCTURAL;
   if (prm.fn_random_phantom!="")
      eval_prm.fn_phantom=fn_phantom;
   else {
      eval_prm.fn_phantom=proj_prm.fn_phantom;
      eval_prm.tell=ONLY_STRUCTURAL;
   }
   eval_prm.fn_recons=fn_recons_root+".vol";
   eval_prm.fn_mask=fn_mask;
   if (fn_mask!="") {
      // Revert the mask for evaluation
      VolumeXmipp aux;
      aux.read(fn_mask);
      aux()*=-1;
      aux()+=1;
      aux.write();
   }
   if (prm.global_radius!=-1) eval_prm.global_radius=prm.global_radius;
   else eval_prm.global_radius=CEIL(proj_prm.proj_Xdim/2);
   eval_prm.tell |= SAVE_HISTOGRAMS | SAVE_MAPS | SHOW_PROCESS;
   cerr << "   Evaluating ...\n";
   ROUT_Evaluate(eval_prm,results);
}
