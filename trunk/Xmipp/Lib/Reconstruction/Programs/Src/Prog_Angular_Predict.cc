/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.uam.es (2002)
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

#include "../Prog_Angular_Predict.hh"
#include "../Prog_FourierFilter.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippHistograms.hh>
#include <XmippData/xmippGeometry.hh>

// Read arguments ==========================================================
void Prog_angular_predict_prm::read(int argc, char **argv) _THROW {
   Prog_parameters::read(argc,argv);
   fn_ref=get_param(argc,argv,"-ref");
   fn_pcaproj=get_param(argc,argv,"-pcaproj");
   fn_ang=get_param(argc,argv,"-ang");
   PCA_list=get_vector_param(argc,argv,"-PCA_list",-1);
   fn_out_ang=get_param(argc,argv,"-oang");
   fn_sym=get_param(argc,argv,"-sym","");
   max_proj_change=AtoF(get_param(argc,argv,"-max_proj_change","-1"));
   max_psi_change=AtoF(get_param(argc,argv,"-max_psi_change","-1"));
   max_shift_change=AtoF(get_param(argc,argv,"-max_shift_change","0"));
   psi_step=AtoF(get_param(argc,argv,"-psi_step","5"));
   shift_step=AtoF(get_param(argc,argv,"-shift_step","1"));
   tell=0;
   if (check_param(argc,argv,"-show_rot_tilt")) tell|=TELL_ROT_TILT;
   if (check_param(argc,argv,"-show_psi_shift")) tell|=TELL_PSI_SHIFT;
   if (check_param(argc,argv,"-show_options")) tell|=TELL_OPTIONS;
   produce_side_info();
}

// Show ====================================================================
void Prog_angular_predict_prm::show() {
   Prog_parameters::show();
   cout << "Reference rootname: " << fn_ref << endl
        << "PCA projection rootname: " << fn_pcaproj << endl
        << "Angle file: " << fn_ang << endl;
   cout << "PCA list:           ";
   if (XSIZE(PCA_list)==0) cout << "All\n";
   else                    cout << PCA_list.transpose() << endl;
   cout << "Ouput angular file: " << fn_out_ang << endl
	<< "Max proj change: " << max_proj_change << endl
	<< "Max psi change: " << max_psi_change << " step: " << psi_step << endl
	<< "Max shift change: " << max_shift_change << " step: " << shift_step << endl
        << "Show level: " << tell << endl
   ;
}

// usage ===================================================================
void Prog_angular_predict_prm::usage() {
   Prog_parameters::usage();
   cerr << "   -ref <ref rootname>      : Rootname passed to Angular_Reference\n"
        << "   -pcaproj <pca rootname>  : Rootname passed to Angular_Project\n"
	<< "  [-PCA_list \"[n,n,n,...]\"] : Number of the PCAs to use for prediction\n"
	<< "                            : by default, all\n"
	<< "   -ang <angle file>        : DocFile with the angles for the reference\n"
	<< "                              produced by xmipp_project\n"
//	<< "  [-prev_ang <angle file>]  : DocFile produced by this same program\n"
//	<< "                              with the predicted angles from previous stages\n"
	<< "   -oang <angle file>       : DocFile with output angles\n"
        << "  [-sym <symmetry file>]    : Symmetry file if any\n"
	<< "  [-max_proj_change <ang=-1>]: Maximum change allowed in rot-tilt\n"
	<< "  [-max_psi_change <ang=-1>]: Maximum change allowed in psi\n"
	<< "  [-max_shift_change <r=0>] : Maximum change allowed in shift\n"
	<< "  [-psi_step <ang=5>]       : Step in psi in degrees\n"
	<< "  [-shift_step <r=1>]       : Step in shift in pixels\n"
	<< "  [-show_rot_tilt]          : Show the rot-tilt process\n"
	<< "  [-show_psi_shift]         : Show the psi-shift process\n"
	<< "  [-show_options]           : Show final options among which\n"
	<< "                              the angles are selected\n"
   ;
}

// Produce side information ================================================
void Prog_angular_predict_prm::produce_side_info() _THROW {
   each_image_produces_an_output=false;
   allow_time_bar=(tell==0);

   // Produce side info of the PCA projector
   PCA_projector.do_not_write=true;
   PCA_projector.fn_ref=fn_ref;
   PCA_projector.PCA_list=PCA_list;
   PCA_projector.fn_in=fn_in;
   PCA_projector.produce_side_info();

   // Produce side info of the angular distance computer
   distance_prm.fn_ang1=distance_prm.fn_ang2="";
   distance_prm.fn_sym=fn_sym;
   distance_prm.produce_side_info();

   // Read the angle file
   DocFile DF;
   DF.read(fn_ang);
   DF.go_first_data_line();
   rot.resize(DF.dataLineNo());
   tilt.resize(DF.dataLineNo());
   int i=0;
   while (!DF.eof()) {
      rot[i]=DF(0);  // Rotational angle
      tilt[i]=DF(1); // Tilting angle
      i++; DF.next_data_line();
   }
   
   // Resize the predicted vectors
   int number_of_images=get_images_to_process();
   current_img=0;
   predicted_rot.resize(number_of_images);
   predicted_tilt.resize(number_of_images);
   predicted_psi.resize(number_of_images);
   predicted_shiftX.resize(number_of_images);
   predicted_shiftY.resize(number_of_images);
   predicted_corr.resize(number_of_images);
   
   // Read the PCA projection of the reference images onto
   // the selected PCAs
   PCA_ref.resize(PCA_projector.PCANo);
   for (int m=0; m<PCA_projector.PCANo; m++) {
      PCA_ref[m]=NULL;
      if (PCA_projector.use_this_PCA[m]) {
         PCA_projector.PCASet(m)->prepare_for_correlation();
      
         PCA_ref[m]=new xmippCTVectors(0,true);
         FileName fn_in=fn_pcaproj+"_PCA_"+ItoA(m,2)+".dat";
	 ifstream fh_in;
	 fh_in.open(fn_in.c_str());
	 if (!fh_in)
	    REPORT_ERROR(1,(string)": Cannot find file"+fn_in);

	 fh_in >> *(PCA_ref[m]);
	 fh_in.close();
      }
   }
}

// Build candidate list ------------------------------------------------------
void Prog_angular_predict_prm::build_ref_candidate_list(const ImageXmipp &I,
   vector<bool> &candidate_list, vector<double> &cumulative_corr,
   vector<double> &sumxy, vector<double> &sumxx, vector<double> &sumyy) {
   int refNo=rot.size();
   candidate_list.resize(refNo);
   cumulative_corr.resize(refNo);
   sumxx.resize(refNo);
   sumyy.resize(refNo);
   sumxy.resize(refNo);
   for (int i=0; i<refNo; i++) {
      candidate_list[i]=true;
      sumxx[i]=sumyy[i]=sumxy[i]=cumulative_corr[i]=0;
      if (max_proj_change!=-1) {
         double dummy_rot=rot[i], dummy_tilt=tilt[i], dummy_psi;
         double ang_distance=distance_prm.check_symmetries(
	    I.rot(),I.tilt(),0,dummy_rot,dummy_tilt,dummy_psi,true);
	 candidate_list[i]=(ang_distance<=max_proj_change);
	 #ifdef DEBUG
	    cout << "(" << I.rot() << "," << I.tilt() << ") and ("
		 << rot[i] << "," << tilt[i] << ") --> " << ang_distance << endl;
	 #endif
      }
   }
}

// Refine candidate list ---------------------------------------------------
void Prog_angular_predict_prm::refine_candidate_list_with_correlation(
   int m,
   xmippVector &PCA, xmippCTVectors &PCA_ref, vector<bool> &candidate_list,
   vector<double> &cumulative_corr,
   vector<double> &sumxy, vector<double> &sumxx, vector<double> &sumyy,
   double &dim, double th) {
   histogram1D hist;
   hist.init(-1,1,201);   

   int dimPCA=PCA.size();
   int dimp=PCA_projector.PCASet(m)->get_eigenDimension();
   dim+=dimp;
   int imax=candidate_list.size();
   for (int i=0; i<imax; i++) {
      if (candidate_list[i]) {
         double sumxyp=PCA_projector.PCASet(m)->prod_mean_mean;
	 double sumxxp=PCA_projector.PCASet(m)->prod_mean_mean;
	 double sumyyp=PCA_projector.PCASet(m)->prod_mean_mean;
	 double x_mean=PCA_projector.PCASet(m)->avg_mean;
	 double y_mean=PCA_projector.PCASet(m)->avg_mean;
	 #define x PCA[d]
	 #define y PCA_ref.itemAt(i)[d]
	 for (int d=0; d<dimPCA; d++) {
	    x_mean+=x*PCA_projector.PCASet(m)->avg_ei[d];
	    y_mean+=y*PCA_projector.PCASet(m)->avg_ei[d];
	    sumxxp+=
	       x*x*PCA_projector.PCASet(m)->prod_ei_ei[d]+
	       2*x*PCA_projector.PCASet(m)->prod_ei_mean[d];
	    sumyyp+=
	       y*y*PCA_projector.PCASet(m)->prod_ei_ei[d]+
	       2*y*PCA_projector.PCASet(m)->prod_ei_mean[d];
	    sumxyp+=
	       x*y*PCA_projector.PCASet(m)->prod_ei_ei[d]+
	       (x+y)*PCA_projector.PCASet(m)->prod_ei_mean[d];
	 }
	 #undef x
	 #undef y
         sumxx[i]+=sumxxp-x_mean*x_mean*dimp;
	 sumyy[i]+=sumyyp-y_mean*y_mean*dimp;
	 sumxy[i]+=sumxyp-x_mean*y_mean*dimp;

	 double corr=sumxy[i]/sqrt(sumxx[i]*sumyy[i]);
	 cumulative_corr[i]/*+*/=corr;
         hist.insert_value(cumulative_corr[i]);

	 if (tell & TELL_ROT_TILT) {
            cout << "Candidate " << i << " corr= " << cumulative_corr[i]
	         << " rot= " << rot[i] << " tilt= " << tilt[i] << endl;
	 }
      }
   }
   double corr_th=hist.percentil(100-th);

   // Remove all those projections below the threshold
   for (int i=0; i<imax; i++)
      if (candidate_list[i] && cumulative_corr[i]<corr_th)
         candidate_list[i]=false;
	 
   // Show the percentil used
   if (tell & TELL_ROT_TILT) {
      cout << "# Percentil " << corr_th << endl << endl;
   }
}

// Predict rot and tilt ----------------------------------------------------
double Prog_angular_predict_prm::predict_rot_tilt_angles(ImageXmipp &I,
   double &assigned_rot, double &assigned_tilt) _THROW {
   // Build initial candidate list
   vector<bool>   candidate_list;
   vector<double> cumulative_corr;
   vector<double> sumxx, sumyy, sumxy;
   build_ref_candidate_list(I,candidate_list,cumulative_corr,sumxy,sumxx,sumyy);
   int imax=candidate_list.size();

   // Project this image onto the PCA space
   vector<xmippVector> Vpca;
   PCA_projector.project_image(I(), Vpca);
   
   // Measure correlation for all possible PCAs
   // These variables are used to compute the correlation at a certain
   // scale
   double dim=0;
   int    PCAs_used=0;
   for (int m=0; m<PCA_projector.PCANo; m++)
      if (PCA_projector.use_this_PCA[m]) {
         // Show image name
         if (tell & TELL_ROT_TILT)
	    cout << "# " << I.name() << " m=" << m
	         << " current rot="  << I.rot()
		 << " current tilt=" << I.tilt() << endl;
         PCAs_used++;
         refine_candidate_list_with_correlation(m,Vpca[m],*(PCA_ref[m]),
	    candidate_list, cumulative_corr, sumxy, sumxx, sumyy, dim, 50);
      }

   // Select the maximum
   int best_i=-1;
   bool first=true;
   int N_max=0;
   for (int i=0; i<imax; i++)
      if (candidate_list[i])
         if (first) {best_i=i; first=false; N_max=1;}
	 else if (cumulative_corr[i]>cumulative_corr[best_i]) best_i=i;
	 else if (cumulative_corr[i]==cumulative_corr[best_i])
	         N_max++;

   if (N_max==0) {
      cerr << "Predict_angles: Empty candidate list for image "
           << I.name() << endl;
      assigned_rot=I.rot();
      assigned_tilt=I.tilt();
      return 0;
   }

   // There are several maxima, choose one randomly
   if (N_max!=1) {
      int selected=FLOOR(rnd_unif(0,3));
      for (int i=0; i<imax && selected>=0; i++)
	 if (candidate_list[i])
            if (cumulative_corr[i]==cumulative_corr[best_i])
	       {best_i=i; selected--;}
   }
   
   assigned_rot    = rot[best_i];
   assigned_tilt   = tilt[best_i];
   if (PCAs_used==0) return 0;
   else return cumulative_corr[best_i]/PCAs_used;
}

// Evaluate candidates by correlation ----------------------------------------
double Prog_angular_predict_prm::evaluate_candidates_by_correlation(
   const vector<double> &vcorr, const vector<int> &candidate_idx,
   vector<double> &candidate_rate, double weight) {
   // Compute maximum and minimum of correlations
   int imax=vcorr.size();
   double min_corr, max_corr;
   min_corr=max_corr=vcorr[0];
   for (int i=1; i<imax; i++) {
      if (vcorr[i]<min_corr) min_corr=vcorr[i];
      if (vcorr[i]>max_corr) max_corr=vcorr[i];
   }
   
   // Divide the correlation segment in as many pieces as candidates
   double corr_step=(max_corr-min_corr)/10;

   int jmax=candidate_idx.size();
   for (int j=0; j<jmax; j++) {
      int i=candidate_idx[j];
      int points;
      if (corr_step!=0) points=FLOOR((vcorr[i]-min_corr)/corr_step);
      else              points=10;
      if (tell & TELL_PSI_SHIFT)
         cout << "Candidate (" << i << ") corr=" << vcorr[i] << " points=" << points << endl;
      candidate_rate[j]+=weight*points;
   }

   return min_corr+7*corr_step;
}

// Group images --------------------------------------------------------------
void Prog_angular_predict_prm::group_views(const vector<double> &vrot,
   const vector<double> &vtilt, const vector<double> &vpsi,
   const vector<int> &best_idx, const vector<int> &candidate_idx, 
   vector< vector<int> > &groups) {
   for (int j=0; j<best_idx.size(); j++) {
      int i=candidate_idx[best_idx[j]];
      double roti=vrot[i];
      double tilti=vtilt[i];
      double psii=vpsi[i];
      // See if there is any suitable group
      bool assigned=false;
      int g;
      for (g=0; g<groups.size(); g++) {
         bool fits=true;
         for (int jp=0; jp<groups[g].size(); jp++) {
	    int ip=candidate_idx[groups[g][jp]];
	    double ang_distance=distance_prm.check_symmetries(
	       vrot[ip],vtilt[ip],vpsi[ip],roti,tilti,psii,false);
	    if (ang_distance>20) {fits=false; break;}
	 }
	 if (fits) {assigned=true; break;}
      }

      if (!assigned) {
         // Create a new group with the first item in the list
         vector<int> group;
	 group.push_back(best_idx[j]);
	 groups.push_back(group);
      } else {
         // Insert the image in the fitting group
         groups[g].push_back(best_idx[j]);
      }
   }

   // Check if there are so many groups as images.
   // If that is the case, then everything is a single group
   if (groups.size()==best_idx.size()) {
      groups.clear();
      vector<int> group;
      for (int j=0; j<best_idx.size(); j++) group.push_back(best_idx[j]);
      groups.push_back(group);
   }
}

// Pick view -----------------------------------------------------------------
int Prog_angular_predict_prm::pick_view(vector< vector<int> >groups,
   const vector<double> &vcorr, 
   const vector<int> &best_idx,
   const vector<int> &candidate_idx, const vector<double> &candidate_rate) {
   // Sum the rates in all groups
   vector<double> group_rate;
   group_rate.reserve(groups.size());
   int best_g;
   double best_group_rate=-1e38;
   for (int g=0; g<groups.size(); g++) {
      double temp_group_rate=0;
      for (int j=0; j<groups[g].size(); j++)
         temp_group_rate+=candidate_rate[groups[g][j]];
      group_rate.push_back(temp_group_rate);
      if (temp_group_rate>best_group_rate) {
         best_g=g;
	 best_group_rate=group_rate[g];
      }
   }
   
   // Check that there are not two groups with the same score
   int groups_with_max_rate=0;
   for (int g=0; g<groups.size(); g++)
      if (group_rate[g]==best_group_rate) groups_with_max_rate++;
   if (groups_with_max_rate>1) {
      cout << "There are two groups with maximum rate\n";
   }

   // Take the best image within that group
   int best_j;
   double best_rate=-1e38;
   for (int j=0; j<groups[best_g].size(); j++) {
      #ifdef NEVER_DEFINED
      // Select the best with the rate
      if (candidate_rate[groups[best_g][j]]>best_rate) {
         best_j=j;
	 best_rate=candidate_rate[groups[best_g][j]];
      }
      #endif
      // Select the best with the correlation
      if (vcorr[candidate_idx[groups[best_g][j]]]>best_rate) {
         best_j=j;
	 best_rate=vcorr[candidate_idx[groups[best_g][j]]];
      }
   }

   // Check that there are not two images with the same rate
   int images_with_max_rate=0;
   for (int j=0; j<groups[best_g].size(); j++)
      #ifdef NEVER_DEFINED
      // Select the best with the rate
      if (candidate_rate[groups[best_g][j]]==best_rate)
         images_with_max_rate++;
      #endif
      // Select the best with correlation
      if (vcorr[candidate_idx[groups[best_g][j]]]==best_rate)
         images_with_max_rate++;
   if (images_with_max_rate>1) {
      // If there are several with the same punctuation take the one
      // with the best correlation
      double best_corr=-1e38;
      for (int j=0; j<groups[best_g].size(); j++) {
         if (vcorr[candidate_idx[groups[best_g][j]]]>best_corr &&
	     candidate_rate[groups[best_g][j]]==best_rate) {
            best_j=j;
	    best_corr=vcorr[candidate_idx[groups[best_g][j]]];
         }
      }
   }
   return groups[best_g][best_j];
}

// Predict shift and psi -----------------------------------------------------
double Prog_angular_predict_prm::predict_angles(ImageXmipp &I,
   double &assigned_shiftX, double &assigned_shiftY,
   double &assigned_rot, double &assigned_tilt, double &assigned_psi) {
   double best_rot, best_tilt, best_psi, best_shiftX, best_shiftY,
          best_corr=0, best_rate;
   ImageXmipp Ip;
   Ip=I;
   matrix1D<double> shift(2);

   // Establish psi limits
   double psi0, psiF;
   if (max_psi_change==-1) {psi0=-180; psiF=180-psi_step;}
   else {psi0=I.psi()-max_psi_change; psiF=I.psi()+max_psi_change;}
//   psi0=psiF=max_psi_change;
   double R2=max_shift_change*max_shift_change;

   // Search in the psi-shift space
   int N_trials=0;
   vector<double> vshiftX, vshiftY, vpsi, vrot, vtilt, vcorr;
   for (double shiftX=-max_shift_change; shiftX<=max_shift_change; shiftX+=shift_step)
      for (double shiftY=-max_shift_change; shiftY<=max_shift_change; shiftY+=shift_step) {
         if (shiftX*shiftX+shiftY*shiftY>R2) continue;
         for (double psi=psi0; psi<=psiF; psi+=psi_step) {
      	    // Shift image if necessary
            if (shiftX==0 && shiftY==0) Ip()=I();
	    else {
	       VECTOR_R2(shift,shiftX,shiftY);
	       I().translate(shift,Ip(),WRAP);
	    }
	    
	    // Rotate image if necessary
	    // Adding 2 is a trick to avoid that the 0, 90, 180 and 270
	    // are treated in a different way
	    Ip().self_rotate(psi+2, WRAP);
	    Ip().self_rotate(-2,WRAP);

      	    // Search for the best tilt, rot angles
            double rotp, tiltp;
            double corrp=predict_rot_tilt_angles(Ip,rotp,tiltp);
      	    if (tell & TELL_PSI_SHIFT)
               cout << "shiftX= " << shiftX << " shiftY= " << shiftY
	            << " psi= " << psi << " rot= " << rotp
		    << " tilt= " << tiltp << " corr= " << corrp
		    /* << " error= " << error_pca */
		    << endl; 
	    
	    vshiftX.push_back(shiftX);
	    vshiftY.push_back(shiftY);
	    vrot.push_back(rotp);
	    vtilt.push_back(tiltp);
	    vpsi.push_back(psi);
	    vcorr.push_back(corrp);
	    
	    N_trials++;
	 }
      }

   // Is the psi range circular?
   bool circular=realWRAP(vpsi[0]-psi_step,-180,180)==
                 realWRAP(vpsi[N_trials-1],-180,180);

   // Smooth the correlation curve
   #ifdef NEVER_DEFINED
   if (tell & TELL_PSI_SHIFT) cout << "Smoothing correlation\n";
   vector<double> aux_corr;
   aux_corr.resize(N_trials);
   for (int i=0; i<N_trials; i++) {
      int il=i-1; if (il==-1       && circular) il=N_trials-1;
      int ir=i+1; if (ir==N_trials) 
                     if (circular) ir=0; else ir=-1;
      aux_corr[i]=vcorr[i]; int N=1;
      if (il!=-1) {aux_corr[i]+=vcorr[il]; N++;}
      if (ir!=-1) {aux_corr[i]+=vcorr[ir]; N++;}
      aux_corr[i]/=N;

      if (tell & TELL_PSI_SHIFT)
         cout << " psi= " << vpsi[i] << " rot= " << vrot[i]
	      << " tilt= " << vtilt[i] << " corr= " << vcorr[i]
	      << " new corr= " << aux_corr[i] << endl; 
   }
   vcorr=aux_corr;
   matrix1D<double> aux_corr;
   aux_corr.resize(N_trials);
   for (int i=0; i<N_trials; i++) aux_corr(i)=vcorr[i];
   FourierMask Filter_corr;
   Filter_corr.FilterShape=RAISED_COSINE;
   Filter_corr.FilterBand=LOWPASS;
   Filter_corr.w1=0.15;
   Filter_corr.raised_w=0.02;
   Filter_corr.apply_mask(aux_corr);
   for (int i=0; i<N_trials; i++) {
      if (tell & TELL_PSI_SHIFT)
         cout << " psi= " << vpsi[i] << " rot= " << vrot[i]
	      << " tilt= " << vtilt[i] << " corr= " << vcorr[i]
	      << " new corr= " << aux_corr(i) << endl; 
      vcorr[i]=aux_corr(i);
   }
   #endif

   // Compute minimum and maximum of the correlation
   double max_corr=vcorr[0], min_corr=vcorr[0];
   double avg_maxima=0;
   vector<int> local_maxima, idx_all;
   idx_all.resize(N_trials);
   if (tell & TELL_PSI_SHIFT) cout << "Local maxima\n";
   for (int i=0; i<N_trials; i++) {
      idx_all[i]=i;

      // Compute maxima
      if (vcorr[i]   <min_corr   ) min_corr   =vcorr[i];
      if (vcorr[i]   >max_corr   ) max_corr   =vcorr[i];

      // Look for the left and right sample
      int il=i-1, ir=i+1;
      if      (i==0 && circular) il=N_trials-1;
      else if (i==N_trials-1)
              if (circular) ir=0; else ir=-1;

      // Check if the correlation is a local maximum
      if (il!=-1)
         if (vcorr[il]>=vcorr[i]) continue;
      if (ir!=-1)
         if (vcorr[ir]>=vcorr[i]) continue;

      // It is a maximum
      local_maxima.push_back(i);
      avg_maxima+=vcorr[i];
      if (tell & TELL_PSI_SHIFT)
         cout << "psi= " << vpsi[i] << " rot= " << vrot[i] << " tilt= "
              << vtilt[i] << " corr= " << vcorr[i] << endl;
   }
   avg_maxima/=local_maxima.size();
   if (tell & TELL_PSI_SHIFT)
      cout << "Avg_maxima=" << avg_maxima << " Max.corr=" << max_corr << " "
           << " Min.corr=" << min_corr << endl;

   // Remove all local maxima below the average
   int jmax=local_maxima.size();
   vector<int> candidate_local_maxima;
   vector<double> candidate_rate;
   if (tell & TELL_PSI_SHIFT) cout << "Keeping ...\n";
   for (int j=0; j<jmax; j++) {
      int i=local_maxima[j];
      if (vcorr[i]>=avg_maxima) {
         candidate_local_maxima.push_back(i);
	 candidate_rate.push_back(0);
	 if (tell & TELL_PSI_SHIFT)
            cout << "psi= " << vpsi[i] << " rot= " << vrot[i] << " tilt= "
        	 << vtilt[i] << " corr= " << vcorr[i] << endl;
      }
   }
   jmax=candidate_local_maxima.size();

   // Evaluate the local maxima according to their correlations
   double th_corr=evaluate_candidates_by_correlation(vcorr,
      candidate_local_maxima,candidate_rate,1);
   if (tell & TELL_PSI_SHIFT)
         cout << "Evaluation on correlation:" << candidate_rate << endl
	      << "Threshold for obtaining a 7 in correlation: " << th_corr << endl;

   // Sort the candidates
   if (tell & TELL_PSI_SHIFT) cout << "\nSelecting image\n";
   matrix1D<double> score(jmax);
   for (int j=0; j<jmax; j++) score(j)=candidate_rate[j];
   matrix1D<int> idx_score=score.index_sort();

   if (tell & (TELL_PSI_SHIFT | TELL_OPTIONS)) {
      cout << I.name() << endl;  cout.flush();
      for (int j=0; j<jmax; j++) {
	 int jp=idx_score(j)-1;
	 double score=candidate_rate[jp];
	 int i=candidate_local_maxima[jp];
         cout << "psi= " << vpsi[i] << " rot= " << vrot[i] << " tilt= "
              << vtilt[i] << " corr= " << vcorr[i]
	      << " rate=" << candidate_rate[jp] << endl;
      }
      cout << endl; cout.flush();
   }
   
   // Consider the top
   int jtop=jmax-1;
   vector<int> best_idx;
   while (jtop>0 &&
          candidate_rate[idx_score(jmax-1)-1]-
	    candidate_rate[idx_score(jtop-1)-1]<=1) {
	  best_idx.push_back(idx_score(jtop)-1);
	  jtop--;
   }
   best_idx.push_back(idx_score(jtop)-1);
   if (tell & TELL_PSI_SHIFT)
      cout << "Best indices: " << best_idx << endl;

   // Pick the best one from the top
   int ibest, jbest;
   if (jtop==jmax-1) {
      // There is only one on the top
      jbest=best_idx[0];
      ibest=candidate_local_maxima[jbest];
   } else if (jtop==jmax-2) {
      // There are two on the top
      // select that with more correlation
      if (vcorr[candidate_local_maxima[best_idx[0]]]>
          vcorr[candidate_local_maxima[best_idx[1]]]) {
	 jbest=best_idx[0];
	 ibest=candidate_local_maxima[jbest];
      } else {
	 jbest=best_idx[1];
	 ibest=candidate_local_maxima[jbest];
      }
   } else {
      // There are more than two in the top
      // Group the different views
      vector< vector<int> > groups;
      group_views(vrot,vtilt,vpsi,best_idx,candidate_local_maxima,groups);
      if (tell & TELL_PSI_SHIFT)
         cout << "Partition: " << groups << endl;
      
      // Pick the best image from the groups
      jbest=pick_view(groups,vcorr,
         best_idx,candidate_local_maxima,candidate_rate);
      ibest=candidate_local_maxima[jbest];
   }
   
   // Take the information of the best image
   best_rot    = vrot[ibest];
   best_tilt   = vtilt[ibest];
   best_psi    = vpsi[ibest];
   best_shiftX = vshiftX[ibest];
   best_shiftY = vshiftY[ibest];
   best_corr   = vcorr[ibest];
   best_rate   = candidate_rate[jbest];
   
   if (tell & (TELL_PSI_SHIFT | TELL_OPTIONS)) {
      cout << "Originally it had, psi=" << I.psi() << " rot=" << I.rot()
           << " tilt=" << I.tilt() << endl;
      cout << "Finally I choose: ";
      if (tell & TELL_PSI_SHIFT) cout << jbest << "\n";
      cout << "psi= " << best_psi << " rot= " << best_rot << " tilt= "
           << best_tilt << " shiftX=" << best_shiftX
	   << " shiftY=" << best_shiftY << " corr= " << best_corr
	   << " rate= " << best_rate << endl << endl;
   }

   // Save results
   assigned_rot    = predicted_rot[current_img]    = best_rot;
   assigned_tilt   = predicted_tilt[current_img]   = best_tilt;
   assigned_psi    = predicted_psi[current_img]    = best_psi;
   assigned_shiftX = predicted_shiftX[current_img] = best_shiftX;
   assigned_shiftY = predicted_shiftY[current_img] = best_shiftY;
                     predicted_corr[current_img]   = best_corr;
   current_img++;
   return best_rate;
}

// Finish processing ---------------------------------------------------------
void Prog_angular_predict_prm::finish_processing() {
   // Save predicted angles
   int p=predicted_rot.size();
   DocFile DF;
   DF.reserve(p+1);
   DF.append_comment("Predicted_Rot Predicted_Tilt Predicted_Psi Predicted_ShiftX Predicted_ShiftY Corr");
   matrix1D<double> v(6);
   for (int i=0; i<p; i++) {
      v(0)=predicted_rot[i];
      v(1)=predicted_tilt[i];
      v(2)=predicted_psi[i];
      v(3)=predicted_shiftX[i];
      v(4)=predicted_shiftY[i];
      v(5)=predicted_corr[i];
      DF.append_data_line(v);
   }
   DF.write(fn_out_ang);
}
