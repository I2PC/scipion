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
#include <XmippData/xmippMasks.hh>
#include <XmippData/xmippSelFiles.hh>
#include "../Prog_evaluate.hh"
#include "../../volume_FOMs.hh"

#include <fstream>

/* Compute special statistics ============================================== */
/* Average of a vector without counting the -1, starting at index i0 */
void avg_without__1(matrix1D<double> &v, double &avg, int i0) {
   int N=0; avg=0;
   for (int i=i0; i<XSIZE(v); i++)
      if (v(i)!=-1) {N++; avg += v(i);}
   if (N!=0) avg /= N; else avg=-1;
}

/* Average and standard deviation  without counting the -1 */
void stats__1(const matrix1D<double> &v, double &avg, double &stddev) {
   int N=0; avg=stddev=0;
   for (int i=0; i<XSIZE(v); i++)
      if (v(i)!=-1) {N++; avg += v(i); stddev += v(i)*v(i);}
   if (N!=0) {avg /= N; stddev=sqrt(stddev/N-avg*avg);}
   else {avg=-1; stddev=0;}
}

/* Default values ========================================================== */
void Prog_Evaluate_Parameters::default_values() {
   back_radius    = 0;
   back_factor    = 1.5;
   back_mode      = ENLARGE_MODE;
   tell           = 0;
   fn_sel         = "";
   fn_phantom     = "";
   fn_recons      = "";
   percent_mass   = 99;
   global_radius  = 0;
   fn_mask        = "";
   fit_gray_scales= FALSE;
   RSrot          = 0;
   RStilt         = 0;
}

/* Read Evaluate parameters from command line ============================== */
void Prog_Evaluate_Parameters::read(int argc, char **argv) {
   int i;

   // By default the enlarge factor is chosen
   default_values();

   // Read from command line
   fn_phantom       = get_param(argc,argv,"-p","");
   if (check_param(argc,argv,"-sel"))
      fn_sel=get_param(argc,argv,"-sel");
   else {
      fn_recons        =      get_param(argc,argv,  "-r"                );
   }
   percent_mass     = AtoF(get_param(argc,argv,  "-mass"       , "99"));
   global_radius    = AtoF(get_param(argc,argv,  "-R"          , "0"));
   fn_mask          = get_param(argc,argv,"-mask","");
   fit_gray_scales  = check_param(argc,argv,"-fit_gray");
   if (check_param(argc,argv,"-back_radius")) {
      back_radius   = AtoF(get_param(argc,argv,  "-back_radius"));
      back_mode=SPHERE_MODE;}
   if (check_param(argc,argv,"-back_factor")) {
      back_factor   = AtoF(get_param(argc,argv,  "-back_factor"));
      back_mode=ENLARGE_MODE;}
   if ((i=position_param(argc,argv,"-dir"))!=-1) {
      if ((++i) < argc) {
         if      (strcmp(argv[i],"X")==0) {RSrot= 0; RStilt=90;}
         else if (strcmp(argv[i],"Y")==0) {RSrot=90; RStilt=90;}
         else if (strcmp(argv[i],"Z")==0) {RSrot= 0; RStilt= 0;}
         else {
            RSrot=AtoF(argv[i]);
            if ((++i) < argc) RStilt = AtoF (argv[i]);
         }
      }
   }
   if (check_param(argc,argv,"-save_maps"))       tell |= SAVE_MAPS;
   if (check_param(argc,argv,"-show_values"))     tell |= SHOW_VALUES;
   if (check_param(argc,argv,"-show_process"))    tell |= SHOW_PROCESS;
   if (check_param(argc,argv,"-save_histograms")) tell |= SAVE_HISTOGRAMS;
   if (check_param(argc,argv,"-only_structural")) tell |= ONLY_STRUCTURAL;
}

/* Evaluate usage ========================================================== */
void Prog_Evaluate_Parameters::usage() {
   printf("Error in the arguments\n");
   printf("Usage: \n");
   printf("       evaluate <options>\n");
   printf("-p <phantom>             : it can be either a description or a Spider volume\n");
   printf("-r <reconstruction> |    : a Spider volume\n");
   printf("-sel <selfile>           : with all reconstruction names\n");
   printf("[-mass <%%>]             : leave out this percentage of mass in the histograms\n");
   printf("                           [99]\n");
   printf("[-R <r>]                 : The global error will be measured within this radius\n");
   printf("[-mask <surface mask>]   : surface mask applied during reconstruction\n");
   printf("[-mask]                  : only for sel files where the mask name is\n"
          "                           automatically computed\n");
   printf("[-fit_gray]              : Fit gray scales before evaluating\n");
   printf("[-dir <rot> <tilt>]      : to perform the directional FOM and the slice histograms\n");
   printf("                           a direction must be specified by two Euler angles\n");
   printf("                           by default the Z axis is taken [0 0]\n"
          "                           other useful axis are [0 90] --> X\n"
          "                           and [90 90] --> Y\n"
          "[-dir X|Y|Z]\n           : To perform directional FOMs from any of these axis\n"
          "[-back_radius <radius> |]: if this option is not given the background\n"
          " -back_factor <factor> ]   of the features is supposed to be the same\n"
          "                           feature enlarged by 1.25. With this option\n"
          "                           the background will be a sphere of the given\n"
          "                           radius, or the enlarging factor may be changed\n");
   printf("[-save_maps]             : Save different volume maps\n"
          "[-show_values]           : Show values in the features\n"
          "[-show_process]          : Show more information during calculations\n"
          "[-save_histograms]       : Save involved histograms\n"
          "[-only_structural]       : Only compute the structural consistency FOMs\n");
}

/* Show parameters ========================================================= */
ostream & operator << (ostream &out, const Prog_Evaluate_Parameters &prm) {
   out << "Evaluating parameters ----------------------\n";
   out << "Phantom         : " << prm.fn_phantom        << endl;
   out << "Reconstruction  : " << prm.fn_recons         << endl;
   out << "Percent mass    : " << prm.percent_mass      << endl;
   out << "RSrot           : " << prm.RSrot             << endl;
   out << "RStilt          : " << prm.RStilt            << endl;
   out << "Global radius   : " << prm.global_radius     << endl;
   out << "Surface mask    : " << prm.fn_mask           << endl;
   out << "Fit gray scales : " << prm.fit_gray_scales   << endl;
   out << "Back mode       : " << prm.back_mode         << endl;
   out << "Back radius     : " << prm.back_radius       << endl;
   out << "Back factor     : " << prm.back_factor       << endl;
   out << "Tell            : " << prm.tell              << endl;
   return out;
}

/* Produce Side information ================================================ */
void EVALUATE_Side_Info::produce_Side_Info(
   const Prog_Evaluate_Parameters &prm) {

   // Set background mode ..................................................
   if (prm.back_mode==SPHERE_MODE) back_param=prm.back_radius;
   else                            back_param=prm.back_factor;

   // Read reconstruction ..................................................
   fn_root=prm.fn_recons.without_extension();
   vol_recons.read(prm.fn_recons);
   vol_recons.move_origin_to_center();
   
   // Read phantom and label ...............................................
   if (Is_VolumeXmipp(prm.fn_phantom)) {
      descr_mode=XMIPP_PHANTOM;
      vol_phantom.read(prm.fn_phantom);
      vol_phantom.move_origin_to_center();
      vol_label().resize(vol_phantom());
      vol_label().init_constant(1);
      num_feat=0;
   } else {
      cerr << "Generating phantom ...\n";
      phantom_descr.read(prm.fn_phantom);
      phantom_descr.draw_in(&vol_phantom);
      phantom_descr.label(&vol_label);
      num_feat=phantom_descr.FeatNo();
      descr_mode=MATH_PHANTOM;
   }
   
   // Check that both dimensions are equal .................................
   if ((vol_phantom().SliNo()!=vol_recons().SliNo()) ||
       (vol_phantom().RowNo()!=vol_recons().RowNo()) ||
       (vol_phantom().ColNo()!=vol_recons().ColNo())) {
       cout << "Be careful!!!, volumes with different sizes\n";
       cout << "Phantom:        " << vol_phantom().SliNo() << " x " <<
          vol_phantom().RowNo() << " x " << vol_phantom().ColNo() << endl;
       cout << "Reconstruction: " << vol_recons().SliNo() << " x " <<
          vol_recons().RowNo() << " x " << vol_recons().ColNo() << endl;
      
       cut_to_common_size(vol_phantom(),vol_recons());
       cout << "Cutting to common size " << vol_phantom().SliNo() << " x " <<
          vol_phantom().RowNo() << " x " << vol_phantom().ColNo() << endl;
       
       cut_to_common_size(vol_label(),vol_recons());
   }

   // Generate global mask .................................................
   cerr << "Generating mask ...\n";
   if (prm.fn_mask!="") {
      vol_mask.read(prm.fn_mask);
      invert_binary_mask(vol_mask());
      vol_mask().set_Xmipp_origin();

      // Find minimum and maximum plane used for the surface
      int ztop=STARTINGZ(vol_mask());
      int zbottom=FINISHINGZ(vol_mask());
      for (int i=STARTINGY(vol_mask()); i<=FINISHINGY(vol_mask()); i++)
         for (int j=STARTINGX(vol_mask()); j<=FINISHINGX(vol_mask()); j++) {
            int state;
            state=0;
            for (int k=STARTINGZ(vol_mask()); k<=FINISHINGZ(vol_mask()); k++) {
                if (state==0 && VOLVOXEL(vol_mask,k,i,j)==1) {
                   ztop=MAX(ztop,k); state=1;
                } else if (state==1 && VOLVOXEL(vol_mask,k,i,j)==0) {
                   zbottom=MIN(zbottom,k); state=2;
                }
            }
         }

      // Now replace all mask values outside these two by 0
      FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_mask())
         if (k<ztop || k>zbottom) VOLVOXEL(vol_mask,k,i,j)=0;
   } else {
      vol_mask().resize(vol_recons());
      vol_mask().init_constant(1);
   }

   if (prm.global_radius!=0) {
      float global_radius2=prm.global_radius*prm.global_radius;
      FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(vol_mask)) {
         float r2=k*k+i*i+j*j;
         if (r2>global_radius2) VOLVOXEL(vol_mask,k,i,j)=0;
      }
   }

   // Fit gray values ......................................................
   if (prm.fit_gray_scales) {
      cerr << "Fitting gray values ...\n";
      range_adjust_within_mask(&(vol_mask()), vol_phantom(), vol_recons());
   }

   // Computing distance map ...............................................
   if (!((prm.tell&ONLY_STRUCTURAL) || descr_mode==XMIPP_PHANTOM)) {
      cerr << "Computing distance map ...\n";
      compute_distance_map(&vol_label, phantom_descr, &vol_mask,
         &vol_distance);
   }
}

/* Compute FOMs ============================================================ */
void compute_FOMs(const Prog_Evaluate_Parameters &prm,
   EVALUATE_Side_Info &side, EVALUATE_results &results) {
   matrix1D<double> feat_voxels;
   histogram1D hist_recons;

/* Structural consistency FOMs --------------------------------------------- */
   // Global measures
   cerr << "Computing global structural consistency ...\n";
   compute_sc_FOMs(&side.vol_phantom,&side.vol_recons,&side.vol_label,
      &side.vol_mask, -1, results.scL2_FOM, results.scL1_FOM,
      results.scmu_FOM, results.scdev_FOM, results.scrange_FOM,
      results.sccorr_FOM, results.scinf_FOM,
      prm.tell&SHOW_PROCESS);
   #define COMPUTE_THROUGH_SINGLE_VALUE
   #ifdef COMPUTE_THROUGH_SINGLE_VALUE
      compute_resolution(side.vol_phantom,side.vol_recons,results.resol_FOM);
   #else
      matrix1D<double> frequency, FSC;
      results.resol_FOM=compute_FSC(side.vol_phantom,side.vol_recons,
         1, frequency, FSC);
   #endif

   // Local measures
   results.scL2_FOMs.init_zeros(side.num_feat+1); // 0, 1, ..., FeatNo()
   results.scL2_FOMs.init_constant(-1);
   results.scL1_FOMs=results.scmu_FOMs=results.scdev_FOMs=
      results.scrange_FOMs=results.sccorr_FOMs=results.scinf_FOMs=
      results.scL2_FOMs;
   if (side.descr_mode==MATH_PHANTOM) {
      cerr << "Computing Local structural consistency ...\n";
      for (int i=0; i<=side.num_feat; i++)
         compute_sc_FOMs(&side.vol_phantom,&side.vol_recons,&side.vol_label,
            &side.vol_mask, i, results.scL2_FOMs(i), results.scL1_FOMs(i),
            results.scmu_FOMs(i), results.scdev_FOMs(i),
            results.scrange_FOMs(i), results.sccorr_FOMs(i),
	    results.scinf_FOMs(i), prm.tell&SHOW_PROCESS);
   }

   // Weighted L2 and L1 measure
   if (side.descr_mode==MATH_PHANTOM) {
      compute_voxels_in_feat(&side.vol_label, feat_voxels);
      double voxels_inside_feat=0;
      for (int i=0; i<XSIZE(feat_voxels); i++)
          voxels_inside_feat += feat_voxels(i);
      double avg_scL2f, avg_scL1f;
      avg_without__1(results.scL2_FOMs,avg_scL2f,1);
      avg_without__1(results.scL1_FOMs,avg_scL1f,1);
      results.scL2w_FOM=0.5*(feat_voxels.sum()*(results.scL2_FOMs(0)/feat_voxels(0)+
         avg_scL2f/voxels_inside_feat));
      results.scL1w_FOM=0.5*(feat_voxels.sum()*(results.scL1_FOMs(0)/feat_voxels(0)+
         avg_scL1f/voxels_inside_feat));
   } else {
      results.scL2w_FOM=results.scL1w_FOM=-1;
   }

   results.hsmu_FOMs.resize(side.num_feat+1); // 0, 1, ..., FeatNo()
   results.hsmu_FOMs.init_constant(-1);
   results.hsvr_FOMs=results.hsdt_FOMs=results.hsbr_FOMs=results.hsmu_FOMs;
   results.hsmu_FOM=results.hsdt_FOM=results.hsbr_FOM=results.hsvr_FOM=-1;
   results.drrt_FOM=results.dsbl_FOM=results.dsad_FOM=-1;

   if (!(prm.tell&ONLY_STRUCTURAL)) {
   /* Histogram based FOMs ------------------------------------------------- */
      // Local measures
      if (side.descr_mode==MATH_PHANTOM) {
         // FOM for each feature
         cerr << "Computing Histogram based FOMs ...\n";
         for (int i=1; i<=side.num_feat; i++)
            compute_hs_FOMs(&side.vol_phantom,&side.vol_recons,&side.vol_label,
               &side.vol_mask, i, side.phantom_descr, prm.back_mode,
               side.back_param, results.hsmu_FOMs(i), results.hsbr_FOMs(i),
               results.hsdt_FOMs(i), results.hsvr_FOMs(i),
               prm.tell&(SHOW_PROCESS | SAVE_HISTOGRAMS),
               side.fn_root+"_eval_histog.plot");

         // Global FOM
         avg_without__1(results.hsmu_FOMs, results.hsmu_FOM, 1);
         avg_without__1(results.hsbr_FOMs, results.hsbr_FOM, 1);
         avg_without__1(results.hsdt_FOMs, results.hsdt_FOM, 1);
         avg_without__1(results.hsvr_FOMs, results.hsvr_FOM, 1);
      }

   /* Directional FOMs ----------------------------------------------------- */
   // Global measures
      cerr << "Computing directional FOMs ...\n";
      compute_hist(side.vol_recons(),hist_recons,200);
      double threshold=hist_recons.percentil(prm.percent_mass);

      compute_dr_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_mask,
         prm.RSrot, prm.RStilt,
         &results.img_histog, 100, threshold,results.drrt_FOM,
         prm.tell&SHOW_PROCESS,side.fn_root+"_eval_radon.plot");
      if (prm.tell & SAVE_HISTOGRAMS)
         results.img_histog.write(side.fn_root+"_eval_slice_histog.xmp");

   /* Distance map based --------------------------------------------------- */
      if (side.descr_mode==MATH_PHANTOM) {
         cerr << "Computing distance based FOMs ...\n";
         compute_ds_FOMs(&side.vol_phantom, &side.vol_recons, &side.vol_label,
            &side.vol_distance, results.dsbl_FOM, results.dsad_FOM);

         if (prm.tell&SHOW_PROCESS)
            show_shape(&side.vol_phantom, &side.vol_recons, &side.vol_label,
               &side.vol_distance, side.fn_root+"_eval_shape.plot");

      }
   }
}

/* Show FOMs =============================================================== */
void show_FOMs(const Prog_Evaluate_Parameters &prm,
   EVALUATE_Side_Info &side, const EVALUATE_results &results) {

   // Show Parameters ......................................................
   cout << endl;
   cout << "PHANTOM FILE      : " << prm.fn_phantom << endl;
   cout << "RECONSTRUCTED FILE: " << prm.fn_recons  << endl;
   if (!(prm.tell&ONLY_STRUCTURAL)) {
      cout << "Direction for dFOM: (rot=" << prm.RSrot << "," << "tilt="
           << prm.RStilt << ")" << endl;
      cout << "Slice histograms  : ignoring initial " << prm.percent_mass <<
           "% mass\n";
   }
   if (prm.global_radius==0)
      cout << "Global FOMs measured over the whole volume\n";
   else
      cout << "Global FOMs measured over a sphere of radius "
           << prm.global_radius << endl;
   if (prm.back_mode==ENLARGE_MODE)
      cout << "Background mode: ENLARGE by " << prm.back_factor << endl;
   else
      cout << "Background mode: SPHERES of radius " << prm.back_radius << endl;

   if (side.descr_mode==MATH_PHANTOM && (prm.tell & SHOW_PROCESS)) {
      cout << "Phantom description ----------------------------------------\n";
      cout << side.phantom_descr;
   }

   cout << "Structural consistency -------------------------------------\n";

   // Show volume statistics ...............................................
   double avg, stddev, min, max;
   side.vol_phantom().compute_stats(avg, stddev, min, max);
   cout << "Phantom Stats: \n";
   cout << "   "; side.vol_phantom().print_stats(); cout << endl;
   cout << "    range=" << max-min << endl;
   side.vol_recons().compute_stats(avg, stddev, min, max);
   cout << "Recons  Stats: \n";
   cout << "   "; side.vol_recons().print_stats(); cout << endl;
   cout << "    range=" << max-min << endl; cout << endl;
   
   // Show Structural consistency ..........................................
   printf("scL2        FOM:%f\n",results.scL2_FOM);
   printf("scL1        FOM:%f\n",results.scL1_FOM);
   printf("scL2w       FOM:%f\n",results.scL2w_FOM);
   printf("scL1w       FOM:%f\n",results.scL1w_FOM);
   printf("scmu        FOM:%f\n",results.scmu_FOM);
   printf("scdev       FOM:%f\n",results.scdev_FOM);
   printf("scrange     FOM:%f\n",results.scrange_FOM);
   printf("sccorr      FOM:%f\n",results.sccorr_FOM);
   printf("scinf       FOM:%f\n",results.scinf_FOM);
   printf("resolution  FOM:%f\n",results.resol_FOM);
   
   printf("\tFEATURE   scL2     scL1     scmu     scdev   scrange  sccorr    scinf \n");
   printf("\t------- -------- -------- -------- -------- -------- -------- --------\n");
   for (int i=0; i<=side.num_feat; i++) {
      printf("\t  %2d     %1.4f   %1.4f   %1.4f   %1.4f   %1.4f   %1.4f   %1.4f\n",
         i, results.scL2_FOMs(i), results.scL1_FOMs(i), results.scmu_FOMs(i),
         results.scdev_FOMs(i), results.scrange_FOMs(i), results.sccorr_FOMs(i),
	 results.scinf_FOMs(i));
   }
   cout << endl;

   if (!(prm.tell&ONLY_STRUCTURAL)) {
   // Show histogram based FOMS ............................................
   cout << "Histogram based --------------------------------------------\n";
   printf("hsin        FOM:%f\n",results.hsmu_FOM);
   printf("hsbr        FOM:%f\n",results.hsbr_FOM);
   printf("hsdt        FOM:%f\n",results.hsdt_FOM);
   printf("hsvr        FOM:%f\n",results.hsvr_FOM);
   printf("\tFEATURE    hsin       hsbr       hsdt       hsvr\n");
   printf("\t------- ---------- ---------- ---------- ----------\n");
   for (int i=0; i<=side.num_feat; i++) {
      printf("\t  %2d     % 7.2f    % 7.2f    % 7.2f    % 7.2f\n",
         i, results.hsmu_FOMs(i),  results.hsbr_FOMs(i),
         results.hsdt_FOMs(i), results.hsvr_FOMs(i));
   }
   cout << endl;

   // Show Directional FOMs ................................................
   cout << "Directional ------------------------------------------------\n";
   printf("scrt        FOM:%f\n",results.drrt_FOM);
   
   // Show Distance FOMs ...................................................
   cout << "Distance based ---------------------------------------------\n";
   printf("scbl        FOM:%f\n",results.dsbl_FOM);
   printf("scad        FOM:%f\n",results.dsad_FOM);
   }
   
   // Save maps ............................................................
   if (prm.tell&SAVE_MAPS) {
      VolumeXmipp save, error;
      
      // Save generated phantom
      cerr << "Saving generated phantom ...\n";
      side.vol_phantom.write(side.fn_root+"_eval_phantom.vol");

      // Save mask
      cerr << "Saving evaluation mask ...\n";
      side.vol_mask.write(side.fn_root+"_eval_mask.vol");

      // Save label map
      cerr << "Saving label map ...\n";
      side.vol_label.write(side.fn_root+"_eval_label.vol");

      // Save a map of differences
      save()=error()=side.vol_phantom()-side.vol_recons();
      if (side.descr_mode==MATH_PHANTOM) side.phantom_descr.sketch_in(&save);
      cerr << "Saving difference map ...\n";
      save.write(side.fn_root+"_eval_difference_map.vol");
      
      // Save a map of quadratic errors
      save()=error()=error()*error();
      if (side.descr_mode==MATH_PHANTOM) side.phantom_descr.sketch_in(&save);
      cerr << "Saving quadratic errors ...\n";
      save.write(side.fn_root+"_eval_quadratic_map.vol");
      
      // Save a absolute difference map
      save()=ABSnD(side.vol_phantom()-side.vol_recons());
      if (side.descr_mode==MATH_PHANTOM) side.phantom_descr.sketch_in(&save);
      cerr << "Saving absolute difference map ...\n";
      save.write(side.fn_root+"_eval_absolute_map.vol");

      if (!(prm.tell&ONLY_STRUCTURAL)) {
      // Save distance map
         cerr << "Saving distance map ...\n";
         side.vol_distance.write(side.fn_root+"_eval_distance_map.vol");

      // Save blurring map
         save().resize(side.vol_distance());
         FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(side.vol_distance))
            if (side.vol_distance(k,i,j)==-1) save(k,i,j)=0;
            else save(k,i,j)=1/side.vol_distance(k,i,j)*error(k,i,j);
         cerr << "Saving blurring map ...\n";
         save.write(side.fn_root+"_eval_blurring_map.vol");

      // Save appearing map
         FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(side.vol_distance))
            if (side.vol_distance(k,i,j)==-1) save(k,i,j)=0;
            else save(k,i,j)=side.vol_distance(k,i,j)*error(k,i,j);
         cerr << "Saving appearence map ...\n";
         save.write(side.fn_root+"_eval_appearing_map.vol");
      }
   }

   // Show values ..........................................................
   if (prm.tell&SHOW_VALUES) {
      int sel_feat=0;
      string fn_out;
      ofstream fh_out;
      
      cout << "Name of the filename to dump values: ";
      cin >> fn_out;
      fh_out.open(fn_out.c_str(), ios::out);
      if (!fh_out)
         REPORT_ERROR(3005,(string)"Evaluate show: Could not open " + fn_out
            + " for output");
      
      while (sel_feat!=1000) {
         cout << "What feature do you want to see (" << -side.num_feat << ","
              << side.num_feat << ") (1000=finish): ";
         cin  >> sel_feat;
         if (ABS(sel_feat)<=side.num_feat)
            switch(side.descr_mode) {
               case XMIPP_PHANTOM:
                  show_voxels_in_feat(&side.vol_phantom, &side.vol_recons,
                     &side.vol_label, sel_feat, fh_out);
                  break;
               case MATH_PHANTOM:
                  show_voxels_in_feat(&side.vol_phantom, &side.vol_recons,
                     &side.vol_label, side.phantom_descr, sel_feat, fh_out);
                  break;
            }
      }
      fh_out.close();
   }
}

/* Single step ============================================================= */
void Evaluate_single_step(const Prog_Evaluate_Parameters &prm,
   EVALUATE_results &results) {
   cout << prm;

// Read volumes, label them and generate global mask
   EVALUATE_Side_Info side;
   side.produce_Side_Info(prm);

// Compute FOMs
   compute_FOMs(prm, side, results);
   
// Show results
   show_FOMs(prm, side, results);
}

/* Main routine ============================================================ */
void ROUT_Evaluate(Prog_Evaluate_Parameters &prm,
   EVALUATE_results &results) {
   if (prm.fn_sel=="") Evaluate_single_step(prm,results);
   else {
      SelFile SF(prm.fn_sel);
      SelFile SFmask;
      if (prm.fn_mask!="") SFmask.read(prm.fn_mask);
      FOMs foms(SF.ImgNo()), foms_mean(1), foms_stddev(1);
      int k=0;
      bool mathematical_phantom;
      while (!SF.eof()) {
         cerr << "Perfoming measure for test number " << k << endl;
         prm.fn_recons=SF.NextImg();
         if (prm.fn_phantom=="") {
            mathematical_phantom=TRUE;
            prm.fn_phantom=prm.fn_recons.without_extension();
            prm.fn_phantom=prm.fn_phantom.without("_wos");
            int i=prm.fn_phantom.find("_idr");
            if (i!=-1) prm.fn_phantom=prm.fn_phantom.replace(i,6,"");
            prm.fn_phantom +=".descr";
         } else mathematical_phantom=FALSE;
         if (!SFmask.eof()) prm.fn_mask=SFmask.NextImg();
         else               prm.fn_mask="";
         Evaluate_single_step(prm,results);
         foms.set_FOMs(k,results);
         k++;
         if (mathematical_phantom) {
            prm.fn_phantom="";
            mathematical_phantom=FALSE;
         }
      }
      compute_FOMs_stats(foms,0,foms_mean,foms_stddev);
      cout << foms;
      cout << "--------------------------------------------------\n";
      cout << "After combining " << SF.ImgNo() << " volumes\n";
      show_stats(cout,0,foms_mean,foms_stddev);
   }
}

/* FOMs constructor ======================================================== */
FOMs::FOMs(int n) {
   scL2.resize(n);
   scL1.resize(n);
   scL2w.resize(n);
   scL1w.resize(n);
   scmu.resize(n);
   scdev.resize(n);
   scrange.resize(n);
   sccorr.resize(n);
   scinf.resize(n);
   scresol.resize(n);
   scL20.resize(n);
   scL10.resize(n);
   scmu0.resize(n);
   scdev0.resize(n);
   scrange0.resize(n);
   scL21.resize(n);
   scL11.resize(n);
   scmu1.resize(n);
   scdev1.resize(n);
   scrange1.resize(n);
   hsvr.resize(n);
   hsmu.resize(n);
   hsbr.resize(n);
   hsdt.resize(n);
   drrt.resize(n);
   dsbl.resize(n);
   dsad.resize(n);
}

/* Set ===================================================================== */
void FOMs::set_FOMs(int k, EVALUATE_results &results) {
   scL2(k)     = results.scL2_FOM;
   scL1(k)     = results.scL1_FOM;
   scL2w(k)    = results.scL2w_FOM;
   scL1w(k)    = results.scL1w_FOM;
   scmu(k)     = results.scmu_FOM;
   scdev(k)    = results.scdev_FOM;
   scrange(k)  = results.scrange_FOM;
   sccorr(k)   = results.sccorr_FOM;
   scinf(k)    = results.scinf_FOM;
   scresol(k)  = results.resol_FOM;
   scL20(k)    = results.scL2_FOMs(0);
   scL10(k)    = results.scL1_FOMs(0);
   scmu0(k)    = results.scmu_FOMs(0);
   scdev0(k)   = results.scdev_FOMs(0);
   scrange0(k) = results.scrange_FOMs(0);
   if (XSIZE(results.scL2_FOMs)!=1)    scL21(k)    = results.scL2_FOMs(1);
   if (XSIZE(results.scL1_FOMs)!=1)    scL11(k)    = results.scL1_FOMs(1);
   if (XSIZE(results.scmu_FOMs)!=1)    scmu1(k)    = results.scmu_FOMs(1);
   if (XSIZE(results.scdev_FOMs)!=1)   scdev1(k)   = results.scdev_FOMs(1);
   if (XSIZE(results.scrange_FOMs)!=1) scrange1(k) = results.scrange_FOMs(1);
   hsvr(k)     = results.hsvr_FOM;
   hsmu(k)     = results.hsmu_FOM;
   hsbr(k)     = results.hsbr_FOM;
   hsdt(k)     = results.hsdt_FOM;
   drrt(k)     = results.drrt_FOM;
   dsbl(k)     = results.dsbl_FOM;
   dsad(k)     = results.dsad_FOM;
}

/* Compute stats =========================================================== */
void compute_FOMs_stats(const FOMs &foms, int i, FOMs &fmean, FOMs &fstddev) {
   double avg, stddev;
   stats__1(foms.scL2,avg,stddev);     fmean.scL2(i)    = avg; fstddev.scL2(i)    = stddev;
   stats__1(foms.scL1,avg,stddev);     fmean.scL1(i)    = avg; fstddev.scL1(i)    = stddev;
   stats__1(foms.scL2w,avg,stddev);    fmean.scL2w(i)   = avg; fstddev.scL2w(i)   = stddev;
   stats__1(foms.scL1w,avg,stddev);    fmean.scL1w(i)   = avg; fstddev.scL1w(i)   = stddev;
   stats__1(foms.scmu,avg,stddev);     fmean.scmu(i)    = avg; fstddev.scmu(i)    = stddev;
   stats__1(foms.scdev,avg,stddev);    fmean.scdev(i)   = avg; fstddev.scdev(i)   = stddev;
   stats__1(foms.scrange,avg,stddev);  fmean.scrange(i) = avg; fstddev.scrange(i) = stddev;
   stats__1(foms.sccorr,avg,stddev);   fmean.sccorr(i)  = avg; fstddev.sccorr(i)  = stddev;
   stats__1(foms.scinf,avg,stddev);    fmean.scinf(i)   = avg; fstddev.scinf(i)   = stddev;
   stats__1(foms.scresol,avg,stddev);  fmean.scresol(i) = avg; fstddev.scresol(i) = stddev;
   stats__1(foms.scL20,avg,stddev);    fmean.scL20(i)   = avg; fstddev.scL20(i)   = stddev;
   stats__1(foms.scL10,avg,stddev);    fmean.scL10(i)   = avg; fstddev.scL10(i)   = stddev;
   stats__1(foms.scmu0,avg,stddev);    fmean.scmu0(i)   = avg; fstddev.scmu0(i)   = stddev;
   stats__1(foms.scdev0,avg,stddev);   fmean.scdev0(i)  = avg; fstddev.scdev0(i)  = stddev;
   stats__1(foms.scrange0,avg,stddev); fmean.scrange0(i)= avg; fstddev.scrange0(i)= stddev;
   stats__1(foms.scL21,avg,stddev);    fmean.scL21(i)   = avg; fstddev.scL21(i)   = stddev;
   stats__1(foms.scL11,avg,stddev);    fmean.scL11(i)   = avg; fstddev.scL11(i)   = stddev;
   stats__1(foms.scmu1,avg,stddev);    fmean.scmu1(i)   = avg; fstddev.scmu1(i)   = stddev;
   stats__1(foms.scdev1,avg,stddev);   fmean.scdev1(i)  = avg; fstddev.scdev1(i)  = stddev;
   stats__1(foms.scrange1,avg,stddev); fmean.scrange1(i)= avg; fstddev.scrange1(i)= stddev;
  
   stats__1(foms.hsvr,avg,stddev);     fmean.hsvr(i)    = avg; fstddev.hsvr(i)    = stddev;
   stats__1(foms.hsmu,avg,stddev);     fmean.hsmu(i)    = avg; fstddev.hsmu(i)    = stddev;
   stats__1(foms.hsbr,avg,stddev);     fmean.hsbr(i)    = avg; fstddev.hsbr(i)    = stddev;
   stats__1(foms.hsdt,avg,stddev);     fmean.hsdt(i)    = avg; fstddev.hsdt(i)    = stddev;
  
   stats__1(foms.drrt,avg,stddev);     fmean.drrt(i)    = avg; fstddev.drrt(i)    = stddev;
  
   stats__1(foms.dsbl,avg,stddev);     fmean.dsbl(i)    = avg; fstddev.dsbl(i)    = stddev;
   stats__1(foms.dsad,avg,stddev);     fmean.dsad(i)    = avg; fstddev.dsad(i)    = stddev;
}

/* Show ==================================================================== */
ostream & operator << (ostream &out, const FOMs &foms) {
   out << "All results in all tests\n";
   out << "--------------------------------------------------\n";   
   out << "scL2\n"     << foms.scL2.transpose()     << endl;
   out << "scL1\n"     << foms.scL1.transpose()     << endl;
   out << "scL2w\n"    << foms.scL2w.transpose()    << endl;
   out << "scL1w\n"    << foms.scL1w.transpose()    << endl;
   out << "scmu\n"     << foms.scmu.transpose()     << endl;
   out << "scdev\n"    << foms.scdev.transpose()    << endl;
   out << "scrange\n"  << foms.scrange.transpose()  << endl;
   out << "sccorr\n"   << foms.sccorr.transpose()   << endl;
   out << "scinf\n"    << foms.scinf.transpose()    << endl;
   out << "scresol\n"  << foms.scresol.transpose()  << endl;
   out << "scL20\n"    << foms.scL20.transpose()    << endl;
   out << "scL10\n"    << foms.scL10.transpose()    << endl;
   out << "scmu0\n"    << foms.scmu0.transpose()    << endl;
   out << "scdev0\n"   << foms.scdev0.transpose()   << endl;
   out << "scrange0\n" << foms.scrange0.transpose() << endl;
   out << "scL21\n"    << foms.scL21.transpose()    << endl;
   out << "scL11\n"    << foms.scL11.transpose()    << endl;
   out << "scmu1\n"    << foms.scmu1.transpose()    << endl;
   out << "scdev1\n"   << foms.scdev1.transpose()   << endl;
   out << "scrange1\n" << foms.scrange1.transpose() << endl;
  
   out << "hsvr\n"    << foms.hsvr.transpose()    << endl;
   out << "hsmu\n"    << foms.hsmu.transpose()    << endl;
   out << "hsbr\n"    << foms.hsbr.transpose()    << endl;
   out << "hsdt\n"    << foms.hsdt.transpose()    << endl;
  
   out << "drrt\n"    << foms.drrt.transpose()    << endl;
  
   out << "dsbl\n"    << foms.dsbl.transpose()    << endl;
   out << "dsad\n"    << foms.dsad.transpose()    << endl;
   out << "--------------------------------------------------\n";
   return out;
}

/* Show stats ============================================================== */
void show_stats(ostream &out, int i, const FOMs &fmean,
   const FOMs &fstddev) {
   out << "    scL2:     " << fmean.scL2(i)     << "+-" << fstddev.scL2(i)    << endl;
   out << "    scL1:     " << fmean.scL1(i)     << "+-" << fstddev.scL1(i)    << endl;
   out << "    scL2w:    " << fmean.scL2w(i)    << "+-" << fstddev.scL2w(i)   << endl;
   out << "    scL1w:    " << fmean.scL1w(i)    << "+-" << fstddev.scL1w(i)   << endl;
   out << "    scmu:     " << fmean.scmu(i)     << "+-" << fstddev.scmu(i)    << endl;
   out << "    scdev:    " << fmean.scdev(i)    << "+-" << fstddev.scdev(i)   << endl;
   out << "    scrange:  " << fmean.scrange(i)  << "+-" << fstddev.scrange(i) << endl;
   out << "    sccorr:   " << fmean.sccorr(i)   << "+-" << fstddev.sccorr(i)  << endl;
   out << "    scinf:    " << fmean.scinf(i)    << "+-" << fstddev.scinf(i)   << endl;
   out << "    scresol:  " << fmean.scresol(i)  << "+-" << fstddev.scresol(i) << endl;
   out << "    scL20:    " << fmean.scL20(i)    << "+-" << fstddev.scL20(i)   << endl;
   out << "    scL10:    " << fmean.scL10(i)    << "+-" << fstddev.scL10(i)   << endl;
   out << "    scmu0:    " << fmean.scmu0(i)    << "+-" << fstddev.scmu0(i)   << endl;
   out << "    scdev0:   " << fmean.scdev0(i)   << "+-" << fstddev.scdev0(i)  << endl;
   out << "    scrange0: " << fmean.scrange0(i) << "+-" << fstddev.scrange0(i)<< endl;
   out << "    scL21:    " << fmean.scL21(i)    << "+-" << fstddev.scL21(i)   << endl;
   out << "    scL11:    " << fmean.scL11(i)    << "+-" << fstddev.scL11(i)   << endl;
   out << "    scmu1:    " << fmean.scmu1(i)    << "+-" << fstddev.scmu1(i)   << endl;
   out << "    scdev1:   " << fmean.scdev1(i)   << "+-" << fstddev.scdev1(i)  << endl;
   out << "    scrange1: " << fmean.scrange1(i) << "+-" << fstddev.scrange1(i)<< endl;
   out << "    scbl:     " << fmean.dsbl(i)     << "+-" << fstddev.dsbl(i)    << endl;
   out << "    scap:     " << fmean.dsad(i)     << "+-" << fstddev.dsad(i)    << endl;
   out << "    scrt:     " << fmean.drrt(i)     << "+-" << fstddev.drrt(i)    << endl;
   out << "    hsin:     " << fmean.hsmu(i)     << "+-" << fstddev.hsmu(i)    << endl;
   out << "    hsbr:     " << fmean.hsbr(i)     << "+-" << fstddev.hsbr(i)    << endl;
   out << "    hsdt:     " << fmean.hsdt(i)     << "+-" << fstddev.hsdt(i)    << endl;
   out << "    hsvr:     " << fmean.hsvr(i)     << "+-" << fstddev.hsvr(i)    << endl;
   out.flush();
}
