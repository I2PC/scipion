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
//Sun Nov 14 22:07:48 EST 1999: added binary option (R. Marabini)

#include "../Prog_random_phantom.hh"
#include <XmippData/xmippHistograms.hh>
#include "../Prog_FourierFilter.hh"

/* Empty constructor ======================================================= */
Prog_Random_Phantom_Parameters::Prog_Random_Phantom_Parameters() {
   fn_random=fn_output=fn_CTF="";
   min_vol=0;
   discrete=FALSE;
   RPP_distance=RPP_radius=-1;
   N_stats=-1;
   target_SNR=-1;
}

/* Read Random Phantom parameters ========================================== */
void Prog_Random_Phantom_Parameters::read(int argc, char **argv) {
    min_vol       = AtoF(get_param(argc, argv, "-min_volume","0"));
    discrete      =    check_param(argc, argv, "-discrete");
    RPP_distance  = AtoF(get_param(argc, argv, "-distance" , "-1.") );
    RPP_radius    = AtoF(get_param(argc, argv, "-radius"   , "-1.") );
    fn_random     =      get_param(argc, argv, "-i");
    fn_output     =      get_param(argc, argv, "-o","");
    N_stats       = AtoI(get_param(argc, argv, "-stats",      "-1") );
    fn_CTF        =      get_param(argc, argv, "-ctf","");
    Xdim          = AtoI(get_param(argc, argv, "-Xdim",       "-1") );
    Ydim          = AtoI(get_param(argc, argv, "-Ydim",       "-1") );
    target_SNR    = AtoF(get_param(argc, argv, "-target_SNR", "-1.") );
}

/* Usage =================================================================== */
void Prog_Random_Phantom_Parameters::usage() {
   cout << "Usage:\n";
   cout << "random_phantom [Options and Parameters]\n";
   cout << "[Options and Parameters]:\n";
   cout << "   -i <Random      phantom description file | Xmipp Volume>\n";
   cout << "                      The Xmipp volume can only be used for stats\n";
   cout << "  [-o <Realization phantom description file>]\n";
   cout << "  [-min_volume <min_vol>]                     minimum volume\n"
        << "                                              for each feature\n";
   cout << "  [-discrete]         density values are forced to be discrete\n";
   cout << "  [-distance d]       distance between feature centers\n";
   cout << "                      greater than d\n";
   cout << "  [-radius r]         feature centers inside an sphere of radius r\n";
   cout << "  [-stats <N=-1>]     Compute family volume statistics with N samples\n";
   cout << "  [-ctf <Xmipp Fourier Image>] For computing the projection statistics\n";
   cout << "  [-Xdim <dim>        For computing the projection statistics\n";
   cout << "  [-Ydim <dim=Xdim>]  For computing the projection statistics\n";
   cout << "  [-target_SNR <SNR>] For computing the noise power needed for this SNR\n";
}

/* Produce Side Information ================================================ */
void Random_Phantom_Side_Info::produce_Side_Info(
   const Prog_Random_Phantom_Parameters &prm) {
// Read phantom file
   if (Is_VolumeXmipp(prm.fn_random)) {
      voxel_mode=TRUE;
      VoxelPhantom.read(prm.fn_random);
      VoxelPhantom().set_Xmipp_origin();
   } else {
      voxel_mode=FALSE;
      Random.read(prm.fn_random);
   
      // Check that it meets the conditions
      if (Random.FeatNo()%2!=0)
         EXIT_ERROR(1,(string)"Random_phantom: The number of features in "+
            prm.fn_random+" is not a multiple of 2");

      if (Random.FeatNo()==0)
         EXIT_ERROR(1,(string)"Random_phantom: There is no phantom in "+
            prm.fn_random);

      for (int i=1; i<=Random.FeatNo(); i+=2) {
         if (Random(i)->Type!=Random(i+1)->Type)
            EXIT_ERROR(1,(string)"Random_phantom: Feature number " + ItoA(i) +
               " is not of the same type as " + ItoA(i+1));
         if (Random(i)->Add_Assign!=Random(i+1)->Add_Assign)
            EXIT_ERROR(1,(string)"Random_phantom: Feature number " + ItoA(i) +
               " is not of the same +/= behaviour as " + ItoA(i+1));
      }
   }
}

/* Generate a realization of the random phantom  --------------------------- */
void generate_realization_of_random_phantom(
   const Prog_Random_Phantom_Parameters &prm,
   Random_Phantom_Side_Info &side, Phantom &Realization) {
   Feature     *feat;
   int j;
   int Valid_Feature;
   int loop_conunter;
// double distance2;
// distance2 *= fabs(prm.RPP_distance);
    
   // Global characteristics ...............................................
   Realization.clear();
   Realization.xdim               = side.Random.xdim;
   Realization.ydim               = side.Random.ydim;
   Realization.zdim               = side.Random.zdim;
   Realization.Background_Density = side.Random.Background_Density;
   
   randomize_random_generator();
   for (int i=1; i<=side.Random.FeatNo(); i+=2) {
       // side.Random sphere ...............................................
       Feature *feat_aux;
       feat_aux=NULL;
       loop_conunter=0;
       do {
          if (feat_aux!=NULL) delete feat_aux;
          if        (side.Random(i)->Type=="sph") {
             Sphere *sph=new Sphere;
             feat_aux=sph;
             Sphere *sph_ptr_i  =(Sphere *) side.Random(i);
             Sphere *sph_ptr_i_1=(Sphere *) side.Random(i+1);
             sph->init_rnd(
                   sph_ptr_i->radius,     sph_ptr_i_1->radius,
                   sph_ptr_i->Density,    sph_ptr_i_1->Density,
                XX(sph_ptr_i->Center), XX(sph_ptr_i_1->Center),
                YY(sph_ptr_i->Center), YY(sph_ptr_i_1->Center),
                ZZ(sph_ptr_i->Center), ZZ(sph_ptr_i_1->Center));
          // side.Random cylinder ..........................................
          } else if (side.Random(i)->Type=="cyl") {
             Cylinder *cyl=new Cylinder;
             feat_aux=cyl;
             Cylinder *cyl_ptr_i  =(Cylinder *) side.Random(i);
             Cylinder *cyl_ptr_i_1=(Cylinder *) side.Random(i+1);
             cyl->init_rnd(
                   cyl_ptr_i->xradius,    cyl_ptr_i_1->xradius,
                   cyl_ptr_i->yradius,    cyl_ptr_i_1->yradius,
                   cyl_ptr_i->height,     cyl_ptr_i_1->height,
                   cyl_ptr_i->Density,    cyl_ptr_i_1->Density,
                XX(cyl_ptr_i->Center), XX(cyl_ptr_i_1->Center),
                YY(cyl_ptr_i->Center), YY(cyl_ptr_i_1->Center),
                ZZ(cyl_ptr_i->Center), ZZ(cyl_ptr_i_1->Center),
                   cyl_ptr_i->rot,        cyl_ptr_i_1->rot,
                   cyl_ptr_i->tilt,       cyl_ptr_i_1->tilt,
                   cyl_ptr_i->psi,        cyl_ptr_i_1->psi);
          // side.Random double cylinder ...................................
          } else if (side.Random(i)->Type=="dcy") {
             DCylinder *dcy=new DCylinder;
             feat_aux=dcy;
             DCylinder *dcy_ptr_i  =(DCylinder *) side.Random(i);
             DCylinder *dcy_ptr_i_1=(DCylinder *) side.Random(i+1);
             dcy->init_rnd(
                   dcy_ptr_i->radius,     dcy_ptr_i_1->radius,
                   dcy_ptr_i->height,     dcy_ptr_i_1->height,
                   dcy_ptr_i->separation, dcy_ptr_i_1->separation,
                   dcy_ptr_i->Density,    dcy_ptr_i_1->Density,
                XX(dcy_ptr_i->Center), XX(dcy_ptr_i_1->Center),
                YY(dcy_ptr_i->Center), YY(dcy_ptr_i_1->Center),
                ZZ(dcy_ptr_i->Center), ZZ(dcy_ptr_i_1->Center),
                   dcy_ptr_i->rot,        dcy_ptr_i_1->rot,
                   dcy_ptr_i->tilt,       dcy_ptr_i_1->tilt,
                   dcy_ptr_i->psi,        dcy_ptr_i_1->psi);
          // side.Random cube ..............................................
          } else if (side.Random(i)->Type=="cub") {
             Cube *cub=new Cube;
             feat_aux=cub;
             Cube *cub_ptr_i  =(Cube *) side.Random(i);
             Cube *cub_ptr_i_1=(Cube *) side.Random(i+1);
             cub->init_rnd(
                   cub_ptr_i->xdim,       cub_ptr_i_1->xdim,
                   cub_ptr_i->ydim,       cub_ptr_i_1->ydim,
                   cub_ptr_i->zdim,       cub_ptr_i_1->zdim,
                   cub_ptr_i->Density,    cub_ptr_i_1->Density,
                XX(cub_ptr_i->Center), XX(cub_ptr_i_1->Center),
                YY(cub_ptr_i->Center), YY(cub_ptr_i_1->Center),
                ZZ(cub_ptr_i->Center), ZZ(cub_ptr_i_1->Center),
                   cub_ptr_i->rot,        cub_ptr_i_1->rot,
                   cub_ptr_i->tilt,       cub_ptr_i_1->tilt,
                   cub_ptr_i->psi,        cub_ptr_i_1->psi);
          // side.Random ellipsoid .........................................
          } else if (side.Random(i)->Type=="ell") {
             Ellipsoid *ell=new Ellipsoid;
             feat_aux=ell;
             Ellipsoid *ell_ptr_i  =(Ellipsoid *) side.Random(i);
             Ellipsoid *ell_ptr_i_1=(Ellipsoid *) side.Random(i+1);
             ell->init_rnd(
                   ell_ptr_i->xradius,    ell_ptr_i_1->xradius,
                   ell_ptr_i->yradius,    ell_ptr_i_1->yradius,
                   ell_ptr_i->zradius,    ell_ptr_i_1->zradius,
                   ell_ptr_i->Density,    ell_ptr_i_1->Density,
                XX(ell_ptr_i->Center), XX(ell_ptr_i_1->Center),
                YY(ell_ptr_i->Center), YY(ell_ptr_i_1->Center),
                ZZ(ell_ptr_i->Center), ZZ(ell_ptr_i_1->Center),
                   ell_ptr_i->rot,        ell_ptr_i_1->rot,
                   ell_ptr_i->tilt,       ell_ptr_i_1->tilt,
                   ell_ptr_i->psi,        ell_ptr_i_1->psi);
          // side.Random cone ..............................................
          } else if (side.Random(i)->Type=="con") {
             Cone *con=new Cone;
             feat_aux=con;
             Cone *con_ptr_i  =(Cone *) side.Random(i);
             Cone *con_ptr_i_1=(Cone *) side.Random(i+1);
             con->init_rnd(
                   con_ptr_i->radius,     con_ptr_i_1->radius,
                   con_ptr_i->height,     con_ptr_i_1->height,
                   con_ptr_i->Density,    con_ptr_i_1->Density,
                XX(con_ptr_i->Center), XX(con_ptr_i_1->Center),
                YY(con_ptr_i->Center), YY(con_ptr_i_1->Center),
                ZZ(con_ptr_i->Center), ZZ(con_ptr_i_1->Center),
                   con_ptr_i->rot,        con_ptr_i_1->rot,
                   con_ptr_i->tilt,       con_ptr_i_1->tilt,
                   con_ptr_i->psi,        con_ptr_i_1->psi);
          }

       Valid_Feature=1; // if 1 then feature OK

       for (j=1; j<=Realization.FeatNo(); j++)
           {
           if(Realization(j)->Density != (-0.)) 
              {
               if( ((feat_aux->Center - Realization(j)->Center).module() 
                                                       < prm.RPP_distance) 
                                                       && prm.RPP_distance > 0)

                 Valid_Feature=0;
               if(feat_aux->Center.module() > prm.RPP_radius &&
                                          prm.RPP_radius > 0) 
                 Valid_Feature=0;
              }
           }
       if(loop_conunter > MAX_LOOP_NUMBER)
          {
          cout << "\nOh My Dear after " << MAX_LOOP_NUMBER  << " iterations"
               << "\nI simply can not get a correct feature" 
               << "\nConsider relaxing the constraints" << endl;
          cout << "The troublesome feature is no: "<<i/2 << 
                  "(first one is 0)"<< feat_aux <<endl;     
          Valid_Feature=1;
          }
       loop_conunter++;   

       } while ((feat_aux->volume()<prm.min_vol && 
                prm.min_vol!=0) || (Valid_Feature==0));
       feat=feat_aux;

       // Copy common characteristics ......................................
       if (prm.discrete) feat->Density=rint(feat->Density);
       feat->Add_Assign = side.Random(i)->Add_Assign;
       
       // Store this feature ...............................................
       Realization.add(feat);
   }
}

/* Main Random Phantom Routine --------------------------------------------- */
//#define DEBUG
void ROUT_random_phantom(const Prog_Random_Phantom_Parameters &prm,
   Phantom &Realization) _THROW {

// Produce Side Information
   Random_Phantom_Side_Info side;
   side.produce_Side_Info(prm);
   
   if (prm.N_stats==-1) {
      if (side.voxel_mode)
         REPORT_ERROR(1,"Random_phantom: Cannot generate a random realization"
            " of a voxel phantom");
      // Generate realization and write to disk
      generate_realization_of_random_phantom(prm, side, Realization);
      if (prm.fn_output!="") Realization.write(prm.fn_output);
   } else {
      FourierMask ctf;
      if (prm.fn_CTF!="") {
         ctf.FilterBand=CTF;
         ctf.ctf.enable_CTFnoise=FALSE;
         ctf.ctf.read(prm.fn_CTF);
         ctf.ctf.Produce_Side_Info();
      }
   
      int Xdim=prm.Xdim, Ydim=prm.Ydim;
      if (Xdim==-1) 
         if (!side.voxel_mode) Xdim=side.Random.xdim;
         else                  Xdim=XSIZE(side.VoxelPhantom());
      if (Ydim==-1) Ydim=Xdim;

      matrix1D<double> volume(prm.N_stats);
      matrix1D<double> proj_power(prm.N_stats);
      matrix1D<double> proj_area(prm.N_stats);
      Projection proj;
      double power_avg, power_stddev, area_avg, area_stddev, avg, stddev, dummy;
      init_progress_bar(prm.N_stats);
      for (int n=0; n<prm.N_stats; n++) {
	 if (!side.voxel_mode)
            generate_realization_of_random_phantom(prm, side, Realization);

         // Compute phantom volume
	 if (!side.voxel_mode) volume(n)=Realization.volume();
         else volume(n)=
                 side.VoxelPhantom().count_threshold("above",0,0);
         
         // Compute projection
         if (!side.voxel_mode)
            Realization.project_to(proj, Ydim, Xdim,
                  rnd_unif(0,360),rnd_unif(0,180),rnd_unif(0,360));
         else
            project_Volume(side.VoxelPhantom(), proj, Ydim, Xdim,
                  rnd_unif(0,360),rnd_unif(0,180),rnd_unif(0,360));

         // Apply CTF
         if (prm.fn_CTF!="") {
            if (n==0) ctf.generate_mask(proj());
            ctf.apply_mask_Space(proj());
         }

         // Compute projection area
         proj_area(n)=0;
         FOR_ALL_ELEMENTS_IN_MATRIX2D(proj())
            if (ABS(proj(i,j))>1e-6) proj_area(n)++;
	 #ifdef DEBUG
	    cout << "Area: " << proj_area(n) << endl;
	    proj.write("inter.xmp");
	    cout << "Press any key\n";
	    char c; cin >> c;
	 #endif

         // Compute projection power
         proj().compute_stats(avg,proj_power(n),dummy,dummy);


         if (n%30==0) progress_bar(n);
      }
      progress_bar(prm.N_stats);
      histogram1D hist_vol, hist_proj, hist_area;
      compute_hist(volume,hist_vol,300);
      compute_hist(proj_power,hist_proj,300);
      compute_hist(proj_area,hist_area,300);
      volume.compute_stats(avg,stddev,dummy,dummy);
      cout << "# Volume average: " << avg << endl
           << "# Volume stddev:  " << stddev << endl
           << "# Volume percentil  2.5%: " << hist_vol.percentil(2.5) << endl
           << "# Volume percentil 97.5%: " << hist_vol.percentil(97.5) << endl
	   << hist_vol << endl << endl;
      proj_power.compute_stats(power_avg,power_stddev,dummy,dummy);
      cout << "# Projection power average: " << power_avg << endl
           << "# Projection power stddev:  " << power_stddev << endl
           << "# Projection percentil  2.5%: " << hist_proj.percentil(2.5) << endl
           << "# Projection percentil 97.5%: " << hist_proj.percentil(97.5) << endl
	   << hist_proj;
      proj_area.compute_stats(area_avg,area_stddev,dummy,dummy);
      cout << "# Projection area average: " << area_avg << endl
           << "# Projection area stddev:  " << area_stddev << endl
           << "# Area percentil  2.5%:    " << hist_area.percentil(2.5) << endl
           << "# Area percentil 97.5%:    " << hist_area.percentil(97.5) << endl
	   << hist_area << endl;
      if (prm.target_SNR!=-1) {
         cout << endl;
         cout << "For an SNR of " << prm.target_SNR
              << ", a total noise with a standard deviation of "
              << sqrt(power_avg*power_avg*Xdim*Ydim/(prm.target_SNR*area_avg))
              << " is needed\n";
      }
   }
}
#undef DEBUG
