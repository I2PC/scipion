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

#include "../Prog_project.hh"
#include <XmippData/xmippArgs.hh>

/* Read from command line ================================================== */
void Prog_Project_Parameters::read(int argc, char **argv) {
   fn_proj_param = get_param(argc, argv, "-i");
   fn_sel_file   = get_param(argc, argv, "-o","");
   fn_crystal    = get_param(argc, argv, "-crystal","");
   only_create_angles=check_param(argc,argv,"-only_create_angles");
   if (check_param(argc,argv,"-show_angles"))
      tell |= TELL_SHOW_ANGLES;
}

/* Usage =================================================================== */
void Prog_Project_Parameters::usage() {
   printf("\nUsage:\n\n");
   printf("project -i <Parameters File> \n"
          "       [-o <sel_file>]\n"
          "       [-show_angles]\n"
          "       [-only_create_angles]\n"
          "       [-crystal <crystal_parameters_file>]\n");
   printf(
      "\tWhere:\n"
      "\t<Parameters File>:  File containing projection parameters\n"
      "\t                    check the manual for a description of the parameters\n"
      "\t<sel_file>:         This is a selection file with all the generated\n"
      "\t                    projections\n");
}

/* Projection parameters from program parameters =========================== */
void Projection_Parameters::from_prog_params(
   const Prog_Project_Parameters &prog_prm) {
   read(prog_prm.fn_proj_param);
   tell=prog_prm.tell;
}

/* Read Projection Parameters ============================================== */
int translate_randomness(char * str) {
   if (str==NULL)                     return ANGLE_RANGE_DETERMINISTIC;
   if (strcmp(str,"random_group")==0) return ANGLE_RANGE_RANDOM_GROUPS;
   if (strcmp(str,"random")==0)       return ANGLE_RANGE_RANDOM;
   if (strcmp(str,"even")==0)         return ANGLE_EVENLY;
   REPORT_ERROR(3007,
      (string)"Prog_Project_Parameters::read: Not recognized randomness: "
         +str);
}

void Projection_Parameters::read(FileName fn_proj_param) _THROW {
   FILE    *fh_param;
   char    line[201];
   int     lineNo=0;
   char    *auxstr;

   if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,
         (string)"Prog_Project_Parameters::read: There is a problem "
         "opening the file "+fn_proj_param);

   while (fgets (line, 200,fh_param ) != NULL) {
      if (line[0]==0)    continue;
      if (line[0]=='#')  continue;
      if (line[0]=='\n') continue;
      switch(lineNo) {
         case 0:
            fn_phantom=first_word(line,3007,
               "Prog_Project_Parameters::read: Phantom name not found");
            lineNo=1;
            break;
         case 1:
            fn_projection_seed=
               first_word(line,3007,
                  "Prog_Project_Parameters::read: Error in Projection seed");
            // Next two parameters are optional
            auxstr=next_token();
            if (auxstr!=NULL) starting=
                  AtoI(auxstr,3007,
                     "Prog_Project_Parameters::read: Error in First "
                     "projection number");
            fn_projection_extension=next_token();
            lineNo=2;
            break;
         case 2:
            proj_Xdim=AtoI(first_token(line),3007,
               "Prog_Project_Parameters::read: Error in X dimension");
            proj_Ydim=AtoI(next_token(),3007,
               "Prog_Project_Parameters::read: Error in Y dimension");
            lineNo=3;
            break;
         case 3:
            // Angle file
            fn_angle=first_word(line);
            if (fn_angle=="NULL") ;
            else if (exists(fn_angle)) {
               // Angle source file
               try {
                  ang1=next_word(); check_angle_descr(ang1);
                  ang2=next_word(); check_angle_descr(ang2);
                  ang3=next_word(); check_angle_descr(ang3);
               } catch (Xmipp_error XE) {
                  REPORT_ERROR(3007,
                     (string)"An angle order description is "
                     +"missing for file: "+fn_angle);
               }
            }  else
               REPORT_ERROR(3007,(string)"Prog_Project_Parameters::read: "
                  "file " + fn_angle + " doesn't exist");
            lineNo=4;
            break;
         case 4:
            // theta init
            auxstr=first_word(line);
            if (strcmp(auxstr,"NULL")!=0) {
               enable_angle_range=1;
               rot_range.ang0=AtoF(auxstr,3007,
                  "Prog_Project_Parameters::read: Error in Rotational Init");
               auxstr=next_token();
               if (auxstr==NULL) {
                  // Fixed mode
                  rot_range.randomness=ANGLE_RANGE_DETERMINISTIC;
                  rot_range.angF=rot_range.ang0;
                  rot_range.samples=1;
               } else {
                  rot_range.angF=AtoF(auxstr,3007,
                     "Prog_Project_Parameters::read: Error in Rotational Final");
                  rot_range.samples=AtoI(next_token(),3007,
                     "Prog_Project_Parameters::read: Error in Rotational "
                     "Samples");
                  rot_range.randomness=translate_randomness(next_token());
               }
               lineNo=5;
            } else {
               enable_angle_range=0;
               lineNo=7;
            }
            break;
         case 5:
            tilt_range.ang0=AtoF(first_token(line),3007,
               "Prog_Project_Parameters::read: Error in Tilting Init");
            auxstr=next_token();
            if (auxstr==NULL) {
               // Fixed mode
               tilt_range.randomness=ANGLE_RANGE_DETERMINISTIC;
               tilt_range.angF=tilt_range.ang0;
               tilt_range.samples=1;
            } else {
               tilt_range.angF=AtoF(auxstr,3007,
                  "Prog_Project_Parameters::read: Error in Tilting Final");
               tilt_range.samples=AtoI(next_token(),3007,
                  "Prog_Project_Parameters::read: Error in Tilting Samples");
               tilt_range.randomness=translate_randomness(next_token());
            }
            lineNo=6;
            break;
         case 6:
            psi_range.ang0=AtoF(first_token(line),3007,
               "Prog_Project_Parameters::read: Error in Psi Init");
            auxstr=next_token();
            if (auxstr==NULL) {
               // Fixed mode
               psi_range.randomness=ANGLE_RANGE_DETERMINISTIC;
               psi_range.angF=psi_range.ang0;
               psi_range.samples=1;
            } else {
               psi_range.angF=AtoF(auxstr,3007,
                  "Prog_Project_Parameters::read: Error in Psi Final");
               psi_range.samples=AtoI(next_token(),3007,
                  "Prog_Project_Parameters::read: Error in Psi Samples");
               psi_range.randomness=translate_randomness(next_token());
            }
            lineNo=7;
            break;
         case 7:
            rot_range.Ndev=AtoF(first_word(line),3007,
               "Prog_Project_Parameters::read: Error in Rotational noise");
            auxstr=next_token();
	    if (auxstr!=NULL)
               rot_range.Navg=AtoF(auxstr,3007,
                  "Prog_Project_Parameters::read: Error in Rotational bias");
	    else rot_range.Navg=0;
            lineNo=8;
            break;
         case 8:
            tilt_range.Ndev=AtoF(first_word(line),3007,
               "Prog_Project_Parameters::read: Error in tilting noise");
            auxstr=next_token();
	    if (auxstr!=NULL)
               tilt_range.Navg=AtoF(auxstr,3007,
                  "Prog_Project_Parameters::read: Error in tilting bias");
	    else tilt_range.Navg=0;
            lineNo=9;
            break;
         case 9:
           psi_range.Ndev=AtoF(first_word(line),3007,
              "Prog_Project_Parameters::read: Error in psi noise");
           auxstr=next_token();
	   if (auxstr!=NULL)
              psi_range.Navg=AtoF(auxstr,3007,
                 "Prog_Project_Parameters::read: Error in psi bias");
	   else psi_range.Navg=0;
           lineNo=10;
           break;
         case 10:
           Npixel_dev=AtoF(first_word(line),3007,
              "Prog_Project_Parameters::read: Error in pixel noise");
           auxstr=next_token();
	   if (auxstr!=NULL)
              Npixel_avg=AtoF(auxstr,3007,
                 "Prog_Project_Parameters::read: Error in pixel bias");
	   else Npixel_avg=0;
           lineNo=11;
           break;
         case 11:
           Ncenter_dev=AtoF(first_word(line),3007,
              "Prog_Project_Parameters::read: Error in center noise");
           auxstr=next_token();
	   if (auxstr!=NULL)
              Ncenter_avg=AtoF(auxstr,3007,
                 "Prog_Project_Parameters::read: Error in center bias");
	   else Ncenter_avg=0;
           lineNo=12;
           break;
      } /* switch end */  
   } /* while end */
   if (lineNo!=12)
      REPORT_ERROR(3007,(string)"Prog_Project_Parameters::read: I "
         "couldn't read all parameters from file " + fn_proj_param);

   fclose(fh_param);
}

/* Write =================================================================== */
void Projection_Parameters::write(FileName fn_proj_param) _THROW {
   FILE *fh_param;

   if ((fh_param = fopen(fn_proj_param.c_str(), "w")) == NULL)
      REPORT_ERROR(3005,
         (string)"Prog_Project_Parameters::write: There is a problem "
         "opening the file "+fn_proj_param+" for output");

   fprintf(fh_param,
      "# Volume and projection files -----------------------------------\n");
   fprintf(fh_param,
      "# volume description file or volume file\n");
   fprintf(fh_param,
      "%s\n",fn_phantom.c_str());
   fprintf(fh_param,
      "# projection seed, first projection number (by default, 1) and extension\n");
   fprintf(fh_param,
      "%s %d %s\n", fn_projection_seed.c_str(), starting,
         fn_projection_extension.c_str());
   fprintf(fh_param,
      "# Y and X projection dimensions\n");
   fprintf(fh_param,
      "%d %d\n",proj_Ydim, proj_Xdim);
   fprintf(fh_param,
      "#\n");
   
   fprintf(fh_param,
      "# Angle Definitions ---------------------------------------------\n");
   if (fn_angle!="") {
      fprintf(fh_param,
         "%s %s %s %s\n",fn_angle.c_str(), ang1.c_str(), ang2.c_str(),
         ang3.c_str());
   } else {
      fprintf(fh_param,"%d ",rot_range.ang0);
      if (rot_range.angF!=rot_range.ang0) {
         fprintf(fh_param,"%d %d ",rot_range.angF, rot_range.samples);
         switch (rot_range.randomness) {
            case (ANGLE_RANGE_RANDOM_GROUPS):
               fprintf(fh_param,"random_group\n"); break;
            case (ANGLE_RANGE_RANDOM):
               fprintf(fh_param,"random\n"); break;
	    case (ANGLE_EVENLY):
	       fprintf(fh_param,"even\n"); break;
            default: fprintf(fh_param,"\n"); break;
         }
      }

      fprintf(fh_param,"%d ",tilt_range.ang0);
      if (tilt_range.angF!=tilt_range.ang0) {
         fprintf(fh_param,"%d %d ",tilt_range.angF, tilt_range.samples);
         switch (tilt_range.randomness) {
            case (ANGLE_RANGE_RANDOM_GROUPS):
               fprintf(fh_param,"random_group\n"); break;
            case (ANGLE_RANGE_RANDOM):
               fprintf(fh_param,"random\n"); break;
	    case (ANGLE_EVENLY):
	       fprintf(fh_param,"even\n"); break;
            default: fprintf(fh_param,"\n"); break;
         }
      }

      fprintf(fh_param,"%d ",psi_range.ang0);
      if (psi_range.angF!=psi_range.ang0) {
         fprintf(fh_param,"%d %d ",psi_range.angF, psi_range.samples);
         switch (psi_range.randomness) {
            case (ANGLE_RANGE_RANDOM_GROUPS):
               fprintf(fh_param,"random_group\n"); break;
            case (ANGLE_RANGE_RANDOM):
               fprintf(fh_param,"random\n"); break;
	    case (ANGLE_EVENLY):
	       fprintf(fh_param,"even\n"); break;
            default: fprintf(fh_param,"\n"); break;
         }
      }
   }
   fprintf(fh_param,
      "# Noise description ----------------------------------------------\n");
   fprintf(fh_param,
      "#     noise (and bias) applied to rotational angle\n");
   fprintf(fh_param,"%f ",rot_range.Ndev);
   if (rot_range.Navg!=0) fprintf(fh_param,"%f \n",rot_range.Navg);
   else fprintf(fh_param,"\n");

   fprintf(fh_param,
      "#     noise (and bias) applied to tilting angle\n");
   fprintf(fh_param,"%f ",tilt_range.Ndev);
   if (tilt_range.Navg!=0) fprintf(fh_param,"%f \n",tilt_range.Navg);
   else fprintf(fh_param,"\n");

   fprintf(fh_param,
      "#     noise (and bias) applied to psi angle\n");
   fprintf(fh_param,"%f ",psi_range.Ndev);
   if (psi_range.Navg!=0) fprintf(fh_param,"%f \n",psi_range.Navg);
   else fprintf(fh_param,"\n");

   fprintf(fh_param,
      "#     Noise (and bias) applied to pixels\n");
   fprintf(fh_param,"%f ",Npixel_dev);
   if (Npixel_avg!=0) fprintf(fh_param,"%f \n",Npixel_avg);
   else fprintf(fh_param,"\n");

   fprintf(fh_param,
      "#     Noise (and bias) applied to particle center coordenates\n");
   fprintf(fh_param,"%f ",Ncenter_dev);
   if (Ncenter_avg!=0) fprintf(fh_param,"%f \n",Ncenter_avg);
   else fprintf(fh_param,"\n");
   
   fclose(fh_param);
}

/* Generate angles ========================================================= */
// This function generates the angles for a given angle ("rot", "tilt"
// or "psi") according to the projection parameters. The output document
// file is supposed to be large enough to hold all angles
// Some aliases
#define Nrot  prm.rot_range.samples
#define Ntilt prm.tilt_range.samples
#define Npsi  prm.psi_range.samples
#define proj_number(base,irot,itilt,ipsi) base+irot*Ntilt*Npsi+itilt*Npsi+ipsi
void generate_angles(int ExtProjs, const Angle_range &range,
   DocFile &DF, char ang_name, const Projection_Parameters &prm) {
   double ang;
   int   N1,N2;
   int   i,j,k;
   int   iproj, idx;
   int   limit;

   // Select loop limit ....................................................
   switch (range.randomness) {
      case ANGLE_RANGE_DETERMINISTIC: limit=range.samples; break;
      case ANGLE_RANGE_RANDOM_GROUPS: limit=range.samples; break;
      case ANGLE_RANGE_RANDOM       : limit=Nrot*Ntilt*Npsi; break;
   }
   
   // Which column to write in the document file ...........................
   switch(ang_name) {
      case 'r': idx=0; break;
      case 't': idx=1; break;
      case 'p': idx=2; break;
   }

   double unif_min=cos(DEG2RAD(range.angF));
   double unif_max=cos(DEG2RAD(range.ang0));
   for (i=0; i<limit; i++) {
       // Select angle .....................................................
       if (range.randomness==ANGLE_RANGE_DETERMINISTIC) {
          if (range.samples>1)
             ang=range.ang0 +
                (range.angF-range.ang0)/(double)(range.samples-1)*i;
          else ang=range.ang0;
       }  else {
	    switch(ang_name) {
	       case 'r':
	       case 'p': ang=rnd_unif(range.ang0,range.angF);            break;
	       case 't': ang=RAD2DEG(acos(rnd_unif(unif_min,unif_max))); break;
	    }
       }
          

       // Copy this angle to those projections belonging to this group .....
       // If there is any group
       if (range.randomness!=ANGLE_RANGE_RANDOM) {
          switch(ang_name) {
             case 'r': N1=Ntilt; N2=Npsi;  break;
             case 't': N1=Nrot;  N2=Npsi;  break;
             case 'p': N1=Nrot;  N2=Ntilt; break;
          }
          for (j=0; j<N1; j++)
              for (k=0; k<N2; k++) {
                  switch(ang_name) {
                     case 'r': iproj=proj_number(ExtProjs,i,j,k)+DF.FirstKey();
                        break;
                     case 't': iproj=proj_number(ExtProjs,j,i,k)+DF.FirstKey();
                        break;
                     case 'p': iproj=proj_number(ExtProjs,j,k,i)+DF.FirstKey();
                        break;
                  }
                  DF.set(iproj,idx,ang);
              }
       } else
          DF.set(ExtProjs+i+DF.FirstKey(),idx,ang);
   }
}

/* Generate evenly distributed angles ====================================== */
void generate_even_angles(int ExtProjs, int Nrottilt, DocFile &DF,
   const Projection_Parameters &prm) {
   // We will run over the tilt angle in a deterministic way
   // then for every tilt angle, a rot_step is computed so that
   // it keeps the same distance in the circle generated by tilt
   // as the sample distance at the equator (tilt=90).
   int N=0;
   int limit=prm.tilt_range.samples;
   double rot_step_at_equator=(prm.rot_range.angF-prm.rot_range.ang0)/
      (double)(Nrot-1);
   for (int i=0; i<limit; i++) {
      // Compute the corresponding deterministic tilt
      double tilt=prm.tilt_range.ang0 +
                (prm.tilt_range.angF-prm.tilt_range.ang0)/
		(double)(Ntilt-1)*i;
      // Now compute the corresponding rotational angles
      double rot_step;
      if (tilt!=0 && tilt!=180) rot_step=rot_step_at_equator/sin(DEG2RAD(tilt));
      else rot_step=prm.rot_range.angF-prm.rot_range.ang0+1;
      for (double rot=prm.rot_range.ang0; rot<=prm.rot_range.angF; rot+=rot_step) {
          // Copy this angle to those projections belonging to this group .....
          // If there is any group
	  for (int k=0; k<Npsi; k++) {
          // Select psi
	     double psi;
	     if (prm.psi_range.randomness==ANGLE_RANGE_DETERMINISTIC) {
        	if (prm.psi_range.samples>1)
        	   psi=prm.psi_range.ang0 +
                      (prm.psi_range.angF-prm.psi_range.ang0)/
		      (double)(prm.psi_range.samples-1)*k;
        	else psi=prm.psi_range.ang0;
	     }  else
        	psi=rnd_unif(prm.psi_range.ang0,prm.psi_range.angF);

	     int iproj=ExtProjs+N+Nrottilt*k+DF.FirstKey();
	     DF.set(iproj,0,rot);
	     DF.set(iproj,1,tilt);
	     DF.set(iproj,2,psi);
	  }
	  N++;
      }
   }
}

// See generate_even_angles for comments
int count_even_angles(const Projection_Parameters &prm) {
   int N=0;
   int limit=prm.tilt_range.samples;
   double rot_step_at_equator=(prm.rot_range.angF-prm.rot_range.ang0)/
      (double)(Nrot-1);
   for (int i=0; i<limit; i++) {
      double tilt=prm.tilt_range.ang0 +
                (prm.tilt_range.angF-prm.tilt_range.ang0)/
		(double)(Ntilt-1)*i;
      double rot_step;
      if (tilt!=0 && tilt!=180) rot_step=rot_step_at_equator/sin(DEG2RAD(tilt));
      else rot_step=prm.rot_range.angF-prm.rot_range.ang0+1;
      for (double rot=prm.rot_range.ang0; rot<=prm.rot_range.angF; rot+=rot_step)
         N++;
   }
   N++; // This shouldn't be necessary but some GCC optimization
        // sometimes doesn't do well its work. For instance if we
        // add cout << N after N++ in the loop, then it works perfectly
   return N;
}

/* Assign angles =========================================================== */
int Assign_angles(DocFile &DF, const Projection_Parameters &prm) {
   int ExtProjs=0, IntProjs=0;        // External and internal projections
   int Nrottilt;                      // Number of evenly distributed
      	             	      	      // projections

   DF.clear();
   DF.FirstKey()=prm.starting;

// External generation mode
   if (prm.fn_angle!="NULL")
      ExtProjs=read_Euler_document_file(prm.fn_angle,
         prm.ang1, prm.ang2, prm.ang3, DF);

// Internal generation mode
   if (prm.enable_angle_range) {
      randomize_random_generator();
      if (prm.rot_range.randomness!=ANGLE_EVENLY)
         IntProjs=Nrot*Ntilt*Npsi;
      else {
         Nrottilt=count_even_angles(prm);
         IntProjs=Nrottilt*Npsi;
      }
      DF.append_data_line(IntProjs);  // Create lines for the angles
      if (prm.rot_range.randomness!=ANGLE_EVENLY) {
	 generate_angles(ExtProjs,prm.rot_range,  DF, 'r', prm);
	 generate_angles(ExtProjs,prm.tilt_range, DF, 't', prm);
         generate_angles(ExtProjs,prm.psi_range,  DF, 'p', prm);
      } else {
         generate_even_angles(ExtProjs,Nrottilt,DF,prm);
         // Check if the last entry is empty
         DF.locate(DF.LineNo());
         if (DF.get_current_line().get_no_components()==0) 
            DF.remove_current();
      }
   }

// Exit
   DF.go_first_data_line();
   return ExtProjs+IntProjs;
}

/* Produce Side Information ================================================ */
void PROJECT_Side_Info::produce_Side_Info(const Projection_Parameters &prm)
   _THROW {
// Generate Projection angles
   Assign_angles(DF,prm);

// Load Phantom and set working mode
   if (Is_VolumeXmipp(prm.fn_phantom)) {
      phantom_vol.read(prm.fn_phantom);
      phantom_vol().set_Xmipp_origin();
      voxel_mode=1;
   } else {
      phantom_descr.read(prm.fn_phantom);
      voxel_mode=0;
   }
}

/* Effectively project ===================================================== */
int PROJECT_Effectively_project(const Projection_Parameters &prm,
    PROJECT_Side_Info &side, const Crystal_Projection_Parameters &prm_crystal,
    Projection &proj, SelFile &SF) {

   int NumProjs=0;
   SF.clear();
   cerr << "Projecting ...\n";
   if (!(prm.tell&TELL_SHOW_ANGLES)) init_progress_bar(side.DF.dataLineNo());
   SF.reserve(side.DF.dataLineNo());
   
   DocFile DF_movements;
   DF_movements.append_comment("True rot, tilt and psi; rot, tilt, psi, X and Y shifts applied");
   matrix1D<double> movements(8);
   while (!side.DF.eof()) {
      double rot, tilt, psi;         // Actual projecting angles
      FileName fn_proj;              // Projection name
      fn_proj.compose(prm.fn_projection_seed,side.DF.get_current_key(),
         prm.fn_projection_extension);

      // Choose angles .....................................................
      if (prm.tell&TELL_SHOW_ANGLES) side.DF.show_line(cout);
      else if ((NumProjs%MAX(1,side.DF.dataLineNo()/60))==0)
              progress_bar(NumProjs);
      movements(0) = rot  = side.DF(0);
      movements(1) = tilt = side.DF(1);
      movements(2) = psi  = side.DF(2);
      
      // Choose Center displacement ........................................
      double shiftX=rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
      double shiftY=rnd_gaus(prm.Ncenter_avg, prm.Ncenter_dev);
      movements(6) = shiftX;
      movements(7) = shiftY;

      // Really project ....................................................
      if (side.voxel_mode) {
         project_Volume(side.phantom_vol(),proj,prm.proj_Ydim,prm.proj_Xdim,
            rot,tilt,psi);
         IMGMATRIX(proj).self_translate(vector_R2(shiftX,shiftY));
      } else {
         Phantom aux;
         aux=side.phantom_descr;
         aux.shift(shiftX,shiftY,0);
         if (prm_crystal.crystal_Xdim==0)
            // Project a single mathematical volume
            aux.project_to(proj,
               prm.proj_Ydim, prm.proj_Xdim,
               rot, tilt, psi);
         else
            // Project mathematical volume as a crystal
            project_crystal(aux,proj,prm,side,prm_crystal,rot,tilt,psi);
      }
            
      // Add noise in angles and voxels ....................................
      rot  += rnd_gaus(prm.rot_range.Navg,  prm.rot_range.Ndev);
      tilt += rnd_gaus(prm.tilt_range.Navg, prm.tilt_range.Ndev);
      psi  += rnd_gaus(prm.psi_range.Navg,  prm.psi_range.Ndev);
      movements(3) = rot-movements(0);
      movements(4) = tilt-movements(1);
      movements(5) = psi-movements(2);
      proj.set_eulerAngles(rot,tilt,psi);
      IMGMATRIX(proj).add_noise(prm.Npixel_avg, prm.Npixel_dev,"gaussian");
      
      // Save ..............................................................
      proj.write(fn_proj);
      NumProjs++;
      SF.insert(fn_proj,SelLine::ACTIVE);
      DF_movements.append_data_line(movements);
      
      side.DF.next_data_line();
   }
   if (!(prm.tell&TELL_SHOW_ANGLES)) progress_bar(side.DF.dataLineNo());
   
   DF_movements.write(prm.fn_projection_seed+"_movements.txt");
   return NumProjs;
}
    
/* ROUT_project ============================================================ */
int ROUT_project(Prog_Project_Parameters &prm, Projection &proj, SelFile &SF) {
   randomize_random_generator();
// Read projection parameters and produce side information
   Projection_Parameters proj_prm;
   PROJECT_Side_Info side;
   proj_prm.from_prog_params(prm);
   side.produce_Side_Info(proj_prm);

   Crystal_Projection_Parameters crystal_proj_prm;
   if (prm.fn_crystal!="") crystal_proj_prm.read(prm.fn_crystal);

   int ProjNo=0;
   if (!prm.only_create_angles) {
      // Really project
      ProjNo=PROJECT_Effectively_project(proj_prm, side, crystal_proj_prm,
         proj, SF);

      // Save SelFile
      if (prm.fn_sel_file!="") SF.write(prm.fn_sel_file);
   } else {
      cout << side.DF;
   }
   return ProjNo;
}
