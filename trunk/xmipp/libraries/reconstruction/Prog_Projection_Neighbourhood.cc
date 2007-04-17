/***************************************************************************
 *
 * Authors:    Sjors Scheres                     scheres@cnb.uam.es
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

#include "../Prog_Projection_Neighbourhood.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippHeader.hh>

// Read arguments ==========================================================
void Prog_projection_neighbourhood_prm::read(int argc, char **argv) {
   int i,maxcol;
   fn_sel=get_param(argc,argv,"-i");
   fn_ref=get_param(argc,argv,"-nbh");
   fn_root_out=get_param(argc,argv,"-oroot","nbhood");
   maxdist=AtoF(get_param(argc,argv,"-dist","10"));
   fn_sym=get_param(argc,argv,"-sym","");
   if (fn_sym!="") SL.read_sym_file(fn_sym);
   DF2.read(fn_ref);
   SF1.read(fn_sel);
   if (DF2.FirstLine_ColNo()<2) {
     REPORT_ERROR(1,"Projection Neighbourhood: Neighbourhoods docfile has less than 2 columns.");
   }
}

// Extract angles ==========================================================
void Prog_projection_neighbourhood_prm::get_angles(SelFile &SF_in, DocFile &DF_out) {
   DF_out.clear();
   int i=0;
   double phi,theta,psi;
   time_config();
   cerr << "Extracting angles ...\n";
   init_progress_bar(SF_in.ImgNo());
   while (!SF_in.eof()) {
     headerXmipp H;
     H.read(SF_in.NextImg());
     H.get_eulerAngles(phi,theta,psi);  
     DF_out.append_angles(phi,theta,psi,"rot","tilt","psi");
     i++;
     if (i%10==0) progress_bar(i);
   }
   progress_bar(SF_in.ImgNo());
}

// Show ====================================================================
void Prog_projection_neighbourhood_prm::show() {
   cerr << "Selfile                      : " << fn_sel        << endl
        << "Neighbourhoods docfile       : " << fn_ref        << endl	
	<< "Output root                  : " << fn_root_out   << endl
	<< "Max. neighbour distance      : " << maxdist       << endl
	<< "Symmetry file                : " << fn_sym        << endl
   ;
}

// usage ===================================================================
void Prog_projection_neighbourhood_prm::usage() {
   cerr << "   -i     <SelFile>       : Selfile containing the images \n"
        << "   -nbh   <DocFile>       : Document file with the defined neighbourhood directions\n"
	<< "  [-oroot <name=nbhood> ] : Rootname for output files \n"
	<< "  [-dist  <d=10>        ] : Maximum neighbourhood distance \n"
	<< "  [-sym <symmetry file> ] : Symmetry file if any\n"
   ;
}

// Check symmetries --------------------------------------------------------
//#define DEBUG
double Prog_projection_neighbourhood_prm::check_symmetries(double rot1, double tilt1, double &rot2, double &tilt2) {

   int imax=SL.SymsNo()+1;
   matrix2D<double>  L(4,4), R(4,4);    // A matrix from the list
   double best_ang_dist=9999;
   double best_rot2,best_tilt2,diff_rot,diff_tilt;
   double ang_dist;
   double psi1=0.;
   double psi2=0.;
   
   for (int i=0; i<imax; i++) {
      double rot2p, tilt2p, psi2p;
      if (i==0) {rot2p=rot2; tilt2p=tilt2; psi2p=psi2;}
      else {
         SL.get_matrices(i-1,L,R);
         L.resize(3,3); // Erase last row and column
         R.resize(3,3); // as only the relative orientation
                        // is useful and not the translation
         Euler_apply_transf(L,R,rot2,tilt2,psi2,rot2p,tilt2p,psi2p);
      }
      diff_rot=rot1-rot2p;
      diff_tilt=tilt1-tilt2p;
      // Some (like virus) symmetries can give equivalent images which are rotated in psi!
      // Neglect these for now!
      diff_rot=ABS(realWRAP(diff_rot,-180,180));
      diff_tilt=ABS(realWRAP(diff_tilt,-180,180));
      ang_dist=sqrt((diff_rot*diff_rot)+(diff_tilt*diff_tilt));
      if (ang_dist<best_ang_dist) {
	best_rot2=rot2p; best_tilt2=tilt2p; 
	best_ang_dist=ang_dist;
      }
      Euler_another_set(rot2p,tilt2p,psi2p,rot2p,tilt2p,psi2p);
      diff_rot=rot1-rot2p;
      diff_tilt=tilt1-tilt2p;
      diff_rot=ABS(realWRAP(diff_rot,-180,180));
      diff_tilt=ABS(realWRAP(diff_tilt,-180,180));
      ang_dist=sqrt((diff_rot*diff_rot)+(diff_tilt*diff_tilt));
      if (ang_dist<best_ang_dist) {
	best_rot2=rot2p; best_tilt2=tilt2p; 
	best_ang_dist=ang_dist; 
      }      
   }
   rot2=best_rot2;
   tilt2=best_tilt2;
   return best_ang_dist;
}

#define DEBUG
// Compute Projection Neighbourhood -----------------------------------------
void Prog_projection_neighbourhood_prm::compute_neighbourhood() {
  double dist=0.;
  double rot1, tilt1;
  double rot2, tilt2;
  double distp;
  int i=0;
  SelFile SF_out;
  SelLine selline;

  get_angles(SF1,DF1);
  cerr << "Calculating ...\n";
  DF2.go_first_data_line();
  while (!DF2.eof()) {
    // Read reference projection direction
    rot1=realWRAP(DF2(0),-180,180);
    tilt1=realWRAP(DF2(1),-180,180);
    SF_out.reserve(DF1.LineNo());
    SF_out.go_beginning();

    DF1.go_first_data_line();
    while (!DF1.eof()) {
      // Read assigned angles from document file
      rot2=realWRAP(DF1(0),-180,180);
      tilt2=realWRAP(DF1(1),-180,180);
      distp=check_symmetries(rot1,tilt1,rot2,tilt2);
      // Fill the output result
      if (distp<=maxdist) {
	SF1.go_beginning();
	SF1.jump(DF1.get_current_key()-1);
	selline=SF1.current();
	SF_out.insert(selline);
      }      
      // Move to next data line
      DF1.next_data_line();
    } 
    // finished reading all particles for this neighbourhood
    i++;
    FileName fn_sel_out;
    fn_sel_out.compose(fn_root_out,i,"sel");
    SF_out.write(fn_sel_out);
    SF_out.clear();
    DF2.next_data_line();
  }
  return;
}

