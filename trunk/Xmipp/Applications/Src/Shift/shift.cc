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
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippGeometry.hh>

class Shift_Scale_parameters: public Prog_parameters {
public:
   matrix1D<double> shift;
   matrix1D<double> scale;
   bool             wrap;
   DocFile          DF_shifts;
   DocFile          DF_scales;
   int              colX_shift;
   int              colX_scale;
   int              Docfile;
   

   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      int i_shift=position_param(argc,argv,"-shift");
      int i_scale=position_param(argc,argv,"-scale");
      
      
      if (i_shift==-1 && i_scale==-1)
         REPORT_ERROR(1,"Shift_Scale:: Cannot find -shift or -scale");
      else if (ABS(i_shift-i_scale)<=1)
         REPORT_ERROR(1,"Shift_cale: Not enough parameters after -shift or -scale");
      Docfile=check_param(argc, argv, "-colX_shift")||
              check_param(argc, argv, "-colX_scale");
      if (Docfile){
         if(i_shift>0){
	    DF_shifts.read(argv[i_shift+1]);
  	    colX_shift=AtoI(get_param(argc,argv,"-colX_shift"));
	    colX_shift -=3;
	    if(colX_shift<0) {cout << "colX_shift must be no less than 3" 
	                            << endl; exit(1);
			      }
	    }
         else colX_shift= -1;
         if(i_scale>0){
            DF_scales.read(argv[i_scale+1]);
  	    colX_scale=AtoI(get_param(argc,argv,"-colX_scale","5"));
	    colX_scale -=3;
	    if(colX_scale<0) {cout << "colX_scale must be no less than 3" 
	                            << endl; exit(1);
			      }
	    }
         else colX_scale= -1;
      }
      else{
         shift=get_vector_param(argc,argv,"-shift",-1);
         scale=get_vector_param(argc,argv,"-scale",-1);
      }	 
      wrap=!check_param(argc,argv,"-dont_wrap");
   }
   
   void show() {
      Prog_parameters::show();
      if (wrap) cout << "Wrapping image/volume\n";
      else      cout << "Not wrapping image/volume\n";
      if (shift.get_dim() > 1 )
         cout << "Shift: " << shift.transpose() << endl;
      else if (DF_shifts.name()!=""){
         cout << "Shift docfile: " << DF_shifts.name() << endl;
	 cout << "colX_shift:  " << colX_shift << endl;
      }
      if (scale.get_dim() > 1 )
         cout << "Scale: " << scale.transpose() << endl;
      else if (DF_scales.name()!=""){
         cout << "Scale: docfile: "       << DF_scales.name() << endl;
	 cout << "colX_scale:  " << colX_scale << endl;
      }
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "   -shift [<x>,<y>[,<z>]]   : Shift by (x,y,z)\n"
           << "   -scale [<x>,[<y>,<z>]]   : Scale by (x,y,z)\n"
           << "   -Docfile                 : Shift and/or Scales are stored in a Docfile\n"
           << "   -shift <DocFile>         : Shifts are stored in a Docfile\n"
           << "   -scale <DocFile>         : Scales are stored in a Docfile (may be the same\n"
	   << "                              Docfile used for shifts\n"
	   << "  [-colXshift <col>]           : Column with  the X shift\n"
	   << "                              First column with data is number 1.\n"
	   << "  [-colXscale <col>]       : Column with the scale information\n"
           << "  [-dont_wrap]              : By default, the image is wrapped\n";
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Shift_Scale_parameters *eprm=(Shift_Scale_parameters *) prm;
   matrix2D<double> A(3,3); A.init_identity(); 

//   if(eprm->shift.get_dim()>1)
//        ;//shift is already filled
//   else 
   if (eprm->DF_shifts.name()!=""){
      eprm->shift.resize(2);
      XX(eprm->shift)=eprm->DF_shifts(eprm->colX_shift);
      YY(eprm->shift)=eprm->DF_shifts(eprm->colX_shift+1);
      eprm->DF_shifts.next_data_line();
   }
   else if (eprm->Docfile!=FALSE){
      eprm->shift.resize(2); 
      eprm->shift.init_constant(1.);
   }
   A(0,2)=XX(eprm->shift);
   A(1,2)=YY(eprm->shift);
//   if(eprm->scale.get_dim()>1)
//      ;//scale already filled
//   else 
   if (eprm->DF_scales.name()!=""){
      eprm->scale.resize(2);
      XX(eprm->scale)=eprm->DF_scales(eprm->colX_scale);
      YY(eprm->scale)=eprm->DF_scales(eprm->colX_scale+1);
      eprm->DF_scales.next_data_line();
   }
   else if (eprm->Docfile!=FALSE){
      eprm->scale.resize(2); 
      eprm->scale.init_constant(1.);
   }
   A(0,0)=XX(eprm->scale);
   A(1,1)=YY(eprm->scale);
   img().self_apply_geom(A,IS_NOT_INV,eprm->wrap);
   return TRUE;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   matrix2D<double> A(4,4); A.init_identity; 
   Shift_Scale_parameters *eprm=(Shift_Scale_parameters *) prm;
//   if (eprm->DF_shifts.name()=="")
//        ;//shift is already filled
//   else 
  if (eprm->DF_shifts.name()!=""){
      eprm->shift.resize(3);
      XX(eprm->shift)=eprm->DF_shifts(eprm->colX_shift);
      YY(eprm->shift)=eprm->DF_shifts(eprm->colX_shift+1);
      ZZ(eprm->shift)=eprm->DF_shifts(eprm->colX_shift+2);
      eprm->DF_shifts.next_data_line();
   }
   else if (eprm->Docfile!=FALSE){
      eprm->shift.resize(3); 
      eprm->shift.init_constant(1.);
   }
   A(0,3)=XX(eprm->shift);
   A(1,3)=YY(eprm->shift);
   A(2,3)=ZZ(eprm->shift);
//   if(eprm->scale.get_dim()>1)
//      ;//scale already filled
//   else 
   if (eprm->DF_scales.name()!=""){
      eprm->scale.resize(3);
      XX(eprm->scale)=eprm->DF_scales(eprm->colX_scale);
      YY(eprm->scale)=eprm->DF_scales(eprm->colX_scale+1);
      ZZ(eprm->scale)=eprm->DF_scales(eprm->colX_scale+2);
      eprm->DF_scales.next_data_line();
   }
   else if (eprm->Docfile!=FALSE){
      eprm->scale.resize(3); 
      eprm->scale.init_constant(1.);
   }
   A(0,0)=XX(eprm->scale);
   A(1,1)=YY(eprm->scale);
   A(2,2)=ZZ(eprm->scale);
   vol().self_apply_geom(A,IS_NOT_INV,eprm->wrap);
   return TRUE;
}

int main (int argc, char **argv) {
   Shift_Scale_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Shift {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Shift/Help/shift.html";
      help="Shift volumes and images";
      OPEN MENU menu_shift;
      COMMAND LINES {
	+ usual: xmipp_shift
               #include "prog_line.mnu"
                -shift "["$X","$Y[","$Z]"]"
               [-dont_wrap]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
            $X  {type=float; label="Shift X ";}
            $Y  {type=float; label="Shift Y ";}
        OPT($Z) {type=float; label="Shift Z ";}
        OPT(-dont_wrap) {label="Do not wrap";}
      }
   }

   MENU menu_shift {
      #include "prog_menu.mnu"
      "Shift parameters"
      $X
      $Y
      OPT($Z)
      OPT(-dont_wrap)
   }
*/
