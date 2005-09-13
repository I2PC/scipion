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

#include "../Prog_symmetrize.hh"
#include <XmippData/xmippArgs.hh>

/* Read parameters --------------------------------------------------------- */
void Symmetrize_Parameters::read(int argc, char **argv) {
    fn_in  = get_param(argc,argv,"-i");
    fn_out = get_param(argc,argv,"-o","");
    fn_sym = get_param(argc,argv,"-sym");
    do_not_generate_subgroup=check_param(argc,argv,"-no_group");
    wrap=!check_param(argc,argv,"-dont_wrap");
}

/* Usage ------------------------------------------------------------------- */
void Symmetrize_Parameters::usage() {
  cout
     << "Usage: symmetrize [Purpose and Parameters]\n"
     << "Purpose: Symmetrize a 3D volume\n"
     << "Parameter Values: (notice space before value)\n"
     << "    -i <file_in>        : input 3D Xmipp file\n"
     << "   [-o <file_out>]      : if no name is given then the input file is\n"
     << "                          rewritten\n"
     << "    -sym <sym_file>     : symmetry file (see the manual)\n"
     << "   [-no_group]          : do not generate symmetry subgroup\n"
     << "   [-dont_wrap]         : by default, the volume is wrapped\n";
}

/* Show -------------------------------------------------------------------- */
ostream & operator << (ostream &out, const Symmetrize_Parameters &prm) {
   out << "File in:       " << prm.fn_in  << endl
       << "File out:      " << prm.fn_out << endl
       << "Symmetry file: " << prm.fn_sym << endl
       << "Generate group:" << !prm.do_not_generate_subgroup << endl
       << "Wrapping:      "; print(out,prm.wrap); out << endl;
   return out;
}

/* Really symmetrize ------------------------------------------------------- */
//#define DEBUG
void symmetrize(const SymList &SL, VolumeXmipp &V_in, VolumeXmipp &V_out,
   bool wrap, bool show_progress) {
   matrix2D<double> L(4,4), R(4,4); // A matrix from the list
   VolumeXmipp V_aux, V_aux2;
   matrix1D<double> sh(3);
   V_out=V_in;
   if (show_progress) {
      cerr << "Symmetrizing ...\n";
      init_progress_bar(SL.SymsNo());
   }

   for (int i=0; i<SL.SymsNo(); i++) {
      SL.get_matrices(i,L,R);

   SL.get_shift(i,sh);
   R(3,0)=sh(0) * V_aux().ColNo();
   R(3,1)=sh(1) * V_aux().RowNo();
   R(3,2)=sh(2) * V_aux().SliNo();

      /* *** CO: I don't know why the compiler doesn't allow me
         to reuse V_in !!!, this is very memory wasting */
      apply_geom(V_aux(),R.transpose(),V_in(),IS_NOT_INV,wrap);
      #ifdef DEBUG
         V_aux.write((string)"PPPsym_"+ItoA(i)+".vol");
      #endif

      /* *** CO: I am not very sure about the reason for this, but it
         seems to work */
      // apply_geom(V_aux2(),L,V_aux(),IS_NOT_INV,prm.wrap);
      array_by_array(V_out(),V_aux(),V_out(),'+');
      if (show_progress) progress_bar(i);
   }
   if (show_progress) progress_bar(SL.SymsNo());
   array_by_scalar(V_out(),SL.SymsNo()+1.0f,V_out(),'/');
}
#undef DEBUG

/* Main program ------------------------------------------------------------ */
void ROUT_symmetrize(const Symmetrize_Parameters &prm) {
   SymList         SL;
   VolumeXmipp     V_in;
   VolumeXmipp     V_out;
   FileName        fn_out;

   double accuracy=(prm.do_not_generate_subgroup)?-1:1e-6;
   SL.read_sym_file(prm.fn_sym,accuracy);
   V_in.read(prm.fn_in);
   
   cerr << prm;
   symmetrize(SL,V_in,V_out,prm.wrap,TRUE);
   if (prm.fn_out=="") fn_out=V_in.name(); else fn_out=prm.fn_out;
   V_out.write(fn_out);
}
