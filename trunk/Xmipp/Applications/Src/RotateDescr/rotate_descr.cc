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
#include <XmippData/xmippGeometry.hh>
#include <Reconstruction/phantom.hh>

class Rotate_descr_parameters {
public:
   FileName fn_in, fn_out;
   bool Euler_mode;
      double rot, tilt, psi;
   bool Align_mode;
   bool Axis_mode;
      matrix1D<double> axis;
      double ang;
   
   matrix2D<double> A3D;

   void read(int argc, char **argv) {
      fn_in=get_param(argc,argv,"-i");
      fn_out=get_param(argc,argv,"-o","");
      if (fn_out=="") fn_out=fn_in;
      Euler_mode=Align_mode=Axis_mode=false;
      if (check_param(argc,argv,"-euler")) {
         Euler_mode=true;
         int i=position_param(argc,argv,"-euler");
         if (i+3>=argc)
            REPORT_ERROR(1,"Not enough parameters after -euler");
         rot  = AtoF(argv[i+1]);
         tilt = AtoF(argv[i+2]);
         psi  = AtoF(argv[i+3]);
         A3D=Euler_rot3D_matrix(rot,tilt,psi);
      } else if (check_param(argc,argv,"-align_with_Z")) {
         Align_mode=true;
         axis=get_vector_param(argc,argv,"-align_with_Z",3);
         A3D=align_with_Z(axis);
      } else {
         Axis_mode=true;
         if (check_param(argc, argv, "-axis"))
            axis=get_vector_param(argc, argv, "-axis", 3);
         else 
            axis=vector_R3(0.,0.,1.);
         ang=AtoF(get_param(argc, argv, "-ang"));
         A3D=rot3D_matrix(ang,axis);
      }
      A3D.resize(3,3);
   }
   
   void show() {
      cout << "Input phantom:  " << fn_in << endl;
      cout << "Output phantom: " << fn_out << endl;
      if (Euler_mode)
         cout << "Euler angles (rot, tilt, psi): " << rot << " " << tilt
              << " " << psi << endl;
      else if (Align_mode)
         cout << "Aligning " << axis.transpose() << " with Z\n";
      else if (Axis_mode)
         cout << "Rotating " << ang << " degrees around " << axis.transpose()
              << endl;
   }

   void usage() {
      cerr << "Usage: rotate_descr\n"
           << "   -i <Phantom.descr>               : Input phantom description file\n"
           << "  [-o <Phantom.descr>]              : By default, the input is rewritten\n"
      ;
      cerr << "  [-euler <rot> <tilt> <psi>        : Rotate with these Euler angles\n"
           << "  [-align_with_Z \"[<x>,<y>,<z>]\"]   : Align (x,y,z) with Z\n"
           << "                                      Notice that brackets for the\n"
           << "                                      vector must be written and do not\n"
           << "                                      represent optional parameters\n"
           << "  [-axis \"[<x>,<y>,<z>]\" -ang <ang>]: Rotate <ang> degrees around (x,y,z),\n"
           << "                                      by default (0,0,1)\n"
      ;
   }
};

int main (int argc, char **argv) {
   Rotate_descr_parameters prm;
   try {
      prm.read(argc,argv);
   } catch (Xmipp_error XE) {cout << XE; prm.usage(); return 1;}
   
   try {
      Phantom phantom;
      phantom.read(prm.fn_in);
      phantom.rotate(prm.A3D);
      phantom.write(prm.fn_out);
   } catch (Xmipp_error XE) {cout << XE; return 2;}
   return 0;
}
