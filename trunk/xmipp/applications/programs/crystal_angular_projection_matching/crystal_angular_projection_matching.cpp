/***************************************************************************
 *
 * Authors: Roberto Marabini  
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

/* INCLUDES ---------------------------------------------------------------- */
#include <reconstruction/crystal_angular_projection_matching.h> 

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {

  double                        sumCC;
  DocFile                       DFlib,DFo;
  SymList                       SL;
  Prog_projection_matching_crystal_prm  prm;
  FileName                      fn_tmp;

  // Get input parameters
  try {
    // Read command line & produce side info
    prm.read(argc,argv);
    prm.show();

    // Project reference volume etc.
    prm.produce_Side_info();

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}
    
  try {

     DFo.clear();
     DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), \
scale(4), Xoff (5), Yoff (6), Refno (7), maxCC (8)");

    // Process all images
    prm.PM_loop_over_all_images(DFo,sumCC);

//    cerr << " Average maxCC = "<<sumCC/prm.SF.ImgNo()<<endl;
    fn_tmp=prm.fn_root+".doc";
    DFo.write(fn_tmp);

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}

}




