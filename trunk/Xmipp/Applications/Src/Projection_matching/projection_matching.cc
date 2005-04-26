/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)   
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
#include <Reconstruction/Programs/Prog_projection_matching.hh> 

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {

  double                        sumCC,sumZ;
  DocFile                       DFo;
  Prog_projection_matching_prm  prm;
  FileName                      fn_tmp;

  // Get input parameters
  try {
    // Read command line & produce side info
    prm.read(argc,argv);
    prm.show();

    prm.produce_Side_info();
    prm.project_reference_volume();

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}
    
  try {

    DFo.clear();
    DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), maxCC (6), Z-score (7)");

    // Process all images
    prm.PM_loop_over_all_images(prm.SF,DFo,sumCC,sumZ);
 
    cerr << " Average maxCC = "<<sumCC/prm.SF.ImgNo()<<" average Z-score = "<<sumZ/prm.SF.ImgNo() <<endl;
    fn_tmp=prm.fn_root+".doc";
    DFo.write(fn_tmp);

  } catch (Xmipp_error XE) {cout << XE; prm.usage(); exit(0);}

}




