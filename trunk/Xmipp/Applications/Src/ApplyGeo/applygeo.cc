/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
     cerr << "Inside Image\n";
     //set shifts to zero
     img.Xoff()=0.; img.Yoff()=0.;
     img.psi()=0.; if(img.tilt()==0) img.rot()=0;
     //set angles to zero
     
     return TRUE;
}//images end

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
     cerr << "Applygeo does not work with volumes\n";
     return TRUE;
}//volume end

int main (int argc, char **argv) {
   Prog_parameters prm;
   prm.apply_geo=TRUE;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
