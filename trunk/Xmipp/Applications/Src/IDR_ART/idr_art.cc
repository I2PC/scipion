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

#include <Reconstruction/Programs/Prog_IDR_art.hh>

/* ------------------------------------------------------------------------- */
/* Main                                                                      */
/* ------------------------------------------------------------------------- */
int main (int argc, char *argv[]) {
// Variables
   Prog_IDR_ART_Parameters   idr_art_prm;
   VolumeXmipp               vol_recons;

// Read Art Parameters
   try {
      idr_art_prm.read(get_param(argc,argv,"-i"));
   } catch (Xmipp_error &XE) {
      cout << XE;
      cerr << "Purpose: Reconstruct and correct CTF\n"
           << "Usage: idr_art -i <parameter file>\n"
	   << "   With the following parameters inside\n"
	   << "-------------------------------------------------------------\n";
      idr_art_prm.Usage();
      exit(1);
   }

// Call main ART routine
   try {
      idr_art_prm.produce_side_info();
      Basic_ROUT_IDR_Art(idr_art_prm,vol_recons);
   } catch (Xmipp_error XE) {cout << XE; exit(1);}
}

