/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

#include <XmippData/xmippTypes.hh>
#include <XmippData/xmippArgs.hh>
#include <Classification/xmippCTVectors.hh>
#include <fstream>
#include <stdio.h>

int main(int argc, char **argv) {
  char *iname, *oname;
  float offX, offY, offZ;
  double vsize;
      
  // Read arguments
  
  try { 
    iname = get_param(argc, argv, "-iname");
    oname = get_param(argc, argv, "-oname");
    vsize = (double) AtoF(get_param(argc, argv, "-vsize", "1"));
    offX = AtoF(get_param(argc, argv, "-offX", "0"));
    offY = AtoF(get_param(argc, argv, "-offY", FtoA(offX).c_str()));
    offZ = AtoF(get_param(argc, argv, "-offZ", FtoA(offX).c_str()));
  } 
  catch (Xmipp_error) {
    cout << "Usage:" << endl;
    cout << "-iname       : Input file name (data file)" << endl;
    cout << "-oname       : Output file name (PDB File)" << endl;
    cout << "-vsize       : Voxel size in Angstrom (Default: 1)" << endl;
    cout << "-offX        : X offset (default = 0)" << endl;
    cout << "-offY        : Y offset (default X)" << endl;
    cout << "-offZ        : Z offset (default X)" << endl;
    exit(1);
   }
    
  cout << endl << "Given parameters are: " << endl;
  cout << "iname = " << iname << endl;
  cout << "oname = " << oname << endl;
  cout << "vsize = " << vsize << endl;
  cout << "offX  = "<< offX << endl;
  cout << "offY  = "<< offY << endl;
  cout << "offZ  = "<< offZ << endl;


  cout << endl << "Reading input file...." << endl; 
  ifstream iStream(iname);
  if (!iStream) {
    cerr << argv[0] << ": can't open file " << iname << endl;
    exit(EXIT_FAILURE);
  }
  xmippCTVectors ts(0, false);
  iStream >> ts;
  
  FILE  *fout;
  fout = fopen(oname, "w");
  if( fout == NULL ) {
    cerr << argv[0] << ": can't open file " << oname << endl;
    exit(EXIT_FAILURE);   
  }

  cout << "Converting to PDB...." << endl; 
  for (int i = 0; i < ts.size(); i++)
      fprintf(fout, "ATOM  %5d XMIP XMIP    1 %11.3f%8.3f%8.3f %5.2f  0.00      XMIP\n",i+1, (ts.theItems[i][0] + offX)*vsize,(ts.theItems[i][1] + offY)*vsize,(ts.theItems[i][2] + offZ)*vsize);
  
  fclose(fout);      // close output file

  cout << "Done!" << endl << endl; 
  exit(0);  
} 


