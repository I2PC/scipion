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

#include <data/args.h>
#include <classification/training_vector.h>

#include <fstream>
#include <cstdio>

int main(int argc, char **argv) {
  char *iname, *oname;
  float offX, offY, offZ;
  int dimX, dimY, dimZ;
  double radius, vsize;

  // Read arguments

  try {
    iname = get_param(argc, argv, "-iname");
    oname = get_param(argc, argv, "-oname");
    radius = (double) AtoF(get_param(argc, argv, "-radius", "0"));
    vsize = (double) AtoF(get_param(argc, argv, "-vsize", "1"));
    dimX = AtoI(get_param(argc, argv, "-dimX", "60"));
    dimY = AtoI(get_param(argc, argv, "-dimY", FtoA(dimX).c_str()));
    dimZ = AtoI(get_param(argc, argv, "-dimZ", FtoA(dimX).c_str()));
    offX = AtoF(get_param(argc, argv, "-offX", "0"));
    offY = AtoF(get_param(argc, argv, "-offY", FtoA(offX).c_str()));
    offZ = AtoF(get_param(argc, argv, "-offZ", FtoA(offX).c_str()));
  }
  catch (Xmipp_error) {
    cout << "Usage:" << endl;
    cout << "-iname       : Input file name (Data File)" << endl;
    cout << "-oname       : Output file name (Phantom Description File)" << endl;
    cout << "-radius      : Radius of the sphere (Default: none)" << endl;
    cout << "-vsize       : Voxel size in Angstrom (Default: 1)" << endl;
    cout << "-dimX        : X dimension (Default 60)" << endl;
    cout << "-dimY        : Y dimension (Default dimX)" << endl;
    cout << "-dimZ        : Z dimension (Default dimX)" << endl;
    cout << "-offX        : X offset (Default 0)" << endl;
    cout << "-offY        : Y offset (Default offX)" << endl;
    cout << "-offZ        : Z offset (Default offX)" << endl;
    exit(1);
   }

  cout << endl << "Given parameters are: " << endl;
  cout << "iname = " << iname << endl;
  cout << "oname = " << oname << endl;
  if (radius > 0)
	  cout << "radius = " << radius << endl;
  cout << "vsize = " << vsize << endl;
  cout << "dimX  = "<< dimX << endl;
  cout << "dimY  = "<< dimY << endl;
  cout << "dimZ  = "<< dimZ << endl;
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

  fprintf(fout, "%8d%8d%8d  0\n", dimX, dimY, dimZ);
  cout << "Converting to Phantom file...." << endl;
  for (int i = 0; i < ts.size(); i++) {
      if (radius > 0)
      	fprintf(fout, "sph = 1 %8.3f%8.3f%8.3f%8.3f\n", (ts.theItems[i][0] + offX)*vsize,(ts.theItems[i][1] + offY)*vsize,(ts.theItems[i][2] + offZ)*vsize, radius);
      else
        fprintf(fout, "sph = 1 %8.3f%8.3f%8.3f  xxx\n", (ts.theItems[i][0] + offX)*vsize,(ts.theItems[i][1] + offY)*vsize,(ts.theItems[i][2] + offZ)*vsize);
  }

  fclose(fout);      // close output file

  cout << "Done!" << endl << endl;
  exit(0);
}


