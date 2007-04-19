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
#include <data/volume.h>
#include <classification/training_vector.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

int main(int argc, char **argv) {
  FILE *fp;
  float T;
  float tmpR;
  char *fname, *iname, *bmname, *imgName, *ext;
  string selname;
  VolumeXmipp mask;
  vector < vector <float> > dataPoints;
  vector < string > labels;
  bool nomask = false;
  bool noBB = true;
  FileName  tmpN;
  int rows, cols, planes;


  // Read arguments

  try {
    selname = get_param(argc, argv, "-sel");
    fname = get_param(argc, argv, "-iname", "out.dat");
    imgName = get_param(argc, argv, "-imgName", "img");
    ext = get_param(argc, argv, "-ext", "spi");
    bmname = get_param(argc, argv, "-mname", "mask.spi");
    if (check_param(argc, argv, "-nomask")) {
      nomask = true;
      rows = AtoI(get_param(argc, argv, "-rows"));
      cols = AtoI(get_param(argc, argv, "-cols"));
      planes = AtoI(get_param(argc, argv, "-planes"));
    }
    if (check_param(argc, argv, "-noBB"))
      noBB = true;
  }
  catch (Xmipp_error) {
    cout << "data2img: Convert a data set into a set of volumes" << endl;
    cout << "Usage:" << endl;
    cout << "-sel           : Output sel file name" << endl;
    cout << "[-iname]       : Input file name (default: out.dat)" << endl;
    cout << "[-imgName]     : first letters of the images' names (default: img)" << endl;
    cout << "[-ext]         : Extension of the output volumes (default: spi)" << endl;
    cout << "[-mname]       : Input Mask file name (default: mask.spi)" << endl;
    cout << "[-nomask]      : set if the mask is not going to be used" << endl;
    cout << "[-rows]        : Number of rows if the mask is not going to be used" << endl;
    cout << "[-cols]        : Number of columns if the mask is not going to be used" << endl;
    cout << "[-planes]      : Number of planes if the mask is not going to be used" << endl;
    exit(1);
   }

  cout << "Given parameters are: " << endl;
  cout << "sel = " << selname << endl;
  if (!nomask)
     cout << "mname = " << bmname << endl;
  else {
      cout << "No mask is going to be used" << endl;
      cout << "Number of rows of the generated volumes: " << rows << endl;
      cout << "Number of columns of the generated volumes: " << cols << endl;
      cout << "Number of planes of the generated volumes: " << planes << endl;
  }
  cout << "iname = " << fname << endl;
  cout << "imgName = " << imgName << endl;

  // Read spider mask
   if (!nomask) {
   	cout << endl << "reading mask " << bmname << "......" << endl << endl;
   	mask.read(bmname);           // Reads the mask
        //Adjust the range to 0-1
        mask().range_adjust(0, 1);   // just in case
	//if (noBB)
           mask().set_Xmipp_origin();   // sets origin at the center of the mask.
   	cout << mask;		     // Output Volumen Information
   }

   cout << endl << "Reading input file...." << endl;

   ifstream iStream(fname);
   if (!iStream) {
     cerr << argv[0] << ": can't open file " << iname << endl;
     exit(EXIT_FAILURE);
   }
   xmippCTVectors ts(0, true);
   iStream >> ts;

   FILE  *fout;
   fout = fopen(selname.c_str(), "w");
   if( fout == NULL ) {
    cerr << argv[0] << ": can't open file " << selname << endl;
    exit(EXIT_FAILURE);
   }

   if (nomask && (planes*rows*cols != ts.theItems[0].size())) {
      cerr << argv[0] << ": Images size doesn't coincide with data file " << endl;
      exit(EXIT_FAILURE);
   }

   cout << "generating volumes......" << endl;

   for (int i = 0; i < ts.size(); i++) {
      VolumeXmipp image;
      if (nomask)
      	  image().resize(planes, rows, cols);   // creates image
      else {
      	 image().resize(mask());         // creates image
      }
      image().set_Xmipp_origin();       // sets origin at the center of the image.
      int counter = 0;
      double minVal = MAXFLOAT;
     for (int z = STARTINGZ(image()); z <= FINISHINGZ(image()); z++)
      for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
        for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++) {
	// Checks if voxel is different from zero (it's inside the binary mask)
	   if (nomask || mask(z,y,x) != 0) {
	        image(z,y,x) = (double) ts.theItems[i][counter];
		if (ts.theItems[i][counter] < minVal)
		   minVal = ts.theItems[i][counter];
		counter++;
	   }
      } // for x
      if (!nomask) {
      for (int z = STARTINGZ(image()); z <= FINISHINGZ(image()); z++)
        for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
          for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++) {
	     // Checks if pixel is zero (it's outside the binary mask)
	     if (nomask || mask(z,y,x) == 0)
	        image(z, y,x) = (double) minVal;
          } // for x
      } // if nomask.

      tmpN = (string) imgName + ItoA(i) + (string) "." + (string) ext;
      image.write(tmpN);
      fprintf(fout, "%s 1 \n", tmpN.c_str());
   }
   fclose(fout);      // close output file
   exit(0);
}


