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

#include <data/volume.h>
#include <data/args.h>
#include <data/selfile.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

int main(int argc, char **argv) {
  FILE *fp;
  float tmpR;
  char *fname, *iname, *bmname;
  string selname;
  VolumeXmipp mask;
  vector < vector <float> > dataPoints;
  vector < string > labels;
  bool nomask = false;
  bool verb = false;


  // Read arguments

  try {
    selname = get_param(argc, argv, "-sel");
    fname = get_param(argc, argv, "-fname", "out.dat");
    bmname = get_param(argc, argv, "-mname", "mask.spi");
    if (check_param(argc, argv, "-nomask"))
     nomask = true;
    if (check_param(argc, argv, "-verb"))
     verb = true;
  }
  catch (Xmipp_error) {
    cout << "img2data: Convert a set of volumes into a set of data vectors" << endl;
    cout << "Usage:" << endl;
    cout << "-sel           : Sel file name" << endl;
    cout << "-mname         : Input Mask file name (default: mask.spi)" << endl;
    cout << "[-nomask]      : set if the mask is not going to be used" << endl;
    cout << "[-fname]       : Output file name (default: out.dat)" << endl;
    cout << "[-verb]        : Verbosity (default: false)" << endl;
    exit(1);
   }


  cout << "Given parameters are: " << endl;
  cout << "sel = " << selname << endl;
  if (!nomask) {
     cout << "mname = " << bmname << endl;
  } else
      cout << "No mask is going to be used" << endl;
  cout << "fname = " << fname << endl;


  // Read spider mask

   if (!nomask) {
   	cout << endl << "reading mask " << bmname << "......" << endl << endl;
   	mask.read(bmname);        // Reads the mask
        //Adjust the range to 0-1
        mask().range_adjust(0, 1); // just in case
   	cout << mask;		  // Output Volumen Information
   }

  cout << "generating data......" << endl;

  SelFile SF((FileName) selname);
  // Read Sel file
  while (!SF.eof()) {
      string image_name = SF.NextImg();
      if (verb)
        cout << "generating points for image " << image_name << "......" << endl;
      VolumeXmipp image(image_name);      // reads image

      // Extract the data

      image().set_Xmipp_origin();  // sets origin at the center of the image.
      mask().set_Xmipp_origin();   // sets origin at the center of the mask.

      // Generates coordinates (data points)

      vector <float> imagePoints;

     for (int z = STARTINGZ(image()); z <= FINISHINGZ(image()); z++)
      for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
        for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++) {
	// Checks if voxel is different from zero (it's inside the binary mask)
           bool cond;
	   if (!nomask)
	      cond = mask(z,y,x) != 0;
	   else
	      cond = true;
	   if (cond) {
	        imagePoints.push_back(image(z,y,x));
	   }
        } // for x
      labels.push_back(image_name);
      dataPoints.push_back(imagePoints);
   } // while

   cout << endl << "Saving points......" << endl;
   fp=fopen(fname,"w");
   fprintf(fp, "%d %d\n", dataPoints[0].size(), dataPoints.size()); // Save dimension
   for (unsigned i=0; i < dataPoints.size(); i++) {
       for (unsigned j=0; j < dataPoints[0].size(); j++)
           fprintf(fp, "%3.3f ", dataPoints[i][j]);
       fprintf(fp, "%s \n", labels[i].c_str());
    }
   fclose(fp);    // close file

    exit(0);
}


