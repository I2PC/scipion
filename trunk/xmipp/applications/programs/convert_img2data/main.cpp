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
  ImageXmipp mask;
  vector < vector <float> > dataPoints;
  vector < string > labels;
  bool nomask = false;
  bool verb = false;
  bool apply_geo = true;
  bool radial_avg = false;

  // Read arguments
  try {
    selname = get_param(argc, argv, "-sel");
    fname = get_param(argc, argv, "-fname", "out.dat");
    bmname = get_param(argc, argv, "-mname", "mask.spi");
    if (check_param(argc, argv, "-nomask"))
     nomask = true;
    if (check_param(argc, argv, "-verb"))
     verb = true;
    apply_geo=!check_param(argc,argv,"-dont_apply_geo");
    radial_avg=check_param(argc,argv,"-radial_avg");
  }
  catch (Xmipp_error) {
    cout << "img2data: Convert a set of images into a set of data vectors" << endl;
    cout << "Usage:" << endl;
    cout << "-sel                 : Sel file name" << endl;
    cout << "-mname               : Input Mask file name (default: mask.spi)" << endl;
    cout << "[-nomask]            : set if the mask is not going to be used" << endl;
    cout << "[-radial_avg]        : set if only the radial avg should be output" << endl;
    cout << "[-fname]             : Output file name (default: out.dat)" << endl;
    cout << "[-verb]              : Verbosity (default: false)" << endl;
    cout << "[-dont_apply_geo]    : Do not apply transformation stored in the header of 2D-images"<< endl;
    exit(1);
   }

try {
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
       cout << mask;             // Output Volumen Information
  }

  cout << "generating data......" << endl;

  SelFile SF((FileName) selname);
  // Read Sel file
  while (!SF.eof()) {
      string image_name = SF.NextImg();
      if (verb)
        cout << "generating points for image " << image_name << "......" << endl;
      ImageXmipp image(image_name,apply_geo);      // reads image

      // Extract the data
      image().set_Xmipp_origin();  // sets origin at the center of the image.
      mask().set_Xmipp_origin();   // sets origin at the center of the mask.

      // Generates coordinates (data points)
      vector<float> imagePoints;

      if (!radial_avg) {
         // If pixel mode
         for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
           for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++) {
	   // Checks if pixel is different from zero (it's inside the binary mask)
              bool cond;
	      if (!nomask) cond = mask(y,x) != 0;
	      else         cond = true;
	      if (cond)
	           imagePoints.push_back(image(y,x));
         } // for x
      } else {
         // If radial average mode
         // Apply the mask
         if (!nomask)
            FOR_ALL_ELEMENTS_IN_MATRIX2D(image())
               if (mask(i,j)==0) image(i,j)=0;

         // Compute the radial average
         matrix1D<int> center_of_rot(2); VECTOR_R2(center_of_rot,0,0);
         matrix1D<int> radial_count;
         matrix1D<double> radial_mean;
         radial_average(image(),center_of_rot,radial_mean,radial_count);

         // Copy radial_mean to vector<float>
         FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
            imagePoints.push_back((float)radial_mean(i));
      }

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
    }catch (Xmipp_error XE) {cout << XE;}
   fclose(fp);    // close file
   exit(0);
}


