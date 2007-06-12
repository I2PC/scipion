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
#include <data/image.h>
#include <classification/training_vector.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

int main(int argc, char **argv)
{
    FILE *fp;
    float T;
    float tmpR;
//  char *fname, *iname, *bmname, *imgName, *ext;
    FileName fname, iname, bmname, imgName, ext;
    string selname, basename;
    ImageXmipp mask;
    vector < vector <float> > dataPoints;
    vector < string > labels;
    bool nomask = false;
    bool noBB = true;
    FileName  tmpN;
    int rows, cols;


    // Read arguments

    try
    {

        fname = getParameter(argc, argv, "-i");
        basename = fname.get_baseName();
        selname = basename + (string) ".sel";
        selname = getParameter(argc, argv, "-o", selname.c_str());
        imgName = getParameter(argc, argv, "-imgName", basename.c_str());
        ext = getParameter(argc, argv, "-ext", "spi");
        bmname = getParameter(argc, argv, "-mask", "mask.spi");
        if (checkParameter(argc, argv, "-nomask"))
        {
            nomask = true;
            rows = textToInteger(getParameter(argc, argv, "-rows"));
            cols = textToInteger(getParameter(argc, argv, "-cols"));
        }
        if (checkParameter(argc, argv, "-noBB"))
            noBB = true;
    }
    catch (Xmipp_error)
    {
        cout << "data2img: Convert a data set into a set of images" << endl;
        cout << "Usage:" << endl;
        cout << "-i             : Input file name (iname)" << endl;
        cout << "[-o]           : Output sel file name (default: iname.sel)" << endl;
        cout << "[-imgName]     : first letters of the images' names (default: iname)" << endl;
        cout << "[-ext]         : Extension of the output images (default: spi)" << endl;
        cout << "[-mask]        : Input Mask file name (default: mask.spi)" << endl;
//    cout << "[-noBB]        : Images will be inside the Mask's bounding box (default: yes)" << endl;
        cout << "[-nomask]      : set if the mask is not going to be used" << endl;
        cout << "[-rows]        : Number of rows if the mask is not going to be used" << endl;
        cout << "[-cols]        : Number of columns if the mask is not going to be used" << endl;
        exit(1);
    }

    cout << "Given parameters are: " << endl;
    cout << "output = " << selname << endl;
    if (!nomask)
        cout << "mask = " << bmname << endl;
    else
    {
        cout << "No mask is going to be used" << endl;
        cout << "Number of rows of the generated images: " << rows << endl;
        cout << "Number of columns of the generated images: " << cols << endl;
    }
//  if (!noBB)
//     cout << "Generated images will be inside the mask's bounding box" << endl;
    cout << "input = " << fname << endl;
    cout << "imgName = " << imgName << endl;

    // Read spider mask
    if (!nomask)
    {
        cout << endl << "reading mask " << bmname << "......" << endl << endl;
        mask.read(bmname);           // Reads the mask
        //Adjust the range to 0-1
        mask().range_adjust(0, 1);   // just in case
        if (noBB)
            mask().setXmippOrigin();   // sets origin at the center of the mask.
        cout << mask;       // Output Volumen Information
    }

    int minXPixel = 32000, maxXPixel = 0; int minYPixel = 32000, maxYPixel = 0;
    int NewXDim, NewYDim;
    if ((!noBB) && (!nomask))
    {
        cout << endl << "Calculating the mask's minimum bounding box...." << endl;
        for (int y = 0; y < mask().rowNumber(); y++)
            for (int x = 0; x < mask().colNumber(); x++)
            {
                // Checks if pixel is zero (it's outside the binary mask)
                if (mask(y, x) != 0)
                {
                    if (y < minYPixel) minYPixel = y;
                    if (x < minXPixel) minXPixel = x;
                    if (y > maxYPixel) maxYPixel = y;
                    if (x > maxXPixel) maxXPixel = x;
                }
            } // for x
        NewXDim = (maxXPixel - minXPixel) +  1;
        NewYDim = (maxYPixel - minYPixel) +  1;
        cout << "minX = " << minXPixel << " maxX = " << maxXPixel << " DimX = " << NewXDim << endl;
        cout << "minY = " << minYPixel << " maxY = " << maxYPixel << " DimY= " << NewYDim << endl;
        mask().moveOriginTo(minYPixel + NewYDim / 2, minXPixel + NewXDim / 2);   // sets origin at the center of the mask.
    }

    cout << endl << "Reading input file...." << endl;

    ifstream iStream(fname.c_str());
    if (!iStream)
    {
        cerr << argv[0] << ": can't open file " << iname << endl;
        exit(EXIT_FAILURE);
    }
    xmippCTVectors ts(0, true);
    iStream >> ts;

    FILE  *fout;
    fout = fopen(selname.c_str(), "w");
    if (fout == NULL)
    {
        cerr << argv[0] << ": can't open file " << selname << endl;
        exit(EXIT_FAILURE);
    }

    if (nomask && (rows*cols != ts.theItems[0].size()))
    {
        cerr << argv[0] << ": Images size doesn't coincide with data file " << endl;
        exit(EXIT_FAILURE);
    }

    cout << "generating images......" << endl;

    for (int i = 0; i < ts.size(); i++)
    {
        ImageXmipp image;
        if (nomask)
            image().resize(rows, cols);   // creates image
        else
        {
            if (noBB)
                image().resize(mask());         // creates image
            else
                image().resize(NewYDim, NewXDim);
        }
        image().setXmippOrigin();       // sets origin at the center of the image.
        int counter = 0;
        double minVal = MAXFLOAT;
        for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
            for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++)
            {
                // Checks if pixel is different from zero (it's inside the binary mask)
                if (!noBB || nomask || mask(y, x) != 0)
                {
                    image(y, x) = (double) ts.theItems[i][counter];
                    if (ts.theItems[i][counter] < minVal)
                        minVal = ts.theItems[i][counter];
                    counter++;
                }
            } // for x
        if (!nomask && noBB)
        {
            for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
                for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++)
                {
                    // Checks if pixel is zero (it's outside the binary mask)
                    if (nomask || mask(y, x) == 0)
                        image(y, x) = (double) minVal;
                } // for x
        } // if nomask.

        tmpN = (string) imgName + integerToString(i) + (string) "." + (string) ext;
        image.write(tmpN);
        fprintf(fout, "%s 1 \n", tmpN.c_str());
    }
    fclose(fout);      // close output file
    exit(0);
}


