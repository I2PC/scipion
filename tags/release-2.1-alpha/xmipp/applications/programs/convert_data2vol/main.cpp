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

int main(int argc, char **argv)
{
    FILE *fp;
    float T;
    float tmpR;
    char *fname, *iname, *bmname, *imgName, *ext;
    std::string selname;
    VolumeXmipp mask;
    std::vector < std::vector <float> > dataPoints;
    std::vector < std::string > labels;
    bool nomask = false;
    bool noBB = true;
    FileName  tmpN;
    int rows, cols, planes;


    // Read arguments

    try
    {
        selname = getParameter(argc, argv, "-sel");
        fname = getParameter(argc, argv, "-iname", "out.dat");
        imgName = getParameter(argc, argv, "-imgName", "img");
        ext = getParameter(argc, argv, "-ext", "spi");
        bmname = getParameter(argc, argv, "-mname", "mask.spi");
        if (checkParameter(argc, argv, "-nomask"))
        {
            nomask = true;
            rows = textToInteger(getParameter(argc, argv, "-rows"));
            cols = textToInteger(getParameter(argc, argv, "-cols"));
            planes = textToInteger(getParameter(argc, argv, "-planes"));
        }
        if (checkParameter(argc, argv, "-noBB"))
            noBB = true;
    }
    catch (Xmipp_error)
    {
        std::cout << "data2img: Convert a data set into a set of volumes" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "-sel           : Output sel file name" << std::endl;
        std::cout << "[-iname]       : Input file name (default: out.dat)" << std::endl;
        std::cout << "[-imgName]     : first letters of the images' names (default: img)" << std::endl;
        std::cout << "[-ext]         : Extension of the output volumes (default: spi)" << std::endl;
        std::cout << "[-mname]       : Input Mask file name (default: mask.spi)" << std::endl;
        std::cout << "[-nomask]      : set if the mask is not going to be used" << std::endl;
        std::cout << "[-rows]        : Number of rows if the mask is not going to be used" << std::endl;
        std::cout << "[-cols]        : Number of columns if the mask is not going to be used" << std::endl;
        std::cout << "[-planes]      : Number of planes if the mask is not going to be used" << std::endl;
        exit(1);
    }

    std::cout << "Given parameters are: " << std::endl;
    std::cout << "sel = " << selname << std::endl;
    if (!nomask)
        std::cout << "mname = " << bmname << std::endl;
    else
    {
        std::cout << "No mask is going to be used" << std::endl;
        std::cout << "Number of rows of the generated volumes: " << rows << std::endl;
        std::cout << "Number of columns of the generated volumes: " << cols << std::endl;
        std::cout << "Number of planes of the generated volumes: " << planes << std::endl;
    }
    std::cout << "iname = " << fname << std::endl;
    std::cout << "imgName = " << imgName << std::endl;

    // Read spider mask
    if (!nomask)
    {
        std::cout << std::endl << "reading mask " << bmname << "......" << std::endl << std::endl;
        mask.read(bmname);           // Reads the mask
        //Adjust the range to 0-1
        mask().rangeAdjust(0, 1);   // just in case
        //if (noBB)
        mask().setXmippOrigin();   // sets origin at the center of the mask.
        std::cout << mask;       // Output Volumen Information
    }

    std::cout << std::endl << "Reading input file...." << std::endl;

    std::ifstream iStream(fname);
    if (!iStream)
    {
        std::cerr << argv[0] << ": can't open file " << iname << std::endl;
        exit(EXIT_FAILURE);
    }
    xmippCTVectors ts(0, true);
    iStream >> ts;

    FILE  *fout;
    fout = fopen(selname.c_str(), "w");
    if (fout == NULL)
    {
        std::cerr << argv[0] << ": can't open file " << selname << std::endl;
        exit(EXIT_FAILURE);
    }

    if (nomask && (planes*rows*cols != ts.theItems[0].size()))
    {
        std::cerr << argv[0] << ": Images size doesn't coincide with data file " << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "generating volumes......" << std::endl;

    for (int i = 0; i < ts.size(); i++)
    {
        VolumeXmipp image;
        if (nomask)
            image().resize(planes, rows, cols);   // creates image
        else
        {
            image().resize(mask());         // creates image
        }
        image().setXmippOrigin();       // sets origin at the center of the image.
        int counter = 0;
        double minVal = MAXFLOAT;
        for (int z = STARTINGZ(image()); z <= FINISHINGZ(image()); z++)
            for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
                for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++)
                {
                    // Checks if voxel is different from zero (it's inside the binary mask)
                    if (nomask || mask(z, y, x) != 0)
                    {
                        image(z, y, x) = (double) ts.theItems[i][counter];
                        if (ts.theItems[i][counter] < minVal)
                            minVal = ts.theItems[i][counter];
                        counter++;
                    }
                } // for x
        if (!nomask)
        {
            for (int z = STARTINGZ(image()); z <= FINISHINGZ(image()); z++)
                for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
                    for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++)
                    {
                        // Checks if pixel is zero (it's outside the binary mask)
                        if (nomask || mask(z, y, x) == 0)
                            image(z, y, x) = (double) minVal;
                    } // for x
        } // if nomask.

        tmpN = (std::string) imgName + integerToString(i) + (std::string) "." + (std::string) ext;
        image.write(tmpN);
        fprintf(fout, "%s 1 \n", tmpN.c_str());
    }
    fclose(fout);      // close output file
    exit(0);
}


