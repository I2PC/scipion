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

int main(int argc, char **argv)
{
    FILE *fp;
    float tmpR;
    char *fname, *iname, *bmname;
    std::string selname;
    ImageXmipp mask;
    std::vector < std::vector <float> > dataPoints;
    std::vector < std::string > labels;
    bool nomask = false;
    bool verb = false;
    bool apply_geo = true;
    bool radial_avg = false;

    // Read arguments
    try
    {
        selname = getParameter(argc, argv, "-i");
        bmname = getParameter(argc, argv, "-mask", "mask.spi");
        fname = getParameter(argc, argv, "-o", "out.dat");
        if (checkParameter(argc, argv, "-nomask"))
            nomask = true;
        if (checkParameter(argc, argv, "-verb"))
            verb = true;
        apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");
        radial_avg = checkParameter(argc, argv, "-radial_avg");
    }
    catch (Xmipp_error)
    {
        std::cout << "img2data: Convert a set of images into a set of data vectors" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "-i                   : Input selfile name" << std::endl;
        std::cout << "-mask                : Input Mask file name (default: mask.spi)" << std::endl;
        std::cout << "[-o]                 : Output file name (default: out.dat)" << std::endl;
        std::cout << "[-nomask]            : set if the mask is not going to be used" << std::endl;
        std::cout << "[-radial_avg]        : set if only the radial avg should be output" << std::endl;
        std::cout << "[-verb]              : Verbosity (default: false)" << std::endl;
        std::cout << "[-dont_apply_geo]    : Do not apply transformation stored in the header of 2D-images" << std::endl;
        exit(1);
    }

    try
    {
        std::cout << "Given parameters are: " << std::endl;
        std::cout << "sel = " << selname << std::endl;
        if (!nomask)
        {
            std::cout << "mname = " << bmname << std::endl;
        }
        else
            std::cout << "No mask is going to be used" << std::endl;
        std::cout << "fname = " << fname << std::endl;

        // Read spider mask
        if (!nomask)
        {
            std::cout << std::endl << "reading mask " << bmname << "......" << std::endl << std::endl;
            mask.read(bmname);        // Reads the mask
            //Adjust the range to 0-1
            mask().range_adjust(0, 1); // just in case
            std::cout << mask;             // Output Volumen Information
        }

        std::cout << "generating data......" << std::endl;

        SelFile SF((FileName) selname);
        // Read Sel file
        while (!SF.eof())
        {
            std::string image_name = SF.NextImg();
            if (verb)
                std::cout << "generating points for image " << image_name << "......" << std::endl;
            ImageXmipp image(image_name, apply_geo);     // reads image

            // Extract the data
            image().setXmippOrigin();  // sets origin at the center of the image.
            mask().setXmippOrigin();   // sets origin at the center of the mask.

            // Generates coordinates (data points)
            std::vector<float> imagePoints;

            if (!radial_avg)
            {
                // If pixel mode
                for (int y = STARTINGY(image()); y <= FINISHINGY(image()); y++)
                    for (int x = STARTINGX(image()); x <= FINISHINGX(image()); x++)
                    {
                        // Checks if pixel is different from zero (it's inside the binary mask)
                        bool cond;
                        if (!nomask) cond = mask(y, x) != 0;
                        else         cond = true;
                        if (cond)
                            imagePoints.push_back(image(y, x));
                    } // for x
            }
            else
            {
                // If radial average mode
                // Apply the mask
                if (!nomask)
                    FOR_ALL_ELEMENTS_IN_MATRIX2D(image())
                    if (mask(i, j) == 0) image(i, j) = 0;

                // Compute the radial average
                Matrix1D<int> center_of_rot(2);
                VECTOR_R2(center_of_rot, 0, 0);
                Matrix1D<int> radial_count;
                Matrix1D<double> radial_mean;
                radialAverage(image(), center_of_rot, radial_mean, radial_count);

                // Copy radial_mean to std::vector<float>
                FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
                imagePoints.push_back((float)radial_mean(i));
            }

            labels.push_back(image_name);
            dataPoints.push_back(imagePoints);
        } // while

        std::cout << std::endl << "Saving points......" << std::endl;
        fp = fopen(fname, "w");
        fprintf(fp, "%d %d\n", dataPoints[0].size(), dataPoints.size()); // Save dimension
        for (unsigned i = 0; i < dataPoints.size(); i++)
        {
            for (unsigned j = 0; j < dataPoints[0].size(); j++)
                fprintf(fp, "%3.3f ", dataPoints[i][j]);
            fprintf(fp, "%s \n", labels[i].c_str());
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    fclose(fp);    // close file
    exit(0);
}


