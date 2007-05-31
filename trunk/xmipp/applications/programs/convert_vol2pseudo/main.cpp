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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

typedef vector< vector< float > > coordinates;

int main(int argc, char **argv)
{
    FILE *fp;
    float T;
    float tmpR;
    char *fname, *vname, *bmname, *vmname;
    float minVoxel, maxVoxel, slope;
    float minCoord, maxCoord;
    int tmpCoord;
    VolumeXmipp mask/*, vol_mask*/;
    coordinates randomCoord;
    bool sampling = true;
    bool nomask = false;
    bool FourD = false;
    int npoints;


    // Read arguments

    try
    {
        vname = getParameter(argc, argv, "-vname");
        T = AtoF(getParameter(argc, argv, "-T", "2"));
        fname = getParameter(argc, argv, "-fname", "out.dat");
        bmname = getParameter(argc, argv, "-bmname", "mask.spi");
        vmname = getParameter(argc, argv, "-vmname", "vol_mask.spi");
        minCoord = AtoF(getParameter(argc, argv, "-minCoord", "0"));
        maxCoord = AtoF(getParameter(argc, argv, "-maxCoord", "10"));
        if (!checkParameter(argc, argv, "-sampling"))
            sampling = false;
        if (checkParameter(argc, argv, "-nomask"))
            nomask = true;
        if (checkParameter(argc, argv, "-4"))
            FourD = true;
        npoints = AtoI(getParameter(argc, argv, "-npoints", "100000"));
    }
    catch (Xmipp_error)
    {
        cout << "Usage:" << endl;
        cout << "-vname         : Input Volume file name" << endl;
        cout << "-bmname        : Input Mask file name (default: mask.spi)" << endl;
        cout << "[-vmname]      : Volume Mask file name (default: vol_mask.spi)" << endl;
        cout << "[-fname]       : Output file name (default: out.dat)" << endl;
        cout << "[-nomask]      : set if the mask is not going to be used" << endl;
        cout << "[-T]           : Threshold (default: 2)" << endl;
        cout << "[-sampling]    : Use to create a dataset based on sampling" << endl;
        cout << "[-npoints]     : Number of points of the dataset (default: 100000)" << endl;
        cout << "[-minCoord]    : 'Density-Coordinates' minimum value (default: 0)" << endl;
        cout << "[-maxCoord]    : 'Density-Coordinates' maximum value (default: 10)" << endl;
        cout << "[-4]           : Set if the dataset will be 4-dimensional (default: false)" << endl;
        exit(1);
    }


    cout << "Given parameters are: " << endl;
    cout << "vname = " << vname << endl;
    if (!nomask)
    {
        cout << "bmname = " << bmname << endl;
        cout << "vmname = " << vmname << endl;
    }
    else
        cout << "No mask is going to be used" << endl;
    cout << "fname = " << fname << endl;
    cout << "T = " << T << endl;
    if (sampling)
    {
        cout << "Sampling method is used to generate the dataset" << endl;
        cout << "number of points = " << npoints << endl;
    }
    else
    {
        cout << "Exhaustive method is used to generate the dataset" << endl;
        cout << "minCoord = " << minCoord << endl;
        cout << "maxCoord = " << maxCoord << endl;
    }
    if (FourD)
        cout << "Dataset will be 4-Dimensional" << endl;


    // Read spider volumen

    cout << endl << "reading volume " << vname << "......" << endl << endl;
    VolumeXmipp V(vname);    // Reads the volumen
    cout << V;      // Output Volumen Information

    // Read spider mask

    if (!nomask)
    {
        cout << endl << "reading mask " << bmname << "......" << endl << endl;
        mask.read(bmname);        // Reads the mask
        cout << mask;    // Output Volumen Information
    }


    // Extract the data

    V().setXmippOrigin();          // sets origin at the center of the volume.
    mask().setXmippOrigin();       // sets origin at the center of the mask.
    VolumeXmipp vol_mask(V);

//  vol_mask().resize(V());       // Resizes volumen_mask.
    vol_mask().setXmippOrigin();   // sets origin at the center of the volumen mask.


    cout << endl << "Finding minimum and maximum......" << endl;

    // Find Minimum and Maximum density values inside the mask.

    minVoxel = MAXFLOAT;
    maxVoxel = -MAXFLOAT;
    for (int z = STARTINGZ(V()); z <= FINISHINGZ(V()); z++)
    {
        for (int y = STARTINGY(V()); y <= FINISHINGY(V()); y++)
            for (int x = STARTINGX(V()); x <= FINISHINGX(V()); x++)
            {
                if (nomask || VOLVOXEL(mask, z, y, x) != 0)
                {
                    if (VOLVOXEL(V, z, y, x) > maxVoxel)
                        maxVoxel = VOLVOXEL(V, z, y, x);
                    if (VOLVOXEL(V, z, y, x) < minVoxel)
                        minVoxel = VOLVOXEL(V, z, y, x);
                } // if
            } // for x
    } // for z

    cout << endl << "minimum: " << minVoxel << endl << "maximum: " << maxVoxel << endl;


    // Generates coordinates (data points)

    cout << endl << "Generating coordinates......" << endl;

    if (!sampling)
        cout << endl << "using exhaustive generation......" << endl;

    for (int z = STARTINGZ(V()); z <= FINISHINGZ(V()); z++)
        for (int y = STARTINGY(V()); y <= FINISHINGY(V()); y++)
            for (int x = STARTINGX(V()); x <= FINISHINGX(V()); x++)
            {

                // Checks if voxel is different from zero (it's inside the binary mask)

                bool cond;
                if (!nomask)
                    cond = VOLVOXEL(mask, z, y, x) != 0;
                else
                    cond = true;
                if (cond && (VOLVOXEL(V, z, y, x) > T))
                {

                    VOLVOXEL(vol_mask, z, y, x) = VOLVOXEL(V, z, y, x);

                    if (!sampling)
                    {
                        // Lineary transform the density into "number of coordinates"

                        if (FourD)
                        {
                            vector <float> v;
                            v.push_back(x);
                            v.push_back(y);
                            v.push_back(z);
                            v.push_back((float)VOLVOXEL(V, z, y, x));
                            randomCoord.push_back(v);
                        }
                        else
                        { // !FourD
                            if (minVoxel != maxVoxel)
                                slope = (double)(maxCoord - minCoord) / (double)(maxVoxel - minVoxel);
                            else
                                slope = 0;

                            tmpCoord = (int)(minCoord + (slope * (double)(VOLVOXEL(V, z, y, x) - minVoxel)));

                            // Saves the coordinates proportional to the density of the voxel.

                            for (unsigned i = 0; i < tmpCoord; i++)
                            {
                                vector <float> v;
                                v.push_back(x);
                                v.push_back(y);
                                v.push_back(z);
                                randomCoord.push_back(v);
                            }
                        } // FourD

                    }
                    else
                    { // sampling
                        // If sampling, then normalize the masked volume (devide by
                        // maxdensity)
                        VOLVOXEL(vol_mask, z, y, x) /= maxVoxel;
                    }

                }
                else
                { // if VOLVOXEL
                    VOLVOXEL(vol_mask, z, y, x) = 0;
                }

            } // for x

    if (sampling)
    { // Generate "npoints" coordinates
        cout << endl << "using sampling generation......" << endl;
        randomize_random_generator();
        for (int i = 0; i < npoints; i++)
        {
            bool found = false;
            int X, Y, Z;
            do
            {
                X = (int) rnd_unif(STARTINGX(vol_mask()), FINISHINGX(vol_mask()));
                Y = (int) rnd_unif(STARTINGY(vol_mask()), FINISHINGY(vol_mask()));
                Z = (int) rnd_unif(STARTINGZ(vol_mask()), FINISHINGZ(vol_mask()));
                if (VOLVOXEL(vol_mask, Z, Y, X) > rnd_unif()) found = true;
            }
            while (!found);
            vector <float> v;
            v.push_back(X);
            v.push_back(Y);
            v.push_back(Z);
            if (FourD)
                v.push_back((float)VOLVOXEL(vol_mask, Z, Y, X));
            randomCoord.push_back(v);
        }
    }

    cout << endl << "Saving masked volumen......" << endl;
    if (!nomask)
        vol_mask.write(vmname);

    cout << endl << "Saving coordinates at random......" << endl;

    fp = fopen(fname, "w");
    if (FourD)
        fprintf(fp, "4 %d\n", randomCoord.size()); // 4-dimensional
    else
        fprintf(fp, "3 %d\n", randomCoord.size()); // 3-dimensional
    if (!sampling)
        random_shuffle(randomCoord.begin(), randomCoord.end());
    for (unsigned i = 0; i < randomCoord.size(); i++)
        if (FourD)
            fprintf(fp, "%d  %d %d %3.2f \n", (int) randomCoord[i][0], (int) randomCoord[i][1], (int) randomCoord[i][2], randomCoord[i][3]);
        else
            fprintf(fp, "%d  %d %d \n", (int) randomCoord[i][0], (int) randomCoord[i][1], (int) randomCoord[i][2]);

    fclose(fp);    // close file

    exit(0);
}


