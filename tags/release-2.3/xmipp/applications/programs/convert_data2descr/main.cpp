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

int main(int argc, char **argv)
{
    char *iname, *oname;
    float offX, offY, offZ;
    int dimX, dimY, dimZ;
    double radius, vsize;

    // Read arguments

    try
    {
        iname = getParameter(argc, argv, "-iname");
        oname = getParameter(argc, argv, "-oname");
        radius = (double) textToFloat(getParameter(argc, argv, "-radius", "0"));
        vsize = (double) textToFloat(getParameter(argc, argv, "-vsize", "1"));
        dimX = textToInteger(getParameter(argc, argv, "-dimX", "60"));
        dimY = textToInteger(getParameter(argc, argv, "-dimY", floatToString(dimX).c_str()));
        dimZ = textToInteger(getParameter(argc, argv, "-dimZ", floatToString(dimX).c_str()));
        offX = textToFloat(getParameter(argc, argv, "-offX", "0"));
        offY = textToFloat(getParameter(argc, argv, "-offY", floatToString(offX).c_str()));
        offZ = textToFloat(getParameter(argc, argv, "-offZ", floatToString(offX).c_str()));
    }
    catch (Xmipp_error)
    {
        std::cout << "Usage:" << std::endl;
        std::cout << "-iname       : Input file name (Data File)" << std::endl;
        std::cout << "-oname       : Output file name (Phantom Description File)" << std::endl;
        std::cout << "-radius      : Radius of the sphere (Default: none)" << std::endl;
        std::cout << "-vsize       : Voxel size in Angstrom (Default: 1)" << std::endl;
        std::cout << "-dimX        : X dimension (Default 60)" << std::endl;
        std::cout << "-dimY        : Y dimension (Default dimX)" << std::endl;
        std::cout << "-dimZ        : Z dimension (Default dimX)" << std::endl;
        std::cout << "-offX        : X offset (Default 0)" << std::endl;
        std::cout << "-offY        : Y offset (Default offX)" << std::endl;
        std::cout << "-offZ        : Z offset (Default offX)" << std::endl;
        exit(1);
    }

    std::cout << std::endl << "Given parameters are: " << std::endl;
    std::cout << "iname = " << iname << std::endl;
    std::cout << "oname = " << oname << std::endl;
    if (radius > 0)
        std::cout << "radius = " << radius << std::endl;
    std::cout << "vsize = " << vsize << std::endl;
    std::cout << "dimX  = " << dimX << std::endl;
    std::cout << "dimY  = " << dimY << std::endl;
    std::cout << "dimZ  = " << dimZ << std::endl;
    std::cout << "offX  = " << offX << std::endl;
    std::cout << "offY  = " << offY << std::endl;
    std::cout << "offZ  = " << offZ << std::endl;


    std::cout << std::endl << "Reading input file...." << std::endl;
    std::ifstream iStream(iname);
    if (!iStream)
    {
        std::cerr << argv[0] << ": can't open file " << iname << std::endl;
        exit(EXIT_FAILURE);
    }
    xmippCTVectors ts(0, false);
    iStream >> ts;

    FILE  *fout;
    fout = fopen(oname, "w");
    if (fout == NULL)
    {
        std::cerr << argv[0] << ": can't open file " << oname << std::endl;
        exit(EXIT_FAILURE);
    }

    fprintf(fout, "%8d%8d%8d  0\n", dimX, dimY, dimZ);
    std::cout << "Converting to Phantom file...." << std::endl;
    for (int i = 0; i < ts.size(); i++)
    {
        if (radius > 0)
            fprintf(fout, "sph = 1 %8.3f%8.3f%8.3f%8.3f\n", (ts.theItems[i][0] + offX)*vsize, (ts.theItems[i][1] + offY)*vsize, (ts.theItems[i][2] + offZ)*vsize, radius);
        else
            fprintf(fout, "sph = 1 %8.3f%8.3f%8.3f  xxx\n", (ts.theItems[i][0] + offX)*vsize, (ts.theItems[i][1] + offY)*vsize, (ts.theItems[i][2] + offZ)*vsize);
    }

    fclose(fout);      // close output file

    std::cout << "Done!" << std::endl << std::endl;
    exit(0);
}


