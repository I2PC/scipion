/***************************************************************************
 *
 * Authors:  	David Maluenda (dmaluenda@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "coordinates_carbon_filter.h"

void ProgCoordinatesCarbonFilter::defineParams()
{
    // Parameters
    addParamsLine(" -c <coordinates> : Input coordinates");
    addParamsLine(" -m <micrograph> : Reference volume");
    addParamsLine(" [--patchSize <n=10>] : Patch size for the variance/mean filter");
    addParamsLine(" --oroot <fnRoot> : Rootname for the output coordinates");
}

void ProgCoordinatesCarbonFilter::readParams()
{
	fnInCoord = getParam("-c");
	fnInMic = getParam("-m");
	fnRoot = getParam("--oroot");
	patchSize = getIntParam("--patchSize");
}

void ProgCoordinatesCarbonFilter::show()
{
    if (verbose)
		std::cout
		<< "Input Coordinates:  " << fnInCoord  << std::endl
		<< "Input Micrograph:   " << fnInMic  << std::endl
		<< "Output rootname:    " << fnRoot << std::endl
		<< "Patch size:         " << patchSize << std::endl
		;
}

//#define DEBUG
void ProgCoordinatesCarbonFilter::run()
{
    show();

    int patchSize_2 = patchSize/2;

    Micrograph mic;
    mic.open_micrograph(fnInMic);
    mic.read_coordinates(0, (String)"particles@"+fnInCoord);
	mic.add_label("");

	int nPart = mic.ParticleNo();

	std::cout << "Number of particles: " << nPart << std::endl; 

	Image<double> im;
	im.read(fnInMic);

    MultidimArray<double> &matrixMic = im();
    
    // Appling filters
    MultidimArray<double> binMask=matrixMic;
    noisyZonesFilter(binMask, patchSize);

    // if (verbose>1)
    // {
        Image<double> imComb(binMask);
        imComb.write(fnRoot+"_mask.mrc");
    // }


    std::vector<Particle_coords> goodCoords;
    std::vector<Particle_coords> badCoords;
    std::vector<Particle_coords> &allCoords=mic.coords;
    for (size_t i=0; i<allCoords.size(); ++i)
    {
    	Particle_coords &coord=allCoords[i];
    	if (DIRECT_A2D_ELEM(binMask,coord.Y,coord.X))
    		goodCoords.push_back(coord);
    	else
    		badCoords.push_back(coord);
    }
    mic.coords = goodCoords;
    mic.write_coordinates(0,0,(String)"particles@"+fnRoot+"_good.pos");
    mic.coords = badCoords;
    mic.write_coordinates(0,0,(String)"particles@"+fnRoot+"_bad.pos");
    std::cout << "        - Good particles: " << goodCoords.size() << std::endl;
    std::cout << "        - Bad particles: " << badCoords.size() << std::endl;
}
#undef DEBUG
