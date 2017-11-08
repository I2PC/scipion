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

#include "coordinates_noisy_zones_filter.h"

void ProgCoordinatesNoisyZonesFilter::defineParams()
{
    // Parameters
    addParamsLine(" --pos <coordinates> : Input coordinates");
    addParamsLine(" --mic <micrograph> : Reference volume");
    addParamsLine(" [--patchSize <n=10>] : Patch size for the variance/mean filter");
    addParamsLine(" [--giniCoeff <g=0.3>] : Gini coefficient threshold to applay the filter");
    addParamsLine(" --oroot <fnRoot> : Rootname for the output coordinates");
}

void ProgCoordinatesNoisyZonesFilter::readParams()
{
	fnInCoord = getParam("--pos");
	fnInMic = getParam("--mic");
	fnRoot = getParam("--oroot");
	patchSize = getIntParam("--patchSize");
    giniCoeffTh = getDoubleParam("--giniCoeff");
}

void ProgCoordinatesNoisyZonesFilter::show()
{
    if (verbose)
		std::cout
		<< "Input Coordinates:     " << fnInCoord  << std::endl
		<< "Input Micrograph:      " << fnInMic  << std::endl
		<< "Output rootname:       " << fnRoot << std::endl
		<< "Patch size:            " << patchSize << std::endl
        << "Gini coeff. threshold: " << giniCoeffTh << std::endl
		;
}

//#define DEBUG
void ProgCoordinatesNoisyZonesFilter::run()
{
    show();

    Micrograph mic;
    mic.open_micrograph(fnInMic);
    mic.read_coordinates(0, (String)"particles@"+fnInCoord);
	mic.add_label("");

	int nPart = mic.ParticleNo();

	Image<double> im;
	im.read(fnInMic);

    MultidimArray<double> &matrixMic = im();
    
    // Appling filters
    MultidimArray<double> binMask=matrixMic;
    double giniV = grainFilter(binMask, patchSize, giniCoeffTh);

    // if (giniV<giniCoeffTh)
    // {
    //     std::cout << "   > > >   " << fnInMic << " : " 
    //               << "   > > >   G I N I   Coef : " 
    //               << giniV << "   > > >   OK!" << std::endl << std::endl;
    // }else if (giniV<(1-giniCoeffTh)){
    //     std::cout << "   > > >   " << fnInMic << " : " 
    //               << "   > > >   G I N I   Coef : " 
    //               << giniV << "   > > >   In the middle!!" << std::endl;
    // }else{
    //     std::cout << "   > > >   " << fnInMic << " : " 
    //               << "   > > >   G I N I   Coef : " 
    //               << giniV << "   > > >   To Process!!" << std::endl;
    // }
    Image<double> imMask(binMask);
    imMask.write(fnRoot+"_mask.mrc");
    

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
    if (goodCoords.size()>0)
    {
        mic.coords = goodCoords;
        mic.write_coordinates(0,0,(String)"particles@"+fnRoot+"_good.pos");
    }
    if (badCoords.size()>0)
    {
        mic.coords = badCoords;
        mic.write_coordinates(0,0,(String)"particles@"+fnRoot+"_bad.pos");
    }

    std::cout << "Discarted " << badCoords.size() << "/" << nPart << 
                 " particles by the noisy zone filter." << std::endl;
}
#undef DEBUG
