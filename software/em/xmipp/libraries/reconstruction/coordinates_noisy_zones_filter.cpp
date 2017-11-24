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
    addParamsLine(" [--patchSize <n=10>] : Patch size for the variance filter");
    addParamsLine(" [-o <coordinates=pos>] : Output coordinates, if not passed the input is sustituted");
}

void ProgCoordinatesNoisyZonesFilter::readParams()
{
	fnInCoord = getParam("--pos");
	fnInMic = getParam("--mic");
    fnOut = getParam("-o");
    if (fnOut=="pos")
        fnOut = fnInCoord;
	patchSize = getIntParam("--patchSize");
}

void ProgCoordinatesNoisyZonesFilter::show()
{
    if (verbose)
		std::cout
		<< "Input Coordinates:     " << fnInCoord  << std::endl
		<< "Input Micrograph:      " << fnInMic  << std::endl
		<< "Output coordinates:    " << fnOut << std::endl
		<< "Patch size:            " << patchSize << std::endl
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
    double giniV = giniCoeff(matrixMic, patchSize);

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
    if (verbose)
    {
        Image<double> imVar(matrixMic);
        imVar.write(fnInCoord.withoutExtension()+"_varianceFilter.mrc");
    }
    

    std::vector<Particle_coords> goodCoords;
    std::vector<Particle_coords> badCoords;
    std::vector<Particle_coords> &allCoords=mic.coords;
    std::vector<double> varValue;
    for (size_t i=0; i<allCoords.size(); ++i)
    {
    	Particle_coords &coord=allCoords[i];
    	// if (DIRECT_A2D_ELEM(binMask,coord.Y,coord.X))
    	// 	goodCoords.push_back(coord);
    	// else
    	// 	badCoords.push_back(coord);
        varValue.push_back(DIRECT_A2D_ELEM(matrixMic,coord.Y,coord.X));
    }
    // if (goodCoords.size()>0)
    // {
    //     mic.coords = goodCoords;
    //     mic.write_coordinates(0,0,(String)"particles@"+fnRoot+"_good.pos");
    // }
    // if (badCoords.size()>0)
    // {
    //     mic.coords = badCoords;
    //     mic.write_coordinates(0,0,(String)"particles@"+fnRoot+"_bad.pos");
    // }
    // std::cout << "Discarted " << badCoords.size() << "/" << nPart << 
    //              " particles by the noisy zone filter." << std::endl;

    std::cout << " >>> Gini Coeff: " << giniV << std::endl;

    MetaData MD;

    int imax = allCoords.size();
    size_t id;
    for (int i = 0; i < imax; i++)
    {
        if (allCoords[i].valid && allCoords[i].cost > 0
                && allCoords[i].label == 0)
        {
            id = MD.addObject();
            MD.setValue(MDL_XCOOR, allCoords[i].X, id);
            MD.setValue(MDL_YCOOR, allCoords[i].Y, id);
            MD.setValue(MDL_SCORE_BY_VAR, varValue[i], id);
            MD.setValue(MDL_SCORE_BY_GINI, giniV, id);            
        }
    }
    
    MD.write((String)"particles@"+fnOut);


}
#undef DEBUG
