/***************************************************************************
 * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
 *
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

#include "eliminate_nonparticles.h"

#include <data/xmipp_funcs.h>
#include <data/mask.h>
#include <data/filters.h>
#include <vector>
#include <numeric>
#include <set>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>


// Read arguments ==========================================================
void ProgEliminateNonParticles::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    threshold = getIntParam("-t");
}

// Show ====================================================================
void ProgEliminateNonParticles::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:           " << fnIn       << std::endl
    << "Output selfile:          " << fnOut      << std::endl
    << "Threshold:               " << threshold  << std::endl
    ;
}

// Usage ===================================================================
void ProgEliminateNonParticles::defineParams()
{
    addUsageLine("Eliminates samples containing only noise");
    addParamsLine("  -i <selfile>           : Selfile containing set of input particles");
    addParamsLine("  -o <selfile>           : Output selfile");
    addParamsLine("  -t <int>               : Threshold used by algorithm");
}

void ProgEliminateNonParticles::run()
{
    // Read the input metadata
    SF.read(fnSel);
    FileName fnImg, fnClass, fnTemp;
    int itemId = 0, countItems = 0;
    MDRow row;
    MetaData MDsummary, MDclass;
    std::vector<double> fv;
    std::vector<Point> points;
    std::vector<Cluster> clusters;
    srand (time(NULL));

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        countItems++;
    	SF.getValue(MDL_IMAGE, fnImg,__iter.objId);
    	Iref.read(fnImg);
    	Iref().setXmippOrigin();
    	CorrelationAux aux;
    	centerImageTranslationally(Iref(), aux);

        if (isParticle(__iter.objId))
        {
            itemId++;
            fv = feature_extraction();
            Point p(countItems, fv);
            points.push_back(p);
        }
    }
    std::vector<double> vars;
    double mean_all = 0.0;
    for (int yy = 1; yy <= 4; yy++)
    {
        int y_max = YSIZE(Iref()) / 4 * yy;
        int y_min = YSIZE(Iref()) / 4 * (yy-1);
        for (int xx = 1; xx <= 4; xx++)
        {
            int x_max = XSIZE(Iref()) / 4 * xx;
            int x_min = XSIZE(Iref()) / 4 * (xx-1);

            double mean = 0.0;
            double var = 0.0;
            int count = 0;
            for (int y = y_min; y < y_max; y++)
            {
                for (int x = x_min; x < x_max; x++)
                {
                    mean += DIRECT_A2D_ELEM(Iref(),y,x);
                    count++;
                }
            }
            mean = mean / count;
            for (int y = y_min; y < y_max; y++)
            {
                for (int x = x_min; x < x_max; x++)
                {
                    var += (DIRECT_A2D_ELEM(Iref(),y,x) - mean) * (DIRECT_A2D_ELEM(Iref(),y,x) - mean);
                }
            }
            vars.push_back(var / count);
            mean_all += mean;
        }
    }

    //double var_all = 0.0;
    //FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iref())
    //    var_all += (DIRECT_A2D_ELEM(Iref(),i,j) - mean_all) * (DIRECT_A2D_ELEM(Iref(),i,j) - mean_all);
    //var_all = var_all / (XSIZE(Iref()) * YSIZE(Iref()));

    //double average = std::accumulate(vars.begin(), vars.end(), 0.0) / vars.size();
    //double variance = 0.0;
    //for (int i = 0; i < vars.size(); i++)
	//    variance += (vars[i] - average) * (vars[i] - average);

    std::cout << variance / vars.size() << std::endl;
    if (variance / vars.size() < threshold) return true;
    else return false;
}