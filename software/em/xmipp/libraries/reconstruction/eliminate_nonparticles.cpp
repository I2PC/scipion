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
#include <data/filters.h>


// Read arguments ==========================================================
void ProgEliminateNonParticles::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    fnElim = getParam("-e");
    threshold = getDoubleParam("-t");
}

// Show ====================================================================
void ProgEliminateNonParticles::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:           " << fnIn       << std::endl
    << "Output selfile:          " << fnOut      << std::endl
    << "Eliminated selfile:      " << fnElim     << std::endl
    << "Threshold:               " << threshold  << std::endl
    ;
}

// Usage ===================================================================
void ProgEliminateNonParticles::defineParams()
{
    addUsageLine("Eliminates samples containing only noise");
    addParamsLine("  -i <selfile>                       : Selfile containing set of input particles");
    addParamsLine("  [-o <selfile=\"output.xmd\">]      : Output selfile");
    addParamsLine("  [-e <selfile=\"eliminated.xmd\">]  : Eliminated particles selfile");
    addParamsLine("  -t <float>                         : Threshold used by algorithm");
}


bool ProgEliminateNonParticles::isParticle(const MultidimArray<double> &I)
{
    double var_i_sum = 0.0;
    double var_o_sum = 0.0;
    for (int yy = 1; yy <= 4; yy++)
    {
        int y_max = YSIZE(I) / 4 * yy;
        int y_min = YSIZE(I) / 4 * (yy-1);
        for (int xx = 1; xx <= 4; xx++)
        {
            int x_max = XSIZE(I) / 4 * xx;
            int x_min = XSIZE(I) / 4 * (xx-1);

            double mean = 0.0;
            int count = 0;
            double var_i = 0.0;
            double var_o = 0.0;

            for (int y = y_min; y < y_max; y++)
            {
                for (int x = x_min; x < x_max; x++)
                {
                    mean += DIRECT_A2D_ELEM(I,y,x);
                    count++;
                }
            }
            mean = mean / count;

            for (int y = y_min; y < y_max; y++)
            {
                for (int x = x_min; x < x_max; x++)
                {
                    if (yy > 1 && yy < 4 && xx > 1 && xx < 4)
                        var_i += (DIRECT_A2D_ELEM(I,y,x) - mean) *
                                 (DIRECT_A2D_ELEM(I,y,x) - mean);
                    else
                        var_o += (DIRECT_A2D_ELEM(I,y,x) - mean) *
                                 (DIRECT_A2D_ELEM(I,y,x) - mean);
                }
            }

            if (yy > 1 && yy < 4 && xx > 1 && xx < 4)
            {
                var_i_sum += var_i / count;
                features.push_back(var_i / count);
            }
            else
            {
                var_o_sum += var_o / count;
                features.push_back(var_o / count);
            }
        }
    }

    features.push_back((var_i_sum / 4) / (var_o_sum / 12));
    if ((var_i_sum / 4) / (var_o_sum / 12) > threshold) return true;
    else return false;
}

void ProgEliminateNonParticles::run()
{
    MetaData SF;
    SF.read(fnIn);
    FileName fnImg;
    Image<double> I;
    int countItems = 0;
    CorrelationAux aux;
    MDRow row;
    MetaData MDclass, MDclass_elim;

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        countItems++;
    	SF.getValue(MDL_IMAGE, fnImg, __iter.objId);
    	SF.getRow(row, countItems);

    	I.read(fnImg);
    	I().setXmippOrigin();
    	centerImageTranslationally(I(), aux);
    	denoiseTVFilter(I(), 50);

        if (isParticle(I()))
        {
            size_t recId = MDclass.addRow(row);
            MDclass.setValue(MDL_SCORE_BY_VARIANCE, features, recId);
            MDclass.write(formatString("@%s", fnOut.c_str()), MD_APPEND);
        }
        else
        {
            size_t recId = MDclass_elim.addRow(row);
            MDclass_elim.setValue(MDL_SCORE_BY_VARIANCE, features, recId);
            MDclass_elim.write(formatString("@%s", fnElim.c_str()), MD_APPEND);
        }
        features.clear();
    }
}