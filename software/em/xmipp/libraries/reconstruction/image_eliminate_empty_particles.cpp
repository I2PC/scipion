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

#include "image_eliminate_empty_particles.h"
#include "classify_extract_features.h"
#include <data/filters.h>
#include <fstream>

// Read arguments ==========================================================
void ProgEliminateEmptyParticles::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    fnElim = getParam("-e");
    threshold = getDoubleParam("-t");
    addFeatures = checkParam("--addFeatures");
    useDenoising = checkParam("--useDenoising");
    denoise = getIntParam("-d");
}

// Show ====================================================================
void ProgEliminateEmptyParticles::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:           " << fnIn         << std::endl
    << "Output selfile:          " << fnOut        << std::endl
    << "Eliminated selfile:      " << fnElim       << std::endl
    << "Threshold:               " << threshold    << std::endl
	<< "Add features:            " << addFeatures  << std::endl
	<< "Turn on denoising:       " << useDenoising << std::endl
	<< "Denosing parameter:      " << denoise      << std::endl
    ;
}

// Usage ===================================================================
void ProgEliminateEmptyParticles::defineParams()
{
    addUsageLine("Eliminates empty particles (false positives from picking)");
    addParamsLine("  -i <selfile>                       : Selfile containing set of input particles");
    addParamsLine("  [-o <selfile=\"output.xmd\">]      : Output selfile");
    addParamsLine("  [-e <selfile=\"eliminated.xmd\">]  : Eliminated particles selfile");
    addParamsLine("  [-t <float=-1>]                    : Threshold used by algorithm. Set to -1 for no elimination.");
    addParamsLine("  [--addFeatures]                    : Add the emptiness features to the input particles");
    addParamsLine("  [--useDenoising]                   : Option for turning on denoising method while computing emptiness feature");
    addParamsLine("  [-d <int=50>]                      : Parameter for denoising, higher value means stronger denoising and slower computation");
}

void ProgEliminateEmptyParticles::run()
{
    MetaData SF;
    SF.read(fnIn);
    FileName fnImg;
    Image<double> I;
    CorrelationAux aux;
    MDRow row;
    MetaData MDclass, MDclassEl, MDclassT, MDclassElT;
    ProgExtractFeatures ef;
    int countItems = 0;
    std::vector<double> fv;

    init_progress_bar(SF.size());
    std::size_t extraPath = fnOut.find_last_of("/");
    // these files are for streaming and will be later removed in Scipion
    FileName fnDone = fnOut.substr(0, extraPath+1) + "outTemp.xmd";
    FileName fnElDone = fnOut.substr(0, extraPath+1) + "elimTemp.xmd";

    std::ifstream ifile1(fnOut.c_str());
    if (ifile1) MDclass.read(fnOut);
    std::ifstream ifile2(fnElim.c_str());
    if (ifile2) MDclassEl.read(fnElim);

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        countItems++;
    	SF.getValue(MDL_IMAGE, fnImg, __iter.objId);
    	SF.getRow(row, countItems);

    	I.read(fnImg);
    	I().setXmippOrigin();
    	centerImageTranslationally(I(), aux);

        if (useDenoising)
    	    denoiseTVFilter(I(), denoise);

        ef.extractVariance(I(), fv);

        double ratio = fv.back();
        row.setValue(MDL_SCORE_BY_EMPTINESS, ratio);
        if (addFeatures)
        	row.setValue(MDL_SCORE_BY_VARIANCE, fv);
        if (threshold<0 || ratio > threshold)
        {
            MDclass.addRow(row);
            MDclassT.addRow(row);
        }
        else
        {
            MDclassEl.addRow(row);
            MDclassElT.addRow(row);
        }

        fv.clear();
        if (countItems%100==0)
        	progress_bar(countItems);
    }
    progress_bar(SF.size());

    if (MDclass.size()>0)
    	MDclass.write(formatString("@%s", fnOut.c_str()), MD_APPEND);
    if (MDclassEl.size()>0)
    	MDclassEl.write(formatString("@%s", fnElim.c_str()), MD_APPEND);
    if (MDclassT.size()>0)
    	MDclassT.write(formatString("@%s", fnDone.c_str()), MD_APPEND);
    if (MDclassElT.size()>0)
    	MDclassElT.write(formatString("@%s", fnElDone.c_str()), MD_APPEND);
}
