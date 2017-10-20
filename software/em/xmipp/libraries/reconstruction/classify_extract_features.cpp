/***************************************************************************
 *
 * Authors:    Tomas Majtner            tmajtner@cnb.csic.es (2017)
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

#include "classify_extract_features.h"

#include <data/xmipp_funcs.h>
#include <data/mask.h>
#include <data/filters.h>
#include <vector>
#include <string>


// Read arguments ==========================================================
void ProgExtractFeatures::readParams()
{
    fnSel = getParam("-i");
    fnOut = getParam("-o");
    useLBP = checkParam("--lbp");
    useEntropy = checkParam("--entropy");
}

// Show ====================================================================
void ProgExtractFeatures::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:             " << fnSel      << std::endl
    << "Output selfile:            " << fnOut      << std::endl
    << "Extract LBP features:      " << useLBP     << std::endl
    << "Extract entropy features:  " << useEntropy << std::endl
    ;
}

// Usage ===================================================================
void ProgExtractFeatures::defineParams()
{
    addUsageLine("Clusters a set of images");
    addParamsLine("  -i <selfile>                  : Selfile containing images to be clustered");
    addParamsLine("  [-o <selfile=\"\">]           : Output selfile");
    addParamsLine("  [--lbp]                       : Extract LBP features");
    addParamsLine("  [--entropy]                   : Extract entropy features");
}


std::vector<double> ProgExtractFeatures::extractLBP(const MultidimArray<double> &I)
{
    std::vector<double> fv;
    std::vector<double> min_idxs, min_idxs_sort;

    unsigned char code;
    double center;
    int lbp_hist[256] = {};

    for (int i = 0; i < 256; i++)
    {
        code = i;
        int code_min = (int) code;
        for (int ii = 0; ii < 7; ii++)
        {
            unsigned char c = code & 1;
            code >>= 1;
            code |= (c << 7);
            if ((int) code < code_min)
                code_min = (int) code;
        }
        min_idxs.push_back(code_min);
    }
    min_idxs_sort = min_idxs;
    std::sort(min_idxs_sort.begin(), min_idxs_sort.end());
    std::unique(min_idxs_sort.begin(), min_idxs_sort.end());

    for (int y = 1; y < (YSIZE(I)-1); y++)
    {
        for (int x = 1; x < (XSIZE(I)-1); x++)
        {
            code = 0;
            center = DIRECT_A2D_ELEM(I,y,x);
            code |= (DIRECT_A2D_ELEM(I,y-1,x-1) > center) << 7;
            code |= (DIRECT_A2D_ELEM(I,y-1,x  ) > center) << 6;
            code |= (DIRECT_A2D_ELEM(I,y-1,x+1) > center) << 5;
            code |= (DIRECT_A2D_ELEM(I,y,  x+1) > center) << 4;
            code |= (DIRECT_A2D_ELEM(I,y+1,x+1) > center) << 3;
            code |= (DIRECT_A2D_ELEM(I,y+1,x  ) > center) << 2;
            code |= (DIRECT_A2D_ELEM(I,y+1,x-1) > center) << 1;
            code |= (DIRECT_A2D_ELEM(I,y  ,x-1) > center) << 0;
            int idx = min_idxs[(int) code];
            lbp_hist[idx]++;
        }
    }

    for (int i = 0; i < 36; i++)
    {
        int idx = min_idxs_sort[i];
        fv.push_back(lbp_hist[idx]);
    }
    return fv;
}


std::vector<double> ProgExtractFeatures::extractEntropy(const MultidimArray<double> &I, MultidimArray<double> &Imasked)
{
    std::vector<double> fv;
    int hist[256] = {};
    int val;

    // entropy of entire image
    double m, M;
    I.computeDoubleMinMax(m,M);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I)
    {
        val = floor(((DIRECT_MULTIDIM_ELEM(I, n) - m) * 255.0) / (M-m));
        hist[val]++;
    }

    double entropy = 0;
    for (int i = 0; i < 256; i++)
        entropy += std::max(hist[i], 1) * log2((std::max(hist[i], 1)));

    fv.push_back(-1*entropy);

    // preparing masks for feature extraction, will be performed only once
    if (XSIZE(masks[0]) < 1)
    {
        for (int i = 0; i < (sizeof(masks)/sizeof(*masks)); i++)
        {
            masks[i].resize(I);
            masks[i].setXmippOrigin();
        }

        int wave_size = XSIZE(I) / 2;
        int wave_size_step = XSIZE(I) / 32;

        // simple inner-circle regions without negative values
//        for (int i=2; i<(sizeof(masks)/sizeof(*masks)); i++)
//        {
//            BinaryCircularMask(masks[0], wave_size);
//            BinaryCircularMask(masks[1], wave_size - wave_size_step);
//
//            // if we want overlapping regions (rings)
//            // BinaryCircularMask(masks[1], wave_size - 2*wave_size_step);
//
//            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(masks[0])
//                DIRECT_MULTIDIM_ELEM(masks[i],n) = DIRECT_MULTIDIM_ELEM(masks[0],n) -
//                                                   DIRECT_MULTIDIM_ELEM(masks[1],n);
//
//            wave_size -= wave_size_step;
//        }

        for (int i = 2; i < (sizeof(masks)/sizeof(*masks)); i++)
        {
            BinaryCircularMask(masks[0], wave_size);
            BinaryCircularMask(masks[1], wave_size - wave_size_step);

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(masks[0])
                DIRECT_MULTIDIM_ELEM(masks[i],n) = 2*DIRECT_MULTIDIM_ELEM(masks[1],n) -
                                                   DIRECT_MULTIDIM_ELEM(masks[0],n);

            BinaryCircularMask(masks[1], wave_size - 2*wave_size_step);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(masks[0])
                DIRECT_MULTIDIM_ELEM(masks[i],n) = DIRECT_MULTIDIM_ELEM(masks[i],n) -
                                                   DIRECT_MULTIDIM_ELEM(masks[1],n);


            wave_size -= wave_size_step;     // overlapping regions (rings)
            //wave_size -= 2*wave_size_step;   // non-overlapping regions (rings)
        }
    }

    // extracting entropy from unmasked regions
    for (int i = 2; i < (sizeof(masks)/sizeof(*masks)); i++)
    {
        apply_binary_mask(masks[i], I, Imasked);
        int hist[256] = {};
        int idx;
        Imasked.computeDoubleMinMax(m,M);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imasked)
        {
            idx = floor(((DIRECT_MULTIDIM_ELEM(Imasked, n) - m) * 255.0) / (M-m));
            hist[idx]++;
        }

        double entropy = 0;
        for (int i = 0; i < 256; i++)
            entropy += std::max(hist[i], 1) * log2((std::max(hist[i], 1)));

        fv.push_back(-1*entropy);
    }

    return fv;
}

void ProgExtractFeatures::run()
{
    MetaData SF;
    SF.read(fnSel);
    Image<double> I, Imasked;
    FileName fnImg;
    MDRow row;
	CorrelationAux aux;

	FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
    	SF.getValue(MDL_IMAGE, fnImg, __iter.objId);
    	I.read(fnImg);
    	I().setXmippOrigin();
    	centerImageTranslationally(I(), aux);
    	denoiseTVFilter(I(), 50);

        if (useLBP)
            SF.setValue(MDL_SCORE_BY_LBP, extractLBP(I()), __iter.objId);

        if (useEntropy)
            SF.setValue(MDL_SCORE_BY_ENTROPY, extractEntropy(I(),Imasked()), __iter.objId);

    }
	if (fnOut=="")
		fnOut=fnSel;
	SF.write(fnOut,MD_APPEND);
}
