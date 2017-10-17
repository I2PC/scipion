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
    noDenoising = checkParam("--noDenoising");
    useLBP = checkParam("--lbp");
    useEntropy = checkParam("--entropy");
    useVariance = checkParam("--variance");
    useZernike = checkParam("--zernike");
}

// Show ====================================================================
void ProgExtractFeatures::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:             " << fnSel        << std::endl
    << "Output selfile:            " << fnOut        << std::endl
    << "Turn off denoising:        " << noDenoising  << std::endl
    << "Extract LBP features:      " << useLBP       << std::endl
    << "Extract entropy features:  " << useEntropy   << std::endl
    << "Extract variance features: " << useVariance  << std::endl
    << "Extract Zernike moments:   " << useZernike   << std::endl
    ;
}

// Usage ===================================================================
void ProgExtractFeatures::defineParams()
{
    addUsageLine("Clusters a set of images");
    addParamsLine("  -i <selfile>                  : Selfile containing images to be clustered");
    addParamsLine("  [-o <selfile=\"\">]           : Output selfile");
    addParamsLine("  [--noDenoising]               : Turn off denoising");
    addParamsLine("  [--lbp]                       : Extract LBP features");
    addParamsLine("  [--entropy]                   : Extract entropy features");
    addParamsLine("  [--variance]                  : Extract variance features");
    addParamsLine("  [--zernike]                   : Extract Zernike moments");
}


int ProgExtractFeatures::facs(int n)
{
    return (n == 1 || n == 0) ? 1 :
           (n == 2) ? 2 :
           (n == 3) ? 6 : 24;
}

void ProgExtractFeatures::extractLBP(const MultidimArray<double> &I,
                                     std::vector<double> &fv)
{
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
}


void ProgExtractFeatures::extractEntropy(const MultidimArray<double> &I,
                                         MultidimArray<double> &Imasked,
                                         std::vector<double> &fv)
{
    int hist[256] = {};
    int val;

    // entropy of entire image
    double m, M;
    I.computeDoubleMinMax(m,M);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I)
    {
        val = floor(((DIRECT_MULTIDIM_ELEM(I, n) - m) * 255.0) / (M - m));
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
}

void ProgExtractFeatures::extractVariance(const MultidimArray<double> &I,
                                          std::vector<double> &fv)
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
                fv.push_back(var_i / count);
            }
            else
            {
                var_o_sum += var_o / count;
                fv.push_back(var_o / count);
            }
        }
    }

    fv.push_back((var_i_sum / 4) / (var_o_sum / 12));
}


void ProgExtractFeatures::extractZernike(const MultidimArray<double> &I,
                                         std::vector<double> &fv)
{
    MultidimArray<double> R, Theta, Rad;
    R.resize(I); R.setXmippOrigin();
    Theta.resize(I); Theta.setXmippOrigin();
    Rad.resize(I); Rad.setXmippOrigin();

    double c;
    int Sy = YSIZE(I);
    int Sx = XSIZE(I);
    const std::complex<double> i(0.0, 1.0);

    for (int y = 0; y < Sy; y++)
    {
        int r2 = 2*(y+1)-Sy-1;

        for (int x = 0; x < Sx; x++)
        {
            int r1 = 2*(x+1)-Sy-1;
            DIRECT_A2D_ELEM(R,y,x) = sqrt(r1*r1 + r2*r2) / Sy;
            if (DIRECT_A2D_ELEM(R,y,x) > 1) DIRECT_A2D_ELEM(R,y,x) = 0;

            DIRECT_A2D_ELEM(Theta,y,x) = atan2((Sy+1-2*(y+1)), (2*(x+1)-Sy-1));
        }
    }

    for (int n = 1; n < 5; n++)
    {
        for (int m = -n; m < 0; m+=2)
        {
            int mn = (n - abs(m)) / 2;
            int nm = (n + abs(m)) / 2;
            std::complex<double> product = 0.0;

            for (int y = 0; y < Sy; y++)
            {
                for (int x = 0; x < Sx; x++)
                {
                    DIRECT_A2D_ELEM(Rad,y,x) = 0;
                    for (int s = 0; s <= mn; s++)
                    {
                        int ns = n - 2*s;
                        c = ((s%2 == 0) ? 1 : -1) *
                            facs(n-s) / (facs(s) *
                            facs(nm-s) *
                            facs(mn-s));

                        DIRECT_A2D_ELEM(Rad,y,x) += c *
                            pow(DIRECT_A2D_ELEM(R,y,x), ns);
                    }

                    product += DIRECT_A2D_ELEM(I,y,x) *
                               DIRECT_A2D_ELEM(Rad,y,x) *
                               exp(-1.0 * i * (double)m *
                                   DIRECT_A2D_ELEM(Theta,y,x));
                }
            }
            fv.push_back(std::abs(product));
        }
    }
}


void ProgExtractFeatures::run()
{
    MetaData SF;
    SF.read(fnSel);
    Image<double> I, Imasked;
    FileName fnImg;
    MDRow row;
	CorrelationAux aux;
	std::vector<double> fv;

	FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
    	SF.getValue(MDL_IMAGE, fnImg, __iter.objId);
    	I.read(fnImg);
    	I().setXmippOrigin();
    	centerImageTranslationally(I(), aux);

    	if (!noDenoising)
    	    denoiseTVFilter(I(), 50);

        if (useLBP)
        {
            extractLBP(I(), fv);
            SF.setValue(MDL_SCORE_BY_LBP, fv, __iter.objId);
            fv.clear();
        }

        if (useEntropy)
        {
            extractEntropy(I(), Imasked(), fv);
            SF.setValue(MDL_SCORE_BY_ENTROPY, fv, __iter.objId);
            fv.clear();
        }

        if (useVariance)
        {
            extractVariance(I(), fv);
            SF.setValue(MDL_SCORE_BY_VARIANCE, fv, __iter.objId);
            fv.clear();
        }

        if (useZernike)
        {
            extractZernike(I(), fv);
            SF.setValue(MDL_SCORE_BY_ZERNIKE, fv, __iter.objId);
            fv.clear();
        }
    }

	if (fnOut == "") fnOut = fnSel;

	SF.write(fnOut, MD_APPEND);
}