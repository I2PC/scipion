/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "threshold.h"

/* Read parameters --------------------------------------------------------- */
void ProgThreshold::readParams()
{
    XmippMetadataProgram::readParams();
    selectionMethod=getParam("--select");
    if (selectionMethod=="abs_below")
        iSelectionMethod=0;
    else if (selectionMethod=="below")
        iSelectionMethod=1;
    else if (selectionMethod=="above")
        iSelectionMethod=2;
    threshold=getDoubleParam("--select",1);
    substitutionMethod=getParam("--substitute");
    if (substitutionMethod=="value")
        newValue=getDoubleParam("--substitute",1);
    else if (substitutionMethod=="noise")
    {
        noiseAvg=getDoubleParam("--substitute",1);
        noiseStddev=getDoubleParam("--substitute",2);
    }
}

/* Usage ------------------------------------------------------------------- */
void ProgThreshold::defineParams()
{
    addUsageLine("Threshold volumes and images ");
    each_image_produces_an_output=true;
    XmippMetadataProgram::defineParams();
    addSeeAlsoLine("transform_mask, transform_morphology");
    addParamsLine("   --select <mode>                        : Select pixels meeting");
    addParamsLine("     where <mode>");
    addParamsLine("           abs_below <th>                 : Absolute value below a threshold");
    addParamsLine("           below <th>                     : Below a threshold");
    addParamsLine("           above <th>                     : Above a threshold");
    addParamsLine("   --substitute <substitutionMode=value>  : Substitute selected pixels by");
    addParamsLine("     where <substitutionMode>");
    addParamsLine("           binarize                       : Selected are set to 0, non-selected to 1");
    addParamsLine("           value <new=0>                  : New value");
    addParamsLine("           noise <avg=0> <stddev=1>       : Gaussian noise");
    addParamsLine("           avg                            : Average of non-selected");
    addExampleLine("Threshold a volume below a threshold",false);
    addExampleLine("xmipp_transform_threshold -i volume.vol -o volumeThresholded.vol --select below 0.01 --substitute value 0");
    addExampleLine("Generate a binary mask based on a threshold and apply it",false);
    addExampleLine("xmipp_transform_threshold -i volume.vol -o mask.vol --select below 0.5 --substitute binarize");
    addExampleLine("xmipp_transform_morphology -i mask.vol --dil");
    addExampleLine("xmipp_transform_mask -i volume.vol -o volumeMasked.vol --mask mask.vol");
}

/* Show ------------------------------------------------------------------- */
void ProgThreshold::show()
{
    if (verbose==0)
        return;
    XmippMetadataProgram::show();
    std::cout
    << "Selection method:    " << selectionMethod    << std::endl
    << "Threshold:           " << threshold          << std::endl
    << "Substitution method: " << substitutionMethod << std::endl;
    if (substitutionMethod=="value")
        std::cout << "New value:           " << newValue << std::endl;
    else if (substitutionMethod=="noise")
    {
        std::cout << "Noise average:       " << noiseAvg << std::endl;
        std::cout << "Noise std.dev.:      " << noiseStddev << std::endl;
    }
}

/* Process image ------------------------------------------------------------- */
void ProgThreshold::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> I;
    I.read(fnImg);
    MultidimArray<double> &mI=I();

    // Compute substitute value
    double substituteValue=0.0;
    if (substitutionMethod=="value")
        substituteValue=newValue;
    else if (substitutionMethod=="noise")
        substituteValue=rnd_gaus(noiseAvg,noiseStddev);
    else if (substitutionMethod=="avg")
    {
        double N=0;
        switch (iSelectionMethod)
        {
        case 0:
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
            {
            	double pixval=DIRECT_MULTIDIM_ELEM(mI,n);
                if (fabs(pixval)>threshold)
                {
                	substituteValue+=pixval;
                	++N;
                }
            }
            break;
        case 1:
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
            {
            	double pixval=DIRECT_MULTIDIM_ELEM(mI,n);
                if (pixval>threshold)
                {
                	substituteValue+=pixval;
                	++N;
                }
            }
            break;
        case 2:
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
            {
            	double pixval=DIRECT_MULTIDIM_ELEM(mI,n);
                if (pixval<threshold)
                {
                	substituteValue+=pixval;
                	++N;
                }
            }
            break;
        }
        substituteValue/=N;
    }

    // Apply threshold
    bool binarize=substitutionMethod=="binarize";
    switch (iSelectionMethod)
    {
    case 0:
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
        {
            if (fabs(DIRECT_MULTIDIM_ELEM(mI,n))<threshold)
            	DIRECT_MULTIDIM_ELEM(mI,n)=substituteValue;
            else if (binarize)
            	DIRECT_MULTIDIM_ELEM(mI,n)=1;
        }
        break;
    case 1:
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
        {
            if (DIRECT_MULTIDIM_ELEM(mI,n)<threshold)
            	DIRECT_MULTIDIM_ELEM(mI,n)=substituteValue;
            else if (binarize)
            	DIRECT_MULTIDIM_ELEM(mI,n)=1;
        }
        break;
    case 2:
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
        {
            if (DIRECT_MULTIDIM_ELEM(mI,n)>threshold)
            	DIRECT_MULTIDIM_ELEM(mI,n)=substituteValue;
            else if (binarize)
            	DIRECT_MULTIDIM_ELEM(mI,n)=1;
        }
        break;
    }

    // Write result
    I.write(fnImgOut);
}
