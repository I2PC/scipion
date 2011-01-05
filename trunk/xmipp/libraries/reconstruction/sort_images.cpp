/***************************************************************************
 *
 * Authors:    Carlos Oscar           coss@cnb.csic.es (2010)
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

#include "sort_images.h"
#include <data/filters.h>
#include <data/mask.h>
#include <data/image_collection.h>

// Read arguments ==========================================================
void ProgSortImages::readParams()
{
    fnSel  = getParam("-i");
    fnRoot = getParam("--oroot");
}

// Show ====================================================================
void ProgSortImages::show()
{
    std::cerr << "Input selfile:    " << fnSel           << std::endl
              << "Output rootname:  " << fnRoot          << std::endl
    ;
}

// usage ===================================================================
void ProgSortImages::defineParams()
{
    addUsageLine("Sort a set of images by local similarity");
    addParamsLine("   -i <selfile>        : selfile of images");
    addParamsLine("   --oroot <rootname>  : output rootname");
}

// Produce side info  ======================================================
//#define DEBUG
void ProgSortImages::produceSideInfo()
{
    fnStack=fnRoot+".stk";
    if (exists(fnStack))
    	unlink(fnStack.c_str());

    // Read input selfile and reference
    ImageCollection SF;
    SF.read(fnSel);
    int idx=0;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        if (idx==0)
        {
            SF.getValue(MDL_IMAGE,fnImg);
            lastImage.read(fnImg);
            centerImage(lastImage());
            lastImage.write(fnStack,idx,true,WRITE_APPEND);
            SFout.addObject();
            FileName fnImageStack;
            fnImageStack.compose(idx,fnStack);
            SFout.setValue(MDL_IMAGE,fnImageStack);
            SFout.setValue(MDL_IMAGE_ORIGINAL,fnImg);
        }
        else
        {
            SF.getValue(MDL_IMAGE,fnImg);
            toClassify.push_back(fnImg);
        }
        idx++;
    }

    // Prepare mask
    mask.resize(lastImage());
    mask.setXmippOrigin();
    BinaryCircularMask(mask,XSIZE(lastImage())/2, INNER_MASK);
}

// Choose next image =======================================================
void ProgSortImages::chooseNextImage()
{
    int imax=toClassify.size();
    Image<double> bestImage, Iaux;
    Matrix2D<double> M;
    double bestCorr=-1;
    int bestIdx=-1;
    for (int i=0; i<imax; i++)
    {
        Iaux.read(toClassify[i]);
        Iaux().setXmippOrigin();
        
        // Choose between this image and its mirror
        MultidimArray<double> I, Imirror;
        I=Iaux();
        Imirror=I;
        Imirror.selfReverseX();
        Imirror.setXmippOrigin();
        
        alignImages(lastImage(),I,M);
        alignImages(lastImage(),Imirror,M);
        double corr=correlation_index(lastImage(),I,&mask);
        double corrMirror=correlation_index(lastImage(),Imirror,&mask);
        if (corr>bestCorr)
        {
            bestCorr=corr;
            bestImage()=I;
            bestIdx=i;
        }
        if (corrMirror>bestCorr)
        {
            bestCorr=corrMirror;
            bestImage()=Imirror;
            bestIdx=i;
        }
    }
    
    int idxStack=SFout.size();
    FileName fnImageStack;
    fnImageStack.compose(idxStack,fnStack);
    bestImage.write(fnStack,idxStack,true,WRITE_APPEND);
    lastImage=bestImage;
    SFout.addObject();
    SFout.setValue(MDL_IMAGE,fnImageStack);
    SFout.setValue(MDL_IMAGE_ORIGINAL,toClassify[bestIdx]);
    toClassify.erase(toClassify.begin()+bestIdx);
}

// Run  ====================================================================
void ProgSortImages::run()
{
	show();
	produceSideInfo();

	std::cout << "Images to go: ";
    while (toClassify.size()>0) 
    {
        std::cout << toClassify.size() << " ";
        std::cout.flush();
        chooseNextImage();
    }
    std::cout << toClassify.size() << std::endl;
    SFout.write(fnRoot+".sel");
}
