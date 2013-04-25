/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#include "transform_dimred.h"
#include <data/mask.h>
#include "diffusionMaps.h"

void ProgTransformDimRed::readParams()
{
	fnIn = getParam("-i");
	fnOut = getParam("-o");
	dimRefMethod = getParam("-m");
	outputDim  = getIntParam("-d");
}


// Show ====================================================================
void ProgTransformDimRed::show()
{
    if (verbose>0)
        std::cerr
        << "Input metadata file:    " 										<< fnIn          << std::endl
        << "Output metadata:        " 										<< fnOut         << std::endl
        << "Dim Red Method:     "     										<< dimRefMethod  << std::endl
        << "Number of dimensions after dimensionality reduction:     "     << dimRefMethod  << std::endl
        ;
}

// usage ===================================================================
void ProgTransformDimRed::defineParams()
{
    addParamsLine("   -i <metadatafile>             : metadata file  with images");
    addParamsLine("  [-o <metadatafile=\"\">]       : output metadata with distances between images");
    addParamsLine("   -m <dimRefMethod>             : Dimensionality Reduction method selected");
    addParamsLine("  [-d <N=2>]            			: Number of dimensions after the dimensionality reduction");
    addExampleLine("xmipp_transform_dimred -i images.xmd -i distances.xmd -m LTSA");
}

// Produce side info  ======================================================
//#define DEBUG
void ProgTransformDimRed::produceSideInfo()
{
    // Read input selfile and reference
    SFin.read(fnIn);
    if (SFin.size()==0)
        return;

    if (SFin.containsLabel(MDL_ENABLED))
        SFin.removeObjects(MDValueEQ(MDL_ENABLED, -1));

    MetaData SFaux;
    SFaux.removeDuplicates(SFin,MDL_IMAGE);
    SFin=SFaux;

    // Design Mask
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSize(SFin,Xdim,Ydim,Zdim,Ndim);

    mask.type = BINARY_CIRCULAR_MASK;
    mask.mode = INNER_MASK;
    mask.R1 = std::sqrt(0.25*Xdim*Xdim+0.25*Ydim*Ydim);
    mask.resize(Ydim,Xdim);
    mask.generate_mask();

    const MultidimArray<int> &mMask=mask.get_binary_mask();
    X.resizeNoCopy(SFin.size(),mMask.sum());

    // Copy images
    Image<double> img;
    size_t index = 0;
    FOR_ALL_OBJECTS_IN_METADATA(SFin)
    {
    	img.readApplyGeo(SFin,__iter.objId);
    	insertImageInDataMatrix(index++,img());
    }
}

void ProgTransformDimRed::insertImageInDataMatrix(size_t index, const MultidimArray<double> &mImg)
{
    size_t index2 = 0;
    const MultidimArray<int> &mMask=mask.get_binary_mask();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mImg)
		if (DIRECT_MULTIDIM_ELEM(mMask,n))
		{
			MAT_ELEM(X,index,index2)=DIRECT_MULTIDIM_ELEM(mImg,n);
			index2++;
		}
}

void ProgTransformDimRed::extractImageFromDataMatrix(size_t index, MultidimArray<double> &mImg)
{
    size_t index2 = 0;
    const MultidimArray<int> &mMask=mask.get_binary_mask();
    mImg.initZeros(mMask);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mImg)
		if (DIRECT_MULTIDIM_ELEM(mMask,n))
		{
			DIRECT_MULTIDIM_ELEM(mImg,n)=MAT_ELEM(X,index,index2);
			index2++;
		}
}

double ProgTransformDimRed::progCorrelationDistance(size_t i1, size_t i2)
{
	extractImageFromDataMatrix(i1,I1);
	extractImageFromDataMatrix(i2,I2);
	double corr=alignImages(I1,I2,M,true,aux,aux2,aux3);
	std::cout << "Correlation " << i1 << " " << i2 << " -> " << corr << std::endl;
	return 1.-corr;
}

ProgTransformDimRed *prog;
double correlationDistance(const Matrix2D<double> &X, size_t i1, size_t i2)
{
	return prog->progCorrelationDistance(i1,i2);
}

// Run  ====================================================================
void ProgTransformDimRed::run()
{
	prog=this;
	show();
    produceSideInfo();

    //double dim1=intrinsicDimensionality(X, "MLE", false, &correlationDistance);
    //std::cout << "dim1=" << dim1 << std::endl;
    //double dim2=intrinsicDimensionality(X, "CorrDim", false, &correlationDistance);
    //std::cout << "dim2=" << dim2 << std::endl;

    DiffusionMaps dimred;
    dimred.setInputData(X);
    dimred.setSpecificParameters();
    dimred.setOutputDimensionality(outputDim);
    dimred.distance=&correlationDistance;
    dimred.reduceDimensionality();

    MetaData SFout;
    SFout=SFin;
    std::vector<double> dimredProj;
    dimredProj.resize(outputDim);
    int i=0;
    const Matrix2D<double> &Y=dimred.getReducedData();
    FOR_ALL_OBJECTS_IN_METADATA(SFout)
    {
    	memcpy(&dimredProj[0],&MAT_ELEM(Y,i,0),outputDim*sizeof(double));
    	SFout.setValue(MDL_DIMRED,dimredProj,__iter.objId);
    }
    if (fnOut=="")
    	fnOut=fnIn;
    SFout.write(fnOut);
}
