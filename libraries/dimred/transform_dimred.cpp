/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
 *              Carlos Oscar Sorzano (coss@cnb.csic.es)
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
#include <data/matrix1d.h>
#include <data/matrix2d.h>
#include <time.h>
#include "diffusionMaps.h"

void ProgTransformDimRed::readParams()
{
	ProgDimRed::readParams();
    if (checkParam("--randomSample"))
    {
    	fnRandomSampling  = getParam("--randomSample",0);
    	numGrids  =  getIntParam("--randomSample",1);
    }
    distance=getParam("--distance");
}

// Show ====================================================================
void ProgTransformDimRed::show()
{
    if (verbose>0)
    {
    	ProgDimRed::show();
    	std::cout << "Distance:      " << distance << std::endl;
    	if (fnRandomSampling!="")
    		std::cout << "Random sample output:     " << fnRandomSampling << std::endl
    				  << "Number of sampling grids: " << numGrids         << std::endl;
    }
}

// usage ===================================================================
void ProgTransformDimRed::defineParams()
{
	addUsageLine("This program takes an input metadata and projects each image onto a lower dimensional space using the selected method");
	setDefaultComment("-i","Input metadata");
	setDefaultComment("-o","Output metadata");
	ProgDimRed::defineParams();
    addParamsLine("  [--distance <d=Correlation>]    : Distance between images");
    addParamsLine("    where <d>");
    addParamsLine("         Euclidean:    Euclidean distance between images, no alignment");
    addParamsLine("         Correlation:  Correlation between images after alignment");
	addParamsLine("  [--randomSample <file> <num=3>] : Generates a random sample of the reduced map with num grids in each direction");
    addExampleLine("xmipp_transform_dimred -i images.xmd --randomSample randomSample.xmd");
}

// Produce side info  ======================================================
//#define DEBUG
double ProgTransformDimRed::progCorrelationDistance(size_t i1, size_t i2)
{
    extractImageFromDataMatrix(i1,I1);
    extractImageFromDataMatrix(i2,I2);
    double corr=alignImages(I1,I2,M,true,aux,aux2,aux3);
    return 1.-corr;
}

ProgTransformDimRed *prog;
double correlationDistance(const Matrix2D<double> &X, size_t i1, size_t i2)
{
    return prog->progCorrelationDistance(i1,i2);
}

void ProgTransformDimRed::produceSideInfo()
{
	ProgDimRed::produceSideInfo();

	// Read input selfile
    SFin.read(fnIn);
    if (SFin.size()==0)
        return;
    SFin.removeDisabled();

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
	algorithm->setInputData(X);

    // Set distance
    prog=this;
    if (distance=="Correlation")
    	algorithm->distance=&correlationDistance;
    else
    	algorithm->distance=NULL;
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

void ProgTransformDimRed::extractRandomProjections()
{
    //number of grids in each axis
    int num = numGrids;

    //Metadata that has the dimrefcoeffs
    MetaData SFin;
    SFin.read(fnIn);
    SFin.removeDisabled();

    size_t numParticles = SFin.size();
    std::vector<double> dimredProj;
    std::vector< Matrix1D< double > > coor;
    Matrix1D<double> dummy;
    dummy.resizeNoCopy(numParticles);
    for (int n=0; n<outputDim; ++n)
    	coor.push_back(dummy);

    // Keep the reduced coefficients
    size_t i=0;
    FOR_ALL_OBJECTS_IN_METADATA(SFin)
    {
        SFin.getValue(MDL_DIMRED,dimredProj,__iter.objId);
        for (int n=0; n<outputDim; ++n)
        	VEC_ELEM(coor[n],i)=dimredProj[n];
        i++;
    }

    // We obtain the max and min coordinate values
    std::vector<int> minCoor, maxCoor;
    for (int n=0; n<outputDim; ++n)
    {
    	double minval, maxval;
    	coor[n].computeMinMax(minval, maxval);
    	minCoor.push_back(minval);
    	maxCoor.push_back(maxval);
    }

    Matrix2D <int> squares;
    squares.resizeNoCopy(numParticles, (int)pow(num+1,outputDim));

    Matrix1D <int> numElems;
    numElems.resizeNoCopy(MAT_XSIZE(squares));

    for (size_t i=0; i<numParticles; i++)
    {
    	int index=0;
    	for (int n=0; n<outputDim; ++n)
    	{
    		index*=num;
    		index+=(int)floor(((VEC_ELEM(coor[n],i) - minCoor[n]) / (maxCoor[n] - minCoor[n]))*num);
    	}
        MAT_ELEM(squares,VEC_ELEM(numElems,index),index)=i;
        VEC_ELEM(numElems,index)++;
    }

    int numElemFirstRow = 0;
    for (size_t k=0; k < MAT_XSIZE(squares); k++)
        if (VEC_ELEM(numElems,k)!=0)
            numElemFirstRow++;

    std::vector<Matrix1D<double> > selectedCoor;
    dummy.resizeNoCopy(numElemFirstRow);
    for (int n=0; n<outputDim; ++n)
    	selectedCoor.push_back(dummy);

    //Metadata with the well sampled projection and random projections assigned
    MetaData SFout;
    int indx = 0;
    randomize_random_generator();
    std::vector<String> fnImg;
    SFin.getColumnValues(MDL_IMAGE,fnImg);
    for (size_t k=0; k <MAT_XSIZE(squares); k++)
    {
        if (MAT_ELEM(squares,0,k)!=0)
        {
        	int randomN = round((rnd_unif(0,1))*VEC_ELEM(numElems,k));
            for (int n=0; n<outputDim; ++n)
            	VEC_ELEM(selectedCoor[n],indx)=VEC_ELEM(coor[n],MAT_ELEM(squares,randomN,k));
            size_t id=SFout.addObject();
            SFout.setValue(MDL_IMAGE,fnImg[MAT_ELEM(squares,randomN,k)],id);
            SFout.setValue(MDL_ANGLE_ROT,(rnd_unif(0,360)-180),id);
            SFout.setValue(MDL_ANGLE_TILT,(rnd_unif(0,180)),id);
            SFout.setValue(MDL_ANGLE_TILT,(rnd_unif(0,180)),id);
            SFout.setValue(MDL_ANGLE_PSI,(rnd_unif(0,360)),id);
            indx++;
        }
    }

    SFout.write(fnOut);
}

// Run  ====================================================================
void ProgTransformDimRed::run()
{
    show();
    produceSideInfo();
    if (outputDim<0)
    	estimateDimension();

    if (!SFin.containsLabel(MDL_DIMRED))
    {
        algorithm->reduceDimensionality();

        std::vector<double> dimredProj;
        dimredProj.resize(outputDim);
        int i=0;
        const Matrix2D<double> &Y=algorithm->getReducedData();
        FOR_ALL_OBJECTS_IN_METADATA(SFin)
        {
            memcpy(&dimredProj[0],&MAT_ELEM(Y,i,0),outputDim*sizeof(double));
            SFin.setValue(MDL_DIMRED,dimredProj,__iter.objId);
            i++;
        }

        SFin.write(fnOut,MD_APPEND);
    }

    if (fnRandomSampling!="")
    	extractRandomProjections();
}
