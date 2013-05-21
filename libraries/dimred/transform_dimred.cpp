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
#include <data/matrix1d.h>
#include <data/matrix2d.h>
#include <time.h>
#include "diffusionMaps.h"

void ProgTransformDimRed::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    dimRefMethod = getParam("-m");
    outputDim  = getIntParam("-d");
    numGrids  =  getIntParam("-n");
}


// Show ====================================================================
void ProgTransformDimRed::show()
{
    if (verbose>0)
        std::cerr
        << "Input metadata file:    "           << fnIn          << std::endl
        << "Output metadata:        "           << fnOut         << std::endl
        << "Dim Red Method:     "               << dimRefMethod  << std::endl
        << "Number of dimensions after dimensionality reduction:     "      << dimRefMethod  << std::endl
        ;
}

// usage ===================================================================
void ProgTransformDimRed::defineParams()
{
    addParamsLine("   -i <metadatafile>             : metadata file  with images");
    addParamsLine("  [-o <metadatafile=\"\">]       : output metadata with distances between images");
    addParamsLine("   -m <dimRefMethod>             : Dimensionality Reduction method selected");
    addParamsLine("  [-d <N=2>]               : Number of dimensions after the dimensionality reduction");
    addParamsLine("  [-n <num=3>]             : Number of grids in one dimension");
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

void ProgTransformDimRed::extractRandomProjections()
{
    //number of grids in each axis
    size_t num = (int)numGrids;

    //Metadata that has the dimrefcoeffs
    MetaData SFin;
    SFin.read(fnIn);

    size_t numParticles = SFin.size();
    std::vector<double> dimredProj;
    size_t i=0;
    Matrix1D< double > xcorr, ycorr;
    xcorr.resizeNoCopy(numParticles);
    ycorr.resizeNoCopy(numParticles);

    //We store the x and y coordinates in xcorr and ycorr matrix1d
    FOR_ALL_OBJECTS_IN_METADATA(SFin)
    {
        SFin.getValue(MDL_DIMRED,dimredProj,__iter.objId);
        VEC_ELEM(xcorr,i) = dimredProj[0];
        VEC_ELEM(ycorr,i) = dimredProj[1];
        i++;
    }

    //We obtain the max and min coordiante values
    int minIndX = xcorr.minIndex();
    int minIndY = ycorr.minIndex();
    int maxIndX = xcorr.maxIndex();
    int maxIndY = ycorr.maxIndex();

    double maxX = VEC_ELEM(xcorr,maxIndX);
    double minX = VEC_ELEM(xcorr,minIndX);
    double maxY = VEC_ELEM(ycorr,maxIndY);
    double minY = VEC_ELEM(ycorr,minIndY);

    Matrix2D <int> squares;
    squares.resizeNoCopy(numParticles, (num+1)*(num+1));

    Matrix1D <int> numElems;
    numElems.resizeNoCopy((num+1)*(num+1));

    for (size_t i=0; i <  numParticles; i++)
    {
        int indexX = floor(((VEC_ELEM(xcorr,i) - minX) / (maxX - minX))*num);
        int indexY = floor(((VEC_ELEM(ycorr,i) - minY) / (maxY - minY))*num);
        int index = (indexY)*num+indexX;
        MAT_ELEM(squares,VEC_ELEM(numElems,index),index)=i;
        VEC_ELEM(numElems,index)++;
    }

    int numElemFirstRow = 0;
    for (size_t k=0; k <  (num+1)*(num+1); k++)
    {
        if (VEC_ELEM(numElems,k))
            numElemFirstRow++;
    }

    Matrix1D< double > xcorrChoiced, ycorrChoiced;
    xcorrChoiced.resizeNoCopy(numElemFirstRow);
    ycorrChoiced.resizeNoCopy(numElemFirstRow);

    int indx = 0, n;

    randomize_random_generator();

    //Metadata with the well sampled projection and random projections assigned
    MetaData SFout;
    FileName fnMic;
    for (size_t k=0; k <  (num+1)*(num+1); k++)
    {
        if (MAT_ELEM(squares,0,k))
        {
            n = round((rnd_unif(0,1))*VEC_ELEM(numElems,k));
            VEC_ELEM(xcorrChoiced,indx)=VEC_ELEM(xcorr,MAT_ELEM(squares,n,k));
            VEC_ELEM(ycorrChoiced,indx)=VEC_ELEM(ycorr,MAT_ELEM(squares,n,k));
            size_t temp = (size_t)(MAT_ELEM(squares,n,k)+1);
            SFin.getValue(MDL_IMAGE,fnMic,temp);
            size_t id=SFout.addObject();
            SFout.setValue(MDL_IMAGE,fnMic,id);
            SFout.setValue(MDL_ANGLE_ROT,(rnd_unif(0,360)-180),id);
            SFout.setValue(MDL_ANGLE_TILT,(rnd_unif(0,180)),id);
            SFout.setValue(MDL_ANGLE_TILT,(rnd_unif(0,180)),id);
            SFout.setValue(MDL_ANGLE_PSI,(rnd_unif(0,360)),id);
            indx++;
        }
    }

    std::cout << " fnOut : " << fnOut << std::endl;
    SFout.write(fnOut);
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

    MetaData SF;
    SF=SFin;
    if (!SF.containsLabel(MDL_DIMRED))
    {

        DiffusionMaps dimred;
        dimred.setInputData(X);
        dimred.setSpecificParameters();
        dimred.setOutputDimensionality(outputDim);
        dimred.distance=&correlationDistance;
        dimred.reduceDimensionality();

        std::vector<double> dimredProj;
        dimredProj.resize(outputDim);
        int i=0;
        const Matrix2D<double> &Y=dimred.getReducedData();
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            memcpy(&dimredProj[0],&MAT_ELEM(Y,i,0),outputDim*sizeof(double));
            SF.setValue(MDL_DIMRED,dimredProj,__iter.objId);
            i++;
        }

        SF.write(fnIn,MD_OVERWRITE);

    }

    extractRandomProjections();
}
