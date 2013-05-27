/***************************************************************************
 *
 * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
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

#include "knn_classifier.h"

KNN::KNN(int k)
{
    setK(k);
    neighborsIndex.resize(1,1,1,K);
}

void KNN::train(MultidimArray<double> &dataset,MultidimArray<double> &dataLabel,
                MultidimArray<double> &labelset)
{
    __dataset=dataset;
    __dataLabel=dataLabel;
    __labelSet=labelset;
}

void KNN::KNearestNeighbors(MultidimArray<double> &sample)
{
    double maximum,distance;
    int maximumIndex;
    maxDist.resize(1,1,1,K);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(maxDist)
    {
        DIRECT_A1D_ELEM(maxDist,i)=euclideanDistance(sample,i,-1.0);
        DIRECT_A1D_ELEM(neighborsIndex,i)=i;
    }
    maximumIndex=findMaxIndex(maxDist);
    maximum=DIRECT_A1D_ELEM(maxDist,maximumIndex);
    for (size_t i=K;i<YSIZE(__dataset);i++)
    {
        distance=euclideanDistance(sample,i,maximum);
        if (distance==-1)
            continue;
        if (distance<maximum)
        {
            DIRECT_A1D_ELEM(maxDist,maximumIndex)=distance;
            DIRECT_A1D_ELEM(neighborsIndex,maximumIndex)=i;
            maximumIndex=findMaxIndex(maxDist);
            maximum=DIRECT_A1D_ELEM(maxDist,maximumIndex);
        }
    }
}

int KNN::predict(MultidimArray<double> &sample,double &score)
{
    MultidimArray<double> voteArray;
    int index;
    voteArray.initZeros(XSIZE(__labelSet));
    KNearestNeighbors(sample);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(neighborsIndex)
    {
        index=DIRECT_A1D_ELEM(neighborsIndex,i);
        for (size_t j=0;j<XSIZE(__labelSet);++j)
            if (DIRECT_A1D_ELEM(__labelSet,j)==DIRECT_A1D_ELEM(__dataLabel,index))
                DIRECT_A1D_ELEM(voteArray,j)+=1;
    }
    index=findMaxIndex(voteArray);
    score=DIRECT_A1D_ELEM(voteArray,index)/double(K);
    if (DIRECT_A1D_ELEM(voteArray,index)>(K*0.5))
        return (int)DIRECT_A1D_ELEM(__labelSet,index);
    index=findMinIndex(maxDist);
    return (int)DIRECT_A1D_ELEM(__dataLabel,DIRECT_A1D_ELEM(neighborsIndex,index));
}

int KNN::findMaxIndex(MultidimArray<double> &inputArray)
{
    double maximum;
    int maximumIndex;
    maximum=DIRECT_A1D_ELEM(inputArray,0);
    maximumIndex=0;
    for (size_t i=1;i<XSIZE(inputArray);i++)
        if (maximum<DIRECT_A1D_ELEM(inputArray,i))
        {
            maximum=DIRECT_A1D_ELEM(inputArray,i);
            maximumIndex=i;
        }
    return maximumIndex;
}

int KNN::findMinIndex(MultidimArray<double> &inputArray)
{
    double minimum;
    int minimumIndex;
    minimum=DIRECT_A1D_ELEM(inputArray,0);
    minimumIndex=0;
    for (size_t i=1;i<XSIZE(inputArray);i++)
        if (minimum>DIRECT_A1D_ELEM(inputArray,i))
        {
            minimum=DIRECT_A1D_ELEM(inputArray,i);
            minimumIndex=i;
        }
    return minimumIndex;
}

double KNN::euclideanDistance(MultidimArray<double> &sample,int index,double maximumDist)
{
    double dist=0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sample)
    {
        double tmp=DIRECT_A1D_ELEM(sample,i)-DIRECT_A2D_ELEM(__dataset,index,i);
        dist+=tmp*tmp;
        if (maximumDist>0 && sqrt(dist)>maximumDist)
            return -1.0;
    }
    return sqrt(dist);
}

double KNN::cityBlockDistance(MultidimArray<double> &sample,int index,double maximumDist)
{
    double dist=0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sample)
    {
        double tmp=fabs(DIRECT_A1D_ELEM(sample,i)-DIRECT_A2D_ELEM(__dataset,index,i));
        dist+=tmp;
        if (maximumDist>0 && dist>maximumDist)
            return -1.0;
    }
    return dist;
}

void KNN::saveModel(const FileName &fn)
{
    std::ofstream fh;
    fh.open(fn.c_str());
    fh<<XSIZE(__labelSet)<<std::endl;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(__labelSet)
    fh<<DIRECT_A1D_ELEM(__labelSet,i)<<" ";
    fh<<std::endl;
    fh<<YSIZE(__dataset)<<" "<<XSIZE(__dataset)<<std::endl;
    for (size_t i=0;i<YSIZE(__dataset);++i)
    {
        fh<<DIRECT_A1D_ELEM(__dataLabel,i)<<" ";
        for (size_t j=0;j<XSIZE(__dataset);++j)
        {
            fh<<DIRECT_A2D_ELEM(__dataset,i,j)<<" ";
        }
        fh<<std::endl;
    }
    fh.close();
}

void KNN::loadModel(const FileName &fn)
{
    std::ifstream fh;
    fh.open(fn.c_str());
    int x,y;
    fh>>x;
    __labelSet.resize(1,1,1,x);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(__labelSet)
    fh>>DIRECT_A1D_ELEM(__labelSet,i);
    fh>>y>>x;
    __dataset.resize(1,1,y,x);
    __dataLabel.resize(1,1,1,y);
    for (size_t i=0;i<YSIZE(__dataset);++i)
    {
        fh>>DIRECT_A1D_ELEM(__dataLabel,i);
        for (size_t j=0;j<XSIZE(__dataset);++j)
        {
            fh>>DIRECT_A2D_ELEM(__dataset,i,j);
        }
    }
    fh.close();
}

void KNN::setK(int k)
{
    K=k;
}


