/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2013)
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

#include "laplacianEigenmaps.h"

void LaplacianEigenmap::setSpecificParameters(double sigma, size_t numberOfNeighbours)
{
	this->sigma=sigma;
	this->numberOfNeighbours=numberOfNeighbours;
}

void LaplacianEigenmap::reduceDimensionality()
{
	Matrix2D<double> G,L;
	size_t xsize=MAT_XSIZE(G);
	//std::cout<<"Constructing neighborhood graph...";
	computeDistanceToNeighbours(*X,numberOfNeighbours,G);
	//G.write("dimred/G2.txt");
	FOR_ALL_ELEMENTS_IN_MATRIX2D(G)
		MAT_ELEM(G,i,j)*=MAT_ELEM(G,i,j);
	G/=G.computeMax();
	//G.write("dimred/G.txt");

	//Function component

	Matrix1D<int> component;
	connectedComponentsOfUndirectedGraph(G,component);

    int max=0;
    int indMAX=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(component)
    {
        if(VEC_ELEM(component,i)>max)
        {
            max=VEC_ELEM(component,i);
            indMAX=i;
        }
    }

    Matrix1D<int> count;
    count.resizeNoCopy(max);

    for(int j=0;j<max;j++)
    {
        int NoElem=0;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(component)
        {
            if(VEC_ELEM(component,i)==i)
                ++NoElem;
        }
        VEC_ELEM(count,j)=NoElem;
    }

    Matrix1D<int> conn_comp;
    //conn_comp.initZeros(VEC_XSIZE(component));
    conn_comp.resizeNoCopy(VEC_XSIZE(component));
    size_t ord=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(component)
    {
        if(VEC_ELEM(component,i)==indMAX)
        {
            VEC_ELEM(conn_comp,ord)=i;
            ++ord;
        }
    }

    Matrix1D<int> conn_comp2;
    for(size_t j=0;j<ord-1;j++)
    		VEC_ELEM(conn_comp2,j)=VEC_ELEM(conn_comp,j);



	//std::cout<<"Computing weigth matrices...";
	double K=-1.0/(2*sigma*sigma);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(G)
	{
		if(MAT_ELEM(G,i,j)!=0)
			MAT_ELEM(G,i,j)=exp(-MAT_ELEM(G,i,j)*K);
	}
	FOR_ALL_ELEMENTS_IN_MATRIX2D(G)
	{
		if(i!=j)
			MAT_ELEM(L,i,j)=-MAT_ELEM(G,i,j);
		else
		{
			size_t RowSum=0;
			for(size_t jj=0;jj<xsize;++jj)
				RowSum+=MAT_ELEM(G,i,jj);
			MAT_ELEM(L,i,j)=RowSum-MAT_ELEM(G,i,j);
		}
	}
	//std::cout<<"Constructing Eigenmaps...";
}
