/***************************************************************************
 *
 * Author:    Itziar Benito Ortega     (2013)
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

#include "probabilisticPCA.h"
#include <data/multidim_array.h>

void ProbabilisticPCA::setSpecificParameters(size_t Niters)
{
    this->Niters=Niters;
}

void ProbabilisticPCA::reduceDimensionality()
{
    size_t N=MAT_YSIZE(*X);  // N= number of rows of X
    size_t D=MAT_XSIZE(*X);  // D= number of columns of X

    bool converged=false;
    size_t iter=0;
    double sigma2=rnd_unif()*2;

    // S=cov(X)
    Matrix2D<double> S, W, inW, invM;
    subtractColumnMeans(*X);
    matrixOperation_AtA(*X,S);
    S/=(double)N;

    W.initRandom(D,outputDim,0,2,RND_UNIFORM);
    matrixOperation_AtA(W,inW);

    Matrix2D <double> Ez,WtX;
    MultidimArray <double> Ezz(outputDim,outputDim,MAT_YSIZE(*X));
    while (!converged && iter<=Niters)
    {
        ++iter;

        // Perform E-step
        // inW+=sigma2*I
        for (size_t i=0; i<D; ++i)
        	MAT_ELEM(inW,i,i)+=sigma2;
        inW.inv(invM);

        //Ez=invM*W^t*X
        matrixOperation_AtB(W,*X,WtX);
        matrixOperation_AB(invM,WtX,Ez);

        for (size_t i=0; i<N; ++i)
        {

        }

        // Transpuesta de W
        Matrix2D <double> transW = W.transpose();

        //Creo la matriz res = invM_por_Wtrans; res=invM*W'
        Matrix2D <double> res;
        res.initZeros(MAT_YSIZE(invM),MAT_XSIZE(transW));

        // Con este bucle voy recorriendo las filas de invM por un lado y las
        // filas de W' por otro. Modificando k recorro las columnas de ambos
        // y accedo a cada elemento de la matriz res guardando invM*W'.
        for (size_t i=0; i< MAT_YSIZE(invM); ++i)
        {
            for (size_t j=0; j< MAT_XSIZE(transW); ++j)
            {
                for (size_t k=0; k< MAT_XSIZE(transW); ++k)
                {
                    MAT_ELEM(res,i,j) += MAT_ELEM(invM,i,k) * MAT_ELEM(transW,k,j);
                }
            }

        }

        //Escribo un fichero por cada iteracion con el contenido de res para probar
        //char nombrefichero[10000];
        //sprintf(nombrefichero,"dimred/res%ld.txt",iter);
        //res.write(nombrefichero);


        // Ahora calculo Ez(:,i) =invM*W'*X(:,i)= res*X(:,i)

        for (size_t index = 0; index < MAT_YSIZE(*X); ++index)
        {
            //Ez(:,index) =invM*W'*X(:,index);
            Matrix1D <double> columnaX;
            X->getCol(index,columnaX);

            Matrix1D <double> col_tmp;
            // Las columnas tienen de tamagno el numero de filas de res
            col_tmp.initZeros(MAT_YSIZE(res));

            // Multiplico res por columnaX el resultado lo guardo en col_tmp
            // Recorro las filas de res y las multiplico por columnaX

            for(size_t i=0; i<MAT_YSIZE(res);++i)
            {
                for(size_t k=0; k<MAT_YSIZE(res); ++k)
                {
                    VEC_ELEM(col_tmp,i) += MAT_ELEM(res,i,k) * VEC_ELEM(columnaX,k);
                }
            }
            // Guardo col_tmp en la columna index de Ez
            Ez.setCol(index,col_tmp);





            //Ezz(:,:,index)=sigma2*invM+Ez(:,index)*Ez(:,index)';

            /*Matrix2D<double> EzEz_(outputDim,outputDim);

            for (size_t i1  = 0; i1 < outputDim; ++i1)
        {
             for (size_t j1 = 0; j1 < outputDim; ++j1)
             {
              MAT_ELEM(EzEz_, i1, j1) = A2D_ELEM(Ez, index, i1) * A2D_ELEM(Ez, index, j1);
             }
        }

            FOR_ALL_ELEMENTS_IN_MATRIX2D(invM)
        {
             A3D_ELEM(Ezz,index,i,j) = (sigma2 * MAT_ELEM(invM, i, j)) + MAT_ELEM(EzEz_, i, j);
        }*/

        }

    }





    //Creo las matrices Wp1 y Wp2 vacias para luego poder iterar sobre ellas.

    Matrix2D <double> Wp1;
    Wp1.initZeros(D,outputDim);
    Matrix2D <double> Wp2;
    Wp2.initZeros(outputDim,outputDim);


    //Para index de 1 a N
    for (size_t index = 0; index<N; ++index)
    {
        //Saco las columnas que voy a multiplicar
        Matrix1D <double> columnaX; //matrix1D = vector
        X->getCol(index,columnaX);
        Matrix1D <double> columnaEz;
        Ez.getCol(index,columnaEz);

        //Aqui guardo el resultado de multiplicar los dos vectores.
        //La multiplicacion de dos vectores da una matriz de 2 dimensiones(3x3)
        //La matriz resultado tendra de tamagno la longitud de un vector
        //por la longitud del otro.Con vdim obtengo el tamagno del vector
        Matrix2D <double> colXcolEz;
        colXcolEz.initZeros(columnaX.vdim, columnaEz.vdim);


        for(size_t i=0; i<columnaX.vdim;++i)
        {
            for(size_t j=0; j<columnaEz.vdim; ++j)
            {
                MAT_ELEM(colXcolEz,i,j) = VEC_ELEM(columnaX,i) * VEC_ELEM(columnaEz,j);
            }
        }

    }
}
