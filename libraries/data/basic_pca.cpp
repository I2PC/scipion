/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#include "basic_pca.h"
#include "matrix2d.h"

/* Subtract average ------------------------------------------------------- */
void PCAMahalanobisAnalyzer::subtractAvg()
{
    int N=v.size();
    if (N==0)
        return;
    // Compute average

    typeCast(v[0],avg);
    for (int n=1; n<N; n++)
    {
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(avg)
        DIRECT_A1D_ELEM(avg,i)+=DIRECT_A1D_ELEM(aux,i);
    }
    avg/=N;

    // Subtract average and compute stddev
    MultidimArray<float> avgF;
    typeCast(avg,avgF);
    for (int n=0; n<N; n++)
        v[n]-=avgF;
}

/* Standardize variables -------------------------------------------------- */
void PCAMahalanobisAnalyzer::standardarizeVariables()
{
    int N=v.size();
    if (N<=1)
        return;
    // Compute average
    MultidimArray<double> avg;
    typeCast(v[0],avg);
    for (int n=1; n<N; n++)
    {
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(avg)
        DIRECT_A1D_ELEM(avg,i)+=DIRECT_A1D_ELEM(aux,i);
    }
    avg/=N;

    // Subtract average and compute stddev
    MultidimArray<float> avgF;
    typeCast(avg,avgF);
    MultidimArray<double> stddev;
    stddev.initZeros(avg);
    for (int n=0; n<N; n++)
    {
        v[n]-=avgF;
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(aux)
        {
            float f=DIRECT_A1D_ELEM(aux,i);
            DIRECT_A1D_ELEM(stddev,i)+=f*f;
        }
    }
    stddev/=N-1;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(stddev)
    DIRECT_A1D_ELEM(stddev,i)=sqrt(DIRECT_A1D_ELEM(stddev,i));

    // Divide by stddev
    MultidimArray<float> istddevF;
    istddevF.resize(stddev);
    size_t NoInformation=0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(stddev)
    if (DIRECT_A1D_ELEM(stddev,i)>XMIPP_EQUAL_ACCURACY)
    	DIRECT_A1D_ELEM(istddevF,i)=(float)(1.0/DIRECT_A1D_ELEM(stddev,i));
    else
    {
    	DIRECT_A1D_ELEM(istddevF,i)=0.0;
    	NoInformation++;
    }
    if (NoInformation==XSIZE(stddev))
    	REPORT_ERROR(ERR_NUMERICAL,"There is no information to perform PCA");
    for (int n=0; n<N; n++)
    {
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(aux)
        DIRECT_A1D_ELEM(aux,i)*=DIRECT_A1D_ELEM(istddevF,i);
    }
}

/* Add vector ------------------------------------------------------------- */
void PCAMahalanobisAnalyzer::projectOnPCABasis(Matrix2D<double> &CtY)
{
    int N=v.size();
    int NPCA=PCAbasis.size();
    CtY.initZeros(NPCA,N);
    for (int ii=0; ii<N; ii++)
    {
        const MultidimArray<float> &Iii=v[ii];
        const size_t unroll=4;
        size_t nmax=unroll*(MULTIDIM_SIZE(Iii)/unroll);
        for (int jj=0; jj<NPCA; jj++)
        {
            const MultidimArray<double> &Ijj=PCAbasis[jj];

            double dotProduct=0;
            size_t n=0;
            const float *ptrii=MULTIDIM_ARRAY(Iii);
            const double *ptrjj=MULTIDIM_ARRAY(Ijj);
            for (size_t n=0; n<nmax; n+=unroll, ptrii+=unroll, ptrjj+=unroll)
            {
                dotProduct += *ptrii * *ptrjj;
                dotProduct += *(ptrii+1) * *(ptrjj+1);
                dotProduct += *(ptrii+2) * *(ptrjj+2);
                dotProduct += *(ptrii+3) * *(ptrjj+3);
            }
            for (n=nmax, ptrii=MULTIDIM_ARRAY(Iii)+nmax, ptrjj=MULTIDIM_ARRAY(Ijj)+nmax;
                 n<MULTIDIM_SIZE(Iii); ++n, ++ptrii, ++ptrjj)
                dotProduct += *ptrii * *ptrjj;
            MAT_ELEM(CtY,jj,ii)=dotProduct;
        }
    }
}

void PCAMahalanobisAnalyzer::learnPCABasis(size_t NPCA, size_t Niter)
{
    // Take the first vectors for the PCA basis
    MultidimArray<double> vPCA;
    NPCA=XMIPP_MIN(NPCA,v.size());
    for (size_t n=0; n<NPCA; n++)
    {
    	size_t nRandom=(size_t)round(v.size()*rnd_unif(0,1));
        typeCast(v[nRandom],vPCA);
        PCAbasis.push_back(vPCA);
    }

    size_t N=v.size();
    for (size_t n=0; n<Niter; n++)
    {
        // E-step ..........................................................
        // Compute C^t*C
        Matrix2D<double> CtC(NPCA,NPCA);
        for (size_t ii=0; ii<NPCA; ii++)
        {
            const MultidimArray<double> &Iii=PCAbasis[ii];
            for (size_t jj=ii; jj<NPCA; jj++)
            {
                const MultidimArray<double> &Ijj=PCAbasis[jj];
                CtC(ii,jj)=0;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Ijj)
                MAT_ELEM(CtC,ii,jj)+=
                    DIRECT_A1D_ELEM(Iii,i)*DIRECT_A1D_ELEM(Ijj,i);

                if (ii!=jj)
                    CtC(jj,ii)=CtC(ii,jj);
            }
        }

        Matrix2D<double> CtY, X;
        projectOnPCABasis(CtY);
        X=CtC.inv()*CtY;

        // M-step ..........................................................
        Matrix2D<double> XtXXtinv=X.transpose()*(X*X.transpose()).inv();
        for (size_t ii=0;ii<NPCA;ii++)
        {
            MultidimArray<double> &Ipca=PCAbasis[ii];
            const size_t unroll=4;
            size_t nmax=unroll*(MULTIDIM_SIZE(Ipca)/unroll);
            Ipca.initZeros();
            for (size_t jj=0; jj<N; jj++)
            {
                const MultidimArray<float> &I=v[jj];
                const float *ptrI=MULTIDIM_ARRAY(I);
                double *ptrPCA=MULTIDIM_ARRAY(Ipca);
                double val=MAT_ELEM(XtXXtinv,jj,ii);
                size_t n;
                for (n=0; n<nmax; n+=unroll, ptrPCA+=unroll, ptrI+=unroll)
                {
                    *ptrPCA += *ptrI * val;
                    *(ptrPCA+1) += *(ptrI+1) * val;
                    *(ptrPCA+2) += *(ptrI+2) * val;
                    *(ptrPCA+3) += *(ptrI+3) * val;
                }
                for (n=nmax, ptrI=MULTIDIM_ARRAY(I)+nmax, ptrPCA=MULTIDIM_ARRAY(Ipca)+nmax;
                     n<MULTIDIM_SIZE(I); ++n, ++ptrI, ++ptrPCA)
                    *ptrPCA += *ptrI * val;
            }
        }
    }

    //Obtain the Orthonormal vectors for the C (According to paper)
    gramSchmidt();

    //Generate a Matrix for true PCA and compute the average C'*data
    Matrix2D<double> data;
    MultidimArray<double> average;
    data.initZeros(NPCA,v.size());
    average.initZeros(NPCA);
    for (size_t ii=0;ii<NPCA;ii++)
    {
        MultidimArray<double> &C=PCAbasis[ii];
        for (size_t jj=0;jj<v.size();jj++)
        {
            MultidimArray<float> &D=v[jj];
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(C)
            MAT_ELEM(data,ii,jj)+=DIRECT_A1D_ELEM(C,i)*DIRECT_A1D_ELEM(D,i);
            DIRECT_A1D_ELEM(average,ii)+=MAT_ELEM(data,ii,jj);
        }
    }
    average/=v.size();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(data)
    MAT_ELEM(data,i,j)-=DIRECT_A1D_ELEM(average,i);

    Matrix2D<double> covarMatrix;
    Matrix2D<double> u;
    Matrix2D<double> v;
    covarMatrix=data*data.transpose();
    svdcmp(covarMatrix,u,w,v);

    /// Do the PCAbasis*v
    size_t xsize=XSIZE(PCAbasis[0]);
    MultidimArray<double> temp;
    for (size_t i=0;i<xsize;i++)
    {
        temp.initZeros(1,NPCA);
        for (size_t j=0;j<NPCA;j++)
            for (size_t k=0;k<NPCA;k++)
            {
                const MultidimArray<double> &Iii=PCAbasis[k];
                DIRECT_A1D_ELEM(temp,j)+=DIRECT_A1D_ELEM(Iii,i)* MAT_ELEM(v,k,j);
            }
        for (size_t k=0;k<NPCA;k++)
        {
            const MultidimArray<double> &Iii=PCAbasis[k];
            DIRECT_A1D_ELEM(Iii,i)=DIRECT_A1D_ELEM(temp,k);
        }
    }

    // Normalize output vectors
    for (size_t ii=0;ii<NPCA;ii++)
    {
        double norm=sqrt(PCAbasis[ii].sum2());
        PCAbasis[ii]/=norm;
    }

}

/** computes the orthonormal basis for a subspace
   *
   *This method computes the orthonormal basis for the vectors of an Euclidean
   *which is formed by the vectors on each column of the matrix. The output
   *of this method will be a matrix in which each column form an orthonormal
   *basis.
   *
   * @code
   * gramSchmidt();
   * @endcode
   */
void PCAMahalanobisAnalyzer::gramSchmidt()
{
    MultidimArray<double> v,orthonormalBasis;
    for (size_t ii=0;ii<PCAbasis.size();ii++)
    {
        /// Take the current vector
        MultidimArray<double> &Iii=PCAbasis[ii];
        if (ii==0)
        {
            double inorm=1.0/sqrt(Iii.sum2());
            for (size_t i=0;i<XSIZE(Iii);i++)
                DIRECT_A1D_ELEM(Iii,i)=DIRECT_A1D_ELEM(Iii,i)*inorm;
            continue;
        }
        /// Take the previous vectors
        v.initZeros(XSIZE(Iii));
        orthonormalBasis.initZeros(XSIZE(Iii));
        for (size_t jj=0;jj<ii;jj++)
        {
            const MultidimArray<double> &Ijj=PCAbasis[jj];
            double dotProduct=Iii.dotProduct(Ijj);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v)
            DIRECT_A1D_ELEM(v,i)-=dotProduct*DIRECT_A1D_ELEM(Ijj,i);
        }
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(orthonormalBasis)
        DIRECT_A1D_ELEM(orthonormalBasis,i)= DIRECT_A1D_ELEM(Iii,i)+ DIRECT_A1D_ELEM(v,i);
        double inorm=1.0/sqrt(orthonormalBasis.sum2());
        for (size_t i=0;i<XSIZE(Iii);i++)
            DIRECT_A1D_ELEM(Iii,i)=DIRECT_A1D_ELEM(orthonormalBasis,i)*inorm;
    }
}

/* Compute Statistics ----------------------------------------------------- */
void PCAMahalanobisAnalyzer::computeStatistics(MultidimArray<double> & avg,
        MultidimArray<double> & stddev)
{
    int N=v.size();
    if (N==0)
    {
        avg.clear();
        stddev.clear();
        return;
    }
    typeCast(v[0],avg);
    MultidimArray<double> aux;
    for (int n=1; n<N; n++)
    {
        typeCast(v[n],aux);
        avg+=aux;
    }
    avg/=N;
    stddev.initZeros(avg);
    for (int n=0; n<N; n++)
    {
        typeCast(v[n],aux);
        aux-=avg;
        aux*=aux;
        stddev+=aux;
    }
    double iN_1=1.0/(N-1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(stddev)
    DIRECT_A1D_ELEM(stddev,i)=sqrt(DIRECT_A1D_ELEM(stddev,i)*iN_1);
}

/* Evaluate score --------------------------------------------------------- */
//#define DEBUG
void PCAMahalanobisAnalyzer::evaluateZScore(int NPCA, int Niter)
{
    int N=v.size();
    if (N==0)
    {
        Zscore.clear();
        idx.clear();
        return;
    }
    else if (N==1)
    {
        Zscore.initZeros(N);
        Zscore.indexSort(idx);
        return;
    }

#ifdef DEBUG
    std::cout << "Input vectors\n";
    for (int n=0; n<N; n++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v[n])
        std::cout << DIRECT_A1D_ELEM(v[n],i) << " ";
        std::cout << std::endl;
    }
#endif
    subtractAvg();
#ifdef DEBUG

    std::cout << "\n\nInput vectors after normalization\n";
    for (int n=0; n<N; n++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v[n])
        std::cout << DIRECT_A1D_ELEM(v[n],i) << " ";
        std::cout << std::endl;
    }
#endif

    learnPCABasis(NPCA, Niter);
#ifdef DEBUG

    std::cout << "\n\nPCA basis\n";
    for (int n=0; n<NPCA; n++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(PCAbasis[n])
        std::cout << DIRECT_A1D_ELEM(PCAbasis[n],i) << " ";
        std::cout << std::endl;
    }
#endif

    Matrix2D<double> proj;
    projectOnPCABasis(proj);

#ifdef DEBUG

    std::cout << "\n\nPCA projected values\n" << proj << std::endl;
#endif

    // Estimate covariance matrix
    Matrix2D<double> cov(NPCA,NPCA);
    for (int ii=0; ii<N; ii++)
    {
        for (int i=0; i<NPCA; i++)
            for (int j=0; j<NPCA; j++)
                MAT_ELEM(cov,i,j)+=MAT_ELEM(proj,i,ii)*MAT_ELEM(proj,j,ii);
    }
    cov/=N;
#ifdef DEBUG

    std::cout << "\n\nCovariance matrix\n" << cov << std::endl;
#endif

    Matrix2D<double> covinv=cov.inv();
    Zscore.initZeros(N);
    for (int ii=0; ii<N; ii++)
    {
        Zscore(ii)=0;
        for (int i=0; i<NPCA; i++)
        {
            double xi=MAT_ELEM(proj,i,ii);
            for (int j=0; j<NPCA; j++)
                DIRECT_A1D_ELEM(Zscore,ii)+=xi*MAT_ELEM(proj,j,ii)*MAT_ELEM(covinv,i,j);
        }
        Zscore(ii)=sqrt(fabs(Zscore(ii)));
    }
#ifdef DEBUG
    std::cout << "\n\nZscores\n";
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Zscore)
    std::cout << DIRECT_A1D_ELEM(Zscore,i) << " ";
    std::cout << std::endl;
#endif

    Zscore.indexSort(idx);
}
#undef DEBUG
