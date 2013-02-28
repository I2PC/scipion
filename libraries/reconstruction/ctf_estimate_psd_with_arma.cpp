/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Javier Angel Velazquez Muriel (javi@cnb.csic.es)
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

#include "ctf_estimate_psd_with_arma.h"

#include <data/args.h>
#include <data/xmipp_fftw.h>
#include <data/filters.h>

/* Read parameters from command line --------------------------------------- */
void ARMA_parameters::readParams(XmippProgram *program)
{
    N_AR         = program->getIntParam("--N_AR");
    M_AR         = program->getIntParam("--M_AR");
    N_MA         = program->getIntParam("--N_MA");
    M_MA         = program->getIntParam("--M_MA");
}

void ARMA_parameters::defineParams(XmippProgram *program)
{
	program->addParamsLine("==++ ARMA models");
	program->addParamsLine("  [--N_AR <N=24>]    : N order of the AR model");
	program->addParamsLine("  [--M_AR <N=24>]    : M order of the AR model");
	program->addParamsLine("  [--N_MA <N=20>]    : N order of the MA model");
	program->addParamsLine("  [--M_MA <N=20>]    : M order of the MA model");
}

// First quadrant neighbours -----------------------------------------------
void First_Quadrant_Neighbors(int N, int M, MultidimArray<double> &Neighbors)
{
    long NumberOfPoints = (N + 1) * M;
    int n;
    Neighbors.resize(NumberOfPoints, 3);

    n = 0; // Number of neighbors found so far.

    for (int p = N; p >= 0; p--)
    {
        for (int q = M; q > 0; q--)
        {
            A2D_ELEM(Neighbors,n, 0) = (double)p;
            A2D_ELEM(Neighbors,n, 1) = (double)q;
            n++;
        }
    }
}

// Second quadrant neighbours ----------------------------------------------
void Second_Quadrant_Neighbors(int N, int M, MultidimArray<double> &Neighbors)
{
    long NumberOfPoints = N * (M + 1);
    int n;
    Neighbors.resize(NumberOfPoints, 3);

    n = 0; // Number of neighbors found so far.

    for (int p = N; p >= 1; p--)
    {
        for (int q = 0; q >= -M; q--)
        {
            A2D_ELEM(Neighbors,n, 0) = (double)p;
            A2D_ELEM(Neighbors,n, 1) = (double)q;
            n++;
        }
    }
}

// Compute ARMA model ------------------------------------------------------
void CausalARMA(MultidimArray<double> &Img, ARMA_parameters &prm)
{
    // Calculate the autocorrelation matrix
    MultidimArray<double> R;
    auto_correlation_matrix(Img, R);
    R.setXmippOrigin();

    /**********************************************************************/
    // Set equation system for AR part of the model
    /**********************************************************************/
    Matrix2D<double> Coeficients;
    Matrix1D<double> Indep_terms, ARcoeficients;
    MultidimArray<double> N3;

    // Assign the support region for the AR part of the model (N1)
    First_Quadrant_Neighbors(prm.N_AR, prm.M_AR, prm.ARParameters);
    // Assign the support region for the MA part of the model (N2)
    Second_Quadrant_Neighbors(prm.N_MA, prm.M_MA, prm.MAParameters);
    // Assign the support region for the AR equations (N3)
    // Here is the same of N1, but it has not to be
    First_Quadrant_Neighbors(prm.N_AR, prm.M_AR, N3);

    long NumberOfARParameters = YSIZE(prm.ARParameters);
    long NumberOfMAParameters = YSIZE(prm.MAParameters);

    Coeficients.resizeNoCopy(NumberOfARParameters, NumberOfARParameters);
    Indep_terms.resizeNoCopy(NumberOfARParameters);
    ARcoeficients.resizeNoCopy(NumberOfARParameters);

    // Generate matrix (eq stands for equation number and co for coeficents)
    for (long eq = 0 ; eq < NumberOfARParameters; eq++)
    {
        // take the independet term from the correlation matrix (or calculate it
        // if it was not calculated before).
        int l = (int)A2D_ELEM(N3,eq, 0);
        int m = (int)A2D_ELEM(N3,eq, 1);
        VEC_ELEM(Indep_terms,eq) = A2D_ELEM(R, l, m);

        // take the coeficients
        for (long co = 0 ; co < NumberOfARParameters; co++)
        {
            // Take the pertinent coeficient form the correlation matrix (or calculate it)
            int alpha1 = (int)(A2D_ELEM(N3,eq, 0) - A2D_ELEM(prm.ARParameters,co, 0));
            int alpha2 = (int)(A2D_ELEM(N3,eq, 1) - A2D_ELEM(prm.ARParameters,co, 1));
            int beta1 = (int)(A2D_ELEM(N3,eq, 0) + A2D_ELEM(prm.ARParameters,co, 0));
            int beta2 = (int)(A2D_ELEM(N3,eq, 1) + A2D_ELEM(prm.ARParameters,co, 1));
            MAT_ELEM(Coeficients,eq, co) = A2D_ELEM(R,alpha1, alpha2) + A2D_ELEM(R,beta1, beta2);
        }
    }

    /**********************************************************************/
    // Solve the equation system to determine the AR model coeficients and sigma.
    /**********************************************************************/
    solveViaCholesky(Coeficients, Indep_terms, ARcoeficients);

    // Assign the results to the matrix given as parameter
    for (long n = 0 ; n < NumberOfARParameters; n++)
        A2D_ELEM(prm.ARParameters, n, 2) = VEC_ELEM(ARcoeficients,n);

    /**********************************************************************/
    // Determine the sigma coeficient from the equation for (p,q)=(0,0)
    /**********************************************************************/
    double dSum = 0;
    for (long n = 0 ; n < NumberOfARParameters; n++)
    {
        int p = (int)A2D_ELEM(prm.ARParameters,n, 0);
        int q = (int)A2D_ELEM(prm.ARParameters,n, 1);
        dSum += A2D_ELEM(prm.ARParameters, n, 2) * A2D_ELEM(R, p, q);
    }

    // And calculate sigma
    prm.dSigma = (A2D_ELEM(R, 0, 0) - 2 * dSum);
    double idSigma = 1.0 / prm.dSigma;

    /**********************************************************************/
    // Determine the MA parameters of the model using the AR parameters and
    // sigma
    /**********************************************************************/
    for (long n = 0 ; n < NumberOfMAParameters; n++)
    {
        dSum = 0;
        double MAn0 = A2D_ELEM(prm.MAParameters, n, 0);
        double MAn1 = A2D_ELEM(prm.MAParameters, n, 1);
        for (long m = 0 ; m < NumberOfARParameters; m++)
        {
            double ARm0 = A2D_ELEM(prm.ARParameters, m, 0);
            double ARm1 = A2D_ELEM(prm.ARParameters, m, 1);
            int alpha1 = (int)(MAn0 - ARm0);
            int alpha2 = (int)(MAn1 - ARm1);
            int beta1 = (int)(MAn0 + ARm0);
            int beta2 = (int)(MAn1 + ARm1);
            dSum += A2D_ELEM(prm.ARParameters, m, 2) * (
                        A2D_ELEM(R, alpha1, alpha2) + A2D_ELEM(R, beta1, beta2));
        }

        int p = (int)MAn0;
        int q = (int)MAn1;
        A2D_ELEM(prm.MAParameters, n, 2) = (A2D_ELEM(R, p, q) - dSum) * idSigma;
    }
}

// Compute the ARMA Filter -------------------------------------------------
void ARMAFilter(MultidimArray<double> &Img, MultidimArray< double > &Filter,
                ARMA_parameters &prm)
{
    bool apply_final_median_filter = false;
    Matrix1D<double> dDigitalFreq(2);

    // Resize the Filter to the image dimensions
    Filter.initZeros(YSIZE(Img), XSIZE(Img));

    // Compute the filter (only half the values are computed)
    // The other half is computed based in symmetry.
    int sizeX = XSIZE(Img);
    int sizeY = YSIZE(Img);
    long NumberOfMAParameters = YSIZE(prm.MAParameters);
    long NumberOfARParameters = YSIZE(prm.ARParameters);
    MultidimArray<int> iMAParameters(NumberOfMAParameters, 2),
    iARParameters(NumberOfARParameters, 2);
    for (long n = 0 ; n < NumberOfMAParameters; n++)
    {
        DIRECT_A2D_ELEM(iMAParameters, n, 0) = (int)DIRECT_A2D_ELEM(prm.MAParameters, n, 0);
        DIRECT_A2D_ELEM(iMAParameters, n, 1) = (int)DIRECT_A2D_ELEM(prm.MAParameters, n, 1);
    }
    for (long n = 0 ; n < NumberOfARParameters; n++)
    {
        DIRECT_A2D_ELEM(iARParameters, n, 0) = (int)DIRECT_A2D_ELEM(prm.ARParameters, n, 0);
        DIRECT_A2D_ELEM(iARParameters, n, 1) = (int)DIRECT_A2D_ELEM(prm.ARParameters, n, 1);
    }
    for (int i = 0;i < sizeY;i++)
        for (int j = 0;j < (sizeX / 2);j++)
        {
            // Compute dDigitalFreq
            XX(dDigitalFreq) = j / (double)sizeX;
            YY(dDigitalFreq) = i / (double)sizeY;

            // Compute B
            double B = 0;
            for (long n = 0 ; n < NumberOfMAParameters; n++)
            {
                int p = DIRECT_A2D_ELEM(iMAParameters, n, 0);
                int q = DIRECT_A2D_ELEM(iMAParameters, n, 1);
                B += DIRECT_A2D_ELEM(prm.MAParameters, n, 2) * 2 *
                     cos((-2) * PI * (p * YY(dDigitalFreq) + q * XX(dDigitalFreq)));
            }
            B = B + 1.0;

            // Compute A
            double A = 0;
            for (long n = 0 ; n < NumberOfARParameters; n++)
            {
                int p = DIRECT_A2D_ELEM(iARParameters, n, 0);
                int q = DIRECT_A2D_ELEM(iARParameters, n, 1);
                A += DIRECT_A2D_ELEM(prm.ARParameters, n, 2) * 2 *
                     cos((-2) * PI * (p * YY(dDigitalFreq) + q * XX(dDigitalFreq)));
            }
            A = 1.0 - A;

            // Check for posible problems calculating the value of A that could lead
            // to ARMA filters with erroneous values due to a value of A near 0.
            if (A < 0)
                apply_final_median_filter = true;

            // This is the filter proposed by original paper
            // As the filter values are symmetrical respect the center (frequency zero),
            // take advantage of this fact.
            double val2 = prm.dSigma * B / A;
            //if   (val2>=0) val=sqrt(val2);
            //else {val=abs(sqrt(std::complex<double>(val2,0)));}
            A2D_ELEM(Filter,sizeY - 1 - i, sizeX - 1 - j) = A2D_ELEM(Filter,i, j) = fabs(val2);
        }

    // Apply final median filter
    if (apply_final_median_filter)
    {
        MultidimArray<double> aux;
        medianFilter3x3(Filter, aux);
        // Copy all but the borders
        for (size_t i = 1; i < YSIZE(Filter) - 1; i++)
            for (size_t j = 1; j < XSIZE(Filter) - 1; j++)
                DIRECT_A2D_ELEM(Filter, i, j) = DIRECT_A2D_ELEM(aux, i, j);
    }
}
