/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Javier �ngel Vel�zquez Muriel (javi@cnb.csic.es)
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

#include <data/filters.h>

#include "spar.h"

/**************************************************************************

   NAME:          ComputeTermA

   DESCRIPTION:   This function computes term called A. (see ARFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
    ARParameters - Parameters of the AR model to perform the
                   calculus.

   OUTPUT:  The value of A

   DATE:        26-1-2001

/**************************************************************************/
double ComputeTermA(Matrix1D<double> &dDigitalFreq, Matrix2D<double> &ARParameters)
{
    double A = 0;


    for (int p = STARTINGY(ARParameters); p <= FINISHINGY(ARParameters); p++)
    {
        for (int q = STARTINGX(ARParameters); q <= FINISHINGX(ARParameters); q++)
        {
            // The term for (p,q)=(0,0) is not part of the AR model. It
            // contains sigma.
            if (!(p == 0 && q == 0))
            {
                A += MAT_ELEM(ARParameters, p, q) *
                     cos((-2) * PI * (p * YY(dDigitalFreq) + q * XX(dDigitalFreq)));
            }
        }
    }
    return A;
}


/**************************************************************************

   NAME:          ComputeTermB

   DESCRIPTION:   This function computes term called B. (see ARFilter function)

    PARAMETERS:   dDigitalFreq - Digital frecuency where to calculate the value
    ARParameters - Parameters of the AR model to perform the
                   calculus.

   OUTPUT:  The value of B

   DATE:        26-1-2001

/**************************************************************************/
double ComputeTermB(Matrix1D<double> &dDigitalFreq, Matrix2D<double> &ARParameters)
{
    double B = 0;

    for (int p = STARTINGY(ARParameters); p <= FINISHINGY(ARParameters); p++)
    {
        for (int q = STARTINGX(ARParameters); q <= FINISHINGX(ARParameters); q++)
        {
            // The term for (p,q)=(0,0) is not part of the AR model. It
            // contains sigma.
            if (!(p == 0 && q == 0))
            {
                B += MAT_ELEM(ARParameters, p, q) *
                     sin((-2) * PI * (p * YY(dDigitalFreq) + q * XX(dDigitalFreq)));
            }
        }
    }
    return B;
}



/**************************************************************************

   NAME:          CausalAR

   DESCRIPTION:   This function determines the coeficients of an 2D - AR model
                                   pf  qf
       Img(y,x)=(-1.0)*sum(sum( AR(p,q)*Img(y-p,x-q)) + sigma^2 * h(y,x)
                    p=0 q=q0
    (except for the case where p=0 and q=0)

        that adjust to the matrix provided as argument. To do
    that work, it solves the Yule-Walker equations for espectral
    estimation, involving the calculus of autocorrelation factors.
    It returns the AR model coeficients in a matrix ARParameters(p,q)
    ARParameters(0,0) contains the sigma^2 value, and it's also
    returned by the function.

    In this program the region of support considered can be the upper NSHP
    (NonSymmetric Half Plane) or lower NSHP.
    {p= 1, ...,pf; q=q0, ..., 0, ...,qF} U {p=0; q=0, ...,qF}
    (See Ultramicroscopy 68 (1997), pp. 276)

    For more details:
    Dudgeon "Multidimensional DSP",
    Prentice Hall, signal proc. series, pp. 325

    PARAMETERS:   Img - The matrix - Here it's supposed that it comes from
              an image
    ordeR, orderC - The order in Rows and Columns directions of the AR
                    model.
    ARParameters - The matrix to store the resulting parameteres
      for the AR model.

   OUTPUT:  The function stores the AR parameters into ARParameters
     value for (0,0) is sigma form the model.
  Sigma is also returned by the function.

   DATE:        19-1-2001

/**************************************************************************/
// #define DEBUG
double CausalAR(Matrix2D<double> &Img,
                int orderR, int orderC, Matrix2D<double> &ARParameters)
{

    int l0;  // initial Rows coeficient for correlations
    int lF;  // final Rows coeficient for correlations
    int m0;  // initial Columns coeficient for correlations
    int mF;  // final Columns coeficient for correlations
    int p0;  // initial Rows coeficient for AR parameters
    int pF;  // final Rows coeficient for AR parameters
    int q0;  // initial Columns coeficient for AR parameters
    int qF;  // final Columns coeficient for AR parameters

    // Choose region of the image that will affect the pixel considered (quadrants)
    // The region considered is the upper NSHP when orderR is greater than zero
    if (orderR >= 0)
    {
        l0 = 0;
        lF = orderR;
        m0 = -orderC;
        mF = orderC;
        p0 = 0;
        pF = orderR;
        q0 = -orderC;
        qF = orderC;
    }
    else
        // The region considered is the lower NSHP when orderR is greater than zero
    {
        l0 = orderR;
        lF = 0;
        m0 = -orderC;
        mF = orderC;
        p0 = orderR;
        pF = 0;
        q0 = -orderC;
        qF = orderC;
    }

    int eq, co; // auxiliary indexes for equation and coeficient

    // Compute correlation matrix (we'll name it R)
    Matrix2D<double> R((lF - p0) - (l0 - pF) + 1, (mF - q0) - (m0 - qF) + 1);
    R.initZeros();
    STARTINGY(R) = l0 - pF;
    STARTINGX(R) = m0 - qF;
    std::cerr << "Generating correlation coefficients ...\n";
    FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
    MAT_ELEM(R, i, j) = correlation(Img, Img, NULL, i, j);

    // Set equation system for AR model
    Matrix2D<double> Coeficients;
    Matrix1D<double> Indep_terms, ARcoeficients;

    Coeficients.resize((lF - l0 + 1)*(mF - m0 + 1) - mF, (lF - l0 + 1)*(mF - m0 + 1) - mF);

    STARTINGX(Coeficients) = 0;
    STARTINGY(Coeficients) = 0;

    Indep_terms.resize((lF - l0 + 1)*(mF - m0 + 1) - mF);
    STARTINGX(Indep_terms) = 0;

    ARcoeficients.resize((lF - l0 + 1)*(mF - m0 + 1) - mF);
    STARTINGX(ARcoeficients) = 0;

    // Generate matrix
    eq = 0; // equation number
    for (int l = lF; l >= l0; l--)
    {
        for (int m = mF; m >= m0; m--)
        {
            // This line is included to avoid points not in the NSHP.
            if (l == 0 && m != 0 && SGN(m) != SGN(orderR))
            {
                continue;
            }
            else
            {

                // take the independet terms from the correlation matrix
                Indep_terms(eq) = (-1.0) * R(l, m);

                co = 0; // coeficient number
                // take the coeficients
                for (int p = pF; p >= p0; p--)
                {
                    for (int q = qF; q >= q0; q--)
                    {
                        // This line is included to avoid points not in the NSHP.
                        if (p == 0 && q != 0 && SGN(q) != SGN(orderR))
                        {
                            continue;
                        }
                        else
                        {
                            // in the site for a(0,0) coeficient we put the sigma coeficient
                            if (p == 0 && q == 0)
                            {
                                // The coeficient for sigma, de std. dev. of the random process
                                // asociated with the AR model is determined here.
                                // It's supposed that the filter asociated to the random process
                                // has a response h(0,0)=1 and h(i,j)=0 elsewhere.
                                if (l == 0 && m == 0)
                                    MAT_ELEM(Coeficients, eq, co) = -1;
                                else
                                    MAT_ELEM(Coeficients, eq, co) = 0;
                            }
                            else
                            {
                                MAT_ELEM(Coeficients, eq, co) = R(l - p, m - q);
                            }
                            // increment the coeficient counter into an equation
                            co++;
                        }
                    }
                }

                // take the next equation
                eq++;
            }
        }
    }

    // Solve the equation system to determine the AR model coeficients and sigma.
    std::cerr << "Solving AR model ...\n";
    /*****************************************/
#ifdef DEBUG
    std::ofstream fichero("coeficientes.txt");
    fichero << Coeficients << std::endl << Indep_terms << std::endl ;
    fichero.close();
#endif
    /******************************************/

    solve(Coeficients, Indep_terms, ARcoeficients);

    // Put the ARcoeficients into the matrix given as parameter
    ARParameters.resize(pF - p0 + 1, qF - q0 + 1);
    STARTINGY(ARParameters) = p0;
    STARTINGX(ARParameters) = q0;

    co = 0;
    for (int p = pF; p >= p0; p--)
        for (int q = qF; q >= q0; q--)
        {
            // This line is to avoid points out the NSHP.
            if (p == 0 && q != 0 && SGN(q) != SGN(orderR))
            {
                continue;
            }
            else
            {
                MAT_ELEM(ARParameters, p, q) = VEC_ELEM(ARcoeficients, co);
                co++;
            }
        }

    // return the sigma coeficient
    return MAT_ELEM(ARParameters, 0, 0);
}
#undef DEBUG

/* Non causal AR ----------------------------------------------------------- */
double NonCausalAR(Matrix2D<double> &Img,
                   int orderR, int orderC, Matrix2D<double> &ARParameters)
{

    int l0;  // initial Rows coeficient for correlations
    int lF;  // final Rows coeficient for correlations
    int m0;  // initial Columns coeficient for correlations
    int mF;  // final Columns coeficient for correlations
    int p0;  // initial Rows coeficient for AR parameters
    int pF;  // final Rows coeficient for AR parameters
    int q0;  // initial Columns coeficient for AR parameters
    int qF;  // final Columns coeficient for AR parameters

    l0 = -orderR;
    lF = orderR;
    m0 = -orderC;
    mF = orderC;
    p0 = -orderR;
    pF = orderR;
    q0 = -orderC;
    qF = orderC;

    int eq, co; // auxiliary indexes for equation and coeficient

    // Compute correlation matrix (we'll name it R)
    Matrix2D<double> R((lF - p0) - (l0 - pF) + 1, (mF - q0) - (m0 - qF) + 1);
    R.initZeros();
    STARTINGY(R) = l0 - pF;
    STARTINGX(R) = m0 - qF;
    std::cerr << "Generating correlation coefficients ...\n";
    FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
    MAT_ELEM(R, i, j) = correlation(Img, Img, NULL, i, j);

    // Set equation system for AR model
    Matrix2D<double> Coeficients;
    Matrix1D<double> Indep_terms, ARcoeficients;

    Coeficients.resize((lF - l0 + 1)*(mF - m0 + 1), (lF - l0 + 1)*(mF - m0 + 1));

    STARTINGX(Coeficients) = 0;
    STARTINGY(Coeficients) = 0;

    Indep_terms.resize((lF - l0 + 1)*(mF - m0 + 1));
    STARTINGX(Indep_terms) = 0;

    ARcoeficients.resize((lF - l0 + 1)*(mF - m0 + 1));
    STARTINGX(ARcoeficients) = 0;

    // Generate matrix
    eq = 0; // equation number
    for (int l = lF; l >= l0; l--)
    {
        for (int m = mF; m >= m0; m--)
        {
            // take the independet terms from the correlation matrix
            Indep_terms(eq) = (-1.0) * R(l, m);

            co = 0; // coeficient number
            // take the coeficients
            for (int p = pF; p >= p0; p--)
            {
                for (int q = qF; q >= q0; q--)
                {
                    // in the site for a(0,0) coeficient we put the sigma coeficient
                    if (p == 0 && q == 0)
                    {
                        // The coeficient for sigma, de std. dev. of the random process
                        // It's supposed that the filter asociated to the random process
                        // has a response h(0,0)=1 and h(i,j)=0 elsewhere.
                        if (l == 0 && m == 0)
                        {
                            MAT_ELEM(Coeficients, eq, co) = -1;
                        }
                        else
                        {
                            MAT_ELEM(Coeficients, eq, co) = 0;
                        }
                    }
                    else
                    {
                        MAT_ELEM(Coeficients, eq, co) = R(l - p, m - q);
                    }
                    // increment the coeficient counter into an equation
                    co++;
                }
            }

            // take the next equation
            eq++;

        }
    }

    // Solve the equation system to determine the AR model coeficients and sigma.
    std::cerr << "Solving AR model ...\n";
    solve(Coeficients, Indep_terms, ARcoeficients);

    // Put the ARcoeficients into the matrix given as parameter
    ARParameters.resize(pF - p0 + 1, qF - q0 + 1);
    STARTINGY(ARParameters) = p0;
    STARTINGX(ARParameters) = q0;

    co = 0;
    for (int p = pF; p >= p0; p--)
        for (int q = qF; q >= q0; q--)
        {
            MAT_ELEM(ARParameters, p, q) = VEC_ELEM(ARcoeficients, co);
            co++;
        }

    // return the sigma coeficient
    return MAT_ELEM(ARParameters, 0, 0);
}

/* AR Filter --------------------------------------------------------------- */
#define DEBUG
void ARFilter(Matrix2D<double> &Img, Matrix2D< std::complex<double> > &Filter,
              Matrix2D<double> &ARParameters)
{

    double A, B;  /* Two Terms involved in calculation
                       of the filter defined by the AR model */

    Filter.resize(Img);

    /* Then, the Fourier Transform of the filter defined by the AR model
       is done */

    // Compute the filter
    Matrix1D<int>    iIndex(2);       // index in the Fourier image
    Matrix1D<double> dDigitalFreq(2); // digital frequency corresponding to and
    // index

#ifdef DEBUG
    std::ofstream filelog("coeficientes.txt");
#endif
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Filter)
    {
        // Compute dDigitalFreq
        XX(dDigitalFreq) = j / XSIZE(Filter);
        YY(dDigitalFreq) = i / YSIZE(Filter);

        // Compute terms A and B for Filter
        A = ComputeTermA(dDigitalFreq, ARParameters);
        B = ComputeTermB(dDigitalFreq, ARParameters);

        double sigma = sqrt(MAT_ELEM(ARParameters, 0, 0));
        Filter(i, j) = std::complex<double> (sigma * (1 + A) / ((1 + A) * (1 + A) + B * B),
                                             sigma * (-B) / ((1 + A) * (1 + A) + B * B));

#ifdef DEBUG
        filelog << "A " << A << " B " << B << " po " << abs(Filter(i, j)) << std::endl ;
#endif
    }

#ifdef DEBUG
    filelog.close();
#endif
}
#undef DEBUG

/* Combine AR filters ------------------------------------------------------ */
void combineARFilters(const Matrix2D< std::complex<double> > &Filter1,
                      const Matrix2D< std::complex<double> > &Filter2,
                      Matrix2D< std::complex<double> > &Filter,
                      const std::string &method)
{
    Filter.resize(Filter1);

    int imethod;
    if (method == "arithmetic_mean")  imethod = 0;
    else if (method == "armonic_mean")     imethod = 1;
    else if (method == "geometric_mean")   imethod = 2;
    else if (method == "armonic_spectrum") imethod = 3;
    else                                 imethod = 0;

    double abs1, abs2;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Filter)
    {
        switch (imethod)
        {
        case 0:
            Filter(i, j) = 0.5 * (Filter1(i, j) + Filter2(i, j));
            break;
        case 1:
            Filter(i, j) = 2.0 / (1.0 / Filter1(i, j) + 1.0 / Filter2(i, j));
            break;
        case 2:
            Filter(i, j) = sqrt(Filter1(i, j) * Filter2(i, j));
            break;
        case 3:
            abs1 = abs(Filter1(i, j));
            abs2 = abs(Filter2(i, j));
            Filter(i, j) = sqrt(2.0 / (1.0 / (abs1 * abs1) + 1.0 / (abs2 * abs2)));
            break;
        }
    }
}
