/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

/* ------------------------------------------------------------------------- */
/* VECTORS                                                                   */
/* ------------------------------------------------------------------------- */
#include "matrix1d.h"

/* Vector R2 and R3 -------------------------------------------------------- */
Matrix1D<double> vectorR2(double x, double y)
{
    Matrix1D<double> result(2);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    return result;
}

Matrix1D<double> vectorR3(double x, double y, double z)
{
    Matrix1D<double> result(3);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    VEC_ELEM(result, 2) = z;
    return result;
}

Matrix1D<int> vectorR3(int x, int y, int z)
{
    Matrix1D<int> result(3);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    VEC_ELEM(result, 2) = z;
    return result;
}

/* Random permutation ------------------------------------------------------ */
void randomPermutation(int N, Matrix1D<int>& result)
{
    Matrix1D<double> aux;
    aux.resize(N);
    aux.initRandom(0,1);
    
    result=aux.indexSort();
    result-=1;
}

/* Powell's optimizer ------------------------------------------------------ */
void powellOptimizer(Matrix1D<double> &p, int i0, int n,
                      double(*f)(double *x, void *), void * prm,
                      double ftol, double &fret,
                      int &iter, const Matrix1D<double> &steps, bool show)
{
    double *xi = NULL;

    // Adapt indexes of p
    double *pptr = p.adaptForNumericalRecipes1D();
    double *auxpptr = pptr + (i0 - 1);

    // Form direction matrix
    ask_Tvector(xi, 1, n*n);
    int ptr;
    for (int i = 1, ptr = 1; i <= n; i++)
        for (int j = 1; j <= n; j++, ptr++)
            xi[ptr] = (i == j) ? steps(i - 1) : 0;

    // Optimize
    xi -= n; // This is because NR works with matrices starting at [1,1]
    powell(auxpptr, xi, n, ftol, iter, fret, f, prm, show);
    xi += n;

    // Exit
    free_Tvector(xi, 1, n*n);
}

/* Gaussian interpolator -------------------------------------------------- */
void GaussianInterpolator::initialize(double _xmax, int N, bool normalize)
{
    xmax=_xmax;
    xstep=xmax/N;
    ixstep=1.0/xstep;
    v.initZeros(N);
    double inorm=1.0/sqrt(2*PI);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
    {
        double x=i*xstep;
        v(i)=exp(-x*x/2);
        if (normalize) v(i)*=inorm;
    }
}
