/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

/* ------------------------------------------------------------------------- */
/* VECTORS                                                                   */
/* ------------------------------------------------------------------------- */
#include "matrix1d.h"

/* ************************************************************************* */
/* IMPLEMENTATIONS                                                           */
/* ************************************************************************* */
#define maT matrix1D<T>
#define ma  matrix1D
#include "multidim_basic.inc"
#undef ma
#undef maT
// Special case for complex numbers
template <>
ostream& operator << (ostream& out, const matrix1D< complex<double> > & v)
{
    if (MULTIDIM_SIZE(v) == 0)
        out << "NULL vector\n";
    else
    {
        FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
        {
            if (v.row)
                out << VEC_ELEM(v, i) << ' ';
            else
                out << VEC_ELEM(v, i) << '\n';
        }
    }
    return out;
}

/* Vector R2 and R3 -------------------------------------------------------- */
matrix1D<double> vector_R2(double x, double y)
{
    matrix1D<double> result(2);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    return result;
}

matrix1D<double> vector_R3(double x, double y, double z)
{
    matrix1D<double> result(3);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    VEC_ELEM(result, 2) = z;
    return result;
}

matrix1D<int> vector_R3(int x, int y, int z)
{
    matrix1D<int> result(3);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    VEC_ELEM(result, 2) = z;
    return result;
}

/* Are orthogonal ---------------------------------------------------------- */
int are_orthogonal(matrix1D<double> &v1, matrix1D<double> &v2,
                   matrix1D<double> &v3)
{
    if (XSIZE(v1) != 3 || XSIZE(v2) != 3 || XSIZE(v3) != 3)
        REPORT_ERROR(1002, "Are orthogonal: Some vector do not belong to R3");
    try
    {
        if (dot_product(v1, v2) != 0)
            return 0;
        if (dot_product(v2, v3) != 0)
            return 0;
        if (dot_product(v1, v3) != 0)
            return 0;
    }
    catch (Xmipp_error)
    {
        REPORT_ERROR(1007, "Are orthogonal: Vectors are not all of the same shape");
    }
    return 1;
}

/* Are system? ------------------------------------------------------------- */
int are_system(matrix1D<double> &v1, matrix1D<double> &v2,
               matrix1D<double> &v3)
{
    matrix1D<double> aux(3);
    if (XSIZE(v1) != 3 || XSIZE(v2) != 3 || XSIZE(v3) != 3)
        REPORT_ERROR(1002, "Are orthogonal: Some vector do not belong to R3");
    aux = vector_product(v1, v2);
    if (aux != v3)
        return 0;
    aux = vector_product(v2, v3);
    if (aux != v1)
        return 0;
    aux = vector_product(v3, v1);
    if (aux != v2)
        return 0;
    return 1;
}

/* Powell's optimizer ------------------------------------------------------ */
void Powell_optimizer(matrix1D<double> &p, int i0, int n,
                      double(*f)(double *x), double ftol, double &fret,
                      int &iter, const matrix1D<double> &steps, bool show)
{
    double *xi = NULL;

    // Adapt indexes of p
    double *pptr = p.adapt_for_numerical_recipes();
    double *auxpptr = pptr + (i0 - 1);

    // Form direction matrix
    ask_Tvector(xi, 1, n*n);
    int ptr;
    for (int i = 1, ptr = 1; i <= n; i++)
        for (int j = 1; j <= n; j++, ptr++)
            xi[ptr] = (i == j) ? steps(i - 1) : 0;

    // Optimize
    xi -= n; // This is because NR works with matrices
    // starting at [1,1]
    powell(auxpptr, xi, n, ftol, iter, fret, f, show);
    xi += n;

    // Exit
    free_Tvector(xi, 1, n*n);
}
