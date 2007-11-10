/***************************************************************************
 *
 * Authors:     Javier Rodrguez Falces (jrodriguez@cnb.uam.es)
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

#ifndef INTEGRATION_H
#define INTEGRATION_H

// TODO remove this, NEVER open namespaces in header files
using namespace std;

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
#ifndef __INTEL_COMPILER
#include <external/cuba/cuba.h>
#endif
#include "matrix1d.h"

/** @defgroup NumericalIntegration Numerical integration
 *  @ingroup DataLibrary
 *
 * This code performs numeric integrations as described in the Numerical Recipes
 * Book, in particular it implements the "trapezoidal" (Trapeze) and the Romberg
 * integration. Both are designed for smoothly variant functions.
 *
 * This module can also perform multidimensional integration.
 */

/** Integrate a function using Newton-Cotes formula
 * @ingroup NumericalIntegration
 *
 * Estimate the integral of a function between a and b using N points.
 */
double integrateNewtonCotes(double(*f)(double), double a, double b, int N);

/** A double-returning function.
 * @ingroup NumericalIntegration
 *
 * Auxiliary class for numerical integration.
 */
class doubleFunction
{
public:
    virtual double operator()()=0; // pure virtual function
    virtual ~doubleFunction()
    {} // virtual destructor
};

/** Fast integration routine.
 * @ingroup NumericalIntegration
 *
 * Interpolations are made with lines.
 *
 * Example of use.
 *
 * 1) Define function to NumericalIntegration as class:
 * @code
 * // Actual function to be NumericalIntegration
 * class Func1: public doubleFunction
 * {
 * // This should be in testinteg
 * public:
 *     double x;
 *     double cte1, cte2;
 *     // overloads pure virtual
 *     virtual double operator()()
 *     {
 *         return sqrt(1 + cte1 * cte1 * sin(cte2 * x) * sin(cte2 * x));
 *     }
 * };
 * @endcode
 *
 * 2) In the main code
 *
 * @code
 * #include <data/integration.h>
 *
 * Func1 cosine; // cosine surface
 * cosine.cte1 = fabs(cteA * cteA * cteB * cteB * vx * vx);
 * cosine.cte2 = fabs(cteB * vx);
 *
 * Trapeze Trap(cosine, cosine.x, inte_low, inte_high);
 * integralt = Trap();
 * @endcode
 */
class Trapeze: public doubleFunction
{
    double s;
    doubleFunction& func; // the generic function to be NumericalIntegration
    double a, b; // integral limits
    double& x; // integration variable
    double EPS; // desired accuracy
    int JMAX; // 2**(JMAX) = max number of func. evaluation

public:
    // Overloads pure virtual function.
    virtual double operator()();

    /** With parameters.
     * Parameter: min Integration lower limit
     * Parameter: max Integration upper limit
     * Parameter: precision Maximum error allowed
     * Parameter: max_iter Maximum number of iterations
     */
    double operator()(double min, double max,
                      double precision = 1.0e-7, int max_iter = 20)
    {
        a = min;
        b = max;
        EPS = precision;
        JMAX = max_iter;
        return (*this)();
    }

    /** Constructor.
     *
     * Parameter: f Pointer to function to be integrated
     * Parameter: var Integration variable
     * Parameter: min Integration lower limit
     * Parameter: max Integration upper limit
     * Parameter: precision Maximum error allowed
     * Parameter: max_iter Maximum number of iterations
     */
    Trapeze(doubleFunction& f, double& Var, double min, double max,
            double precision = 1.0e-7, int max_iter = 20) : func(f), x(Var)
    {
        a = min;
        b = max;
        EPS = precision;
        JMAX = max_iter;
    }

    /** Workhorse that doublely does the integral.
     */
    double Trap(int n);
};

/** More accurate integration.
 * @ingroup NumericalIntegration
 *
 * More accurate integration than Trapeze with smaller truncation error
 * (interpolation is made with polynomials)
 *
 * Example of use:
 *
 * 1) Define function to NumericalIntegration as class:
 *
 * @code
 * // Actual function to be NumericalIntegration
 * class Func1: public doubleFunction
 * {
 * // This should be in testinteg
 * public:
 *     double x;
 *     double cte1, cte2;
 *     // Overloads pure virtual
 *     virtual double operator()()
 *     {
 *         return sqrt(1 + cte1 * cte1 * sin(cte2 * x) * sin(cte2 * x));
 *     }
 * };
 * @endcode
 *
 * 2) In the main code
 *
 * @code
 * #include <data/integration.h>
 *
 * Func1 cosine;   // cosine surface
 * cosine.cte1 = fabs(cteA * cteA * cteB * cteB * vx * vx);
 * cosine.cte2 = fabs(cteB * vx);
 *
 * Romberg Rom(cosine, cosine.x, inte_low, inte_high);
 * integralt = Rom();
 * @endcode
 *
 */
class Romberg : public doubleFunction
{
    double s;
    doubleFunction& func; // the function to be Numerical_interationd
    double a, b; // integral limits
    double& x; // integration variable
    double EPS; // desired accuracy

public:
    // TODO Document
    virtual double operator()();

    /** With parameters.
     *
     * Parameter: min Integration lower limit
     * Parameter: max Integration upper limit
     * Parameter: precision Maximum error allowed
     * Parameter: max_iter Maximum number of iterations
     */
    double operator()(double min, double max, double precision = 1.0e-7)
    {
        a = min;
        b = max;
        EPS = precision;
        return (*this)();
    }

    /** Constructor.
     *
     * Parameter: f Pointer to function to be integrated
     * Parameter: var Integration variable
     * Parameter: min Integration lower limit
     * Parameter: max Integration upper limit
     * Parameter: precision Maximum error allowed
     * Parameter: max_iter Maximum number of iterations
     */
    Romberg(doubleFunction& f, double& Var, double min, double max,
            double precision = 1.0e-7) : func(f), x(Var)
    {
        a = min;
        b = max;
        EPS = precision;
    }

    /** Workhorse that doublely does the integral
     */
    double midpnt(int n);
};

#ifndef __INTEL_COMPILER
/** Multidimensional integration
    @ingroup Numerical_interation
    Compute the double, triple, ... integral of a function. For instance

    @code
    integral_x0^xF  integral_y0^yF f(x,y) dx dy
    @endcode

    Our computation is based on the CUBA library (http://www.feynarts.de/cuba/).

    Example of use.

  1) Define the function to integrate:
@code
void integrand(const int *ndim, const double xx[],
               const int *ncomp, double ff[]) {
   #define f ff[0]
   #define x xx[0]
   #define y xx[1]
   f=x+y;
}
@endcode

    ndim is the number of dimensions of the integral. It is equal
    to the number of variables of f. ncomp is the number of components
    of f. The CUBA library is
    prepared for functions returning more than 1 value, but our
    function is not.

 2) In the main code

@code
#include <data/integration.hh>
cout << multidimensionalIntegral(vectorR2(0,0),vector_R2(10,10),
            &integrand) << endl;
@endcode
 */
double multidimensionalIntegral(const Matrix1D<double> &x0,
                                const Matrix1D<double> &xF, integrand_t integrand);
#endif

#endif
