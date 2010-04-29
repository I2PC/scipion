/*

CONDOR 1.06 - COnstrained, Non-linear, Direct, parallel Optimization
              using trust Region method for high-computing load,
              noisy functions
Copyright (C) 2004 Frank Vanden Berghen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation version 2
of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

If you want to include this tools in any commercial product,
you can contact the author at fvandenb@iridia.ulb.ac.be

*/

#ifndef _INCLUDE_POLY_H_
#define _INCLUDE_POLY_H_

#include "MultInd.h"
#include "Vector.h"
//#include "tools.h"
#include "Vector.h"
#include "Matrix.h"

// ORDER BY DEGREE !
class Polynomial
{
protected:

    typedef struct PolynomialDataTag
    {
        double *coeff;  // Coefficients
        unsigned n,     // size of vector of Coefficients
                 dim,   // Dimensions
                 deg;   // Degree
        int ref_count;
    } PolynomialData;
    PolynomialData *d;
    void init(int _dim, int _deg, double *data=NULL);
    void destroyCurrentBuffer();

public:
    Polynomial()
    {
        init(0,0);
    };
    Polynomial( unsigned Dim, unsigned deg=0, double *data=0 );
    Polynomial( unsigned Dim, double val ); // Constant polynomial
    Polynomial( MultInd& );  // Monomials
    Polynomial(char *name);


    // Accessor
    inline unsigned dim()
    {
        return d->dim;
    };
    inline unsigned deg()
    {
        return d->deg;
    };
    inline unsigned sz()
    {
        return d->n;
    };
    inline operator double*() const
    {
        return d->coeff;
    };

    // allow shallow copy:
    ~Polynomial();
    Polynomial(const Polynomial &A);
    Polynomial& operator=( const Polynomial& A );
    Polynomial clone();
    void copyFrom(Polynomial a);

    // Arithmetic operations

    // friend Polynomial operator*( const double&, const Polynomial& );
    Polynomial operator*( const double );
    Polynomial operator/( const double );
    Polynomial operator+( Polynomial );
    Polynomial operator-( Polynomial );

    // Unary
    Polynomial operator-( void ); // the opposite (negative of)
    Polynomial operator+( void )
    {
        return *this;
    }

    // Assignment+Arithmetics

    Polynomial operator+=( Polynomial );
    Polynomial operator-=( Polynomial );
    Polynomial operator*=( const double );
    Polynomial operator/=( const double );

    // simple math tools

//    double simpleEval( Vector P);
    double shiftedEval( Vector Point, double minusVal);
    double operator()( Vector );
    Polynomial derivate(int i);
    void gradient(Vector P, Vector G);
    void gradientHessian(Vector P, Vector G, Matrix H);
    void translate(Vector translation);

    // Comparison

    inline int operator==( const Polynomial q)
    {
        return d==q.d;
    };
    int equals( Polynomial q );

    // Output

    void print();
    void save(char *name);

    //ostream& PrintToStream( ostream& ) const;

    //behaviour
    static const unsigned int NicePrint;
    static const unsigned int Warning;
    static const unsigned int Normalized;  // Use normalized monomials

    static unsigned int flags;
    void       setFlag( unsigned int val )
    {
        flags |= val;
    }
    void     unsetFlag( unsigned int val )
    {
        flags &= ~val;
    }
    unsigned queryFlag( unsigned int val )
    {
        return flags & val;
    }

    static Polynomial emptyPolynomial;
};

unsigned long choose( unsigned n, unsigned k );

// operator * defined on double:
inline Polynomial operator*( const double& dou, Polynomial& p  )
{
    // we can use operator * defined on Polynomial because of commutativity
    return p * dou;
}

#endif  /* _MPI_POLY_H_ */
