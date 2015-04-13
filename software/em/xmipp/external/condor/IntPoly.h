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

#include "Poly.h"
#include "Vector.h"
#include "ObjectiveFunction.h"

#ifndef _INTPOLY_H_
#define _INTPOLY_H_


class InterPolynomial : public Polynomial 
{
public:
    
    double M;
    unsigned nPtsUsed, nUpdateOfM;

    // (*this) = sum_i newBasis[i]*NewtonCoefPoly[i]
	Polynomial *NewtonBasis;
    // double *NewtonCoefPoly;

    // data:
    Vector *NewtonPoints,vBase;
    double *valuesF;
    int kbest;

    double *NewtonCoefficient(double *);
    void ComputeLagrangeBasis(double *, unsigned nPtsTotal);
    void GenerateBasis(double rho,double rhosmall, Matrix data,ObjectiveFunction *of);

/*
    InterPolynomial() : Polynomial() {}
    InterPolynomial( const Polynomial& p ) : Polynomial( p ) {};
    InterPolynomial( const InterPolynomial& p ) : Polynomial( p ) {};
*/
    
//    InterPolynomial( unsigned _deg, unsigned nPtsTotal, Vector *_Pp, double *_Yp );
    InterPolynomial( unsigned _deg, double rho, Vector vBase, Matrix data, ObjectiveFunction *of);

    int findAGoodPointToReplace(int excludeFromT,double rho, 
                                                    Vector pointToAdd, double *modelStep=NULL);
    void replace(int t, Vector pointToAdd, double valueF);
    int maybeAdd(Vector pointToAdd, unsigned k, double rho, double valueF);

    void updateM(Vector newPoint, double valueF);
    int checkIfValidityIsInBound(Vector dd, unsigned k, double bound, double rho);
    int getGoodInterPolationSites(Matrix d, int k, double rho, Vector *v=NULL);
    double interpError(Vector Point);

    void translate(int k);
    void translate(Vector translation);

//    void test();
//    void check(Vector Base, double (*f)(  Vector ) );

        // allow shallow copy:
    ~InterPolynomial();
    InterPolynomial(const InterPolynomial &A);
    InterPolynomial& operator=( const InterPolynomial& A );
    InterPolynomial clone();
    void copyFrom(InterPolynomial a);
    void copyFrom(Polynomial a);
    InterPolynomial(unsigned dim, unsigned deg); 

protected:
    void destroyCurrentBuffer();
    
};

#ifndef NOOBJECTIVEFUNCTION

#endif

#endif 	/* _MPI_INTPOLY_H_ */
