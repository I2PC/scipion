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

#include <stdio.h>
#include <memory.h>

//#include <crtdbg.h>

#include "ObjectiveFunction.h"
#include "Matrix.h"
#include "IntPoly.h"
#include "tools.h"
#include "VectorChar.h"

Vector FullLambda;

// from QPsolver:
void simpleQPSolve(Matrix mH, Vector vG, Matrix mA, Vector vB,   // in
                   Vector vP, Vector vLambda, int *info);        // out
void restartSimpleQPSolve(Vector vBO,  // in
                          Vector vP);  // out

// from TRSSolver:
Vector L2NormMinimizer(Polynomial q, double delta, 
                       int *infoOut=NULL, int maxIter=1000, double *lambda1=NULL);
Vector L2NormMinimizer(Polynomial q, Vector pointXk, double delta, 
                       int *infoOut=NULL, int maxIter=1000, double *lambda1=NULL);
Vector L2NormMinimizer(Polynomial q, Vector pointXk, double delta, 
                       int *infoOut, int maxIter, double *lambda1, Vector minusG, Matrix H);


// from CTRSSolver: 
char checkForTermination(Vector d, Vector Base, double rhoEnd){return 0;}
void initConstrainedStep(ObjectiveFunction *of){ FullLambda.setSize(0); }

Vector ConstrainedL2NormMinimizer(InterPolynomial poly, int k, double delta, 
                                   int *info, int iterMax, double *lambda1, Vector vOBase, ObjectiveFunction *of)
{
    int dim=poly.dim();
    Matrix mH(dim,dim);
    Vector vG(dim);
    poly.gradientHessian(poly.NewtonPoints[k],vG,mH);

    if (of->isConstrained) 
	    printf("Limited version! Ignoring constraints !\n");
    return L2NormMinimizer(poly, poly.NewtonPoints[k], delta, info, iterMax, lambda1, vG, mH);

//    return ConstrainedL2NormMinimizer(mH,vG,delta,info,iterMax,lambda1,vOBase+poly.NewtonPoints[k],of);
}

void projectionIntoFeasibleSpace(Vector vFrom, Vector vBase,ObjectiveFunction *of)
{
    vBase.copyFrom(vFrom);
}















