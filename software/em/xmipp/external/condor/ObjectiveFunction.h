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

#ifndef OBJECTIVEFUNCTION_INCLUDE
#define OBJECTIVEFUNCTION_INCLUDE

#include <stdio.h>
#include "Vector.h"
#include "Matrix.h"

#define INF 1.7E+308 

class ObjectiveFunction
{
  public:
    friend class CorrectScaleOF;

    char name[9], startPointIsGiven;
    Vector xStart, xBest, xOptimal;
    // xOptimal is the theoretical exact solution of the optimization problem 
    //          if this value is set, it is used by printstat. It's an optional value to set.
    // xBest is the solution given by the optimization algorithm.
    double valueOptimal, valueBest, noiseAbsolute, noiseRelative, objectiveConst;
    // valueOptimal is the value of the obj.funct. at the theoretical exact solution of the optimization problem.
    // valueBest is the value of the obj.funct. at the solution given by the optimization algorithm.
    // objectiveConst is used inside method "printStats" to give a correction to the evaluation of the obj.funct.
    // noiseAbsolute, noiseRelative are the absolute and relative noise on the evaluation of the obj.funct.
    //                              These are optional values to set. If not given, the optimizer will assumea value
    //                              of zero.
    Matrix data;
    int t,nNLConstraints, isConstrained;
    // nfe: total number of function evalution
    // nfe2: number of function evalution needed to reach xBest without knowing it's the optimal solution
    // t: type of the OF (used inside "getObjectiveFunction")
    // nNLConstraints: number of non-linear constraints

    // CONSTRAINTS:
    // for lower/upper bounds (box constraints)
    Vector bl, bu;
    
	// for linear constraints Ax>=b
	Matrix A;
	Vector b;
    
    // for non-linear constraints
    virtual double evalNLConstraint(int j, Vector v, int *nerror)=0;
    virtual Vector evalGradNLConstraint(int j, Vector v, int *nerror);
    virtual void evalGradNLConstraint(int j, Vector v, Vector result, int *nerror)=0;

    // tolerances for constraints (relative and absolute)
    double tolRelFeasibilityForNLC, tolNLC;
    double tolRelFeasibilityForLC, tolLC;

    ObjectiveFunction() : startPointIsGiven(0), valueOptimal(INF), valueBest(INF), noiseAbsolute(0.0), 
           noiseRelative(0.0), objectiveConst(0.0), nNLConstraints(0), 
           isConstrained(1), tolRelFeasibilityForNLC(1e-9), tolNLC(1e-6), 
           tolRelFeasibilityForLC(1e-6), tolLC(1e-8), 
           saveFileName(NULL), dfold(INF), maxNormLC(0.0), maxNormNLC(0.0), nfe(0), nfe2(0) { };
    virtual ~ObjectiveFunction(){ if (saveFileName) free(saveFileName); };
    virtual double eval(Vector v, int *nerror)=0;
    int dim();
    void initData();
    virtual void saveValue(Vector tmp,double valueOF, int nerror);
    virtual void printStats(char cc=1);
    void saveStats(char *filename, Vector vG, Matrix mH, Vector vLambda);
    virtual void finalize(Vector vG, Matrix mH, Vector vLambda){};
    void setName(char *s);
    void setSaveFile(char *b=NULL);
    void updateCounter(double df, Vector vX, int nerror=0);
    char isFeasible(Vector vx, double *d=NULL);
    void initBounds();
    void endInit();
    void initTolLC(Vector vX);
    void initTolNLC(Vector c, double delta);
    virtual int getNFE() { return nfe; }
    virtual int getNFE2() { return nfe2; }

private:
    char *saveFileName;
    double dfold, dfref, maxNormLC, maxNormNLC;
    void addClosestFeasiblePointInData(Vector vX);
protected:
    int nfe,nfe2;
};

class UnconstrainedObjectiveFunction : public ObjectiveFunction
{
  public:
    UnconstrainedObjectiveFunction(): ObjectiveFunction(){ isConstrained=0; }
    ~UnconstrainedObjectiveFunction() {};

    virtual double evalNLConstraint(int j, Vector v, int *nerror=NULL){ return 0; };
    virtual Vector evalGradNLConstraint(int j, Vector v, int *nerror=NULL){ return Vector::emptyVector; };
    virtual void evalGradNLConstraint(int j, Vector v, Vector result, int *nerror=NULL) { result=Vector::emptyVector; };
};


class CorrectScaleOF: public ObjectiveFunction
{
  public:
    Vector rescaling;
    CorrectScaleOF(int _t, ObjectiveFunction *_of, Vector _rescaling);
    CorrectScaleOF(int _t, ObjectiveFunction *_of);
    ~CorrectScaleOF(){};
    void saveValue(Vector X,double valueOF,int nerror);
    double eval(Vector v, int *nerror=NULL);
    virtual double evalNLConstraint(int j, Vector v, int *nerror=NULL);
    virtual void evalGradNLConstraint(int j, Vector v, Vector result, int *nerror=NULL);
    virtual void finalize(Vector vG, Matrix mH, Vector vLambda);
    virtual int getNFE() { return of->nfe; }
    virtual int getNFE2() { return of->nfe2; }
  private:
    void init();
    ObjectiveFunction *of;
    Vector xTemp;
};

#endif

