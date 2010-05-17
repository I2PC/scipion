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

#ifdef WIN32
#include <crtdbg.h>
#endif

#include "ObjectiveFunction.h"
#include "Solver.h"
#include "Matrix.h"
#include "tools.h"
#include "KeepBests.h"
#include "IntPoly.h"
#include "parallel.h"
#include "MultInd.h"
#include "VectorChar.h"


#ifdef WIN32
void (__cdecl *quickHack)(InterPolynomial,int)=NULL;
#else
void (*quickHack)(InterPolynomial,int)=NULL;
#endif

// from CTRSSolver: 
Vector ConstrainedL2NormMinimizer(InterPolynomial poly, int k, 
                                  double delta, int *info, int iterMax, double *lambda1, 
                                  Vector vOBase, ObjectiveFunction *of);
char checkForTermination(Vector d, Vector Base, double rhoEnd);
void initConstrainedStep(ObjectiveFunction *of);

extern Vector FullLambda;

void fullUpdateOfM(double rho,Vector Base,Matrix data,InterPolynomial poly)
{
    int dim=data.nColumn()-2,i;
    Vector r(dim);
    for (i=0; i<data.nLine(); i++)
    {
        data.getLine(i,r,dim);
        if (r.euclidianDistance(Base)<1.5*rho)
        {
            r-=Base;
            poly.updateM(r, ((double**)data)[i][dim]);
        }
    }
}

void CONDOR( double rhoStart, double rhoEnd, int niter, 
              ObjectiveFunction *of, int nnode)
{
    rhoStart=mmax(rhoStart,rhoEnd);
    int dim=of->dim(), info, k, t, nerror;
    double rho=rhoStart, delta=rhoStart, rhoNew,
           lambda1, normD=rhoEnd+1.0, modelStep, reduction, r, valueOF, valueFk, bound, noise;
    Vector d, tmp;
    bool improvement, forceTRStep=true, evalNeeded;
    
    initConstrainedStep(of);

    // pre-create the MultInd indexes to prevent multi-thread problems:
    cacheMultInd.get(dim,1); cacheMultInd.get(dim,2);    

    parallelInit(nnode, dim, of);

    of->initData();
//    points=getFirstPoints(&ValuesF, &nPtsTotal, rhoStart,of);

    InterPolynomial poly(2, rho, Vector::emptyVector, of->data, of);
    if (poly.NewtonBasis==NULL)
    {
        printf("cannot construct lagrange basis.\n");
        exit(255);
    }

    
    
    //Base=poly.vBase;
    k=poly.kbest;
    valueFk=poly.valuesF[k];

    //fprintf(stderr,"init part 1 finished.\n");

/*    InterPolynomial poly(dim,2,of);
    for (;;)
    {
        // find index k of the best (lowest) value of the function 
        k=findK(ValuesF, nPtsTotal, of, points);
        Base=points[k].clone();
    //    Base=of->xStart.clone();
        valueFk=ValuesF[k];
        // translation:
        t=nPtsTotal; while (t--) points[t]-=Base;

        // exchange index 0 and index k (to be sure best point is inside poly):
        tmp=points[k];
        points[k]=points[0];
        points[0]=tmp;
        ValuesF[k]=ValuesF[0];
        ValuesF[0]=valueFk;
        k=0;

        poly=InterPolynomial(2, nPtsTotal, points, ValuesF);
        if (poly.NewtonBasis!=NULL) break;

        // the construction of the first polynomial has failed
        // delete[] points;
        free(ValuesF);

        int kbest=findBest(of->data, of); // return 0 if startpoint is given
        Vector BaseStart=of->data.getLine(kbest,dim);
        double vBaseStart=((double**)of->data)[kbest][dim];
        points=GenerateData(rho, BaseStart, vBaseStart, of, &ValuesF);
        nPtsTotal=n;
    }
*/
    // update M:
    fullUpdateOfM(rho,poly.vBase,of->data,poly);
    
    //fprintf(stderr,"init part 2 finished.\n");
    //fprintf(stderr,"init finished.\n");

    // first of init all variables:
    parallelImprove(&poly, &k, rho, &valueFk, poly.vBase);

    // really start in parallel:
    startParallelThread();

    while (true)
    {
//        fprintf(stderr,"rho=%e; fo=%e; NF=%i\n", rho,valueFk,QP_NF);
        while (true)
        {
            // trust region step
            while (true)
            {
//                poly.print();
                parallelImprove(&poly, &k, rho, &valueFk, poly.vBase);

                niter--;
                if ((niter==0)
                   ||(of->isConstrained&&checkForTermination(poly.NewtonPoints[k], poly.vBase, rhoEnd)))
                {
                    poly.vBase+=poly.NewtonPoints[k];
                    //fprintf(stderr,"rho=%e; fo=%e; NF=%i\n", rho,valueFk,of->getNFE());
                    of->valueBest=valueFk;
                    of->xBest=poly.vBase;
                    // to do : compute H and Lambda

                    Vector vG(dim);
                    Matrix mH(dim,dim);
                    poly.gradientHessian(poly.NewtonPoints[k],vG,mH);
                    of->finalize(vG,mH,FullLambda.clone());
                    return;
                }
   
                // to debug:
                //fprintf(stderr,"Best Value Objective=%e (nfe=%i)\n", valueFk, of->getNFE());
                
                d=ConstrainedL2NormMinimizer(poly,k,delta,&info,1000,&lambda1,poly.vBase,of);

//                if (d.euclidianNorm()>delta)
//                {
//                    printf("Warning d to long: (%e > %e)\n", d.euclidianNorm(), delta);
//                }

                normD=mmin(d.euclidianNorm(), delta);
                d+=poly.NewtonPoints[k];

//              next line is equivalent to reduction=valueFk-poly(d); 
//              BUT is more precise (no rounding error)
                reduction=-poly.shiftedEval(d,valueFk);

                //if (normD<0.5*rho) { evalNeeded=true; break; }
                if ((normD<0.5*rho)&&(!forceTRStep)) { evalNeeded=true; break; }

                //  IF THE MODEL REDUCTION IS SMALL, THEN WE DO NOT SAMPLE FUNCTION
                //  AT THE NEW POINT. WE THEN WILL TRY TO IMPROVE THE MODEL.

                noise=0.5*mmax(of->noiseAbsolute*(1+of->noiseRelative), condorAbs(valueFk)*of->noiseRelative);
                if ((reduction<noise)&&(!forceTRStep)) { evalNeeded=true; break; }
                forceTRStep=false; evalNeeded=false;

                if (quickHack) (*quickHack)(poly,k);
                tmp=poly.vBase+d; nerror=0; valueOF=of->eval(tmp, &nerror); 
                of->saveValue(tmp,valueOF,nerror);
                if (nerror)
                {
                    // evaluation failed
                    delta*=0.5;
                    if (normD>=2*rho) continue;
                    break;
                }
                if (!of->isFeasible(tmp, &r))
                {
                    printf("violation: %e\n",r);
                }

                // update of delta:
                r=(valueFk-valueOF)/reduction;
                if (r<=0.1) delta=0.5*normD;
                else if (r<0.7) delta=mmax(0.5*delta, normD);
                     else delta=mmax(rho+ normD, mmax(1.25*normD, delta));
            // powell's heuristics:
                if (delta<1.5*rho) delta=rho;

                if (valueOF<valueFk) 
                {
                    t=poly.findAGoodPointToReplace(-1, rho, d,&modelStep);
                    k=t; valueFk=valueOF; 
                    improvement=true;
//                    fprintf(stderr,"Value Objective=%e\n", valueOF);
                } else
                {
                    t=poly.findAGoodPointToReplace(k, rho, d,&modelStep);
                    improvement=false;
//                    fprintf(stderr,".");
                };

                if (t<0) { poly.updateM(d, valueOF); break; }

                // If we are along constraints, it's more important to update
                // the polynomial with points which increase its quality.
                // Thus, we will skip this update to use only points coming
                // from checkIfValidityIsInBound

                if ((!of->isConstrained)||(improvement)||(reduction>0.0)||(normD<rho)) poly.replace(t, d, valueOF); 

                if (improvement) continue;
//                if (modelStep>4*rho*rho) continue;
                if (modelStep>2*rho) continue;
                if (normD>=2*rho) continue;
                break;
            }
            // model improvement step
            forceTRStep=true;

//            fprintf(stderr,"improvement step\n");
            bound=0.0;
            if (normD<0.5*rho) 
            {
                bound=0.5*sqr(rho)*lambda1;
                if (poly.nUpdateOfM<10) bound=0.0;
            }

            parallelImprove(&poly, &k, rho, &valueFk, poly.vBase);

            // !! change d (if needed):
            t=poly.checkIfValidityIsInBound(d, k, bound, rho );
            if (t>=0)
            {
                if (quickHack) (*quickHack)(poly,k);
                tmp=poly.vBase+d; nerror=0; valueOF=of->eval(tmp, &nerror); 
                if (nerror) 
                {
                    Vector GXk(dim);
                    Matrix H(dim,dim);
                    poly.NewtonBasis[t].gradientHessian(poly.NewtonPoints[k],GXk,H);
                    double rhot=rho,vmax;

                    while (nerror)
                    {
                        rhot*=.5;
			            d=LAGMAXModified(GXk,H,rhot,vmax);
	                    d+=poly.NewtonPoints[k];
                        tmp=poly.vBase+d; nerror=0; valueOF=of->eval(tmp, &nerror); 
                        of->saveValue(tmp,valueOF,nerror);
                    }
                }
                poly.replace(t, d, valueOF); 
                if ((valueOF<valueFk)&&
                    (of->isFeasible(tmp))) { k=t; valueFk=valueOF; };
                continue;
            }

            // the model is perfect for this value of rho:
            // OR
            // we have crossed a non_linear constraint which prevent us to advance
            if ((normD<=rho)||(reduction<0.0)) break;
        }

        
        
        // change rho because no improvement can now be made:
        if (rho<=rhoEnd) break;

        //fprintf(stderr,"rho=%e; fo=%e; NF=%i\n", rho,valueFk,of->getNFE());

        if (rho<16*rhoEnd) rhoNew=rhoEnd;
        else if (rho<250*rhoEnd) rhoNew=sqrt(rho*rhoEnd);
             else rhoNew=0.1*rho;
        delta=mmax(0.5*rho,rhoNew);
        rho=rhoNew;

        
        // update of the polynomial: translation of x[k].
        // replace BASE by BASE+x[k]
        poly.translate(poly.NewtonPoints[k]);
    }
    parallelFinish();
    
    Vector vG(dim);
    Matrix mH(dim,dim);

    
    
    if (evalNeeded)
    {
        tmp=poly.vBase+d; nerror=0; valueOF=of->eval(tmp,&nerror); 
        of->saveValue(tmp,valueOF,nerror);
        if ((nerror)||(valueOF<valueFk))
        {
            poly.vBase+=poly.NewtonPoints[k];
            poly.gradientHessian(poly.NewtonPoints[k],vG,mH);
        }
        else
        {
            valueFk=valueOF; poly.vBase=tmp;
            poly.gradientHessian(d,vG,mH);
        }
    } else 
    {
        poly.vBase+=poly.NewtonPoints[k];
        poly.gradientHessian(poly.NewtonPoints[k],vG,mH);
    }


//    delete[] points; :not necessary: done in destructor of poly which is called automatically:
    //fprintf(stderr,"rho=%e; fo=%e; NF=%i\n", rho,valueFk,of->getNFE());
    
    of->valueBest=valueFk;
    of->xBest=poly.vBase;
    of->finalize(vG,mH,FullLambda.clone());
    
}

