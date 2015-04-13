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

#include "Matrix.h"
#include "tools.h"
#include "VectorChar.h"


// #define POWELLQP
#ifdef POWELLQP

int ql0001_(int *m,int *me,int *mmax,int *n,int *nmax,int *mnn,
            double *c,double *d,double *a,double *b,double *xl,
            double *xu,double *x,double *u,int *iout,int *ifail,
            int *iprint,double *war,int *lwar,int *iwar,int *liwar,
            double *eps1);

void simpleQPSolve(Matrix G, Vector vd, Matrix Astar, Vector vBstar,   // in
                   Vector vXtmp, Vector vLambda)                       // out
{
    Astar.setNLine(Astar.nLine()+1);
    Matrix thisG, At=Astar.transpose();
    Astar.setNLine(Astar.nLine()-1);
    Vector minusvd=vd.clone();
    minusvd.multiply(-1.0); vBstar.multiply(-1.0);
    int m=Astar.nLine(), me=0, mmax=Astar.nLine()+1, n=vd.sz(), nmax=n,
        mnn=m+n+n, iout=0, ifail, iprint=0, lwar=3*nmax*nmax/2 + 10*nmax+ 2*mmax+1,
        liwar=n;
    if (G==Matrix::emptyMatrix) { 
        thisG.setSize(n,n); thisG.diagonal(1.0); }
    else thisG=G;
    double *c=*((double**)thisG), *d=minusvd, *a=*((double**)At), *b=vBstar; 
    Vector vxl(n),vxu(n), temp(lwar);
    VectorInt itemp(liwar);
    vLambda.setSize(mnn);
    double *xl=vxl, *xu=vxu, *x=vXtmp, *u=vLambda, *war=temp, eps1=1e-20;
    int *iwar=itemp;

    int dim=n; while (dim--) { xl[dim]=-INF; xu[dim]=INF; }
    iwar[0]=0;

    ql0001_(&m,&me,&mmax,&n,&nmax,&mnn,c,d,a,b,xl,xu,x,u,&iout,&ifail,
            &iprint,war,&lwar,iwar,&liwar,&eps1);

    vLambda.setSize(m);
}

#else

Matrix mQ_QP;
MatrixTriangle mR_QP;
VectorInt vi_QP, viPermut_QP;
Vector vLastLambda_QP;
char bQRFailed_QP;

void QPReconstructLambda(Vector vLambda, Vector vLambdaScale)
{
    int n=vi_QP.sz()-1, nc=vLambda.sz();
    int *ii=vi_QP;
    double *l=vLambda, *s=vLambdaScale;

    if (n>=0)
        while (nc--)
        {
            if (ii[n]==nc) 
            { 
//                l[nc]=mmax(0.0,l[n]/s[n]); 
                l[nc]=l[n]/s[n]; 
                n--; 
                if (n<0) break; 
            }
            else l[nc]=0.0;
        }
    if (nc) memset(l,0,nc*sizeof(double));
}

void simpleQPSolve(Matrix mH, Vector vG, Matrix mA, Vector vB,   // in
                   Vector vP, Vector vLambda, int *info)         // out
{
    const double tolRelFeasibility=1e-6;
//     int *info=NULL;
    int dim=mA.nColumn(), nc=mA.nLine(), ncr, i,j,k, lastJ=-2, *ii;
    MatrixTriangle M(dim);
    Matrix mAR, mZ(dim,dim-1), mHZ(dim,dim-1), mZHZ(dim-1,dim-1);
    Vector vTmp(mmax(nc,dim)), vYB(dim), vD(dim), vTmp2(dim), vTmp3(nc), vLast(dim), 
           vBR_QP, vLambdaScale(nc);
    VectorChar vLB(nc);
    double *lambda=vLambda, *br, *b=vB, violationMax, violationMax2,
           *rlambda,ax,al,dviolation, **a=mA, **ar=mAR, mymin, mymax, maxb, *llambda,
           **r,t, *scaleLambda=vLambdaScale;
    char finished=0, feasible=0, *lb=vLB;

    if (info) *info=0;

    // remove lines which are null
    k=0; vi_QP.setSize(nc); maxb=0.0;
    for (i=0; i<nc; i++)
    {
        mymin=INF; mymax=-INF;
        for (j=0; j<dim; j++)
        {
            mymin=mmin(mymin,a[i][j]);
            mymax=mmax(mymax,a[i][j]);
            if ((mymin!=mymax)||(mymin!=0.0)||(mymax!=0.0)) break;
        }
        if ((mymin!=mymax)||(mymin!=0.0)||(mymax!=0.0))
        {
            lambda[k]=lambda[i];
            maxb=mmax(maxb,condorAbs(b[i]));
            scaleLambda[k]=mA.euclidianNorm(i);
            vi_QP[k]=i;
            k++;
        }
    }
    nc=k; vi_QP.setSize(nc); ii=vi_QP; maxb=(1.0+maxb)*tolRelFeasibility;
    vLast.zero(); 
    
    for (i=0; i<nc; i++) if (lambda[i]!=0.0) lb[i]=2; else lb[i]=1;

    while (!finished)
    {
        finished=1;
        mAR.setSize(dim,dim); mAR.zero(); ar=mAR;
        vBR_QP.setSize(mmin(nc,dim)); br=vBR_QP;
        ncr=0;
        for (i=0; i<nc; i++)
            if (lambda[i]!=0.0)
            {
//                mAR.setLines(ncr,mA,ii[i],1);
                k=ii[i];
                t=scaleLambda[ncr];
                for (j=0; j<dim; j++) ar[ncr][j]=a[k][j]*t;
                br[ncr]=b[ii[i]]*t;
                ncr++;
            }
        mAR.setSize(ncr,dim);
        vBR_QP.setSize(ncr);
        vLastLambda_QP.copyFrom(vLambda); llambda=vLastLambda_QP;

        if (ncr==0)
        {
            // compute step
            vYB.copyFrom(vG);
            vYB.multiply(-1.0);
            if (mH.cholesky(M))
            {
                M.solveInPlace(vYB);
                M.solveTransposInPlace(vYB);
            } else
            {
                printf("warning: cholesky factorisation failed.\n");
                if (info) *info=2;
            }
            vLambda.zero();
        } else
        {
            Matrix mAR2=mAR.clone(), mQ2;
            MatrixTriangle mR2;
            mAR2.QR(mQ2,mR2);

            mAR.QR(mQ_QP,mR_QP, viPermut_QP); // content of mAR is destroyed here !
            bQRFailed_QP=0;

            r=mR_QP;
            for (i=0; i<ncr; i++)
                if (r[i][i]==0.0)
                {
                    // one constraint has been wrongly added.
                    bQRFailed_QP=1;
                    QPReconstructLambda(vLambda,vLambdaScale); vP.copyFrom(vLast); return; 
                }

            for (i=0; i<ncr; i++)
                if (viPermut_QP[i]!=i)
                {
                  //  printf("whoups.\n");
                }
            if (ncr<dim)
            {
                mQ_QP.getSubMatrix(mZ,0,ncr);

                // Z^t H Z
                mH.multiply(mHZ,mZ);
                mZ.transposeAndMultiply(mZHZ,mHZ);
                mQ_QP.setSize(dim,ncr);            
            }

            // form Yb
            vBR_QP.permutIn(vTmp,viPermut_QP);
            mR_QP.solveInPlace(vTmp);
            mQ_QP.multiply(vYB,vTmp);

            if (ncr<dim)
            {

                // ( vG + H vYB )^t Z : result in vD
                
                mH.multiply(vTmp,vYB);
                vTmp+=vG;
                vTmp.transposeAndMultiply(vD,mZ);

                // calculate current delta (result in vD)
                vD.multiply(-1.0);
                if (mZHZ.cholesky(M))
                {
                    M.solveInPlace(vD);
                    M.solveTransposInPlace(vD);
                }
                else
                {
                    printf("warning: cholesky factorisation failed.\n");
                    if (info) *info=2;
                };

                // evaluate vX* (result in vYB):
                mZ.multiply(vTmp, vD);
                vYB+=vTmp;
            }

            // evaluate vG* (result in vTmp2)
            mH.multiply(vTmp2,vYB);
            vTmp2+=vG;

            // evaluate lambda star (result in vTmp):
            mQ2.transposeAndMultiply(vTmp,vTmp2);
            mR2.solveTransposInPlace(vTmp);

            // evaluate lambda star (result in vTmp):
            mQ_QP.transposeAndMultiply(vTmp3,vTmp2);
            mR_QP.solveTransposInPlace(vTmp3);
            vTmp3.permutOut(vTmp,viPermut_QP);
            rlambda=vTmp;

            ncr=0;
            for (i=0; i<nc; i++)
                if (lambda[i]!=0.0)
                {
                    lambda[i]=rlambda[ncr];
                    ncr++;
                }
        } // end of test on ncr==0

        // find the most violated constraint j among non-active Linear constraints:
        j=-1;
        if (nc>0)
        {
            k=-1; violationMax=-INF; violationMax2=-INF;
            for (i=0; i<nc; i++)
            {
                if (lambda[i]<=0.0) // test to see if this constraint is not active 
                {
                    ax=mA.scalarProduct(ii[i],vYB);
                    dviolation=b[ii[i]]-ax;
                    if (llambda[i]==0.0)
                    {
                        // the constraint was not active this round
                        // thus, it can enter the competition for the next 
                        // active constraint

                        if (dviolation>maxb)
                        {
                            // the constraint should be activated
                            if (dviolation>violationMax2) 
                                { k=i; violationMax2=dviolation; }
                            al=mA.scalarProduct(ii[i],vLast)-ax; 
                            if (al>0.0) // test to see if we are going closer
                            {  
                                dviolation/=al;
                                if (dviolation>violationMax ) 
                                    { j=i; violationMax =dviolation; }
                            }
                        }
                    } else
                    {
                        lb[i]--;
                        if (feasible) 
                        {
                            if (lb[i]==0) 
                            {
                                vLambda.copyFrom(vLastLambda_QP);
                                if (lastJ>=0) lambda[lastJ]=0.0;
                                QPReconstructLambda(vLambda,vLambdaScale); 
                                vP.copyFrom(vYB); 
                                return; 
                            }
                        } else
                        {
                            if (lb[i]==0)
                            {
                                if (info) *info=1;
                                QPReconstructLambda(vLambda,vLambdaScale); 
                                vP.zero(); 
                                return;
                            }
                        }
                        finished=0;  // this constraint was wrongly activated.
                        lambda[i]=0.0;
                    }
                }
            }

            // !!! the order the tests is important here !!!
            if ((j==-1)&&(!feasible)) 
            { 
                feasible=1; 
                for (i=0; i<nc; i++) if (llambda[i]!=0.0) lb[i]=2; else lb[i]=1;
            }
            if (j==-1) { j=k; violationMax=violationMax2; } // change j to k after feasible is set
            if (ncr==mmin(dim,nc)) { 
                if (feasible) { // feasible must have been checked before
                    QPReconstructLambda(vLambda,vLambdaScale); vP.copyFrom(vYB); return; 
                } else
                {
                    if (info) *info=1;
                    QPReconstructLambda(vLambda,vLambdaScale); vP.zero(); return; 
                }
            }
            // activation of constraint only if ncr<mmin(dim,nc)
            if (j>=0) { lambda[j]=1e-5; finished=0; } // we need to activate a new constraint 
//            else if (ncr==dim){ 
//                QPReconstructLambda(vLambda); vP.copyFrom(vYB); return; }
        }

        // to prevent rounding error
        if (j==lastJ) { 
//        if (0) {
            QPReconstructLambda(vLambda,vLambdaScale); vP.copyFrom(vYB); return; }
        lastJ=j;

        vLast.copyFrom(vYB);
    }
    QPReconstructLambda(vLambda,vLambdaScale); vP.copyFrom(vYB); return;
}

void restartSimpleQPSolve(Vector vBO,   // in
                          Vector vP)                    // out
{
    if (bQRFailed_QP) { vP.zero(); return; }
    int i,k=0, *ii=vi_QP, nc2=vi_QP.sz();
    double *lambda=vLastLambda_QP, *b=vBO;
    Vector vTmp(nc2);
    for (i=0; i<nc2; i++)
    {
        if (lambda[i]!=0.0)
        {
            b[k]=b[ii[i]];
            k++;
        }
    }
    vBO.setSize(k);
    vBO.permutIn(vTmp,viPermut_QP);
    mR_QP.solveInPlace(vTmp);
    mQ_QP.multiply(vP,vTmp);
}

#endif
