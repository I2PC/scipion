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
// trust region step solver

#include <stdio.h>
#include <memory.h>

#include "Matrix.h"
#include "tools.h"
#include "Poly.h"

double findAlpha(Vector s,Vector u, double delta, Polynomial &q, Vector pointXk,Vector &output, Vector minusG, Matrix H)
// find root (apha*) of equation L2norm(s+alpha u)=delta
// which makes q(s)=<g,s>+.5*<s,Hs> smallest
// output is (s+alpha* u)
{
    static Vector v1, v2, tmp;
    double a=0,b=0,c=-sqr(delta), *sp=s, *up=u;
    int n=s.sz();
    while (n--)
    {
        a+=sqr(*up); b+=*up * *sp; c+=sqr(*sp);
        sp++; up++;
    }
    double tmp1=-b/a, tmp2= sqrt(b*b-a*c)/a, q1, q2;

    n=s.sz();
    v1.setSize(n); v2.setSize(n); tmp.setSize(n);
    if (!(q==Polynomial::emptyPolynomial))
    {
        v1.copyFrom(u);
        v1.multiply(tmp1+tmp2);
        v1+=s;
        tmp.copyFrom(v1);
        tmp+=pointXk;
        q1=q(tmp);

    // !!! don't do this:
    //    output=v1;
    //    return tmp1+tmp2;

        v2.copyFrom(u);
        v2.multiply(tmp1-tmp2);
        v2+=s;
        tmp.copyFrom(v2);
        tmp+=pointXk;
        q2=q(tmp);

    } else
    {
        v1.copyFrom(u);
        v1.multiply(tmp1+tmp2);
        v1+=s;
        H.multiply(tmp,v1);
        q1=-minusG.scalarProduct(v1)+0.5*v1.scalarProduct(tmp);

        v2.copyFrom(u);
        v2.multiply(tmp1-tmp2);
        v2+=s;
        H.multiply(tmp,v2);
        q2=-minusG.scalarProduct(v2)+0.5*v2.scalarProduct(tmp);
    }
    if (q1>q2) { output=v1; return tmp1+tmp2; }
    output=v2; return tmp1-tmp2;
}

double initLambdaL(double normG,double delta, Matrix H)
{
    int n=H.nLine(),i,j;
    double **h=H, sum,l,a=INF;

    for (i=0; i<n; i++) a=mmin(a,h[i][i]);
    l=mmax(0.0,-a);

    a=0;
    for (i=0; i<n; i++) 
    {
        sum=h[i][i];
        for (j=0; j<n; j++) if (j!=i) sum+=condorAbs(h[i][j]);
        a=mmax(a,sum);
    }
    a=mmin(a,H.frobeniusNorm());
    a=mmin(a,H.LnftyNorm());

    l=mmax(l,normG/delta-a);
    return l;
}

double initLambdaU(double normG,double delta, Matrix H)
{
    int n=H.nLine(),i,j;
    double **h=H, sum,l,a=-INF;

    for (i=0; i<n; i++) 
    {
        sum=-h[i][i];
        for (j=0; j<n; j++) if (j!=i) sum+=condorAbs(h[i][j]);
        a=mmax(a,sum);
    }
    a=mmin(a,H.frobeniusNorm());
    a=mmin(a,H.LnftyNorm());

    l=mmax(0.0,normG/delta+a);
    return l;
}

double initLambdaU2(Matrix H)
{
    int n=H.nLine(),i,j;
    double **h=H, sum,a=-INF;

    for (i=0; i<n; i++) 
    {
        sum=h[i][i];
        for (j=0; j<n; j++) if (j!=i) sum+=condorAbs(h[i][j]);
        a=mmax(a,sum);
    }
    a=mmin(a,H.frobeniusNorm());
    return mmin(a,H.LnftyNorm());
}

// #define POWEL_TERMINATION 1

Vector L2NormMinimizer(Polynomial q, Vector pointXk, double delta, 
                                        int *infoOut, int maxIter, double *lambda1, Vector minusG, Matrix H)
{
// lambda1>0.0 if interior convergence

    const double theta=0.01;
//    const double kappaEasy=0.1, kappaHard=0.2;
    const double kappaEasy=0.01, kappaHard=0.02;

    double normG,lambda,lambdaCorrection, lambdaPlus, lambdaL, lambdaU, 
           uHu, alpha, normS;
    int info=0, n=minusG.sz();
    Matrix HLambda(n,n);
    MatrixTriangle L(n);
    Vector s(n), omega(n), u(n), sFinal;
    bool gIsNull, choleskyFactorAlreadyComputed=false;

//    printf("\nG= "); minusG.print();
//    printf("\nH=\n"); H.print();

    gIsNull=minusG.isNull();
    normG=minusG.euclidianNorm();
    lambda=normG/delta;
    minusG.multiply(-1.0);

    lambdaL=initLambdaL(normG,delta,H);
    lambdaU=initLambdaU(normG,delta,H);

    //    Special case: parl = paru.
    lambdaU= mmax(lambdaU,(1+kappaEasy)*lambdaL);

    lambda=mmax(lambda, lambdaL);
    lambda=mmin(lambda, lambdaU);

    while (maxIter--)
    {
        if (!choleskyFactorAlreadyComputed)
        {
            if (!H.cholesky(L, lambda, &lambdaCorrection))
            {
       //         lambdaL=mmax(mmax(lambdaL,lambda),lambdaCorrection);
                lambdaL=mmax(lambdaL,lambda+lambdaCorrection);
                lambda=mmax(sqrt(lambdaL*lambdaU), lambdaL+theta*(lambdaU-lambdaL));
                continue;
            }
        } else choleskyFactorAlreadyComputed=false;


        // cholesky factorization successful : solve Hlambda * s = -G
        s.copyFrom(minusG);
        L.solveInPlace(s);
        L.solveTransposInPlace(s);
        normS=s.euclidianNorm();

        // check for termination
#ifndef POWEL_TERMINATION
        if (condorAbs(normS-delta)<kappaEasy*delta) 
        { 
            s.multiply(delta/normS);
            info=1; 
            break; 
        }
#else
//      powell check !!!
        HLambda.copyFrom(H);
        HLambda.addUnityInPlace(lambda);
        double sHs=s.scalarProduct(HLambda.multiply(s));
        if (sqr(delta/normS-1)<kappaEasy*(1+lambda*delta*delta/sHs)) 
        {
            s.multiply(delta/normS);
            info=1; 
            break;
        }
#endif

        if (normS<delta)
        {
            // check for termination
            // interior convergence; maybe break;
            if (lambda==0) { info=1; break; }
            lambdaU=mmin(lambdaU,lambda); 
        } else lambdaL=mmax(lambdaL,lambda);

//        if (lambdaU-lambdaL<kappaEasy*(2-kappaEasy)*lambdaL) { info=3; break; };

        omega.copyFrom(s);
        L.solveInPlace(omega);
        lambdaPlus=lambda+(normS-delta)/delta*sqr(normS)/sqr(omega.euclidianNorm());
        lambdaPlus=mmax(lambdaPlus, lambdaL);
        lambdaPlus=mmin(lambdaPlus, lambdaU);

        if (normS<delta)
        {
            L.LINPACK(u);
#ifndef POWEL_TERMINATION
            HLambda.copyFrom(H);
            HLambda.addUnityInPlace(lambda);
#endif
            uHu=u.scalarProduct(HLambda.multiply(u));
            lambdaL=mmax(lambdaL,lambda-uHu);

            alpha=findAlpha(s,u,delta,q,pointXk,sFinal,minusG,H);
            // check for termination
#ifndef POWEL_TERMINATION
            if (sqr(alpha)*uHu<
                kappaHard*(s.scalarProduct(HLambda.multiply(s))   ))//  +lambda*sqr(delta)))
#else
            if (sqr(alpha)*uHu+sHs<
                kappaHard*(sHs+lambda*sqr(delta)))
#endif
            {
                s=sFinal; info=2; break;
            }
        }
        if ((normS>delta)&&(!gIsNull)) { lambda=lambdaPlus; continue; };

        if (H.cholesky(L, lambdaPlus, &lambdaCorrection)) 
        { 
            lambda=lambdaPlus; 
            choleskyFactorAlreadyComputed=true; 
            continue; 
        }

        lambdaL=mmax(lambdaL,lambdaPlus);
        // check lambdaL for interior convergence
//      if (lambdaL==0) return s;
        lambda=mmax(sqrt(lambdaL*lambdaU), lambdaL+theta*(lambdaU-lambdaL));
    }

    if (infoOut) *infoOut=info;
    if (lambda1)
    {
        if (lambda==0.0)
        {
            // calculate the value of the lowest eigenvalue of H
            // to check
            lambdaL=0; lambdaU=initLambdaU2(H);
            while (lambdaL<0.99*lambdaU)
            {
                lambda=0.5*(lambdaL+lambdaU);
                if (H.cholesky(L,-lambda)) lambdaL=lambda;
//                if (H.cholesky(L,-lambda,&lambdaCorrection)) lambdaL=lambda+lambdaCorrection;
                else lambdaU=lambda;
            }
            *lambda1=lambdaL;
        } else *lambda1=0.0;
    }
    return s;
}

Vector L2NormMinimizer(Polynomial q, Vector pointXk, double delta, 
                                        int *infoOut, int maxIter, double *lambda1)
{
    int n=q.dim();
    Matrix H(n,n);
    Vector vG(n);
    q.gradientHessian(pointXk,vG,H);
    return L2NormMinimizer(q,pointXk,delta,infoOut,maxIter,lambda1,vG,H);
}

Vector L2NormMinimizer(Polynomial q, double delta, 
                                        int *infoOut, int maxIter, double *lambda1)
{
    return L2NormMinimizer(q, Vector::emptyVector, delta, infoOut, maxIter, lambda1);
}
