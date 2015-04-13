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

// model step solver

#include <stdio.h>
#include <memory.h>

#include "Matrix.h"
#include "tools.h"
#include "Poly.h"


//int lagmax_(int *n, double *g, double *h__, double *rho, 
//                    double *d__, double *v, double *vmax);

Vector LAGMAXModified(Vector G, Matrix H, double rho,double &VMAX)
{
    //not tested in depth but running already quite good

//    SUBROUTINE LAGMAX (N,G,H,RHO,D,V,VMAX)
//    IMPLICIT REAL*8 (A-H,O-Z)
//    DIMENSION G(*),H(N,*),D(*),V(*)
//
//    N is the number of variables of a quadratic objective function, Q say.
//    G is the gradient of Q at the origin.
//    H is the symmetric Hessian matrix of Q. Only the upper triangular and
//      diagonal parts need be set.
//    RHO is the trust region radius, and has to be positive.
//    D will be set to the calculated vector of variables.
//    The array V will be used for working space.
//    VMAX will be set to |Q(0)-Q(D)|.
//
//     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .EQ. RHO
//    requires of order N**3 operations, but sometimes it is adequate if
//    |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
//    subroutine provides such a solution in only of order N**2 operations,
//    where the claim of accuracy has been tested by numerical experiments. 
/*
    int n=G.sz();
    Vector D(n), V(n);
    lagmax_(&n, (double *)G, *((double**)H), &rho, 
                    (double*)D, (double*)V, &VMAX);
    return D;
*/
    int i,n=G.sz();
    Vector D(n);

    Vector V=H.getMaxColumn();
    D=H.multiply(V);
    double vv=V.square(),
           dd=D.square(),
           vd=V.scalarProduct(D), 
           dhd=D.scalarProduct(H.multiply(D)),
           *d=D, *v=V, *g=G;

//
//     Set D to a vector in the subspace spanned by V and HV that maximizes
//     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
//
    if (sqr(vd)<0.9999*vv*dd)
    {
        double a=dhd*vd-dd*dd,
               b=.5*(dhd*vv-dd*vd),
               c=dd*vv-vd*vd,
               tmp1=-b/a;
        if (b*b>a*c)
        {
            double tmp2=sqrt(b*b-a*c)/a, dd1, dd2, dhd1, dhd2;
            Vector D1=D.clone();
            D1.multiply(tmp1+tmp2);
            D1+=V;
            dd1=D1.square();
            dhd1=D1.scalarProduct(H.multiply(D1));

            Vector D2=D.clone();
            D2.multiply(tmp1-tmp2);
            D2+=V;
            dd2=D2.square();
            dhd2=D2.scalarProduct(H.multiply(D2));

            if (condorAbs(dhd1/dd1)>condorAbs(dhd2/dd2)) { D=D1; dd=dd1; dhd=dhd1; } 
            else { D=D2; dd=dd2; dhd=dhd2; }
            d=(double*)D;
        }
    };


//
//     We now turn our attention to the subspace spanned by G and D. A multiple
//     of the current D is returned if that choice seems to be adequate.
//
    double gg=G.square(),
           normG=sqrt(gg),
           gd=G.scalarProduct(D), 
           temp=gd/gg,
           scale=sign(rho/sqrt(dd), gd*dhd);

    i=n; while (i--) v[i]=d[i]-temp*g[i];
    vv=V.square();

    if ((normG*dd)<(0.5-2*rho*condorAbs(dhd))||(vv/dd<1e-4))
    {
        D.multiply(scale);
        VMAX=condorAbs(scale*(gd+0.5*scale*dhd));
        return D;
    }

//
//     G and V are now orthogonal in the subspace spanned by G and D. Hence
//     we generate an orthonormal basis of this subspace such that (D,HV) is
//     negligible or zero, where D and V will be the basis vectors.
//
    H.multiply(D,G); //  D=HG;
    double ghg=G.scalarProduct(D),
           vhg=V.scalarProduct(D),
           vhv=V.scalarProduct(H.multiply(V));
    double theta, cosTheta, sinTheta;

    if (condorAbs(vhg)<0.01*mmax(condorAbs(vhv),condorAbs(ghg)))
    {
        cosTheta=1.0;
        sinTheta=0.0;
    } else
    {
        theta=0.5*atan(0.5*vhg/(vhv-ghg));
        cosTheta=cos(theta);
        sinTheta=sin(theta);
    }
    i=n;
    while(i--)
    {
        d[i]= cosTheta*g[i]+ sinTheta*v[i];
        v[i]=-sinTheta*g[i]+ cosTheta*v[i];
    };

//
//     The final D is a multiple of the current D, V, D+V or D-V. We make the
//     choice from these possibilities that is optimal.
//

    double norm=rho/D.euclidianNorm();
    D.multiply(norm);
    dhd=(ghg*sqr(cosTheta)+vhv*sqr(sinTheta))*sqr(norm);

    norm=rho/V.euclidianNorm();
    V.multiply(norm);
    vhv=(ghg*sqr(sinTheta)+vhv*sqr(cosTheta)*sqr(norm));

    double halfRootTwo=sqrt(0.5),   // =sqrt(2)/2=cos(PI/4)
           t1=normG*cosTheta*rho,   // t1=condorAbs(D.scalarProduct(G));
           t2=normG*sinTheta*rho,   // t2=condorAbs(V.scalarProduct(G));
           at1=condorAbs(t1),
           at2=condorAbs(t2),
           t3=0.25*(dhd+vhv),
           q1=condorAbs(at1+0.5*dhd),
           q2=condorAbs(at2+0.5*vhv),
           q3=condorAbs(halfRootTwo*(at1+at2)+t3),
           q4=condorAbs(halfRootTwo*(at1-at2)+t3);
    if ((q4>q3)&&(q4>q2)&&(q4>q1))
    {
        double st1=sign(t1*t3), st2=sign(t2*t3);
        i=n; while (i--) d[i]=halfRootTwo*(st1*d[i]-st2*v[i]);
        VMAX=q4;
        return D;
    }
    if ((q3>q2)&&(q3>q1))
    {
        double st1=sign(t1*t3), st2=sign(t2*t3);
        i=n; while (i--) d[i]=halfRootTwo*(st1*d[i]+st2*v[i]);
        VMAX=q3;
        return D;
    }
    if (q2>q1)
    {
        if (t2*vhv<0) V.multiply(-1);
        VMAX=q2;
        return V;
    }
    if (t1*dhd<0) D.multiply(-1);
    VMAX=q1;
    return D;
}

Vector LAGMAXModified(Polynomial q, Vector pointXk, double rho,double &VMAX)
{
    int n=q.dim();
    Matrix H(n,n);
    Vector G(n), D(n);
    q.gradientHessian(pointXk,G,H);
    return LAGMAXModified(G,H,rho,VMAX);
}

Vector LAGMAXModified(Polynomial q, double rho,double &VMAX)
{
    return LAGMAXModified(q,Vector::emptyVector,rho,VMAX);
}
