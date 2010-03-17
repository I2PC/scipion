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
#include "Poly.h"
#include "Vector.h"
#include "tools.h"
#include "KeepBests.h"
#include "IntPoly.h"

Vector LAGMAXModified(Vector G, Matrix H, double rho,double &VMAX);

/*
Vector *GenerateData(double **valuesF, double rho, 
                     Vector Base, double vBase, ObjectiveFunction *of )
// generate points to allow start of fitting a polynomial of second degree 
// around point Base
{
    int j,k,dim=Base.sz(), N=(dim+1)*(dim+2)/2;
    double *vf=(double*)malloc(2*N*sizeof(double));
    int *failed=(int*)(vf+N); // value objective function
    *valuesF=vf;
    Vector *ap=new Vector[ N-1 ], *cp=ap, cur; // ap: allPoints
                                               // cp: current Point
    double *sigma=(double*)malloc(dim*sizeof(double));

    for (j=0; j<dim;j++)
    {
        cur=Base.clone();
        cur[j]+=rho;
        
        *(cp++)=cur;
    }

	calculateNParallelJob(dim,vf,ap,of,failed);

    for (j=0; j<dim; j++)
	{        
        cur=Base.clone();
        if (*(vf++)<vBase) { cur[j]+=2*rho; sigma[j]=rho; }
        else { cur[j]-=rho; sigma[j]=-rho; }
        *(cp++)=cur;
	}

    for (j=0; j<dim; j++)
    {
        for (k=0; k<j; k++)
        {
            cur=Base.clone();
            cur[j]+=sigma[j];
            cur[k]+=sigma[k];
            *(cp++)=cur;
        }
    }

    free(sigma);

	// parallelize here !
	calculateNParallelJob(N-dim-1,vf,ap+dim,of,failed);
	
	return ap;
}
*/

// Note:
// Vectors do come from outside. Newton Basis and associated permutation
// vector are generated internally and can be deleted.

/*
double *InterPolynomial::NewtonCoefficient(double *yy)
{
   // Initialize to local variables
   unsigned N=nPtsUsed,i;
   Polynomial *pp=NewtonBasis;
   Vector *xx=NewtonPoints;
   double *ts=(double*)calloc(N,sizeof(double)), *tt=ts;

   if (!ts)
   {
      printf("NewtonCoefficient : No mem\n");
      getchar(); exit(251);
   }

   for (i=0; i<N; i++)    // 0th difference everywhere
    *(tt++) = yy[i];
  
   unsigned deg=d->deg, Dim=d->dim, Nfrom, Nto, j, curDeg;
   double *ti, *tj;

   for (curDeg=0, Nfrom=0; curDeg<deg; Nfrom=Nto )
   {
      Nto=Nfrom+choose( curDeg+Dim-1,Dim-1 );
      for (ti=ts+Nto, i=Nto; i<N; i++,ti++)
      {
         for (tj=ts+Nfrom, j=Nfrom; j<Nto; j++, tj++) 
            *ti -= *tj ? // Evaluation takes time 
                         *tj * (pp[j])( xx[i] ) : 0;
      }
      curDeg++;
   }
   return ts;
}
*/

void InterPolynomial::ComputeLagrangeBasis(double *yy, unsigned nPtsTotal)
{
#ifdef NOOBJECTIVEFUNCTION
   const double eps  =  1e-6;
#else
   const double eps  =  1e-10;
#endif
   const double good =  1 ;
   unsigned dim=d->dim,i;
  
   Vector *xx=NewtonPoints, xtemp;
   Polynomial *pp= new Polynomial[nPtsUsed], *qq=pp;
   NewtonBasis=pp;

   if (!pp)
   {
      printf("ComputeNewtonBasis( ... ) : Alloc for polynomials failed\n");
      getchar(); exit(251);
   }
  
   MultInd I(dim);
   for (i=0; i<nPtsUsed; i++)
   {
      *(qq++)=Polynomial(I);
      I++;
   }
    
   unsigned k, kmax;
   double v, vmax, vabs;
#ifdef VERBOSE
   printf("Constructing first quadratic ... (N=%i)\n",nPtsUsed);
#endif
   for (i=0; i<nPtsUsed; i++)
   {
#ifdef VERBOSE
      printf(".");
#endif
      vmax = vabs = 0;
      kmax = i;
      if (i==0)
      {
          // to be sure point 0 is always taken:
          vmax=(pp[0])( xx[0] );
      } else
      for (k=i; k<nPtsTotal; k++)    // Pivoting
      {
          v=(pp[i])( xx[k] );
          if (fabs(v) > vabs)
          {
              vmax = v;
              vabs = abs(v);
              kmax = k;
          }
          if (fabs(v) > good ) break;
      }

      // Now, check ...
      if (fabs(vmax) < eps)
      {
          printf("Cannot construct Lagrange basis");
          delete[] NewtonBasis;
          NewtonBasis=NULL;
          return;
      }

      // exchange component i and k of NewtonPoints
      // fast because of shallow copy
      xtemp   =xx[kmax];
      xx[kmax]=xx[i];
      xx[i]   =xtemp;

      // exchange component i and k of newtonData
      v       =yy[kmax];
      yy[kmax]=yy[i];
      yy[i]   =v;

      pp[i]/=vmax;
      for (k=0;  k<i;        k++) pp[k] -= (pp[k])( xx[i] ) * pp[i];
      for (k=i+1; k<nPtsUsed; k++) pp[k] -= (pp[k])( xx[i] ) * pp[i];

      // Next polynomial, break if necessary
   }
#ifdef VERBOSE
   printf("\n");
#endif
}

int findK(double *ValuesF, int n, ObjectiveFunction *of, Vector *points)
{
    // if (of->isConstrained) return 0;
    // find index k of the best value of the function 
    double minimumValueF=INF;
    int i,k=-1; 
    for (i=0; i<n; i++)
        if ((ValuesF[i]<minimumValueF)
             &&((of==NULL)||(of->isFeasible(points[i])))) 
        { k=i; minimumValueF=ValuesF[i]; }
    if (k==-1) k=0;
    return k;
}

InterPolynomial::InterPolynomial( unsigned _deg, double rho,Vector _vBase, Matrix data, ObjectiveFunction *of)
//    unsigned _nPtsTotal, Vector *_Pp, double *_Yp ) 
    : Polynomial(data.nColumn()-2, _deg), M(0.0), nUpdateOfM(0), NewtonBasis(NULL), NewtonPoints(NULL), kbest(-1)
{
    nPtsUsed=choose( _deg+d->dim,d->dim );
    if (data.nLine()==0)
    {
        printf( "InterPolynomial::InterPolynomial( double *) : No data\n");
        getchar(); exit(-1);
    }

    // Generate basis 
    if (_vBase==Vector::emptyVector)
    {
        double rhosmall=rho*ROOT2*.5;
        GenerateBasis(rho,rhosmall, data,NULL);
        if (!NewtonBasis) 
        {
            GenerateBasis(rho,rho*.5, data,of);
        }
    } else
    {
        int i,ii,j=0,k,mdim=dim();
        double norm, *p, *pb=_vBase, rhomin=0;
        KeepBests kb(nPtsUsed*2,mdim);
//        fprintf(stderr,"Value Objective=%e\n", vBase);
        
        ii=data.lineIndex(_vBase);
        if (ii!=-1)
        {
            p=data[ii];
            kb.add(0.0,p[mdim],p); j++; 
            rhomin=rho;
        }

        i=data.nLine();
        while (i--)
        {
            if (i==ii) continue;
            p=data[i];
            if (p[mdim+1]) continue;
            norm=0; k=mdim;
            while (k--) norm+=sqr(p[k]-pb[k]);
            if ((rho==0.0)||
                ((norm<=4.00*rho*rho))//&&(norm>.25*rhomin*rhomin))
               ) 
            {
                kb.add(norm,p[mdim], p); j++;
                if (norm<.25*rho*rho) rhomin=rho;
            }
        }
        if (j<(int)nPtsUsed)
        {
            printf("Not Enough points in DB to build interpolation polynomial\n");
            return;
        }
        // we have retained only the 2*nPtsUsed best points:
        j=mmin(j,(int)(2*nPtsUsed));
        NewtonPoints=new Vector[j];
        valuesF=(double*)malloc(j*sizeof(double));
        for (i=0; i<j; i++)
        {
            valuesF[i]=kb.getValue(i);
            NewtonPoints[i]=Vector(mdim, kb.getOptValue(i));
        }
        //kbest=findK(valuesF,j,NULL,NewtonPoints);
        //vBase=NewtonPoints[kbest].clone();
        vBase=_vBase;
        for (i=0; i<j; i++) NewtonPoints[i]-=vBase;
        ComputeLagrangeBasis(valuesF, j);
    }
    if (NewtonBasis==NULL) return;

//    test();

    // Compute Interpolant
//    double *NewtonCoefPoly=NewtonCoefficient(_Yp);
    double *NewtonCoefPoly=valuesF;

    double *coc= NewtonCoefPoly+nPtsUsed-1;
    Polynomial *ppc= NewtonBasis+nPtsUsed-1;
    this->copyFrom((Polynomial)(*coc * *ppc));        //take highest degree

    int i=nPtsUsed-1;
    if (i)
      while(i--)
        (*this) += *(--coc) * *(--ppc);    
                // No reallocation here because of the order of
                // the summation
    //free(NewtonCoefPoly);
}

#ifdef NOOBJECTIVEFUNCTION
int calculateNParallelJob(int n,double *vf,Vector *cp, ObjectiveFunction *of, int *notsuccess){return 0;}
#else
int calculateNParallelJob(int n,double *vf,Vector *cp, ObjectiveFunction *of, int *notsuccess);
#endif

void InterPolynomial::GenerateBasis(double rho,double rhosmall,Matrix data,ObjectiveFunction *of)
// generate points to allow start of fitting a polynomial of second degree 
// around point data[0]
{
    if (deg()!=2)
    {
        printf("The GenerateBasis function only works for interpolation polynomials of degree 2.\n");
        exit(244);
    }
    int i,j,k,l,dim=data.nColumn()-2, N=(dim+1)*(dim+2)/2,
        *failed=(int*)malloc(N*sizeof(int)),nToCompute, nl=data.nLine(),best;
    Vector Base=data.getLine(0,dim);
    VectorInt vi(dim);
    valuesF=(double*)malloc((N+dim)*sizeof(double));
    double dBase=((double**)data)[0][dim],*vrho=valuesF+N,
           *base=Base,*dcur,*p,*pb, norm, normbest, r; // value objective function
    
    NewtonPoints=new Vector[N];
    Vector *cp, cur; // cp: current Point
    NewtonPoints[N-1]=Base.clone(); valuesF[N-1]=dBase; 

    // generate perfect sampling site
    for (i=0; i<dim; i++) 
    { 
        cur=Base.clone();
        cur[i]+=rho;
        NewtonPoints[i]=cur;
        vi[i]=i;
    }

    nToCompute=0; cp=NewtonPoints;
    // try to find good points in DB, which are near perfect sampling sites
    for (i=0; i<dim; i++)
    {
        normbest=INF; pb=cp[i];
        for (j=1; j<nl; j++)
        {
            p=data[j];
            if (p[dim+1]) continue;
            norm=0; k=dim;
            while (k--) norm+=sqr(p[k]-pb[k]);
            if (normbest>norm) { normbest=norm; best=j; }
        }
        if (normbest<rhosmall*rhosmall)
        {
            cp[i]=data.getLine(best,dim);
            valuesF[i]=((double**)data)[best][dim];
        } else
        {
            cur=cp[nToCompute]; cp[nToCompute]=cp[i]; cp[i]=cur;
            l=vi[nToCompute];   vi[nToCompute]=vi[i]; vi[i]=l;
            valuesF[i]=valuesF[nToCompute];
            nToCompute++;
        }
    }

    if ((!of)&&(nToCompute)) { free(failed); free(valuesF); delete[] NewtonPoints; return; }

    // if some points have not been found in DB, start evaluation
    while (nToCompute>0)
    {
        // parallelize here !
	    calculateNParallelJob(nToCompute,valuesF,NewtonPoints,of,failed);

        k=0;
        for (j=0; j<nToCompute; j++)
        {
            of->saveValue(cp[j],valuesF[j],failed[j]);
            if (failed[j])
	        {
                dcur=cp[j];
                for (i=0; i<dim; i++)
                {
                    dcur[i]=base[i]+(base[i]-dcur[i])*.7;
                }

                // group all missing values at the beginning (around NewtonPoints+dim)
                cur=cp[k]; cp[k]=cp[j]; cp[j]=cur;
                l=vi[k];   vi[k]=vi[j]; vi[j]=l;
                valuesF[j]=valuesF[k];
                k++;
	        }
        }
        nToCompute=k;
    }

    // re-order vectors (based on vi)
    for (j=0; j<dim-1; j++)
	{
        k=j; while (vi[k]!=j) k++;
        if (k==j) continue;
        cur=cp[k];    cp[k]=cp[j];           cp[j]=cur;
        l=vi[k];      vi[k]=vi[j];           vi[j]=l;
        r=valuesF[k]; valuesF[k]=valuesF[j]; valuesF[j]=r;
    }

    // select again some good sampling sites
    cp=NewtonPoints+dim;
    for (j=0; j<dim; j++)
	{
        cur=Base.clone();
        dcur=NewtonPoints[j];
        r=dcur[j]-base[j];
        if (valuesF[j]<dBase) { cur[j]+=r*2.0; vrho[j]=r; }
        else { cur[j]-=r; vrho[j]=-r; }
        *(cp++)=cur;
	}
    for (j=0; j<dim; j++)
    {
        for (k=0; k<j; k++)
        {
            cur=Base.clone();
            cur[j]+=vrho[j];
            cur[k]+=vrho[k];
            *(cp++)=cur;
        }
    }

    nToCompute=0;
    cp=NewtonPoints+dim;
    // try to find good points in DB, which are near perfect sampling sites
    for (i=0; i<N-dim-1; i++)
    {
        normbest=INF; pb=cp[i];
        for (j=1; j<nl; j++)
        {
            p=data[j];
            if (p[dim+1]) continue;
            norm=0; k=dim;
            while (k--) norm+=sqr(p[k]-pb[k]);
            if (normbest>norm) { normbest=norm; best=j; }
        }
        if (normbest<rhosmall*rhosmall)
        {
            cp[i]=data.getLine(best,dim);
            valuesF[i+dim]=((double**)data)[best][dim];
        } else
        {
            cur=cp[nToCompute];
            cp[nToCompute]=cp[i];
            cp[i]=cur;
            valuesF[i+dim]=valuesF[nToCompute+dim];
            nToCompute++;
        }
    }

    if ((!of)&&(nToCompute)) { free(failed); free(valuesF); delete[] NewtonPoints; return; }

    while (nToCompute>0)
    {
	    // parallelize here !
	    calculateNParallelJob(nToCompute,valuesF+dim,cp,of,failed);

        k=0;
        for (j=0; j<nToCompute; j++)
        {
            of->saveValue(cp[j],valuesF[j+dim],failed[j]);
            if (failed[j])
	        {
                dcur=cp[j];
                for (i=0; i<dim; i++)
                {
                    dcur[i]=base[i]+(base[i]-dcur[i])*.65;
                }

                // group all missing values at the beginning (around NewtonPoints+dim)
                cur=cp[k];
                cp[k]=cp[j];
                cp[j]=cur;
                valuesF[j+dim]=valuesF[k+dim];
                k++;
	        }
        }
        nToCompute=k;

    }
    free(failed); 

    // "insertion sort" to always place best points first
    for (i=0; i<N-1; i++)
    {
        kbest=findK(valuesF+i,N-i,of,NewtonPoints+i)+i; 
        
        // to accept infeasible points for i>=1 :
        of=NULL;

        cur=NewtonPoints[i]; NewtonPoints[i]=NewtonPoints[kbest]; NewtonPoints[kbest]=cur;
        r=valuesF[i];        valuesF[i]=valuesF[kbest];           valuesF[kbest]=r;
    }
    kbest=0;
    vBase=NewtonPoints[0].clone();
    for (i=0; i<N; i++) NewtonPoints[i]-=vBase;
    ComputeLagrangeBasis(valuesF, N);
}

void InterPolynomial::updateM(Vector newPoint, double valueF)
{
    //not tested
    unsigned i=nPtsUsed;
    double sum=0,a;
    Polynomial *pp=NewtonBasis;
    Vector *xx=NewtonPoints;

    while (i--)
    {
        a=newPoint.euclidianDistance(xx[i]);
        sum+=abs(pp[i]( newPoint ))*a*a*a;
    }
    M=mmax(M, abs((*this)(newPoint)-valueF)/sum);
    nUpdateOfM++;
}

double InterPolynomial::interpError(Vector Point)
{
    unsigned i=nPtsUsed;
    double sum=0,a;
    Polynomial *pp=NewtonBasis;
    Vector *xx=NewtonPoints;

    while (i--)
    {
        a=Point.euclidianDistance(xx[i]);
        sum+=abs(pp[i]( Point ))*a*a*a;
    }
    return M*sum;
}

int InterPolynomial::findAGoodPointToReplace(int excludeFromT, 
                                  double rho, Vector pointToAdd, double *maxd)
{
    //not tested

    // excludeFromT is set to k if not sucess from optim and we want
    //        to be sure that we keep the best point
    // excludeFromT is set to -1 if we can replace the point x_k by
    //        pointToAdd(=x_k+d) because of the success of optim.
    
    // choosen t: the index of the point inside the newtonPoints
    //            which will be replaced.

    Vector *xx=NewtonPoints;
    Vector XkHat;
    if (excludeFromT>=0) XkHat=xx[excludeFromT];
    else XkHat=pointToAdd;

    int t=-1, i, N=nPtsUsed;
    double a, aa, maxa=-1.0, maxdd=0;
    Polynomial *pp=NewtonBasis;

  //  if (excludeFromT>=0) maxa=1.0;

    for (i=0; i<N; i++)
    {
        if (i==excludeFromT) continue;
        aa=XkHat.euclidianDistance(xx[i]);
        if (aa==0.0) return -1;
        a=aa/rho;
        // because of the next line, rho is important:
        a=mmax(a*a*a,1.0);
        a*=abs(pp[i] (pointToAdd));

        if (a>maxa)
        {
            t=i; maxa=a; maxdd=aa;
        }
    }
    if (maxd) *maxd=maxdd;
    return t;
}
/*
void InterPolynomial::check(Vector Base, double (*f)(  Vector ) )
{
    int i,n=sz();
    double r, bound;

    for (i=0; i<n; i++)
    {
        r=(*f)(NewtonPoints[i]+Base);
        bound=(*this)(NewtonPoints[i]);
        if ((abs(bound-r)>1e-15)&&(abs(bound-r)>1e-3*abs(bound))) 
        {
            printf("error\n"); 
            test();
        }
//                    for (j=0; j<n; j++)
//                        r=poly.NewtonBasis[j](poly.NewtonPoints[i]);
    }
}

void InterPolynomial::test()
{
    unsigned i,j,n=d->n; Matrix M(n,n); double **m=M;
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            m[i][j]=NewtonBasis[i](NewtonPoints[j]);
    M.print();
};
*/
void InterPolynomial::replace(int t, Vector pointToAdd, double valueF)
{
    //not tested
    updateM(pointToAdd, valueF);
    if (t<0) return;

    Vector *xx=NewtonPoints;
    Polynomial *pp=NewtonBasis, t1;
    int i, N=nPtsUsed;
    double t2=(pp[t]( pointToAdd ));

    if (t2==0) return;

    t1=pp[t]/=t2;
  
    for (i=0; i<t; i++)   pp[i]-= pp[i]( pointToAdd )*t1;
    for (i=t+1; i<N; i++) pp[i]-= pp[i]( pointToAdd )*t1;
    xx[t].copyFrom(pointToAdd);

    // update the coefficents of general poly.

    valueF-=(*this)(pointToAdd);
    if (abs(valueF)>1e-11) (*this)+=valueF*pp[t];
    
//    test();
}

int InterPolynomial::maybeAdd(Vector pointToAdd, unsigned k, double rho, double valueF)
// input: pointToAdd, k, rho, valueF
// output: updated polynomial
{

    unsigned i,N=nPtsUsed;
    int j;
    Vector *xx=NewtonPoints, xk=xx[k];
    double distMax=-1.0,dd;
 /*
    Polynomial *pp=NewtonBasis;
    Vector *xx=NewtonPoints, xk=xx[k],vd;
    Matrix H(n,n);
    Vector GXk(n); //,D(n);
*/
    // find worst point/newton poly

    for (i=0; i<N; i++)
    {
        dd=xk.euclidianDistance(xx[i]);
        if (dd>distMax) { j=i; distMax=dd; };
    }
    dd=xk.euclidianDistance(pointToAdd);

    // no tested:

    if (abs(NewtonBasis[j](pointToAdd))*distMax*distMax*distMax/(dd*dd*dd)>1.0) 
    {
        printf("good point found.\n");
        replace(j, pointToAdd, valueF);
        return 1;
    }
    return 0;
}

int InterPolynomial::checkIfValidityIsInBound(Vector ddv, unsigned k, double bound, double rho)
// input: k,bound,rho
// output: j,ddv
{

    // check validity around x_k
    // bound is epsilon in the paper
    // return index of the worst point of J
    // if (j==-1) then everything OK : next : trust region step
    // else model step: replace x_j by x_k+d where d
    // is calculated with LAGMAX

    unsigned i,N=nPtsUsed, n=dim();
    int j;
    Polynomial *pp=NewtonBasis;
    Vector *xx=NewtonPoints, xk=xx[k],vd;
    Vector Distance(N);
    double *dist=Distance, *dd=dist, distMax, vmax, tmp;
    Matrix H(n,n);
    Vector GXk(n); //,D(n);

    for (i=0; i<N; i++) *(dd++)=xk.euclidianDistance(xx[i]);

    while (true)
    {
        dd=dist; j=-1; distMax=2*rho;
        for (i=0; i<N; i++)
        {
            if (*dd>distMax) { j=i; distMax=*dd; };
            dd++;
        }
        if (j<0) return -1;

        // to prevent to choose the same point once again:
        dist[j]=0;
        
        pp[j].gradientHessian(xk,GXk,H);
//        d=H.multiply(xk);
//        d.add(G);

        tmp=M*distMax*distMax*distMax;

        if (tmp*rho*(GXk.euclidianNorm()+0.5*rho*H.frobeniusNorm())>=bound)
        {
/*			vd=L2NormMinimizer(pp[j], xk, rho);
            vd+=xk;
            vmax=abs(pp[j](vd));

			Vector vd2=L2NormMinimizer(pp[j], xk, rho);
            vd2+=xk;
            double vmax2=abs(pp[j](vd));
			
			if (vmax<vmax2) { vmax=vmax2; vd=vd2; }
*/
			vd=LAGMAXModified(GXk,H,rho,vmax);
//            tmp=vd.euclidianNorm();
	        vd+=xk;
			vmax=abs(pp[j](vd));
            
            if (tmp*vmax>=bound) break;
        }
    }
    if (j>=0) ddv.copyFrom(vd);
    return j;
}

int InterPolynomial::getGoodInterPolationSites(Matrix d, int k, double rho, Vector *v)
// input: k,rho
// output: d,r
{
    //not tested

    unsigned i, N=nPtsUsed, n=dim();
    int ii, j,r=0;
    Polynomial *pp=NewtonBasis;
    Vector *xx=NewtonPoints, xk,vd;
    Vector Distance(N);
    double *dist=Distance, *dd=dist, distMax, vmax;
    Matrix H(n,n);
    Vector GXk(n);
    if (k>=0) xk=xx[k]; else xk=*v;

    for (i=0; i<N; i++) *(dd++)=xk.euclidianDistance(xx[i]);

    for (ii=0; ii<d.nLine(); ii++)
    {
        dd=dist; j=-1; 
        distMax=-1.0;
        for (i=0; i<N; i++)
        {
            if (*dd>distMax) { j=i; distMax=*dd; };
            dd++;
        }
        // to prevent to choose the same point once again:
        dist[j]=-1.0;
        
        if (distMax>2*rho) r++;
        pp[j].gradientHessian(xk,GXk,H);
		vd=LAGMAXModified(GXk,H,rho,vmax);
	    vd+=xk;

        d.setLine(ii,vd);
    }
    return r;
}

void InterPolynomial::translate(Vector translation)
{
    if (translation==Vector::emptyVector) return;
    vBase+=translation;
    Polynomial::translate(translation);
    int i=nPtsUsed;
    while (i--) NewtonBasis[i].translate(translation);
    i=nPtsUsed;
    while (i--) if (NewtonPoints[i]==translation) NewtonPoints[i].zero();
                else NewtonPoints[i]-=translation;
}


// to allow shallow copy:

void InterPolynomial::destroyCurrentBuffer()
{
    if (!d) return;
	if (d->ref_count==1)
    {
        if (NewtonBasis)
        {
            delete[] NewtonBasis;
            delete[] NewtonPoints;
            free(valuesF);
        }
    }
}

InterPolynomial::~InterPolynomial()
{
    destroyCurrentBuffer();
}

InterPolynomial::InterPolynomial(const InterPolynomial &A)
{
    // shallow copy for inter poly.
    d=A.d;
    NewtonBasis=A.NewtonBasis;
    NewtonPoints=A.NewtonPoints;
//    ValuesF=A.ValuesF;
	M=A.M;
    nPtsUsed=A.nPtsUsed; 
    nUpdateOfM=A.nUpdateOfM;
	(d->ref_count)++;
    
}

InterPolynomial& InterPolynomial::operator=( const InterPolynomial& A )
{
    // shallow copy
    if (this != &A)
	{
        destroyCurrentBuffer();

        d=A.d;
        NewtonBasis=A.NewtonBasis;
        NewtonPoints=A.NewtonPoints;
//        ValuesF=A.ValuesF;
	    M=A.M;
        nPtsUsed=A.nPtsUsed; 
        nUpdateOfM=A.nUpdateOfM;
        kbest=A.kbest;
        vBase=A.vBase;
        valuesF=A.valuesF;

		(d->ref_count) ++ ;
	}
	return *this;
}

InterPolynomial::InterPolynomial(unsigned _dim, unsigned _deg): Polynomial(_dim, _deg), M(0.0), nUpdateOfM(0)
{ 	
    nPtsUsed=choose( _deg+_dim,_dim );
    NewtonBasis= new Polynomial[nPtsUsed];
	NewtonPoints= new Vector[nPtsUsed];
}

InterPolynomial InterPolynomial::clone()
{
    // a deep copy
    InterPolynomial m(d->dim, d->deg);
    m.copyFrom(*this);
    return m;
}

void InterPolynomial::copyFrom(InterPolynomial m)
{
    if (m.d->dim!=d->dim)
    {
        printf("poly: copyFrom: dim do not agree\n");
        getchar(); exit(254);
    }
    if (m.d->deg!=d->deg)
    {
        printf("poly: copyFrom: degree do not agree\n");
        getchar(); exit(254);
    }
    Polynomial::copyFrom(m);
	M=m.M;
//    nPtsUsed=m.nPtsUsed; // not usefull because dim and degree already agree.
    nUpdateOfM=m.nUpdateOfM;
// ValuesF        

	int i=nPtsUsed;
	while (i--) 
	{ 
//		NewtonBasis[i]=m.NewtonBasis[i];
//		NewtonPoints[i]=m.NewtonPoints[i]; 
		NewtonBasis[i]=m.NewtonBasis[i].clone();
		NewtonPoints[i]=m.NewtonPoints[i].clone(); 
	} 
}

void InterPolynomial::copyFrom(Polynomial m)
{
    Polynomial::copyFrom(m);
}
