/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../xmippFuncs.hh"

/* NUMERICAL UTILITIES ----------------------------------------------------- */
void nrerror(char error_text[]) {
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* RANDOM NUMBERS ---------------------------------------------------------- */
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/* Chapter 7 Section 1: UNIFORM RANDOM NUMBERS */
double ran1(int *idum) {
	static long ix1,ix2,ix3;
	static double r[98];
	double temp;
	static int iff=0;
	int j;

	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) nrerror("RAN1: This cannot happen.");
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

/* Chapter 7 Section 3: GAUSSIAN RANDOM NUMBERS */
double gasdev(int *idum) {
	static int iset=0;
	static double gset;
	double fac,r,v1,v2;
	
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

/* MATRICES ---------------------------------------------------------------- */
#define TINY 1.0e-20;
/* Chapter 2 Section 3: LU DECOMPOSITION */
template <class T>
void ludcmp(T *a, int n, int *indx, T *d)
{
	int i,imax,j,k;
	T big,dum,sum,temp;
	T *vv;

	ask_Tvector(vv,1,n);
	*d=(T)1.0;
	for (i=1;i<=n;i++) {
		big=(T)0.0;
		for (j=1;j<=n;j++)
			if ((temp=(T)fabs((double)a[i*n+j])) > big) big=temp;
		if (big == (T)0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=(T)1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i*n+j];
			for (k=1;k<i;k++) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
		}
		big=(T)0.0;
		for (i=j;i<=n;i++) {
			sum=a[i*n+j];
			for (k=1;k<j;k++)
				sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
			if ( (dum=vv[i]*(T)fabs((double)sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax*n+k];
				a[imax*n+k]=a[j*n+k];
				a[j*n+k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j*n+j] == 0.0) a[j*n+j]=(T) TINY;
		if (j != n) {
			dum=(T)1.0/(a[j*n+j]);
			for (i=j+1;i<=n;i++) a[i*n+j] *= dum;
		}
	}
	free_Tvector(vv,1,n);
}

#undef TINY

/* Chapter 2 Section 3: LU BACKWARD-FORWARD SUBSTITUTION */
template <class T> void lubksb(T *a, int n, int *indx,T b[])
{
	int i,ii=0,ip,j;
	T sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i*n+j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i*n+j]*b[j];
		b[i]=sum/a[i*n+i];
	}
}

/* Chapter 2, Section 1. Gauss-Jordan equation system resolution ----------- */
template <class T>
void gaussj(T *a, int n, T *b, int m)
{
    T temp;
    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l,ll;
    T big,dum;
    double pivinv;

    ask_Tvector(indxc,1,n);
    ask_Tvector(indxr,1,n);
    ask_Tvector(ipiv,1,n);
    for (j=1;j<=n;j++) ipiv[j]=0;
    for (i=1;i<=n;i++) {
	big=(T)0;
	for (j=1;j<=n;j++)
	    if (ipiv[j] != 1)
		for (k=1;k<=n;k++) {
		    if (ipiv[k] == 0) {
			if (fabs((double)a[j*n+k]) >= (double) big) {
			    big=ABS(a[j*n+k]);
			    irow=j;
			    icol=k;
			}
		    } else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
		}
	++(ipiv[icol]);
	if (irow != icol) {
	    for (l=1;l<=n;l++) SWAP(a[irow*n+l],a[icol*n+l],temp)
	    for (l=1;l<=m;l++) SWAP(b[irow*n+l],b[icol*n+l],temp)
	}
	indxr[i]=irow;
	indxc[i]=icol;
	if (a[icol*n+icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
	pivinv=1.0f/a[icol*n+icol];
	a[icol*n+icol]=(T)1;
	for (l=1;l<=n;l++) a[icol*n+l] = (T)(pivinv * a[icol*n+l]);
	for (l=1;l<=m;l++) b[icol*n+l] = (T)(pivinv * b[icol*n+l]);
	for (ll=1;ll<=n;ll++)
	    if (ll != icol) {
		dum=a[ll*n+icol];
		a[ll*n+icol]=(T)0;
		for (l=1;l<=n;l++) a[ll*n+l] -= a[icol*n+l]*dum;
		for (l=1;l<=m;l++) b[ll*n+l] -= b[icol*n+l]*dum;
	    }
    }
    for (l=n;l>=1;l--) {
	if (indxr[l] != indxc[l])
	    for (k=1;k<=n;k++)
		SWAP(a[k*n+indxr[l]],a[k*n+indxc[l]],temp);
    }
    free_Tvector(ipiv,1,n);
    free_Tvector(indxr,1,n);
    free_Tvector(indxc,1,n);
}

/* SORTING ----------------------------------------------------------------- */
/* Chapter 8, Section 3: Indexing */
void indexx(int n, double arrin[], int indx[])
{
	int l,j,ir,indxt,i;
	double q;

	for (j=1;j<=n;j++) indx[j]=j;
	l=(n >> 1) + 1;
	ir=n;
	for (;;) {
		if (l > 1)
			q=arrin[(indxt=indx[--l])];
		else {
			q=arrin[(indxt=indx[ir])];
			indx[ir]=indx[1];
			if (--ir == 1) {
				indx[1]=indxt;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q < arrin[indx[j]]) {
				indx[i]=indx[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		indx[i]=indxt;
	}
}

/* Chapter 8, Section 4: Quicksort */
#define M 7
#define NSTACK 50
#define FM 7875
#define FA 211
#define FC 1663

void qcksrt(int n, double arr[])
{
	int l=1,jstack=0,j,ir,iq,i;
	int istack[NSTACK+1];
	long int fx=0L;
	double a;

	ir=n;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;arr[i]>a && i>0;i--) arr[i+1]=arr[i];
				arr[i+1]=a;
			}
			if (jstack == 0) return;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			i=l;
			j=ir;
			fx=(fx*FA+FC) % FM;
			iq=l+((ir-l+1)*fx)/FM;
			a=arr[iq];
			arr[iq]=arr[l];
			for (;;) {
				while (j > 0 && a < arr[j]) j--;
				if (j <= i) {
					arr[i]=a;
					break;
				}
				arr[i++]=arr[j];
				while (a > arr[i] && i <= n) i++;
				if (j <= i) {
					arr[(i=j)]=a;
					break;
				}
				arr[j--]=arr[i];
			}
			if (ir-i >= i-l) {
				istack[++jstack]=i+1;
				istack[++jstack]=ir;
				ir=i-1;
			} else {
				istack[++jstack]=l;
				istack[++jstack]=i-1;
				l=i+1;
			}
			if (jstack > NSTACK) nrerror("NSTACK too small in QCKSRT");
		}
	}
}

#undef M
#undef NSTACK
#undef FM
#undef FA
#undef FC

/* BESSEL FUNCTIONS -------------------------------------------------------- */
/* CO: They may not come in the numerical recipes but it is not a bad
   idea to put them here, in fact they come from Gabor's group in Feb'84     */
double bessj0(double x) {
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=57568490574.0+y*(-13362590354.0+
			  y*(651619640.7
			     +y*(-11214424.18+
				 y*(77392.33017+
				    y*(-184.9052456)))));
    ans2=57568490411.0+y*(1029532985.0+
			  y*(9494680.718
			     +y*(59272.64853+
				 y*(267.8532712+
				    y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			       +y*(-0.6911147651e-5+y*(0.7621095161e-6
						       -y*0.934935152e-7)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

/*............................................................................*/
double bessi0(double x)
{
double y, ax, ans;
if ((ax=fabs(x)) < 3.75) {
   y=x/3.75;
   y*=y;
   ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
      +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
} else {
   y=3.75/ax;
   ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
      +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
      +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
      +y*0.392377e-2))))))));
}
return ans;
}

/*............................................................................*/
double bessi1(double x)
{
double ax, ans;
double y;
if ((ax=fabs(x)) < 3.75) {
   y=x/3.75;
   y*=y;
   ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
      +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
} else {
   y=3.75/ax;
   ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
      -y*0.420059e-2));
   ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
   +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
   ans *= (exp(ax)/sqrt(ax));
}
return x < 0.0 ? -ans : ans;
}

/* General Bessel functions ------------------------------------------------ */
double chebev(double a, double b, double c[], int m, double x)
{
        double d=0.0,dd=0.0,sv,y,y2;
        int j;

        if ((x-a)*(x-b) > 0.0) nrerror("x not in range in routine chebev");
        y2=2.0*(y=(2.0*x-a-b)/(b-a));
        for (j=m-1;j>=1;j--) {
                sv=d;
                d=y2*d-dd+c[j];
                dd=sv;
        }
        return y*d-dd+0.5*c[0];
}
#define NUSE1 5
#define NUSE2 5

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
        double xx;
        static double c1[] = {
                -1.142022680371172e0,6.516511267076e-3,
                3.08709017308e-4,-3.470626964e-6,6.943764e-9,
                3.6780e-11,-1.36e-13};
        static double c2[] = {
                1.843740587300906e0,-0.076852840844786e0,
                1.271927136655e-3,-4.971736704e-6,-3.3126120e-8,
                2.42310e-10,-1.70e-13,-1.0e-15};

        xx=8.0*x*x-1.0;
        *gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
        *gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
        *gampl= *gam2-x*(*gam1);
        *gammi= *gam2+x*(*gam1);
}

#undef NUSE1
#undef NUSE2

#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
{
        int i,isign,l,nl;
        double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,
                fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,rjl,
                rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
                temp,w,x2,xi,xi2,xmu,xmu2;

        if (x <= 0.0 || xnu < 0.0) nrerror("bad arguments in bessjy");
        nl=(x < XMIN ? (int)(xnu+0.5) : MAX(0,(int)(xnu-x+1.5)));
        xmu=xnu-nl;
        xmu2=xmu*xmu;
        xi=1.0/x;
        xi2=2.0*xi;
        w=xi2/PI;
        isign=1;
        h=xnu*xi;
        if (h < FPMIN) h=FPMIN;
        b=xi2*xnu;
        d=0.0;
        c=h;
        for (i=1;i<=MAXIT;i++) {
                b += xi2;
                d=b-d;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=b-1.0/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=c*d;
                h=del*h;
                if (d < 0.0) isign = -isign;
                if (fabs(del-1.0) < EPS) break;
        }
        if (i > MAXIT) nrerror("x too large in bessjy; try asymptotic expansion");
        rjl=isign*FPMIN;
        rjpl=h*rjl;
        rjl1=rjl;
        rjp1=rjpl;
        fact=xnu*xi;
        for (l=nl;l>=1;l--) {
                rjtemp=fact*rjl+rjpl;
                fact -= xi;
                rjpl=fact*rjtemp-rjl;
                rjl=rjtemp;
        }
        if (rjl == 0.0) rjl=EPS;
        f=rjpl/rjl;
        if (x < XMIN) {
                x2=0.5*x;
                pimu=PI*xmu;
                fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
                d = -log(x2);
                e=xmu*d;
                fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
                beschb(xmu,&gam1,&gam2,&gampl,&gammi);
                ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
                e=exp(e);
                p=e/(gampl*PI);
                q=1.0/(e*PI*gammi);
                pimu2=0.5*pimu;
                fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
                r=PI*pimu2*fact3*fact3;
                c=1.0;
                d = -x2*x2;
                sum=ff+r*q;
                sum1=p;
                for (i=1;i<=MAXIT;i++) {
                        ff=(i*ff+p+q)/(i*i-xmu2);
                        c *= (d/i);
                        p /= (i-xmu);
                        q /= (i+xmu);
                        del=c*(ff+r*q);
                        sum += del;
                        del1=c*p-i*del;
                        sum1 += del1;
                        if (fabs(del) < (1.0+fabs(sum))*EPS) break;
                }
                if (i > MAXIT) nrerror("bessy series failed to converge");
                rymu = -sum;
                ry1 = -sum1*xi2;
                rymup=xmu*xi*rymu-ry1;
                rjmu=w/(rymup-f*rymu);
        } else {
                a=0.25-xmu2;
                p = -0.5*xi;
                q=1.0;
                br=2.0*x;
                bi=2.0;
                fact=a*xi/(p*p+q*q);
                cr=br+q*fact;
                ci=bi+p*fact;
                den=br*br+bi*bi;
                dr=br/den;
                di = -bi/den;
                dlr=cr*dr-ci*di;
                dli=cr*di+ci*dr;
                temp=p*dlr-q*dli;
                q=p*dli+q*dlr;
                p=temp;
                for (i=2;i<=MAXIT;i++) {
                        a += 2*(i-1);
                        bi += 2.0;
                        dr=a*dr+br;
                        di=a*di+bi;
                        if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
                        fact=a/(cr*cr+ci*ci);
                        cr=br+cr*fact;
                        ci=bi-ci*fact;
                        if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
                        den=dr*dr+di*di;
                        dr /= den;
                        di /= -den;
                        dlr=cr*dr-ci*di;
                        dli=cr*di+ci*dr;
                        temp=p*dlr-q*dli;
                        q=p*dli+q*dlr;
                        p=temp;
                        if (fabs(dlr-1.0)+fabs(dli) < EPS) break;
                }
                if (i > MAXIT) nrerror("cf2 failed in bessjy");
                gam=(p-f)/q;
                rjmu=sqrt(w/((p-f)*gam+q));
                rjmu=NRSIGN(rjmu,rjl);
                rymu=rjmu*gam;
                rymup=rymu*(p+q/gam);
                ry1=xmu*xi*rymu-rymup;
        }
        fact=rjmu/rjl;
        *rj=rjl1*fact;
        *rjp=rjp1*fact;
        for (i=1;i<=nl;i++) {
                rytemp=(xmu+i)*xi2*ry1-rymu;
                rymu=ry1;
                ry1=rytemp;
        }
        *ry=rymu;
        *ryp=xnu*xi*rymu-ry1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN

/*............................................................................*/
double bessi0_5(double x) {return (x==0)? 0:sqrt(2/(PI*x))*sinh(x);}
double bessi1_5(double x) {return (x==0)? 0:sqrt(2/(PI*x))*(cosh(x)-sinh(x)/x);}
double bessi2  (double x) {return (x==0)? 0:bessi0  (x)-(2*1  )/x*bessi1  (x);}
double bessi2_5(double x) {return (x==0)? 0:bessi0_5(x)-(2*1.5)/x*bessi1_5(x);}
double bessi3_5(double x) {return (x==0)? 0:bessi1_5(x)-(2*2.5)/x*bessi2_5(x);}
double bessj3_5(double x) {
   double rj, ry, rjp, ryp;
   bessjy(x,3.5, &rj, &ry, &rjp, &ryp);
   return rj;
}

/* Special functions ------------------------------------------------------- */
double gammln(double xx) {
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
	-1.231739516,0.120858003e-2,-0.536382e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
	x += 1.0;
	ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}


double betai(double a, double b, double x) {
    double bt;
    if (x < 0.0 || x > 1.0) nrerror("Bad x in routine BETAI");
    if (x == 0.0 || x == 1.0) bt=0.0;
    else
	bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0))
	return bt*betacf(a,b,x)/a;
    else
	return 1.0-bt*betacf(b,a,1.0-x)/b;

}

#define ITMAX 100
#define EPS 3.0e-7
double betacf(double a, double b, double x) {
    double qap,qam,qab,em,tem,d;
    double bz,bm=1.0,bp,bpp;
    double az=1.0,am=1.0,ap,app,aold;
    int m;

    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    bz=1.0-qab*x/qap;
    for (m=1;m<=ITMAX;m++) {
	em=(double) m;
	tem=em+em;
	d=em*(b-em)*x/((qam+tem)*(a+tem));
	ap=az+d*am;
	bp=bz+d*bm;
	d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
	app=ap+d*az;
	bpp=bp+d*bz;
	aold=az;
	am=ap/bpp;
	bm=bp/bpp;
	az=app/bpp;
	bz=1.0;
	if (fabs(az-aold) < (EPS*fabs(az))) return az;
    }
    nrerror("a or b too big, or ITMAX too small in BETACF");
}
#undef ITMAX
#undef EPS

void instantiate_recipes() {
   double **DD1;
   double *D1;

   double **FF1;
   double *F1;

   int **II1;   
   int *I1;   
   int i1;

   char *C1;

   ask_Tvector(D1,i1,i1); free_Tvector(D1,i1,i1);
   ask_Tvector(F1,i1,i1); free_Tvector(F1,i1,i1);
   ask_Tvector(I1,i1,i1); free_Tvector(I1,i1,i1);
   ask_Tvector(C1,i1,i1); free_Tvector(C1,i1,i1);

   ask_Tmatrix(DD1,i1,i1,i1,i1); free_Tmatrix(DD1,i1,i1,i1,i1);
   ask_Tmatrix(FF1,i1,i1,i1,i1); free_Tmatrix(FF1,i1,i1,i1,i1);
   ask_Tmatrix(II1,i1,i1,i1,i1); free_Tmatrix(II1,i1,i1,i1,i1);
}

/* Optimization ------------------------------------------------------------ */
#define TOL 2.0e-4

int ncom=0;	/* defining declarations */
double *pcom=NULL, *xicom=NULL;
double (*nrfunc)(double *)=NULL;

double f1dim(double x)
{
	int j;
	double f,*xt;

	ask_Tvector(xt,1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_Tvector(xt,1,ncom);
	return f;
}

#undef MAX
#undef SIGN

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax, double *bx, double *cx,
   double *fa, double *fb, double *fc, double (*func)(double))
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SHFT

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
   double *xmin)
{
   int iter;
   double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;

   a=(ax < cx ? ax : cx);
   b=(ax > cx ? ax : cx);
   x=w=v=bx;
   fw=fv=fx=(*f)(x);
   for (iter=1;iter<=ITMAX;iter++) { 
      xm=0.5*(a+b); tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) {*xmin=x; return fx;}
      if (fabs(e) > tol1) {
	 r=(x-w)*(fx-fv); q=(x-v)*(fx-fw);
	 p=(x-v)*q-(x-w)*r; q=2.0*(q-r);
	 if (q > 0.0) p = -p;
	 q=fabs(q);
	 etemp=e;
	 e=d;
	 if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=CGOLD*(e=(x >= xm ? a-x : b-x));
	 else {
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
	 }
      } else { d=CGOLD*(e=(x >= xm ? a-x : b-x));}
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=(*f)(u);
      if (fu <= fx) {
	 if (u >= x) a=x; else b=x;
	 SHFT(v,w,x,u) 
	 SHFT(fv,fw,fx,fu)
      } else {
	 if (u < x) a=u; else b=u;
	 if (fu <= fw || w == x) { v=w; w=u; fv=fw; fw=fu; }
	 else if (fu <= fv || v == x || v == w) { v=u; fv=fu; }
      } 
   }
   nrerror("Too many iterations in brent");
   *xmin=x;
   return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

void linmin(double *p, double *xi, int n, double &fret,
   double (*func)(double *))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	nrfunc=func;
	ask_Tvector(pcom,1,n);
	ask_Tvector(xicom,1,n);
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_Tvector(xicom,1,n);
	free_Tvector(pcom,1,n);
}

#undef TOL

#define ITMAX 200
static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

void powell(double *p, double *xi, int n, double ftol, int &iter,
   double &fret, double (*func)(double *), bool show)
{
	int i,ibig,j;
	double t,fptt,fp,del;
	double *pt,*ptt,*xit;
	bool   different_from_0;

	ask_Tvector(pt,1,n);
	ask_Tvector(ptt,1,n);
	ask_Tvector(xit,1,n);
	fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	
	for (iter=1;;(iter)++) {
	        /* By coss ----- */
	        if (show) {
	           cout << iter << " (" << p[1];
	           for (int co=2; co<=n; co++) cout << "," << p[co];
	           cout << ")--->" << fret << endl;
                }
	        /* ------------- */
	
		fp=fret;
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
      	             	different_from_0=FALSE; // CO
			for (j=1;j<=n;j++) {
			    xit[j]=xi[j*n+i];
			    if (xit[j]!=0) different_from_0=TRUE;
			}
      	             	if (different_from_0) {
			   fptt=fret;
			   linmin(p,xit,n,fret,func);
			   if (fabs(fptt-fret) > del) {
				   del=fabs(fptt-fret);
				   ibig=i;
			   }
	        	   /* By coss ----- */
	        	   if (show) {
	        	      cout << "	  (";
			      if (i==1) cout << "***";
			      cout << p[1];
	        	      for (int co=2; co<=n; co++) {
			          cout << ",";
			          if (co==i) cout << "***";
			          cout << p[co];
			      }
	        	      cout << ")--->" << fret << endl;
                	   }
	        	   /* ------------- */
      	             	}
		}
		if (2.0*fabs(fp-fret) <= ftol*(fabs(fp)+fabs(fret))) {
			free_Tvector(xit,1,n);
			free_Tvector(ptt,1,n);
			free_Tvector(pt,1,n);
			return;
		}
		if (iter == ITMAX) nrerror("Too many iterations in routine POWELL");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) xi[j*n+ibig]=xit[j];
			}
		}
	}
}

#undef ITMAX
#undef SQR

/* Non linear least squares ------------------------------------------------ */
// These routines have been taken from
// http://users.utu.fi/vesoik/userdocs/programs/libpet
// and they implement an algorithm of Lawson-Hanson of
// nonnegative least squares

/* Example of use:
   double a[]={ 5, 0, -2,
               0, 3,  0,
               1, 1, -1,
              -1, 1, -1,
               9, 9, -9};
   double b[]={1, 9, -1};
   double x[5];
   double rnorm;
   int i;
   
   int success=nnls(a,3,5,b,x,&rnorm,NULL,NULL,NULL);
   printf("success=%d\n",success);
   printf("rnorm=%d\n",rnorm);
   for (i=0; i<5; i++)
      printf("%f\n",x[i]);
   
   In this case: x=0 2.666 0 0 0.111
   
   This program resolves A^t*x=b subject to x>=0.
   In terms of basis vectors, the rows of A are the basis axes, b is
   the vector we want to represent in the subspace spanned by the rows of A
   and x are the nonnegative coordinates of the representation of b in A.
*/

/*****************************************************************************
 *
 *  Compute orthogonal rotation matrix:
 *    (C, S) so that (C, S)(A) = (sqrt(A**2+B**2))
 *    (-S,C)         (-S,C)(B)   (   0          )
 *  Compute sig = sqrt(A**2+B**2):
 *    sig is computed last to allow for the possibility that sig may be in
 *    the same location as A or B.
 */
void _nnls_g1(double a, double b, double *cterm, double *sterm, double *sig)
{
  double d1, xr, yr;

  if(fabs(a)>fabs(b)) {
    xr=b/a; d1=xr; yr=sqrt(d1*d1 + 1.); d1=1./yr;
    *cterm=(a>=0.0 ? fabs(d1) : -fabs(d1));
    *sterm=(*cterm)*xr; *sig=fabs(a)*yr;
  } else if(b!=0.) {
    xr=a/b; d1=xr; yr=sqrt(d1*d1 + 1.); d1=1./yr;
    *sterm=(b>=0.0 ? fabs(d1) : -fabs(d1));
    *cterm=(*sterm)*xr; *sig=fabs(b)*yr;
  } else {
    *sig=0.; *cterm=0.; *sterm=1.;
  }
} /* _nnls_g1 */
/****************************************************************************/

/*****************************************************************************
 *
 *  Construction and/or application of a single Householder transformation:
 *           Q = I + U*(U**T)/B
 *
 *  Function returns 0 if succesful, or >0 in case of erroneous parameters.
 *
 */
int _nnls_h12(
  int mode,
  /* mode=1 to construct and apply a Householder transformation, or
     mode=2 to apply a previously constructed transformation */
  int lpivot,     /* Index of the pivot element */
  int l1, int m,
  /* Transformation is constructed to zero elements indexed from l1 to M */
  double *u, int u_dim1, double *up,
  /* With mode=1: On entry, u[] must contain the pivot vector.
     On exit, u[] and up contain quantities defining the vector u[] of
     the Householder transformation. */
  /* With mode=2: On entry, u[] and up should contain quantities previously
     computed with mode=1. These will not be modified. */
  /* u_dim1 is the storage increment between elements. */
  double *cm,
  /* On entry, cm[] must contain the matrix (set of vectors) to which the
     Householder transformation is to be applied. On exit, cm[] will contain
     the set of transformed vectors */
  int ice,        /* Storage increment between elements of vectors in cm[] */
  int icv,        /* Storage increment between vectors in cm[] */
  int ncv         /* Nr of vectors in cm[] to be transformed;
                     if ncv<=0, then no operations will be done on cm[] */
) {
  double d1, d2, b, clinv, cl, sm;
  int incr, k, j, i2, i3, i4;

  /* Check parameters */
  if(mode!=1 && mode!=2) return(1);
  if(m<1 || u==NULL || u_dim1<1 || cm==NULL) return(2);
  if(lpivot<0 || lpivot>=l1 || l1>=m) return(0);
  /* Function Body */
  cl= (d1 = u[lpivot*u_dim1], fabs(d1));
  if(mode==2) { /* Apply transformation I+U*(U**T)/B to cm[] */
    if(cl<=0.) return(0);
  } else { /* Construct the transformation */
    for(j=l1; j<m; j++) { /* Computing MAX */
      d2=(d1=u[j*u_dim1], fabs(d1)); if(d2>cl) cl=d2;}
    if(cl<=0.) return(0);
    clinv=1.0/cl;
    /* Computing 2nd power */
    d1=u[lpivot*u_dim1]*clinv; sm=d1*d1;
    for(j=l1; j<m; j++) {d1=u[j*u_dim1]*clinv; sm+=d1*d1;}
    cl*=sqrt(sm); if(u[lpivot*u_dim1]>0.) cl=-cl;
    *up=u[lpivot*u_dim1]-cl; u[lpivot*u_dim1]=cl;
  }
  if(ncv<=0) return(0);
  b=(*up)*u[lpivot*u_dim1];
  /* b must be nonpositive here; if b>=0., then return */
  if(b>=0.) return(0);
  b=1.0/b; i2=1-icv+ice*lpivot; incr=ice*(l1-lpivot);
  for(j=0; j<ncv; j++) {
    i2+=icv; i3=i2+incr; i4=i3; sm=cm[i2-1]*(*up);
    for(k=l1; k<m; k++) {sm+=cm[i3-1]*u[k*u_dim1]; i3+=ice;}
    if(sm!=0.0) {
      sm*=b; cm[i2-1]+=sm*(*up);
      for(k=l1; k<m; k++) {cm[i4-1]+=sm*u[k*u_dim1]; i4+=ice;}
    }
  }
  return(0);
} /* _nnls_h12 */

/*****************************************************************************
 *  Algorithm NNLS (Non-negative least-squares)
 *
 *  Given an m by n matrix A, and an m-vector B, computes an n-vector X,
 *  that solves the least squares problem
 *      A * X = B   , subject to X>=0
 *
 *  Function returns 0 if succesful, 1, if iteration count exceeded 3*N,
 *  or 2 in case of invalid problem dimensions or memory allocation error.
 *
 *  Instead of pointers for working space, NULL can be given to let this
 *  function to allocate and free the required memory.
 */
int nnls(
  double *a, int m, int n,
  /* On entry, a[n][m] contains the m by n matrix A. On exit, a[][] contains 
     the product matrix Q*A, where Q is an m by n orthogonal matrix generated
     implicitly by this function.*/
  double *b,
  /* On entry, b[] must contain the m-vector B.
     On exit, b[] contains Q*B */
  double *x,
  /* On exit, x[] will contain the solution vector */
  double *rnorm,
  /* On exit, rnorm contains the Euclidean norm of the residual vector */
  double *wp,  /* An n-array of working space, w[]. */
  /* On exit, w[] will contain the dual solution vector.
     w[i]=0.0 for all i in set p and w[i]<=0.0 for all i in set z. */
  double *zzp, /* An m-array of working space, zz[]. */
  int *indexp  /* An n-array of working space, index[]. */
) {
  int pfeas, ret=0, iz, jz, iz1, iz2, npp1, *index;
  double d1, d2, sm, up, ss, *w, *zz;
  int iter, k, j=0, l, itmax, izmax=0, nsetp, ii, jj=0, ip;
  double temp, wmax, t, alpha, asave, dummy, unorm, ztest, cc;


  /* Check the parameters and data */
  if(m<=0 || n<=0 || a==NULL || b==NULL || x==NULL) return(2);
  /* Allocate memory for working space, if required */
  if(wp!=NULL) w=wp; else w=(double*)calloc(n, sizeof(double));
  if(zzp!=NULL) zz=zzp; else zz=(double*)calloc(m, sizeof(double));
  if(indexp!=NULL) index=indexp; else index=(int*)calloc(n, sizeof(int));
  if(w==NULL || zz==NULL || index==NULL) return(2);

  /* Initialize the arrays INDEX[] and X[] */
  for(k=0; k<n; k++) {x[k]=0.; index[k]=k;}
  iz2=n-1; iz1=0; nsetp=0; npp1=0;

  /* Main loop; quit if all coeffs are already in the solution or */
  /* if M cols of A have been triangularized */
  iter=0; itmax=n*3;
  while(iz1<=iz2 && nsetp<m) {
    /* Compute components of the dual (negative gradient) vector W[] */
    for(iz=iz1; iz<=iz2; iz++) {
      j=index[iz]; sm=0.; for(l=npp1; l<m; l++) sm+=a[j*m+l]*b[l];
      w[j]=sm;
    }

    while(1) {
      /* Find largest positive W[j] */
      for(wmax=0., iz=iz1; iz<=iz2; iz++) {
        j=index[iz]; if(w[j]>wmax) {wmax=w[j]; izmax=iz;}}

      /* Terminate if wmax<=0.; */
      /* it indicates satisfaction of the Kuhn-Tucker conditions */
      if(wmax<=0.0) break;
      iz=izmax; j=index[iz];

      /* The sign of W[j] is ok for j to be moved to set P. */
      /* Begin the transformation and check new diagonal element to avoid */
      /* near linear dependence. */
      asave=a[j*m+npp1];
      _nnls_h12(1, npp1, npp1+1, m, &a[j*m+0], 1, &up, &dummy, 1, 1, 0);
      unorm=0.;
      if(nsetp!=0) for(l=0; l<nsetp; l++) {d1=a[j*m+l]; unorm+=d1*d1;}
      unorm=sqrt(unorm);
      d2=unorm+(d1=a[j*m+npp1], fabs(d1)) * 0.01;
      if((d2-unorm)>0.) {
        /* Col j is sufficiently independent. Copy B into ZZ, update ZZ */
        /* and solve for ztest ( = proposed new value for X[j] ) */
        for(l=0; l<m; l++) zz[l]=b[l];
        _nnls_h12(2, npp1, npp1+1, m, &a[j*m+0], 1, &up, zz, 1, 1, 1);
        ztest=zz[npp1]/a[j*m+npp1];
        /* See if ztest is positive */
        if(ztest>0.) break;
      }

      /* Reject j as a candidate to be moved from set Z to set P. Restore */
      /* A[npp1,j], set W[j]=0., and loop back to test dual coeffs again */
      a[j*m+npp1]=asave; w[j]=0.;
    } /* while(1) */
    if(wmax<=0.0) break;

    /* Index j=INDEX[iz] has been selected to be moved from set Z to set P. */
    /* Update B and indices, apply householder transformations to cols in */
    /* new set Z, zero subdiagonal elts in col j, set W[j]=0. */
    for(l=0; l<m; ++l) b[l]=zz[l];
    index[iz]=index[iz1]; index[iz1]=j; iz1++; nsetp=npp1+1; npp1++;
    if(iz1<=iz2) for(jz=iz1; jz<=iz2; jz++) {
      jj=index[jz];
      _nnls_h12(2, nsetp-1, npp1, m, &a[j*m+0], 1, &up,
           &a[jj*m+0], 1, m, 1);
    }
    if(nsetp!=m) for(l=npp1; l<m; l++) a[j*m+l]=0.;
    w[j]=0.;
    /* Solve the triangular system; store the solution temporarily in Z[] */
    for(l=0; l<nsetp; l++) {
      ip=nsetp-(l+1);
      if(l!=0) for(ii=0; ii<=ip; ii++) zz[ii]-=a[jj*m+ii]*zz[ip+1];
      jj=index[ip]; zz[ip]/=a[jj*m+ip];
    }

    /* Secondary loop begins here */
    while(++iter<itmax) {
      /* See if all new constrained coeffs are feasible; if not, compute alpha */
      for(alpha=2.0, ip=0; ip<nsetp; ip++) {
        l=index[ip];
        if(zz[ip]<=0.) {t=-x[l]/(zz[ip]-x[l]); if(alpha>t) {alpha=t; jj=ip-1;}}
      }

      /* If all new constrained coeffs are feasible then still alpha==2. */
      /* If so, then exit from the secondary loop to main loop */
      if(alpha==2.0) break;
      /* Use alpha (0.<alpha<1.) to interpolate between old X and new ZZ */
      for(ip=0; ip<nsetp; ip++) {l=index[ip]; x[l]+=alpha*(zz[ip]-x[l]);}

      /* Modify A and B and the INDEX arrays to move coefficient i */
      /* from set P to set Z. */
      k=index[jj+1]; pfeas=1;
      do {
        x[k]=0.;
        if(jj!=(nsetp-1)) {
          jj++;
          for(j=jj+1; j<nsetp; j++) {
            ii=index[j]; index[j-1]=ii;
            _nnls_g1(a[ii*m+j-1], a[ii*m+j], &cc, &ss, &a[ii*m+j-1]);
            for(a[ii*m+j]=0., l=0; l<n; l++) if(l!=ii) {
              /* Apply procedure G2 (CC,SS,A(J-1,L),A(J,L)) */
              temp=a[l*m+j-1];
              a[l*m+j-1]=cc*temp+ss*a[l*m+j];
              a[l*m+j]=-ss*temp+cc*a[l*m+j];
            }
            /* Apply procedure G2 (CC,SS,B(J-1),B(J)) */
            temp=b[j-1]; b[j-1]=cc*temp+ss*b[j]; b[j]=-ss*temp+cc*b[j];
          }
        }
        npp1=nsetp-1; nsetp--; iz1--; index[iz1]=k;

        /* See if the remaining coeffs in set P are feasible; they should */
        /* be because of the way alpha was determined. If any are */
        /* infeasible it is due to round-off error. Any that are */
        /* nonpositive will be set to zero and moved from set P to set Z */
        for(jj=0; jj<nsetp; jj++) {k=index[jj]; if(x[k]<=0.) {pfeas=0; break;}}
      } while(pfeas==0);

      /* Copy B[] into zz[], then solve again and loop back */
      for(k=0; k<m; k++) zz[k]=b[k];
      for(l=0; l<nsetp; l++) {
        ip=nsetp-(l+1);
        if(l!=0) for(ii=0; ii<=ip; ii++) zz[ii]-=a[jj*m+ii]*zz[ip+1];
        jj=index[ip]; zz[ip]/=a[jj*m+ip];
      }
    } /* end of secondary loop */
    if(iter>itmax) {ret=1; break;}
    for(ip=0; ip<nsetp; ip++) {k=index[ip]; x[k]=zz[ip];}
  } /* end of main loop */
  /* Compute the norm of the final residual vector */
  sm=0.;
  if(npp1<m) for(k=npp1; k<m; k++) sm+=(b[k]*b[k]);
  else for(j=0; j<n; j++) w[j]=0.;
  *rnorm=sqrt(sm);
  /* Free working space, if it was allocated here */
  if(wp==NULL) free(w); if(zzp==NULL) free(zz); if(indexp==NULL) free(index);
  return(ret);
} /* nnls_ */
/****************************************************************************/
/****************************************************************************/
/*
  nnlsWght()

  Algorithm for weighting the problem that is given to nnls-algorithm.
  Square roots of weights are used because in nnls the difference
  w*A-w*b is squared.
  Algorithm returns zero if successful, 1 if arguments are inappropriate.

*/
int nnlsWght(int N, int M, double *A, double *b, double *weight)
{
  int n, m;
  double *w;

  /* Check the arguments */
  if(N<1 || M<1 || A==NULL || b==NULL || weight==NULL) return(1);

  /* Allocate memory */
  w=(double*)malloc(M*sizeof(double)); if(w==NULL) return(2);

  /* Check that weights are not zero and get the square roots of them to w[] */
  for(m=0; m<M; m++) {
    if(weight[m]<=1.0e-20) w[m]=0.0;
    else w[m]=sqrt(weight[m]);
  }
 
  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      A[n*M+m]*=w[m];
    }
    b[m]*=w[m];
  }

  free(w);
  return(0);
}
/****************************************************************************/

/* Singular value descomposition ------------------------------------------- */
static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

#define EPS 1e-40
void svdcmp(double *a,int m,int n,double *w, double *v)
{
	int flag,i,its,j,jj,k,l,nm;
	double c,f,h,s,x,y,z;
	double anorm=0.0,g=0.0,scale=0.0;
	double *rv1=NULL;

	if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
	ask_Tvector(rv1,1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k*n+i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k*n+i] /= scale;
					s += a[k*n+i]*a[k*n+i];
				}
				f=a[i*n+i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*n+i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k*n+i]*a[k*n+j];
						f=s/h;
						for (k=i;k<=m;k++) a[k*n+j] += f*a[k*n+i];
					}
				}
				for (k=i;k<=m;k++) a[k*n+i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i*n+k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i*n+k] /= scale;
					s += a[i*n+k]*a[i*n+k];
				}
				f=a[i*n+l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*n+l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i*n+k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j*n+k]*a[i*n+k];
						for (k=l;k<=n;k++) a[j*n+k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i*n+k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) {
                                        double den=a[i*n+l]*g;
                                        if (ABS(den)>EPS) v[j*n+i]=a[i*n+j]/den;
                                        else              v[j*n+i]=0.0;
					//COSS: v[j][i]=(a[i][j]/a[i][l])/g;
                                 }
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i*n+k]*v[k*n+j];
					for (k=l;k<=n;k++) v[k*n+j] += s*v[k*n+i];
				}
			}
			for (j=l;j<=n;j++) v[i*n+j]=v[j*n+i]=0.0;
		}
		v[i*n+i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i*n+j]=0.0;
		if (ABS(g)>EPS) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k*n+i]*a[k*n+j];
                                        if (ABS(a[i*n+i])>EPS) f=(s/a[i*n+i])*g;
                                        else                  f=0.0;
                                        // COSS: f=(s/a[i*n+i])*g;
					for (k=i;k<=m;k++) a[k*n+j] += f*a[k*n+i];
				}
			}
			for (j=i;j<=m;j++) a[j*n+i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j*n+i]=0.0;
		}
		++a[i*n+i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=PYTHAG(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j*n+nm];
							z=a[j*n+i];
							a[j*n+nm]=y*c+z*s;
							a[j*n+i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j*n+k]=(-v[j*n+k]);
				}
				break;
			}
			if (its == 60) nrerror("No convergence in 60 SVDCMP iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj*n+j];
					z=v[jj*n+i];
					v[jj*n+j]=x*c+z*s;
					v[jj*n+i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj*n+j];
					z=a[jj*n+i];
					a[jj*n+j]=y*c+z*s;
					a[jj*n+i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_Tvector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG
#undef EPS

void svbksb(double *u,double *w,double *v, int m,int n,double *b,double *x)
{
	int jj,j,i;
	double s,*tmp;

	ask_Tvector(tmp,1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i*n+j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j*n+jj]*tmp[jj];
		x[j]=s;
	}
	free_Tvector(tmp,1,n);
}

/* Savitzky-Golay filter coefficients. ------------------------------------- */
// This routine is used in Xmipp to perform Numerical Derivatives of equally spaced
// data
void savgol(double *c, int np, int nl, int nr, int ld, int m)
/* Returns in c[1..np], in wrap-around order a set of Savitzky-Golay
 filter coeficients. nl is the number of leftward (past) data points used,
 while nr is the number of rightward (future) data points,
 making the total number of data points used nl +nr +1.
 ld is the order of the derivative desired (e.g., ld = 0 for smoothed function).
 m is the order of the smoothing polynomial, also equal to the highest conserved
 moment; usual values are m = 2or m = 4. */
{
	int imj,ipj,j,k,kk,mm,*indx;
	double d,fac,sum,*a,*b;
	
	if (np < nl+nr+1 || nl < 0 || nr < 0 || ld > m||nl+nr < m)
		nrerror("SAVGOL: bad arguments");
	ask_Tvector(indx,1,m+1);
	ask_Tvector(a,1,(m+1)*(m+1));
	ask_Tvector(b,1,m+1);
	// Set up the normal equations of the desired least-squares fit. 
	for (ipj=0;ipj<=(m << 1);ipj++)
	{
		sum=(ipj ? 0.0 : 1.0);
		for (k=1;k<=nr;k++) sum += pow((double)k,(double)ipj);
		for (k=1;k<=nl;k++) sum += pow((double)-k,(double)ipj);
		mm=MIN(ipj,2*m-ipj);
   		for (imj = -mm;imj<=mm;imj+=2)
                   a[(1+(ipj+imj)/2)*(m+1)+1+(ipj-imj)/2]=sum;
   	}
   	// Solve them: LU decomposition.
   	ludcmp(a,m+1,indx,&d);
	
   	for (j=1;j<=m+1;j++) b[j]=0.0;
   	b[ld+1]=1.0;
   	// Right-hand side vector is unit vector, depending on which derivative we want.
	lubksb(a,m+1,indx,b); // Get one row of the inverse matrix.
	for (kk=1;kk<=np;kk++) c[kk]=0.0; // Zero the output array (it may be bigger than number of coeficients).
	for (k = -nl;k<=nr;k++)
	{
		sum=b[1]; // Each Savitzky-Golay coeficient is the dot product of powers of an integer with the inverse matrix row.
   		fac=1.0;
   		for (mm=1;mm<=m;mm++) sum += b[mm+1]*(fac *= k);
   		kk=((np-k) % np)+1; // Store in wrap-around order.
		c[kk]=sum;   
   	}
	
    free_Tvector(b,1,m+1);   
    free_Tvector(a,1,(m+1)*(m+1));
    free_Tvector(indx,1,m+1);
}

// Wavelets ----------------------------------------------------------------
void wt1(double a[], unsigned long n, int isign,
	void (*wtstep)(double [], unsigned long, int))
{
	unsigned long nn;

	if (n < 4) return;
	if (isign >= 0) {
		for (nn=n;nn>=4;nn>>=1) (*wtstep)(a,nn,isign);
	} else {
		for (nn=4;nn<=n;nn<<=1) (*wtstep)(a,nn,isign);
	}
}

void wtn(double a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(double [], unsigned long, int))
{
	unsigned long i1,i2,i3,k,n,nnew,nprev=1,nt,ntot=1;
	int idim;
	double *wksp;

	for (idim=1;idim<=ndim;idim++) ntot *= nn[idim];
	ask_Tvector(wksp,1,ntot);
	for (idim=1;idim<=ndim;idim++) {
		n=nn[idim];
		nnew=n*nprev;
		if (n > 4) {
			for (i2=0;i2<ntot;i2+=nnew) {
				for (i1=1;i1<=nprev;i1++) {
					for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) wksp[k]=a[i3];
					if (isign >= 0) {
						for(nt=n;nt>=4;nt >>= 1)
							(*wtstep)(wksp,nt,isign);
					} else {
						for(nt=4;nt<=n;nt <<= 1)
							(*wtstep)(wksp,nt,isign);
					}

					for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) a[i3]=wksp[k];
				}
			}
		}
		nprev=nnew;
	}
	free_Tvector(wksp,1,ntot);
}

typedef struct {
	int ncof,ioff,joff;
	double *cc,*cr;
} wavefilt;

wavefilt wfilt;

void pwtset(int n)
{
	void nrerror(char error_text[]);
	int k;
	float sig = -1.0;
	static double c4[5]={0.0,0.4829629131445341,0.8365163037378079,
			0.2241438680420134,-0.1294095225512604};
	static double c12[13]={0.0,0.111540743350, 0.494623890398, 0.751133908021,
		0.315250351709,-0.226264693965,-0.129766867567,
		0.097501605587, 0.027522865530,-0.031582039318,
		0.000553842201, 0.004777257511,-0.001077301085};
	static double c20[21]={0.0,0.026670057901, 0.188176800078, 0.527201188932,
		0.688459039454, 0.281172343661,-0.249846424327,
		-0.195946274377, 0.127369340336, 0.093057364604,
		-0.071394147166,-0.029457536822, 0.033212674059,
		0.003606553567,-0.010733175483, 0.001395351747,
		0.001992405295,-0.000685856695,-0.000116466855,
		0.000093588670,-0.000013264203};
	static double c4r[5],c12r[13],c20r[21];

	wfilt.ncof=n;
	if (n == 4) {
		wfilt.cc=c4;
		wfilt.cr=c4r;
	}
	else if (n == 12) {
		wfilt.cc=c12;
		wfilt.cr=c12r;
	}
	else if (n == 20) {
		wfilt.cc=c20;
		wfilt.cr=c20r;
	}
	else nrerror("unimplemented value n in pwtset");
	for (k=1;k<=n;k++) {
		wfilt.cr[wfilt.ncof+1-k]=sig*wfilt.cc[k];
		sig = -sig;
	}
	wfilt.ioff = wfilt.joff = -(n >> 1);
}

void pwt(double a[], unsigned long n, int isign)
{
	double ai,ai1,*wksp;
	unsigned long i,ii,j,jf,jr,k,n1,ni,nj,nh,nmod;

	if (n < 4) return;
	ask_Tvector(wksp,1,n);
	nmod=wfilt.ncof*n;
	n1=n-1;
	nh=n >> 1;
	for (j=1;j<=n;j++) wksp[j]=0.0;
	if (isign >= 0) {
		for (ii=1,i=1;i<=n;i+=2,ii++) {
			ni=i+nmod+wfilt.ioff;
			nj=i+nmod+wfilt.joff;
			for (k=1;k<=wfilt.ncof;k++) {
				jf=n1 & (ni+k);
				jr=n1 & (nj+k);
				wksp[ii] += wfilt.cc[k]*a[jf+1];
				wksp[ii+nh] += wfilt.cr[k]*a[jr+1];
			}
		}
	} else {
		for (ii=1,i=1;i<=n;i+=2,ii++) {
			ai=a[ii];
			ai1=a[ii+nh];
			ni=i+nmod+wfilt.ioff;
			nj=i+nmod+wfilt.joff;
			for (k=1;k<=wfilt.ncof;k++) {
				jf=(n1 & (ni+k))+1;
				jr=(n1 & (nj+k))+1;
				wksp[jf] += wfilt.cc[k]*ai;
				wksp[jr] += wfilt.cr[k]*ai1;
			}
		}
	}
	for (j=1;j<=n;j++) a[j]=wksp[j];
	free_Tvector(wksp,1,n);
}

/* Instantiantion ---------------------------------------------------------- */
template <class T>
   void instantiate_Numerical_T(T t) {
      T  *a, *d;
      int n, *indx;
      ludcmp(a, n, indx, d);
      lubksb(a, n, indx, d);
      gaussj(a, n, a, n);
}

void instantiate_Numerical() {
   short          s; instantiate_Numerical_T(s);
   char           h; instantiate_Numerical_T(h);
   int            i; instantiate_Numerical_T(i);
   float          f; instantiate_Numerical_T(f);
   double         d; instantiate_Numerical_T(d);
}
