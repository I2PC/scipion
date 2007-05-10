/***************************************************************************
 *
 * Authors:     Javier Rodríguez Falces (jrodriguez@cnb.uam.es)
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

#include "integration.h"

#define JMAXP 14
#define K 5

//**********************************************************
// Implementation of the integral using the Trapeze method
//**********************************************************

double Trapeze::operator()() { 		//adapted from qtrap
	double s,olds;
	int j;
	olds = -1.0e30;
		for (j=1;j<=JMAX;j++){
		s=Trap(j); 					//changed; Trap is integrating fcn
		if (fabs(s-olds) <= EPS*fabs(olds)) return s;
		if(s==0.0 && olds == 0.0 && j>6 ) 	return s;
		olds=s;
		}
	printf("Too many steps in routine qtrap_y\n");
	exit(1);
	return 0.0;
}


double Trapeze::Trap(int n) { //adapted from trapzd
	double tnm,sum,del;
	int j,it;
	if (n == 1){
		it=1;
		x=a; s=func(); 			//changed
		x=b; s+=func(); 		//changed
		return (s*=0.5*(b-a));
	}else{
		for (it=1,j=1;j<n-1;j++) it <<=1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum+=func(); //changed
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}



//**********************************************************
// Implementation of the integral using the Romberg method
//**********************************************************

double Romberg::operator()(){ 	//adapted from qromb
	int j;
	double ss,dss,h[JMAXP+2],s[JMAXP+2];
	h[1]=1.0;
	for(j=1;j<=JMAXP;j++){
		s[j]=midpnt(j); 		//changed; midpnt is integrating
		if(j>=K){ 				//function
			polint(&h[j-K],&s[j-K],K,0.0,ss,dss);
			if(fabs(dss)<=EPS*fabs(ss)) return ss;
		}
	s[j+1]=s[j];
	h[j+1]=h[j]/9.0;
	}
	cout << "Too many steps in routine Romberg" << endl; exit(1);
	return 0.0;
}

//*
// The midpnt function is used in the Romberg::operator only
//*
double Romberg::midpnt(int n){ 		//adapted from midpoint
	double tnm,sum,del,ddel;
	int it,j;
	if(n==1){
		x=0.5*(a+b); 				//changed; set x
		return (s=(b-a)*func()); 	//changed; evaluate func
	}else{
		for(it=1,j=1;j<n-1;j++) it *=3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for(j=1;j<=it;j++){
		sum+=func(); 				//changed; evaluate func
		x+=ddel; 					//changed; set x
		sum+=func(); 				//changed; evaluate func
		x+=del; 					//changed; set x
		}
	s=(s+(b-a)*sum/tnm)/3.0;
	return s;
	}
}

//*
// The polint function is used in the Romberg::operator only
//*
	double *my_vector(int nl,int nh);
	
	void polint(double xa[],double ya[],int n,double x,double &y,double &dy)
	{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;
	dif=fabs(x-xa[1]);
	c=my_vector(1,n);
	d=my_vector(1,n);
	for(i=1;i<=n;i++){
		if((dift=fabs(x-xa[i])) < dif){
		ns=i;
		dif=dift;
		}
	c[i]=ya[i];
	d[i]=ya[i];
	}
	y=ya[ns--];
	for(m=1;m<n;m++){
		for(i=1;i<=n-m;i++){
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if((den=ho-hp)==0.0){
			printf("error in routine polint\n");
			exit(1);
			}
		den=w/den;
		d[i]=hp*den;
		c[i]=ho*den;
		}
	y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
delete [] d; //replaces free_my_vector
delete [] c; //replaces free_my_vector
}


void error(const char *s)
{
	cerr << endl << s << endl; //1st "endl" is in case program is printing
	cout << endl << s << endl; //something out when the error occurs
	cout.flush(); //write the output buffer to the screen
	//or wherever the output of cout goes.
	abort();
}

#define NR_END 1

double *my_vector(int nl,int nh)
// allocate a double my_vector with subscript range v[nl .. nh]
{
	double *v;
	v=new double[nh-nl+1+NR_END];
	if(!v) error("allocation failure in my_vector1()");
	return v-nl+NR_END;
}


// Multidimensional integration --------------------------------------------
matrix1D<double>  cuhreX0;
matrix1D<double>  cuhreXX;
matrix1D<double>  cuhreRange;
integrand_t       cuhreIntegrand=NULL;

void scaledIntegrand(const int *ndim, const double xx[],
   const int *ncomp, double ff[]) {
   FOR_ALL_ELEMENTS_IN_MATRIX1D(cuhreXX)
     cuhreXX(i)=cuhreX0(i)+cuhreRange(i)*xx[i];
   (*cuhreIntegrand)(ndim,MULTIDIM_ARRAY(cuhreXX),ncomp,ff);
}

double multidimensionalIntegral(const matrix1D<double> &x0,
   const matrix1D<double> &xF, integrand_t integrand) {
   // Set some global variables
   cuhreX0=x0;
   cuhreRange=xF-x0;
   cuhreXX.init_zeros(x0);
   cuhreIntegrand=integrand;

   // Set Cuhre specific variables
   const double EPSREL  = 1e-3;
   const double EPSABS  = 1e-12;
   const int    VERBOSE = 0;
   const int    MINEVAL = 0;
   const int    MAXEVAL = 50000;
   const int    NCOMP   = 1;
   const int    KEY     = 0;
   double integral[NCOMP], error[NCOMP], prob[NCOMP];
   int NDIM=XSIZE(x0);
   
   // Compute integral
   int nregions, neval, fail;
   Cuhre(NDIM, NCOMP, scaledIntegrand,
      EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
      KEY,
      &nregions, &neval, &fail, integral, error, prob);

   // Take into account that Cuhre routine computes the integral
   // between 0 and 1. Thus, the result must be multiplied by the Jacobian
   // of a scaling
   double jacobian=1;
   FOR_ALL_ELEMENTS_IN_MATRIX1D(x0) jacobian*=(xF(i)-x0(i));
   return integral[0]*jacobian;
}
