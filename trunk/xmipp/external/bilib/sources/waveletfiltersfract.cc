/* ---------------------------------------------------------------------------

Filename     : WaveletFilters_fract.c    

Author       : Olivier Saudan

Organization : EPFL, Biomedical Imaging Group

Date         : June 1999                 

---------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>

#include "../configs.h"
#include "../headers/error.h"
#include "../headers/waveletfiltersfract.h"
#include "../headers/messagedisplay.h"

#define EPSILON 0.000001

#define ERR_MEM -1

static void Convolve121( double R[] , int l );
static int FFT( double R[] , double I[] , int m , int* errcode );
static double sinc( double x );
static void SampleTau( double alpha , double dzeta , int N , double R[] , int nbsamples );
static double EstimateSum( double f , double beta , int N );
static int getFilterSize(double alpha);
static int getSumTerms(double alpha);

/* ----------------------------------------------------------------------------
	Function:	WaveletFiltersGetSize_Fract
	
				empirical relation
---------------------------------------------------------------------------- */
extern int WaveletFiltersGetSize_Fract(double Alpha, long *nh, long *ng)
{

/*	*ng=64L; 
	*nh=64L;

*/
	*ng=*nh=getFilterSize(fabs(Alpha))+1;	
	return !ERROR;
}


/* ----------------------------------------------------------------------------
	Function:	WaveletFiltersGetCoef_Fract

	Purpose: 	Extern function to get the coefficient of the filter,
				Special routine for Fractional Orthonormal Spline
 				
	Parameters:	
	 			h:		output value, array of coefficient lowpass filter
	 			g:		output value, array of coefficient highpass filter
 			

---------------------------------------------------------------------------- */

extern int WaveletFiltersGetCoef_Fract(double Alpha, double *h, double *g) 
{
static double *R, *I;
double dzeta,mulfact;
int i,n,m,sumTerms,size;
unsigned long t;

	Alpha=fabs(Alpha);
        size=getFilterSize(Alpha);
	m=(int)ceil(log((long double)size)/log(2.0))+1;
    n=1<<m;
    while(2*n<size)
    { /* So that we're sure we have enough points ! */
    	m++;
    	n<<=1;
    }

    R=(double*)malloc(t=sizeof(*R)*n);
    if(!R)
    {
    	char st[128];
    	
    	sprintf( st , "fractfilters.c: getOrthoCoef(): cannot allocate %d bytes of"
    			 "memory for the R buffer!" , (int)t );
    	MessageDisplay( st );
    	return ERROR;
    }

    I=(double*)malloc(t=sizeof(*I)*n);
    if(!I)
    {
    	char st[128];
    	
    	sprintf( st , "fractfilters.c: getOrthoCoef(): cannot allocate %d bytes of"
    			 "memory for the I buffer!" , (int)t );
    	MessageDisplay( st );
    	free(R);
    	return ERROR;
    }

    sumTerms=getSumTerms(Alpha);

	if(Alpha>=1.)
		dzeta=Alpha-1.;
	else 
		dzeta=Alpha+1.;

		
    SampleTau( Alpha , dzeta , sumTerms , R , n );

    for(i=0;i<n;i++)
    	I[i]=0.;

    FFT(R,I,m,NULL);
   
   		
    if(Alpha>dzeta)
    {	/* Must convolve with (1+z)^2 */
    	R[size]=0;
    	Convolve121( R , size+1 );
    }	

    mulfact=1./sqrt(2.0)/(double)n;

    for(i=0;i<size+1;i++)
    {
		h[i]=R[i]*mulfact;
		if(i&1) /* odd */
			g[i]=h[i];
		else /* even */
			g[i]=-h[i];
    }

   	free(I);
   	free(R);
   	
	return !ERROR;

}

static int getFilterSize(double alpha)
{
	if(alpha<1) return 64;
    return (int)ceil(8+9*alpha);  /* 10^-6 */
}

static int getSumTerms(double alpha)
{
	alpha;
	return 20;
}


/******************************************************************************

  fft.c -  Fast Fourier Transforms


  Author: Olivier Saudan [OSA] <o.saudan@epfl.ch> - BIG/IOA/DMT/EPFL

  History:                                                 
  		March 23, 1999 - Ver 1.0 : [OSA] Created, FFT()

******************************************************************************/

/*-----------------------------------------------------------------------------

  FFT - Computes the FFT of a vector with 2^m points
        Source: "Signaux et systèmes", F.Pellandini, EPFL

  Author: Olivier Saudan <o.saudan@epfl.ch>

  Parameters:
 		R[] : real part of input and output vector
 		I[] : imaginary part of input and output vector
 		m : 2^m=number of points in vector
 		errcode : if not null, the error code is put there

  Returns : error code

  Error codes:
  		!ERROR (0): no error
        ERR_MEM (-1): caanot allocate memory
*/
int FFT( double R[] , double I[] , int m , int* errcode )
{
	double* WnR, *WnI; /* WnR[k]+i*WnI[k] = exp(i*2*Pi/n)^k */
    double tr,ti; /* temporary complex value */
    int n; /* number of points (2^m) */
    int i,j,k,stepsize,shifter;
    int i_j,i_j_s; /* i+j, i+j+stepsize */

    n=1<<m; /* n=2^m */

    /* We allocate both WnR and WnI tables in one step */
    if( (WnR=(double*) malloc( 2*n*sizeof(double) ))==0 )
    {
		if(errcode) *errcode=ERR_MEM;
        return ERR_MEM;
    }
    WnI=WnR+n;

    /* Now we compute the omegan coeffs */
    /* WnR[0]=1.; WnI[0]=0.;*/ /* The first is easy, but is never used */
    for(i=1;i<n;i++)
    {
    	/* We prefer to compute each coeff using explicit cos and sin
        instead of calculating  Wn[1]^i, so that we don't loose precision.
        In fact, that's why we store these coeffs in a table */
	    double arg;

        arg=2.*PI*(double)i / (double)n;
    	WnR[i] =  cos( arg );
        WnI[i] = -sin( arg );
    }

	/* First, we must swap the points so that the order matches our algorithm.
    e.g: for 8 points, the order becomes 0 4 2 6 1 5 3 7 */
    /* Source: Maple V R5 */
    for(i=j=0;i<n-1;i++)
    {
    	if(i<j)
        {   /* Swap points i and j */
			tr=R[i]; ti=I[i];
            R[i]=R[j]; I[i]=I[j];
            R[j]=tr; I[j]=ti;
        }
        k=n>>1;
        while(k<=j)
        {
        	j-=k;
            k/=2;
        }
        j+=k;
    }

    /* Perform the FFT */
    for(stepsize=1,shifter=m-1;stepsize<n;stepsize<<=1,--shifter)
    {
    	/* now we have n/2/stepsize partial FFT steps */
		for(j=0;j<n;j+=stepsize<<1)
        {   /* The j index is the start of the partial FFT */
        	for(i=0;i<stepsize;i++)
            {
            	i_j=i+j;
                i_j_s=i_j+stepsize;
                /* Must multiply the term j+i+stepsize with omegan^(i*2^shifter)*/
                /* (omegan[0])=1, so we test i to avoid useless mutliplications */
            	if(i)
                {
                	tr = WnR[i<<shifter]*R[i_j_s] -
                    	 WnI[i<<shifter]*I[i_j_s] ;
                    I[i_j_s] = WnR[i<<shifter]*I[i_j_s] +
                    		   WnI[i<<shifter]*R[i_j_s] ;
					R[i_j_s] = tr;
                }
                /* Now the "Butterfly" operator */
                tr = R[i_j] - R[i_j_s];
                ti = I[i_j] - I[i_j_s];
                R[i_j] += R[i_j_s];
                I[i_j] += I[i_j_s];
                R[i_j_s] = tr;
                I[i_j_s] = ti;
            }
        }
    }
    free(WnR);
    if(errcode) *errcode=!ERROR;
    return !ERROR;
}


/*
void main()
{
	double tr[16] = { 4, -2 , 3 , 5 , 10 , 3.1 , -7 , -3.14 , 2 , 7 , -4 , 3 , 1 , -1.806 , 0.707 , -0.707 };
    double ti[16] = { 2 , 1 , -4.31 , 0 , 0 , 8 , 12 , -12 , -1 , -0.345 , 78.2 , 2.4 , 3.14159 , 3.14159 , -0.707 , -0.707 };
    int i;

    FFT(tr,ti,4,NULL);

    for(i=0;i<16;i++)
    	printf( "%f " , (double)tr[i] );
    printf( "\n" );
    for(i=0;i<16;i++)
    	printf( "%f " , (double)ti[i] );
    printf( "\n" );
} */

/******************************************************************************

  sample.c -  Sample of the Tau function


  Author: Olivier Saudan [OSA] <o.saudan@epfl.ch> - BIG/IOA/DMT/EPFL

  History:
  		March 23, 1999 - Ver 1.0 : [OSA] Created, EstimateSum(), sinc(),
        	SampleTau()

******************************************************************************/
/*-----------------------------------------------------------------------------

  EstimateSum - Estimates the infinite sum

	                 inf
  				  -------
        	       \            1             1
            	    \       ---------  +  ---------
	                /            beta          beta
    	           /        (k+f)         (k-f)
        	      -------
            	    k=1

  Author: Olivier Saudan <o.saudan@epfl.ch>

  Parameters:
        f: (-1,1)
        beta: 2*alpha+2
        N: number of terms used to compute the sum. Time = O(N)

  Returns : the estimated value of the sum

  Note: Mathematically, the sum is divided into two sums. Each sum is
  		estimated using a correcting factor to speed up convergence.
        The two correcting factors have been choosen so that the error
        is minimal for f=[0,0.5]

*/
static double EstimateSum( double f , double beta , int N )
{
	double sum;
    int k;

    /* These are the two correcting factors */
    sum = 1./( (beta-1.)*pow( (double)N+f , beta-1. ) ) +
    	  1./( (beta-1.)*pow( (double)(N+1)-f , beta-1. ) ) ;

    for(k=1;k<=N;k++)
    	sum += 1./pow((double)k+f,beta) + 1./pow((double)k-f,beta);

    return sum;
}

/*-----------------------------------------------------------------------------

  sinc - Computes sin(x)/x, calculating the limit when x is small

  Author: Olivier Saudan <o.saudan@epfl.ch>

  Note: When x is small, sin(x)/x -> 1 - 1/6 x^2 + O(x^4)
		For x<EPSILON, x^2 is very small compared to 1, so that in this case,
        no calculus is made; 1 is returned.
*/
static double sinc( double x )
{
	if(x<EPSILON)
    	return 1.;
    else
    	return sin(x)/x;
}


/*-----------------------------------------------------------------------------

  SampleTau - Computes the samples of Tau

  Author: Olivier Saudan <o.saudan@epfl.ch>

  Parameters:
  		alpha : order of the fractional B-spline
        dzeta : order of the h*(w) filter (usually alpha-1)
        N : number of terms in sum (usually 20)
 		R[] : table where samples are stored (length of nbsamples)
		nbsamples : number of samples

  Returns : -

  Notes: Tau(f) = Omega(alpha,f) * 2 * Khi(f)^dzeta

                __________     |       -iw | dzeta
               / APhi(f)  |    | 1 + e     |
  Tau(f) = _  / ---------   2  |-----------|
            \/   APhi(2f)      |     2     |

          \______  _______/   \_____  _____/
                 \/                 \/
              Omega(f)            Khi(f)

  Note that if dzeta=alpha+1, Tau is the Fourier transform of the wavelet
  We usually "remove" Khi^2 from it (this means dzeta=alpha-1) and convolve
  the inverse fourier transform of Tau with {1/4;1/2;1/4}
*/
void SampleTau( double alpha , double dzeta , int N , double R[] , int nbsamples )
{
	double beta,PiFact,f;
	int i;

    /* !!! We use the fact that:
       - APhi(0)=1
       - APhi(f+k)=APhi(f) (1-Periodic)
       - APhi(f)=APhi(-f) (Even function)
       !!! */

    /* These will be used often */
    beta=2.*alpha+2.;
    PiFact=pow(PI,-beta);

    /* We first fill R[i] with APhi(f=i/nbsamples) */
    R[0] = 1.; /* The first is easy... */

    /* We use symmetry, that is APhi(1-f)=APhi(f), that is
    R[nbsamples-i] = R[i] for i=[1,nbsamples-1] */
    for(i=1;i<=nbsamples/2;i++)
    {
    	f=(double)i/(double)nbsamples;
    	R[nbsamples-i] = R[i] =
        	pow( fabs(sinc(PI*f)) , beta ) +
        	pow( fabs(sin(PI*f)) , beta ) * PiFact *
        	EstimateSum(f,beta,N) ;
    }

	/* Now we calculate sqrt(APhi(f)/APhi(2f)) for f=[0,0.5), that is
    sqrt(R[i]/R[2*i]) for i=0..nbsamples/2-1. NB: For i=0, the result is 1 ! */
    for(i=1;i<nbsamples/2;i++)
    	R[i] = sqrt( R[i] / R[i<<1] );

    /* For f=0.5, Omega=sqrt(APhi(f)) since APhi(2f)=1 */
    R[i] = sqrt(R[i]);

    /* Now we multiply with 2*khi^dzeta for f=[0,0.5]*/
    for(i=0;i<=nbsamples/2;i++)
	{
    	f = (double)i/(double)nbsamples;
        R[i] *= 2.*pow( .5*sqrt(2.+2.*cos(2.*PI*f)) , dzeta );
    }

    /*f We have calculated Tau for f=[0,0.5]. Since Tau is symmetric, we
    copy the values for f=(0.5,1) */
    for(i=1;i<nbsamples/2;i++)
    	R[nbsamples-i] = R[i] ;

    /* Done ! */
}


/* Convolution of the filter with {1/4 1/2 1/4}

   Author: Olivier Saudan <o.saudan@epfl.ch>

   Parameters:
   		R[] : real part of filters (input and output)
        l : length of convolution = (n/2)+1 if n=size of FFT result

   Note:
   		The convolution has been written to operate on the half of the filter
        If you want to get the whole filter (symmetry), you have to copy, using
        for example this code:

        for(i=1;i<n/2;i++)
        	R[n-i]=R[i];
*/
void Convolve121( double R[] , int l )
{
	double old,res;
    int i;

    if(l<2) return;

    /* "R[-1]"=R[1] */
    old=R[1];
    for(i=0;i<l-1;i++)
    {
    	res=(old+R[i+1])/4+R[i]/2;  /* Convolution with {1/4;1/2;1/4} */
        old=R[i];
        R[i]=res;
    }

    /* The last one (i=l-1), "R[l]"=0 */
    R[i]=(old)/4+R[i]/2;
}

