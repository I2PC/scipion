/***************************************************************************
 *
 * Authors:     Irene Martinez
 *              Roberto Marabini
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

/***************************************************************************/
/* Program for computing the FFT on complex data. The input arrays are x   */
/* (real part) and y (imaginary part). They are overwritten by the result  */
/* of the FFT. Split-radix, decimation in frequency algorithm. For more    */
/* details, see Sorensen et al., "On computing the split radix FFT". IEEE  */
/* Trans. on Acoustics, Speeech and Signal Processing, Vol. ASSP-34, No. 1 */
/* pp  152-156.                                                            */
/* If: sentido == 0 ===> direct FFT                                        */
/* else             ===> inverse "                                         */
/***************************************************************************/
/* This version is intended to be used with column arrays obtained in the  */
/* 2D FFT process. The parameters x & y are two arrays of pointers to      */
/* float, so that x[i][8] is the 8th column of the array.            */
/* The second index remains constant, while the 1st one is varied, thus    */
/* performing a "vertical" FFT on the data.                                */
/* Related files: fft.c is the "normal" fft program                        */
/***************************************************************************/
/* Juan P. Secilla (IBM-MSC)   Oct 86                                      */
/***************************************************************************/


#define TURBO 0         /* Set to 1 for Turbo-C */
#define MSC   1         /* Set to 1 for MSC     */

#define TWOPI 6.2831853071719586
#define SQRT2 1.41421356
#define MAX_LEN 1024

#include "math.h"
#include "stdio.h"
#include "malloc.h"
#include "memory.h"
#include "groe.h"

#ifdef __STDC__
  void fftcol(float **,int,float **,int,int,int,int);
  int image_fft (float **, int, int, int);
  int no_correct (int , int *);
  void rfft ( float *,int,int,int);
  void cfft (float *, float *, int,int,int);
#else
  void fftcol();
  int image_fft ();
  int no_correct ();
  void rfft ();
  void cfft ();
#endif



void fftcol(float **x,int col1,float **y,int col2,int n,int m,int sentido)

{
      /* x:  Real part array and column position */
       /* y:  Imaginary part array and column position */
    register int i0, i1;
    int n2, n4, i2, i3, id, is, i, j, k, size;
    float e, a, r1, r2, s1, s2, s3, xt, cc1, cc3, ss1, ss3;
    static float *ccc1[10], *sss1[10], *ccc3[10], *sss3[10];
    static int last_size = 0, last_m = 0;


           /************************************************************/
           /* As this program was originally written in Fortran, it's  */
    x--;   /* necessary to use this trick because in C the array ind-  */
    y--;   /* exes begin always at 0. Now the program will reference   */
           /* the elements correctly.                                  */
           /************************************************************/

    if (n != last_size)    /* Transform size has changed */
    {   if (last_size != 0)
            for (k = 1; k < last_m; k++)
            {   free ((char *) sss1[k]);
                free ((char *) sss3[k]);   /* Free old tables */
                free ((char *) ccc1[k]);
                free ((char *) ccc3[k]);
            }
        size = n>>2;
        for (k = 1; k < m; k++, size>>=1) /* Get space for new ones */
        {   sss1[k] = (float *) calloc (size + 1, sizeof (float));
            sss3[k] = (float *) calloc (size + 1, sizeof (float));
            ccc1[k] = (float *) calloc (size + 1, sizeof (float));
            ccc3[k] = (float *) calloc (size + 1, sizeof (float));
            if (sss1[k] == NULL || sss3[k] == NULL || ccc1 [k] == NULL ||
                ccc3[k] == NULL)
            {   puts ("Error: no hay memoria (fft)");
                exit (1);
            }
        }
        last_size = n;     /* Update last size */
        last_m = m;
        n2 = n<<1;
        for (k = 1; k < m; k++)
        {   n2 >>=1;
            n4 = n2>>2;
            e = 6.283185307179586/n2;
            a = 0.;
            for (j = 1; j <= n4; j++)
            {   ccc1[k][j] = cos (a);
                sss1[k][j] = sin (a);  /* Compute new tables */
                ccc3[k][j] = cos (3.*a);
                sss3[k][j] = sin (3.*a);
                a = j * e;
            }
        }
    }

    if (sentido != 0)
        for (i0 = 1; i0 <= n; i0 ++)
            y[i0][col2] = -y[i0][col2];             /* If inverse, conjugate input */

    n2 = n<<1;
    for (k = 1; k < m; k++)
    {   n2 >>=1;
        n4 = n2>>2;
        for (j = 1; j <= n4; j++)
        {   is = j;
            id = n2<<1;
            cc1 = ccc1[k][j];
            cc3 = ccc3[k][j];
            ss3 = sss3[k][j];
            ss1 = sss1[k][j];
            do {
                for (i0 = is; i0 < n; i0 += id)
                {   i1 = i0 + n4;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    r1 = x[i0][col1] - x[i2][col1];
                    x[i0][col1] += x[i2][col1];
                    r2 = x[i1][col1] - x[i3][col1];
                    x[i1][col1] += x[i3][col1];
                    s1 = y[i0][col2] - y[i2][col2];
                    y[i0][col2] += y[i2][col2];
                    s2 = y[i1][col2] - y[i3][col2];
                    y[i1][col2] += y[i3][col2];
                    s3 = r1 - s2;
                    r1 += s2;
                    s2 = r2 - s1;
                    r2 += s1;
                    x[i2][col1] = r1*cc1 - s2*ss1;
                    y[i2][col2] = -s2*cc1 - r1*ss1;
                    x[i3][col1] = s3*cc3 + r2*ss3;
                    y[i3][col2] = r2*cc3 - s3*ss3;
                }
                is = (id<<1) - n2 + j;
                id <<=2;
            } while (is < n);
        }
    }
/***************** last stage, length-2 butterfly ****************************/
    is = 1;
    id = 4;
    do {
        for (i0 = is; i0 <= n; i0 += id)
        {   i1 = i0 + 1;
            r1 = x[i0][col1];
            x[i0][col1] = r1 + x[i1][col1];
            x[i1][col1] = r1 - x[i1][col1];
            r1 = y[i0][col2];
            y[i0][col2] = r1 + y[i1][col2];
            y[i1][col2] = r1 - y[i1][col2];
         }
         is = (id<<1) - 1;
         id <<=2;
    } while (is < n);
/*************************** bit reverse counter ****************************/
    j = 1;
/*    n1 = n - 1;    */
    for (i = 1; i < n; i++)
    {   if (i < j)
        {   xt = x[j][col1];
            x[j][col1] = x[i][col1];
            x[i][col1] = xt;
            xt = y[j][col2];
            y[j][col2] = y[i][col2];
            y[i][col2] = xt;
        }
        k = n>>1;
        while (k < j)
        {   j -= k;
            k >>=1;
        }
        j += k;
    }

    if (sentido != 0)
        for (i0 = 1; i0 <= n; i0 ++)
            y[i0][col2] = -y[i0][col2];             /* If inverse, conjugate output */
}

/*************************************************************************/
/* This program performs the FFT on an image. It uses rfft.c to perform  */
/* row transforms and fft.c (J. P. Secilla) to perform column transforms.*/
/* The original image must be at 1 byte/pixel, be real and its dimensions*/
/* have to be an integral power of two.                                  */
/*************************************************************************/
/*        Juan P. Secilla        (IBM/MSC)     Nov. 86                   */
/*************************************************************************/


/****************************************************************************/
/* This routine performs a direct FFT algorithm on the array im_char (1 byte*/
/* /pixel image) and leaves the result in the x array (Fourier transform    */
/* array). Intended only for real images.                                   */
/* Returns ERROR if the dimensions are not an integral power of two,        */
/* OK if everything has worked as planned                                   */
/****************************************************************************/

int image_fft (float **x, int row, int col, int kind)

                    /* x: Fourier transform array */
                    /* row, col: Image dimensions */
                    /* kind: of transform to be performed */

{
int power;                  /* 2**power = nrow or ncol */
int i,j;

/******************* Check image dimensions *********************************/

if (no_correct (row, &power))
    return ERROR;
if (no_correct (col, &power))
    return ERROR;

/*=========================================================================*/
/************************ Perform FFT now **********************************/
/*=========================================================================*/

switch (kind) {
    case DIRECT:            /* direct FFT */

/************************ Perform FFT row by row ***************************/
        for (i = 0; i < row; i++)
            rfft (x[i], col, power, 0);

/*********************** Perform FFT column by column **********************/

        for (i = 0; i <= col/2; i++)
            fftcol (x, 2*i, x, 2*i+1, row, power, 0);
        /* Only half of the columns are required */
        break;

    case INVERSE:       /******* inverse FFT ********/

/*************************** Perform column transform *********************/

        for (i = 0; i <= col/2; i++)
            fftcol (x, 2*i, x, 2*i+1, row, power, 1);

/************************** Perform row transform and rescale ************/

        for (i = 0; i < row; i++)
        {   rfft (x[i], col, power, 1);
            for (j=0; j < col; j++)
                x[i][j] /= row;
        }

/************************** The FFT process has ended at last ***********/
}
return OK;
}

/****************************** End of the FFT routine *********************/

/***************************************************************************/
/* NO_CORRECT (VALUE): returns T if value is not an integral power of 2, F */
/* if it is. (A macro would have also been possible)                       */
/* The power number is returned in the parameter "power".                  */
/***************************************************************************/

int no_correct (int value,int *power)
{

return (2 << ((*power = ((int) (log ((double) value)/log (2.) + .5))) - 1))
        == value? FALSE: TRUE;

/**************************************************************************/
/* Note of the programmer:                                                */
/*       Standard disclaimer: yes, I know it isn't a good programming     */
/*       practice to program this way. It may be a bit difficult to       */
/*       follow but, isn't it beautiful? I promise, anyway, not to use    */
/*       this kind of code any more (at least in this program). J.P.S.    */
/**************************************************************************/
}

/****************************************************************************/
/* FUNCTION:    RFFT, a real-valued, in place, split-radix FFT              */
/* PARAMETERS:  X = Hermitian symmetric or real input (length N+2)          */
/*              N = 2**M, dimensions of the transform (power of 2)          */
/*              KIND   = direction of the transform (0=direct,else inverse) */
/* RETURNS:     Real output in the same array X                             */
/* DESCRIPTION: Decimation in frequency, cos/sin in second loop.            */
/*              Length of input is a power of 2.                            */
/*              Input order:                                                */
/*                 [Re(0),Re(1), ..., Re(N-1)] (direct transform)           */
/*                 [Re(0),Im(0),Re(1),Im(1),...,Re(N-1),Im(N-1)] (inverse)  */
/*              Internal work order:                                        *
/*                 [Re(0)Re(1)...Re(N/2)Im(N/2-1)...Im(1)                   */
/****************************************************************************/
/*              Written by: Jos‚ Acu¤a CBM/CSIC/UAM                         */
/*              modified by Juan Pedro Secilla ===== IBM/MSC                */
/****************************************************************************/

void rfft (float *y, int n, int m,int kind)

{
int j=1,is,id,k,n1,i,i0,n2,n4,n5,n8,i1,i2,i3,i4,i5,i6,i7,i8;
float e,xt,a,r1,t1,t2,t3,t4,t5,t6,cc1,ss1,cc3,ss3,a3,fn;
register float *x;
static float x1[MAX_LEN+2];
static float *ccc1[11], *sss1[11], *ccc3[11], *sss3[11];
static int last_n = 0, last_m = 0;

/******* Check transform length *********************************************/
if (n > MAX_LEN)
{   puts ("Error: imagen demasiado grande (fft)");
    exit (1);
}

/******* Compute table of sines/cosines *************************************/
if (n != last_n)    /* Transform size has changed */
{   if (last_n != 0)     /* Correct , it does not work !!! */
        for (k = 3; k <= last_m; k++)
        {    free ((char *) sss1[k]);
             free ((char *) sss3[k]);
             free ((char *) ccc1[k]);
             free ((char *) ccc3[k]);
        }
    n2 = 4;
    for (k = 3; k <= m; k++) /* Get space for new ones */
    {   n2 <<= 1;
        n8 = n2 >> 3;
        e = a = TWOPI/n2;
        sss1[k] = (float *) malloc ((n8 + 1)*sizeof (float));
        sss3[k] = (float *) malloc ((n8 + 1)*sizeof (float));
        ccc1[k] = (float *) malloc ((n8 + 1)*sizeof (float));
        ccc3[k] = (float *) malloc ((n8 + 1)*sizeof (float));
        if (sss1[k] == NULL || sss3[k] == NULL || ccc1 [k] == NULL ||
            ccc3[k] == NULL)
        {   puts ("Error: no hay memoria (fft)");
            exit (1);
        }
        for (j = 2; j <= n8; j++)
        {   a3 = 3.*a;
            ccc1[k][j] = cos (a);
            sss1[k][j] = sin (a);
            ccc3[k][j] = cos (a3);
            sss3[k][j] = sin (a3);
            a = j*e;
        }
    }
    last_n = n;     /* Update last size */
    last_m = m;
}
/********************* Copy items to temp. array ****************************/
#if TURBO
    movmem ((char *) y, (char *) x1, (n+2)*sizeof(float));
#elif MSC
    memcpy ((char *) x1, (char *) y, (n+2)*sizeof(float));
#endif

if (kind == 0)  /*** Direct transform, operate with original data    ***/
    x = x1 -1;  /*** Dirty trick to simulate that indices begin at 1 ***/
else            /*** Alter to internal format for inverse transform  ***/
{   j = n/2;
    y[0]   = x1[0];
    y[j] =  x1[n];
    for (i = 1, k = 2; i < j; i++, k+=2)
    {   y[i] =   x1[k];
        y[n-i] = x1[k+1];
    }
    x = y-1;    /*** Operate with reversed data ***/
}

if (kind == 0)  /* Direct transform */
/*******************  Digit reverse  counter  *******************************/
{   j=1;
    n1=n-1;
    for(i=1;i<=n1;i++)
    {   if(i<j)
        {   xt=x[j];
            x[j]=x[i];
            x[i]=xt;
        }
        k=n>>1;  /* n/2; */
        while(k<j)
        {   j-=k;
            k >>=1; /* k/=2 */
        }
        j+=k;
    }

/******************  Length two butterflies  *********************************/

    is=1;
    id=4;
    do {
        for(i0=is;i0<=n;i0+=id)
        {   i1=i0+1;
            r1=x[i0];
            x[i0]=r1+x[i1];
            x[i1]=r1-x[i1];
        }
        is=(id<<1)-1;
        id <<=2;        /* id *= 4; */
    } while(is<n);

/*************** L-shaped butterflies  **************************************/

    n2=2;
    for(k=2; k<=m; k++)
    {   n2 <<= 1;   /* n2 *= 2; */
        n4 = n2 >> 2; /* n2/4; */
        n8=n2 >> 3;   /* n2/8; */
        is=0;
        id=n2<<1;
        do {
            for(i=is;i<n;i+=id)
            {   i1=i+1;
                i2=i1+n4;
                i3=i2+n4;
                i4=i3+n4;
                t1=x[i4]+x[i3];
                x[i4]-=x[i3];
                x[i3]=x[i1]-t1;
                x[i1]+=t1;
                if(n4!=1)
                {   i1+=n8;
                    i2+=n8;
                    i3+=n8;
                    i4+=n8;
                    t1=(x[i3]+x[i4])/SQRT2;
                    t2=(x[i3]-x[i4])/SQRT2;
                    x[i4]=x[i2]-t1;
                    x[i3]=-x[i2]-t1;
                    x[i2]=x[i1]-t2;
                    x[i1]+=t2;
                }
            }
            is=(id<<1)-n2;
            id <<=2;  /* id *= 4; */
        } while(is<n);
        for(j=2;j<=n8;j++)
        {   cc1=ccc1[k][j];
            ss1=sss1[k][j];
            cc3=ccc3[k][j];
            ss3=sss3[k][j];
            is=0;
            id=n2<<1;
            do {
                for(i=is;i<n;i+=id)
                {   i1=i+j;
                    i2=i1+n4;
                    i3=i2+n4;
                    i4=i3+n4;
                    i5=i+n4-j+2;
                    i6=i5+n4;
                    i7=i6+n4;
                    i8=i7+n4;
                    t1=x[i3]*cc1+x[i7]*ss1;
                    t2=x[i7]*cc1-x[i3]*ss1;
                    t3=x[i4]*cc3+x[i8]*ss3;
                    t4=x[i8]*cc3-x[i4]*ss3;
                    t5=t1+t3;
                    t6=t2+t4;
                    t3=t1-t3;
                    t4=t2-t4;
                    t2=x[i6]+t6;
                    x[i3]=t6-x[i6];
                    x[i8]=t2;
                    t2=x[i2]-t3;
                    x[i7]=-x[i2]-t3;
                    x[i4]=t2;
                    t1=x[i1]+t5;
                    x[i6]=x[i1]-t5;
                    x[i1]=t1;
                    t1=x[i5]+t4;
                    x[i5]-=t4;
                    x[i2]=t1;
                }
                is=(id<<1)-n2;
                id <<= 2;  /* id *=4;*/
            } while(is<n);
        }
    }
/******* Rerrange data as in program header, copy to y **********************/
    j = n/2;
    y[0] = x1 [0];
    y[1] = 0.;
    y[n] = x1 [j];
    y[n+1] = 0.;
    for (i = 1, k = 2; i < j; i++, k+=2)
    {   y[k] =   x1[i];
        y[k+1] = x1[n-i];
    }
}
else    /*** Inverse transform ***/
{   n2=n<<1;
    for(k=1;k<m;k++)
    {   is=0;
        id=n2;
        n2 >>=1;
        n4=n2>>2;
        n8=n4>>1;
        do {
            for(i=is;i<n;i+=id)
            {   i1=i+1;
                i2=i1+n4;
                i3=i2+n4;
                i4=i3+n4;
                t1=x[i1]-x[i3];
                x[i1]+=x[i3];
                x[i2]*=2.;
                x[i3]=t1-2.*x[i4];
                x[i4]=t1+2.*x[i4];
                if(n4!=1)
                {   i1+=n8;
                    i2+=n8;
                    i3+=n8;
                    i4+=n8;
                    t1=(x[i2]-x[i1])/SQRT2;
                    t2=(x[i4]+x[i3])/SQRT2;
                    x[i1]+=x[i2];
                    x[i2]=x[i4]-x[i3];
                    x[i3]=2.*(-t2-t1);
                    x[i4]=2.*(-t2+t1);
                }
            }
            is=(id<<1)-n2;
            id <<=2;
        } while(is<n-1);
        for(j=2;j<=n8;j++)
        {   cc1=ccc1[m-k+1][j];  /*** Note that cosine table is inverted]] ***/
            ss1=sss1[m-k+1][j];
            cc3=ccc3[m-k+1][j];
            ss3=sss3[m-k+1][j];
            is=0;
            id=n2<<1;
            do {
                for(i=is;i<n;i+=id)
                {   i1=i+j;
                    i2=i1+n4;
                    i3=i2+n4;
                    i4=i3+n4;
                    i5=i+n4-j+2;
                    i6=i5+n4;
                    i7=i6+n4;
                    i8=i7+n4;
                    t1=x[i1]-x[i6];
                    x[i1]+=x[i6];
                    t2=x[i5]-x[i2];
                    x[i5]+=x[i2];
                    t3=x[i8]+x[i3];
                    x[i6]=x[i8]-x[i3];
                    t4=x[i4]+x[i7];
                    x[i2]=x[i4]-x[i7];
                    t5=t1-t4;
                    t1+=t4;
                    t4=t2-t3;
                    t2+=t3;
                    x[i3]=t5*cc1+t4*ss1;
                    x[i7]=-t4*cc1+t5*ss1;
                    x[i4]=t1*cc3-t2*ss3;
                    x[i8]=t2*cc3+t1*ss3;
                }
                is=(id<<1)-n2;
                id<<=2;
            } while(is<n-1);
        }
    }

/*******************  Length two butterflies  *******************************/

    is=1;
    id=4;
    do {
        for(i0=is;i0<=n;i0+=id)
        {   i1=i0+1;
            r1=x[i0];
            x[i0]=r1+x[i1];
            x[i1]=r1-x[i1];
        }
        is=(id<<1)-1;
        id<<=2;
    } while(is<n);

/*******************  Digit reverse counter  ********************************/

    j=1;
    n1=n-1;
    for(i=1;i<=n1;i++)
    {   if(i<j)
        {   xt=x[j];
            x[j]=x[i];
            x[i]=xt;
        }
        k=n>>1;
        while(k<j)
        {   j-=k;
            k >>=1;
        }
        j+=k;
    }

/*********** Re-scale output ************************************************/
    fn = n;
    for(i=1;i<=n;i++)
        x[i]/=fn;
}
}

/********** Fourier transform of complex data (uses rfft twice) *************/

void cfft (float *x, float *y, int n,int m,int sentido)
      /*  x : Real part array and column position */
      /*  y : Imaginary part array and column position */

{   register int i0, i1;
    int n2, n4, i2, i3, id, is, i, j, k, size;
    float e, a, r1, r2, s1, s2, s3, xt, cc1, cc3, ss1, ss3;
    static float *ccc1[10], *sss1[10], *ccc3[10], *sss3[10];
    static int last_size = 0, last_m = 0;


           /************************************************************/
           /* As this program was originally written in Fortran, it's  */
    x--;   /* necessary to use this trick because in C the array ind-  */
    y--;   /* exes begin always at 0. Now the program will reference   */
           /* the elements correctly.                                  */
           /************************************************************/

    if (n != last_size)    /* Transform size has changed */
    {   if (last_size != 0)
            for (k = 1; k < last_m; k++)
            {   free ((char *) sss1[k]);
                free ((char *) sss3[k]);   /* Free old tables */
                free ((char *) ccc1[k]);
                free ((char *) ccc3[k]);
            }
        size = n>>2;
        for (k = 1; k < m; k++, size>>=1) /* Get space for new ones */
        {   sss1[k] = (float *) calloc (size + 1, sizeof (float));
            sss3[k] = (float *) calloc (size + 1, sizeof (float));
            ccc1[k] = (float *) calloc (size + 1, sizeof (float));
            ccc3[k] = (float *) calloc (size + 1, sizeof (float));
            if (sss1[k] == NULL || sss3[k] == NULL || ccc1 [k] == NULL ||
                ccc3[k] == NULL)
            {   puts ("Error: no hay memoria (fft)");
                exit (1);
            }
        }
        last_size = n;     /* Update last size */
        last_m = m;
        n2 = n<<1;
        for (k = 1; k < m; k++)
        {   n2 >>=1;
            n4 = n2>>2;
            e = 6.283185307179586/n2;
            a = 0.;
            for (j = 1; j <= n4; j++)
            {   ccc1[k][j] = cos (a);
                sss1[k][j] = sin (a);  /* Compute new tables */
                ccc3[k][j] = cos (3.*a);
                sss3[k][j] = sin (3.*a);
                a = j * e;
            }
        }
    }

    if (sentido != 0)
        for (i0 = 1; i0 <= n; i0 ++)
            y[i0] = -y[i0];             /* If inverse, conjugate input */

    n2 = n<<1;
    for (k = 1; k < m; k++)
    {   n2 >>=1;
        n4 = n2>>2;
        for (j = 1; j <= n4; j++)
        {   is = j;
            id = n2<<1;
            cc1 = ccc1[k][j];
            cc3 = ccc3[k][j];
            ss3 = sss3[k][j];
            ss1 = sss1[k][j];
            do {
                for (i0 = is; i0 < n; i0 += id)
                {   i1 = i0 + n4;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    r1 = x[i0] - x[i2];
                    x[i0] += x[i2];
                    r2 = x[i1] - x[i3];
                    x[i1] += x[i3];
                    s1 = y[i0] - y[i2];
                    y[i0] += y[i2];
                    s2 = y[i1] - y[i3];
                    y[i1] += y[i3];
                    s3 = r1 - s2;
                    r1 += s2;
                    s2 = r2 - s1;
                    r2 += s1;
                    x[i2] = r1*cc1 - s2*ss1;
                    y[i2] = -s2*cc1 - r1*ss1;
                    x[i3] = s3*cc3 + r2*ss3;
                    y[i3] = r2*cc3 - s3*ss3;
                }
                is = (id<<1) - n2 + j;
                id <<=2;
            } while (is < n);
        }
    }
/***************** last stage, length-2 butterfly ****************************/
    is = 1;
    id = 4;
    do {
        for (i0 = is; i0 <= n; i0 += id)
        {   i1 = i0 + 1;
            r1 = x[i0];
            x[i0] = r1 + x[i1];
            x[i1] = r1 - x[i1];
            r1 = y[i0];
            y[i0] = r1 + y[i1];
            y[i1] = r1 - y[i1];
         }
         is = (id<<1) - 1;
         id <<=2;
    } while (is < n);
/*************************** bit reverse counter ****************************/
    j = 1;
/*    n1 = n - 1;    */
    for (i = 1; i < n; i++)
    {   if (i < j)
        {   xt = x[j];
            x[j] = x[i];
            x[i] = xt;
            xt = y[j];
            y[j] = y[i];
            y[i] = xt;
        }
        k = n>>1;
        while (k < j)
        {   j -= k;
            k >>=1;
        }
        j += k;
    }

    if (sentido != 0)
        for (i0 = 1; i0 <= n; i0 ++)
            y[i0] = -y[i0];             /* If inverse, conjugate output */
}

