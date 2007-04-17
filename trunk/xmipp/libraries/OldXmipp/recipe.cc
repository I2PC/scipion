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

#include <math.h>

#define SIGN(a, b) ((b) < 0? -fabs(a): fabs (a))
#define MAX_ITER 200


#ifdef __STDC__
  void tred2 (float **, int, float *, float *);
  void tqli (float *, float *, int, float **);
#else
  void tred2 ();
  void tqli ();
#endif


/***************************************************************************/
/* Tridiagonalizacion de una matriz simetrica por el metodo de Householder */
/*        (see "Numerical Recipes in C")                                   */
/***************************************************************************/

void tred2 (float **a, int n, float *d, float *e)

{   int l,k,j,i;
    float scale, hh, h, g, f;


/**** Truco para que los ¡ndices empiecen en 1 ***/
d--;
e--;
for (i=0; i < n; i++)
    a[i]--;
a--;

for (i=n; i >= 2; i--)
{   l = i-1;
    h = scale = 0.;
    if (l > 1)
    {   for (k=1; k <= l; k++)
            scale += fabs (a[i][k]);
        if (scale == 0.)
            e[i] = a[i][l];
        else
        {   for (k=1; k <= l; k++)
            {   a[i][k] /= scale;
                h += a[i][k]*a[i][k];
            }
            f = a[i][l];
            g = f > 0? -sqrt (h): sqrt (h);
            e[i] = scale*g;
            h -= f*g;
            a[i][l] = f-g;
            f = 0.;
            for (j=1; j <= l; j++)
            {   a[j][i] = a[i][j]/h;
                g = 0.;
                for (k=1; k <= j; k++)
                    g += a[j][k]*a[i][k];
                for (k=j+1; k <= l; k++)
                    g += a[k][j]*a[i][k];
                e[j] = g/h;
                f += e[j]*a[i][j];
            }
            hh = f/(h+h);
            for (j=1; j <= l; j++)
            {   f = a[i][j];
                e[j] = g = e[j]-hh*f;
                for (k=1; k <= j; k++)
                    a[j][k] -= (f*e[k]+g*a[i][k]);
            }
        }
    }
    else
        e[i] = a[i][l];
    d[i] = h;
}
d[1] = 0.;
e[1] = 0.;
for (i=1; i <= n; i++)
{   l = i-1;
    if (d[i])
    {   for (j=1; j <= l; j++)
        {   g = 0.;
            for (k=1; k <= l; k++)
                g += a[i][k]*a[k][j];
            for (k=1; k <= l; k++)
                a[k][j] -= g*a[k][i];
        }
    }
    d[i] = a[i][i];
    a[i][i] = 1.;
    for (j=1; j <= l; j++)
        a[j][i] = a[i][j] = 0.;
}
a++;
for (i=0; i < n; i++)
    a[i]++;
}

/***************************************************************************/
/* Diagonalizacion de una matriz tridiagonal por el metodo QL con despla-  */
/* zamientos implicitos       ( see "Numerical Recipes in C")              */
/***************************************************************************/

void tqli (float *d, float *e, int n, float **z)

{   int m, l, iter, i, k;
    float s, r, p, g, f, dd, c, b;

/**** Truco para que los ¡ndices empiecen en 1 ***/
d--;
e--;
for (i=0; i < n; i++)
    z[i]--;
z--;

for (i=2; i <= n; i++)
   e[i-1] = e[i];
e[n] = 0.;
for (l=1; l <= n; l++)
{   iter = 0;
    do {
        for (m=l; m <= n-1; m++)
        {   dd = fabs (d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd)
                break;
        }
        if (m != l)
        {   if (iter++ == MAX_ITER)
                puts ("Error: demasiadas iteraciones (tqli)");
            g = (d[l+1]-d[l])/(2.*e[l]);
            r = sqrt ((g*g)+1.);
            g = d[m]-d[l]+e[l]/(g+SIGN(r, g));
            s = c = 1.;
            p = 0.;
            for (i=m-1; i >= l; i--)
            {   f = s*e[i];
                b = c*e[i];
                if (fabs (f) >= fabs (g))
                {   c = g/f;
                    r = sqrt ((c*c)+1.);
                    e[i+1] = f*r;
                    c *= (s = 1./r);
                }
                else
                {   s = f/g;
                    r = sqrt ((s*s)+1.);
                    e[i+1] = g*r;
                    s *= (c = 1./r);
                }
                g = d[i+1]-p;
                r = (d[i]-g)*s+2.*c*b;
                p = s*r;
                d[i+1] = g+p;
                g = c*r-b;
                for (k=1; k <= n; k++)
                {   f = z[k][i+1];
                    z[k][i+1] = s*z[k][i]+c*f;
                    z[k][i] = c*z[k][i]-s*f;
                }
            }
            d[l] = d[l]-p;
            e[l] = g;
            e[m] = 0.;
        }
    } while (m != l);
}
z++;
for (i=0; i < n; i++)
    z[i]++;
}
