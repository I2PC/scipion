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

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spider.h"
#include "groe.h"


/***************************************************************
 This routine normalizes the corelation of two images in such
 way that the result values between -1 and 1. The bigger image
 is "imag1" and its dimensions are "nf1 x nc1".
 The corelations is "corr", and "sum_masc" is the addition of
 the pixels's squares of the second image . The dimensions of
 the second image are "nf2 x nc2".
 ****************************************************************/
void cornor(float **imag1, int nf1, int nc1, float **corr, int nf2, int nc2, double sum_masc)

{
int i,j,k,l;                     /* contadores */
double  suma;                    /* variable auxiliar */
int nc3,nf3;
float *sum1, *sum2, aux;

nc3 = nc2 - nc1 +1;
nf3 = nf2 - nf1 +1;

/******** Miramos condiciones de error ********/
if (nf2 < nf1 || nc2 < nc1)
{   puts ("Dimensiones incorrectas para la normalizaci¢n (cornor)");
    exit (1);
}

/******** Pedimos un poco de memoria *********/

if((sum1 = (float *) malloc(nc3*sizeof(float)))==NULL)
{   puts ("Error: no hay memoria (cornor)");
    exit (1);
}
if((sum2 = (float *) malloc(nc3*sizeof(float)))==NULL)
{   puts ("Error: no hay memoria (cornor)");
    exit (1);
}

/******** Comenzamos a normalizar la correlaci¢n ********/

for(j = 0; j < nc3; j++)
{   suma = 0;
    for (k = 0; k < nf1; k++)
        for (l = j; l < j+nc1; l++)
            suma += imag1[k][l]*imag1[k][l];
    sum1 [j] = suma;
    corr[0][j] /= (float) sqrt (suma*sum_masc);
}

for(i = 1; i < nf3; i++)
{   for(l = 0; l < nc3; l++)
        sum2 [l] = sum1 [l];
    for(l = 0; l < nc1; l++)
    {   suma = imag1[i+nf1-1][l];
        aux  = imag1[i-1][l];
        sum1[0] += suma*suma-aux*aux;
    }
    corr[i][0] /= (float) sqrt (sum1[0]*sum_masc);
    for(j = 1; j < nc3;j++)
    {   suma = sum2[j] - sum2[j-1] + sum1[j-1];
        aux = imag1[i-1][j-1];
        suma += aux*aux;
        aux = imag1[i-1][j+nc1-1];
        suma -= aux*aux;
        aux = imag1[i+nf1-1][j-1];
        suma -= aux*aux;
        aux = imag1[i+nf1-1][j+nc1-1];
        suma += aux*aux;
        sum1[j] = suma;
        corr[i][j] /= (float) sqrt (suma*sum_masc);
    }
}
free(sum1);
free(sum2);
}

/***************************************************************
 This routine transforms an image using 3x3 matrix. The input
 matrix is "im_in" and the transformated image will be "im_out".
 All the interpolations are bilinears.
 ***************************************************************/
void transforma (float **im_in, float **im_out, int fil, int col, double matriz[3][3])

/**** im_in, im_out: input and output images            ***/
/**** fil, col: their dimensions                        ***/
/**** matriz: matrix for transform the input image      ***/
/**** [x' y' 1] = | a d 0 | [x y 1]  ecuation           ***/
/****             | b e 0 |                             ***/
/****             | c f 1 |                             ***/

{    int i, j, k, l, fil1, col1;
     float xprima, yprima, cen_fil, cen_col;
     float valor1, valor2, scalex, scaley, y, x, scale1;
     float mat00, mat01;
/* Si es la identidad no hagas nada */
if(matriz[0][0]==1.0 &&
   matriz[1][0]==0.0 &&
   matriz[2][0]==0.0 &&
   matriz[0][1]==0.0 &&
   matriz[1][1]==1.0 &&
   matriz[2][1]==0.0 
   ) 
{
for (i = 0; i < fil; i++)    
    for (j = 0; j < col; j++)
        im_out[i][j] = im_in[i][j]; 
return;

}

cen_fil = fil/2 - 0.5;  /**** este es el centro ****/
cen_col = col/2 - 0.5;

fil1 = fil-1;
col1 = col-1;
mat00 = matriz[0][0];
mat01 = matriz[0][1];

for (i = 0; i < fil; i++)       /**** Imagen de salida a cero ***/
    for (j = 0; j < col; j++)
        im_out[i][j] = 0.;

for (i = 0; i < fil1; i++)  /*** Vamos al rev‚s (del final al origen) ***/
{   y = i - cen_fil;
    x = -cen_col;
    xprima = x*matriz[0][0] + y*matriz[1][0] + matriz[2][0] + cen_col;
    yprima = x*matriz[0][1] + y*matriz[1][1] + matriz[2][1] + cen_fil;
    /**** este es el punto inicial ****/
    for (j = 0; j < col1; j++)
    {   /**** El nuevo punto es xprima, yprima ****/
        /***** Ahora, a interpolar           ****/
        k = (int)yprima;
        l = (int)xprima;
        if (k < 0 || l < 0 || k >= fil1 || l >= col1)
     /*        im_out[i][j] = 0;*/
             im_out[i][j] = im_in[0][0];
        else
        {   scalex = xprima - l;
            scaley = yprima - k;
            scale1 = 1. - scalex;
            valor1 = scalex*im_in[k][l+1]   + scale1*im_in[k][l];
            valor2 = scalex*im_in[k+1][l+1] + scale1*im_in[k+1][l];
            im_out[i][j] = scaley * valor2 + (1.-scaley) * valor1;
        }
        xprima += mat00; /**** Actualiza el punto *****/
        yprima += mat01;
    }
}

for(i=0; i<fil1; i++)
   im_out[i][col1]=im_out[i][col1-1];

for(j=0; j<col; j++)
   im_out[fil1][j]=im_out[fil1-1][j];

/*** El truco del fmod es para que k, l siempre queden entre 0 y N-1 ***/
/*** es decir, la imagen tiene topolog¡a toroidal. De esta forma, lo ***/
/*** que sale por un lado entra por el otro.                         ***/
}

/*****************************************************************************/
void compon (double m1[3][3], double m2[3][3], double m3[3][3])

{   int i,j,k;
    double aux;

for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
    {   aux = 0.;
        for (k = 0; k < 3; k++)
            aux += m1[i][k] * m2[k][j];
        m3[i][j] = aux;
     }
}
