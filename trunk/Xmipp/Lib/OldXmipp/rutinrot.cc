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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <groe.h>
#include <spider.h>


/* *************************************************************************/
/* *************************************************************************/

/***************************************************************
 This function calculates the angular corelation of two images,
 called "im_1" and "im_2". Both images must have the same
 dimension. The parameter "ang" is the rotational angle of the
 second image.
 ***************************************************************/
float cor_ang(double sum11,float **im_1,float **im_2,int dim, int ang, int r1, int r2)

{
    float sum=0.;
    float sum2=0.;
    float sum22;
    float dim2;
    int in1, in2, fin1, fin2;
    int i, j, k, l, m;
    float xprima1, yprima1, xprima2, yprima2;
    float valor1, valor2, scalex, scaley, y, x1, x2, scale1, val_imagen;
    float mat00, mat01, mat10, mat11;

dim2 = dim/2 - 0.5;             /*** The centre           ***/

mat00 =  cos (PI/180.*ang);     /*** Matriz de giro       ***/
mat01 =  sin (PI/180.*ang);
mat10 = -mat01;
mat11 =  mat00;

for (i=dim/2-r1; i < dim/2+r1; i++)   /*** Corelation only in one circle  ***/
{   bordes (dim2, r1, r2, i, &in1, &fin1, &in2, &fin2);
    y = i - dim2;
    x1 = in1-dim2;
    xprima1 = x1*mat00 + y*mat10 + dim2;
    yprima1 = x1*mat01 + y*mat11 + dim2;
    x2 = in2-dim2;
    xprima2 = x2*mat00 + y*mat10 + dim2;
    yprima2 = x2*mat01 + y*mat11 + dim2;
    /**** este es el punto inicial ****/
    for (j = in1; j < fin1; j++)
    {   k = yprima1;
        l = xprima1;
        scalex = xprima1 - l;
        scaley = yprima1 - k;
        scale1 = 1. - scalex;
        valor1 = scalex*im_2[k][l+1]   + scale1*im_2[k][l];
        valor2 = scalex*im_2[k+1][l+1] + scale1*im_2[k+1][l];
        val_imagen = scaley * valor2 + (1.-scaley) * valor1;
        sum2 += val_imagen*val_imagen;
        sum  += im_1[i][j]*val_imagen;
        xprima1 += mat00;
        yprima1 += mat01;
    }
    for (m = in2; m < fin2; m++)
    {   k = yprima2;
        l = xprima2;
        scalex = xprima2 - l;
        scaley = yprima2 - k;
        scale1 = 1. - scalex;
        valor1 = scalex*im_2[k][l+1]   + scale1*im_2[k][l];
        valor2 = scalex*im_2[k+1][l+1] + scale1*im_2[k+1][l];
        val_imagen = scaley * valor2 + (1.-scaley) * valor1;
        sum2 += val_imagen*val_imagen;
        sum  += im_1[i][m]*val_imagen;
        xprima2 += mat00;
        yprima2 += mat01;
    }
}
sum22 = sqrt(sum2);

if (sum11*sum22==0)
return (0);
return (sum/(sum11*sum22));

}

/***************************************************************
  This routine finds the local maximums in the result of
  the angular corelation. 
 * **************************************************************/
void max_cor (float *array, float *maxi, int *imax, int nmax, int largo)

{   int i,k,l,im;
    float maximo, aux;
    int continua;
    int ini, fin;

maximo = -1e38;
for (i=0; i < largo; i++)
    if (array[i] > maximo)
    {   maximo = array[i];
        im = i;
    }
maxi[0] = maximo;
imax[0] = im;

for (k=1; k < nmax; k++)
{   maximo = -1e38;
    for (i=0; i < largo; i++)
        if ((aux = array[i]) > maximo && aux < maxi[k-1])
            /******* No a los m ximos anteriores ********/
        {   continua = TRUE;
            for (l=-1; l <= 1; l++)
                if (l != 0 )
                   if (aux <= array[(i+l+largo)%largo])  /*** Se enrolla ***/
                   {   continua = FALSE;
                       break;
                   }
            if (continua)
            {   maximo = aux;
                im = i;
            }
        }
    imax[k] = im;
    maxi[k] = maximo;
}
}

void maximin (float **mat, int dim, int r1)

{
float maxi = -1e36, mini = 1e36, media = 0, sigma = 0, raiz;
int i, j,inicial,final;
float frad2,dim2;
long num_puntos = 0;

frad2 = r1*r1;
dim2 = dim/2 - 0.5;
media = num_puntos = 0;
for (i=dim/2-r1; i < dim/2+r1; i++)
{   raiz =  sqrt( fabs (frad2 - (i-dim2)*(i-dim2)) );
    inicial = dim2 - raiz + 1;
    final   = dim2 + raiz;
    for (j=inicial; j < final; j++, num_puntos++)
    {   media += mat[i][j];
        sigma += mat[i][j]*mat[i][j];
        if (mat[i][j] > maxi)
            maxi = mat[i][j];
        if (mat[i][j] < mini)
            mini = mat[i][j];
    }
}
media /= num_puntos;
sigma = sqrt(fabs (sigma/num_puntos-media*media) );

printf ("\n%f %f %f %f\n", maxi, mini, media, sigma);
}

/***************************************************************/
/*  Low pass filter of an image.                               */
/* *************************************************************/
void filtra_pb (float **im_ent, float **im_sal, int fil, int col)

{   int i, j;

/***** Copia la orla a la salida *****/

for (i=0; i < fil; i++)
{    im_sal[i][0] = im_ent[i][0];
     im_sal[i][col-1] = im_ent[i][col-1];
}
for (j=0; j < col; j++)
{    im_sal[0][j] = im_ent[0][j];
     im_sal[fil-1][j] = im_ent[fil-1][j];
}
/***** Promedia salvajemente *********/

for (i=1; i < fil-1; i++)
    for (j=1; j < col-1; j++)
        im_sal[i][j] = (im_ent[i-1][j-1] + im_ent[i-1][j] + im_ent[i-1][j+1] +
                        im_ent[i][j-1]   + im_ent[i][j]   + im_ent[i][j+1]   +
                        im_ent[i+1][j-1] + im_ent[i+1][j] + im_ent[i+1][j+1])/9.;

/****** Copiamos al original *********/

for (i=0; i < fil; i++)
    for (j=0; j < col; j++)
        im_ent[i][j] = im_sal[i][j];
}

/***************************************************************
  This routine calculates the circumference's borders.
 ***************************************************************/
void bordes (float dim2, int r1, int r2, int punto, int *x1, int *x2, int *x3, int *x4)

{
float fr1, fr2;
float raiz1, raiz2;

fr1 = r1* (float) r1;
fr2 = r2* (float) r2;


raiz1 = sqrt(fabs( (fr1 - (punto-dim2)*(punto-dim2)) ) );
        /** external circumference **/
*x1 = dim2 - raiz1 +1;
*x4 = dim2 + raiz1;

if ((fr2 - (punto-dim2)*(punto-dim2)) < 0. )
   *x2 = *x3 = dim2;
else
{   raiz2 = sqrt(fabs( (fr2 - (punto-dim2)*(punto-dim2)) ));
    *x2 = dim2 - raiz2 +1;
    *x3 = dim2 + raiz2;
}

}

/* ************************************************************************/
/* *************************************************************************/
/* *************************************************************************/
/* *************************************************************************/
