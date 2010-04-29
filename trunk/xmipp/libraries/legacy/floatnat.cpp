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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/


/* **********************************************************************

      This file contains many useful routines to convert data from char
to float or from float to char.
When "form_entrada" is NATURAL then the input data is char and
the output data will be float; so, if "form_entrada" is FLOATFMT
the input data is float and it'll be convert to char.

 ************************************************************************/


#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "spider.h"
#include "groe.h"

/***************************************************************

This routine converts the a char matrix into a float matrix.

   INPUT :
           - a "fil x col" matrix of floats or chars
           - the input data type : "form_entrada" (f/n)

   OUTPUT:
           - a "fil x col" matrix of floats or chars, depending
             on "form_entrada".
 ***************************************************************/
void oraCharoraFloat(BYTE **img_nat, float **img_flo, int form_entrada, int fil, int col)
{
    BYTE byte_max, byte_min;
    long aux;
    double suma1 = 0, suma2 = 0;
    register int i, j;
    long size;
    char form;
    float float_max, float_min;
    double escala;
    double fmedia, fsigma;

    if (form_entrada == NATURAL)
        form = 'n';
    else
        form = 'f';

    switch (form)
    {
    case 'f':

        float_max = -1e36;
        float_min = 1e36;

        for (i = 0; i < fil; i++)
            for (j = 0; j < col; j++)
            {
                if (img_flo[i][j] > float_max) float_max = img_flo[i][j];
                if (img_flo[i][j] < float_min) float_min = img_flo[i][j];
            }

        if (float_max == float_min)
        {
            escala = 0.;
        }
        else
            escala = 255. / (float_max - float_min);
        for (i = 0; i < fil; i++)
            for (j = 0; j < col; j++)
                img_nat[i][j] = (unsigned char)(escala * (img_flo[i][j] - float_min));

        break;

    case 'n':

        size = fil * col;
        for (i = 0; i < fil; i++)    /***** Calculo de media y sigma *****/
            for (j = 0; j < col; j++)
            {
                aux = img_nat[i][j];
                suma1 += aux;
                suma2 += aux * aux;
            }

        fmedia = suma1 / (float) size;
        fsigma = sqrt(fabs(suma2 / (float)size - (fmedia * fmedia)));

        /***** Pasamos a flotante con media 0 y sigma 1 *****/

        for (i = 0; i < fil; i++)
        {
            for (j = 0; j < col; j++)
                img_flo[i][j] = (float)((img_nat[i][j] - fmedia) / fsigma);
        }
        break;

    }
}


/***************************************************************
  This routine performs the Bilinear Interpolation of an image.
  The input matrix is "iimage" and its dimensions are "nrow x ncol".
  The output matrix is loaded in "img" and its dimensions will be
  "nrow2 x ncol2"
 ***************************************************************/
void InterpBilinear(BYTE **iimage, BYTE **img, int nrow, int ncol, int nrow2,                           int ncol2)
{
    register int m, i, j, ii, jj;
    BYTE *aux;
    float escala1, escala2;
    float isx, isy;
    int ncol1, nrow1;
    int int1, int2;
    float startx, starty, incx, incy;
    float scalex, scaley;

    scalex = ((float) ncol2) / ncol;
    scaley = ((float) nrow2) / nrow;

    ncol1 = ncol - 1;
    nrow1 = nrow - 1;
    starty = 0;
    incx = 1. / scalex;
    incy = 1. / scaley;

    for (i = 0; i < nrow2; i++)
    {
        isy = starty;
        ii = (int)isy;
        aux = img [i];
        escala1 = isy - ii;
        startx = 0;
        for (j = 0; j < ncol2; j++)
        {
            isx = startx;
            jj = (int) isx;
            escala2 = isx - jj;
            if (ii < nrow1 && jj < ncol1)
            {
                int1 = iimage[ii][jj]  +
                       (int)(escala2 * ((int) iimage[ii][jj+1] - (int) iimage[ii][jj]));
                int2 = iimage[ii+1][jj] +
                       (int)(escala2 * ((int) iimage[ii+1][jj+1] - (int) iimage[ii+1][jj])) ;
                *aux ++ = (int1 + (int)(escala1 * (int2 - int1)));
            }
            else if (ii >= nrow1 && jj >= ncol1)
                *aux++ = iimage [ii][jj];
            else if (ii >= nrow1)
                *aux++  = iimage [ii][jj] +
                          (int)(escala2 * ((int) iimage[ii][jj+1] - (int) iimage [ii][jj]));
            else if (jj >= ncol1)
                *aux++  = iimage [ii][jj] +
                          (int)(escala1 * ((int) iimage[ii+1][jj] - (int) iimage [ii][jj]));
            startx += incx;
        }
        starty += incy;
    }
}


/***************************************************************
 This routine converts the data type when we have an array of
 chars (called "tira") and a matrix of floats (called "img_flo").
 Then, depending on the value of "form_entrada" the input and
 output will be as:

   if "form_entrada" is FLOATFMT  (floats)

     INPUT :  a "fil x col" matrix of floats.

     OUTPUT:  a "fil x col" long array of chars.

   if "form_entrada" is NATURAL ( chars )

     INPUT :  a "fil x col" matrix of floats.

     OUTPUT:  a "fil x col" long array of chars.

 If you're not using a table with the pixel value in a ColorMap,
 that is "xcolors" and "invxcolors", then you must call this
 routine with NULL values in the last two arguments.

 ***************************************************************/
void mFloat_tira(float **img_flo, unsigned char *tira, int fil, int col,                            int form_entrada, unsigned long *xcolors,
                 unsigned long *invxcolors)
{

    int i, j, k;
    unsigned char pixel;
    int aux;
    int char_max = 0, char_min = 255;
    double suma1 = 0, suma2 = 0;
    float float_max, float_min;
    float escala;
    float fsize;
    double fmedia, fsigma;

    switch (form_entrada)
    {
    case FLOATFMT:
        if (((int)cabecero.fImami) != 1.)
        {
            float_max = -1e36;
            float_min = 1e36;

            for (i = 0; i < fil; i++)
                for (j = 0; j < col; j++)
                {
                    if (img_flo[i][j] > float_max) float_max = img_flo[i][j];
                    if (img_flo[i][j] < float_min) float_min = img_flo[i][j];
                }
        }
        else
        {
            float_max = cabecero.fFmax;
            float_min = cabecero.fFmin;
        }
        if (float_max == float_min)
        {
            escala = 0.;
        }
        else
            escala = MAX_GREY / (float_max - float_min);
        for (i = 0, k = 0; i < fil; i++)
            for (j = 0; j < col; j++)
            {
                tira[k] = pixel = (unsigned char)(escala * (img_flo[i][j] - float_min)) + 5;
                if (xcolors != NULL)
                {
                    tira[k] = (unsigned char) xcolors[ tira[k] ];
                    if (invxcolors != NULL) invxcolors[ tira[k] ] = pixel;
                }
                k++;
            }

        break;

    case NATURAL:
        fsize = (float)(fil * col);
        for (i = 0, k = 0; i < fil; i++)    /***** Calculo de media y sigma *****/
            for (j = 0; j < col; j++)
            {
                aux = (int) tira[k];
                suma1 += aux;
                suma2 += aux * aux;
                if (aux > char_max) char_max = aux;
                if (aux < char_min) char_min = aux;
                k++;
            }

        fmedia = (double)(suma1 / fsize);
        fsigma = (double)(sqrt(fabs(suma2 / fsize - (fmedia * fmedia))));

        /***** Pasamos a flotante con media 0 y sigma 1 *****/

        for (i = 0, k = 0; i < fil; i++)
            for (j = 0; j < col; j++)
            {
                img_flo[i][j] = (float)(((double)tira[k] - fmedia) / fsigma);
                k++;
            }
        cabecero.fFmax = (float)(((double)char_max - fmedia) / fsigma);
        cabecero.fFmin = (float)(((double)char_min - fmedia) / fsigma);
        cabecero.fImami = 1.;
        cabecero.fAv = 0.;
        cabecero.fSig = 1.;
        break;

    }
}

/***************************************************************
 This fantastic routine loads an array of char into a matrix
 of chars. Easy!
 ***************************************************************/
void tira_matriz(unsigned char *tira, unsigned char **matriz, int nrow, int ncol)
{
    register int i, f, c;

    for (f = 0, i = 0; f < nrow; f++)
        for (c = 0; c < ncol; c++)
        {
            matriz[f][c] = tira[i];
            i++;
        }
}

/***************************************************************
 And this one loads a "nrow x ncol" matrix into an array of
 chars using a table of pixel values in a ColorMap.
 If you don't have such table ( you're lucky ) you can use
 this routine calling it with NULL value in the xcolors and invxcolors
 arguments.
 ***************************************************************/
void matriz_tira(unsigned char **matriz, unsigned char *tira,
                 int nrow, int ncol, unsigned long *xcolors,
                 unsigned long *invxcolors, int flag, int kte)
{
    register int i, f, c;
    unsigned char maxC, minC;
    unsigned char pixel;
    float escala, mean, dev, tmp;

    maxC = 0;
    minC = 255;
    /* Find max and  min  for change the scale from 0-255 to 5-255  */
    if (flag == 0)
    {
        for (f = 0, i = 0; f < nrow; f++)
            for (c = 0; c < ncol; c++)
            {
                if (matriz[f][c] < minC)     minC = matriz[f][c];
                else if (matriz[f][c] > maxC)  maxC = matriz[f][c];
            }

    }

    else
    {
        /* Calculates mean and standard deviation */
        mean = 0.0;
        dev = 0.0;
        tmp = 0.0;

        for (f = 0; f < nrow; f++)
            for (c = 0; c < ncol; c++)
            {
                mean += matriz[f][c];
                tmp += matriz[f][c] * matriz[f][c];
            }

        mean /= (nrow * ncol);
        dev = (float)sqrt((double)((tmp - nrow * ncol * mean * mean) / (nrow * ncol)));

        minC = (unsigned char)(mean - kte * dev);
        maxC = (unsigned char)(mean + kte * dev);

        for (f = 0; f < nrow; f++)
            for (c = 0; c < ncol; c++)
            {
                if (matriz[f][c] < minC) matriz[f][c] = minC;
                if (matriz[f][c] > maxC) matriz[f][c] = maxC;
            }

    }
    if (maxC == minC)
    {
        escala = 0.;
    }
    else
        escala = (float)(MAX_GREY / (float)(maxC - minC));

    for (f = 0, i = 0; f < nrow; f++)
        for (c = 0; c < ncol; c++)
        {
            tira[i] = pixel = ((unsigned char)((matriz[f][c] - minC) * escala))
                              + (unsigned char)5;
            /*
                      if (tira[i]<5) tira[i]=5;
                      if (tira[i]>255) tira[i]=255;
            */
            if (xcolors != NULL)
            {
                tira[i] = (unsigned char) xcolors[ tira[i] ];
                if (invxcolors != NULL) invxcolors[ tira[i] ] = pixel;
            }
            i++;
        }
}

/***************************************************************
 This routine convert the data when the input and output are
 volumes.

   INPUT :
           - a "slice x fil x col" volume of floats or chars
           - the input data type : "form_entrada"

   OUTPUT:
           - a "slice x fil x col" volume of floats or chars,
             depending on "form_entrada".
 ***************************************************************/
void volCharvolFloat(BYTE ***vol_nat, float ***vol_flo, int form_entrada,
                     int slice , int fil, int col)
{

    long aux;
    float float_max, float_min;
    float escala, f;
    double fmedia, fsigma;
    double suma1 = 0, suma2 = 0;
    double numpuntos;
    register int i, j, k;

    float_max = -1e36;
    float_min = 1e36;

    switch (form_entrada)
    {
    case NATURAL :

        for (i = 0; i < slice; i++)
            for (j = 0; j < fil ; j++)
                for (k = 0; k < col; k++)
                {
                    aux = vol_nat[i][j][k];
                    suma1 += aux;
                    suma2 += aux * aux;
                }
        numpuntos = (double)slice * (double) fil * (double) col;
        fmedia = (double)(suma1 / numpuntos);
        fsigma =  sqrt(fabs((suma2 / numpuntos) - (fmedia * fmedia)));

        /**** float volume with average=0 and sigma=1 ****/
        for (i = 0; i < slice; i++)
            for (j = 0; j < fil ; j++)
                for (k = 0; k < col; k++)
                {
                    f = vol_flo[i][j][k] = (float)((vol_nat[i][j][k] - fmedia) / fsigma);
                    if (f > float_max) float_max = f;
                    else if (f < float_min) float_min = f;
                }

        cabecero.fImami = 1.;
        cabecero.fSig = 1.;
        cabecero.fAv  = 0.;
        cabecero.fFmax = float_max;
        cabecero.fFmin = float_min;

        break;

    case FLOATFMT:
        if (((int)cabecero.fImami) != 1.)
        {
            for (i = 0; i < slice; i++)
                for (j = 0; j < fil ; j++)
                    for (k = 0; k < col; k++)
                    {
                        if (vol_flo[i][j][k] > float_max) float_max = vol_flo[i][j][k];
                        if (vol_flo[i][j][k] < float_min) float_min = vol_flo[i][j][k];
                    }
        }
        else
        {
            float_max = cabecero.fFmax;
            float_min = cabecero.fFmin;
        }
        if (float_max == float_min)
            escala = 0.;
        else
            escala = (float)(255. / (float_max - float_min));
        for (i = 0; i < slice; i++)
            for (j = 0; j < fil ; j++)
                for (k = 0; k < col; k++)
                    vol_nat[i][j][k] = (BYTE)(escala * (vol_flo[i][j][k] - float_min));

        break;
    }
}
/* *************************************************************************/
/* *************************************************************************/

