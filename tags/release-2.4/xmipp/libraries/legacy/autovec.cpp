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

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "groe.h"
#include "sacaf.h"

#ifdef __STDC__
int gen_covar(float **, int, int , float **, float *);
void fill_covar(float **, int);
void med_0_var_1(float *, int);
void proyecImAutov(float *, float **, int, int, float *, float *);
void traspMatrix(float **, int);
void sortAutov(float *, float **, int);
#else
int gen_covar();
void fill_covar();
void med_0_var_1();
void proyecImAutov();
void traspMatrix();
void sortAutov();
#endif



/*****************************************************************************/
/* This kind function creates a covariance matrix of "num_ima" images.       */
/*****************************************************************************/
int gen_covar(float **list_ima, int num_ima, int num_pix,
              float **mat_cov, float *ima_med)

{
    int i, j, k;
    float *lin_cov;    /*** Line in the covariance matrix  ***/
    size_t size;


    if ((lin_cov = (float *) malloc(num_pix * sizeof(float))) == NULL)
    {
        fprintf(stdout, "\n > Sorry, no memory ");
        return(0);
        /*** Desagradecido ***/
    }

    /**** calculate the average image (Ver Gonz lez & Wintz, p g. 122) ****/
    for (j = 0; j < num_pix; j++)
        ima_med[j] = 0;

    for (i = 0; i < num_ima; i++)
        for (j = 0; j < num_pix; j++)
            ima_med[j] += list_ima[i][j];

    for (j = 0; j < num_pix; j++)
        ima_med[j] /= num_ima;

    /**** Calculate the covariance matrix ****/

    for (i = 0; i < num_pix; i++)
    {
        for (j = 0; j <= i; j++)         /*** Only a half ! ***/
            lin_cov[j] = 0;
        for (k = 0; k < num_ima; k++)
            for (j = 0; j <= i; j++)
                lin_cov[j] += list_ima[k][i] * list_ima[k][j];
        for (j = 0; j <= i; j++)
            lin_cov[j] = lin_cov[j] / num_ima - ima_med[i] * ima_med[j];

        size = (i + 1) * sizeof(float);
        memcpy(mat_cov[i], lin_cov, size);
    }
    free(lin_cov);
    return (1);
}


/*****************************************************************************/
/* This procedure copy the covariance matrix in the other half of it .       */
/*****************************************************************************/
void fill_covar(float **mat, int puntos)
{
    int i, j;

    for (i = 0; i < puntos; i++)
        for (j = i + 1; j < puntos; j++)
            mat[i][j] = mat[j][i];
}

/*****************************************************************************/
/***  Put an array of floats with average=0 and variance=1 ....            ***/
/*****************************************************************************/
void med_0_var_1(float *ristra, int longitud)
{
    int i;
    float media = 0., sigma = 0.;

    for (i = 0; i < longitud; i++)
    {
        media += ristra[i];
        sigma += ristra[i] * ristra[i];
    }

    media /= longitud;
    sigma = sqrt(fabs(sigma / longitud - media * media));

    for (i = 0; i < longitud; i++)
        ristra[i] = (ristra[i] - media) / sigma;
}

/***** This function projects an image in an autovector , and calculates  ****/
/***** the scalar product of both.                                        ****/
void proyecImAutov(float *imagen, float **autovec, int nvec, int n,
                   float *proyector, float *ima_med)
{
    int i, j;
    float proyec;

    for (i = 0; i < nvec; i++)
    {
        proyec = 0.;
        for (j = 0; j < n; j++)
            proyec += autovec[i][j] * (imagen[j] - ima_med[j]);
        ;
        proyector[i] = proyec;
    }
}

/**************************** Transpose a matrix *****************************/
void traspMatrix(float **a, int n)
{
    int i, j;
    float tempo;

    for (i = 0; i < n; i++)
        for (j = 0; j < i; j++)
        {
            tempo = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tempo;
        }
}

/*****     Sort a matrix of autovalues and a its associated matrix of     ****/
/*****     autovectors .                                                  ****/

void sortAutov(float *autoval, float **autovec, int n)

{
    int hay_cambios = TRUE;
    int i;
    float tempo;
    float *ptempo;

    while (hay_cambios)
    {
        hay_cambios = FALSE;
        for (i = 1; i < n; i++)
        {
            if (autoval[i-1] < autoval [i])
            {
                hay_cambios = TRUE;
                tempo = autoval[i-1];
                autoval[i-1] = autoval[i];
                autoval[i] = tempo;
                ptempo = autovec[i-1];
                autovec[i-1] = autovec[i];
                autovec[i] = ptempo;
            }
        }
    }
}
/* ***************************************************************************/

