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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <malloc.h>

#include "groe.h"

#ifdef __STDC__
void centrado_t(float **, float **, int, float **, float **);
#else
void centrado_t();
#endif

/*****************************************************************************/
/*****************************************************************************/

void centrado_t(float **fuIm2, float **fimTil, int dim,
                float **fuIm1, float **lectura)

{
    int i, j, imax[4], jmax[4], k, l;
    int dim2;
    float  maxi[4], real1, real2, imag1, imag2;
    float sum_patr;
    float xx_max, yy_max, sum_corr, x_nue, y_nue;
    int n_max;                   /**** Para el centro de gravedad ****/
    int para_gravedad;           /**** Este es un buleano         ****/
    int i_actual, j_actual;
    float scale_x, scale_y, valor1, valor2;
    double temporal[3][3];


    dim2 = 2 * dim;

    for (i = 0; i < dim2; i ++)
        for (j = 0; j < dim2; j ++)
            fuIm1[i][j] = lectura[i][j] = 0.;

    sum_patr = dim * dim;  /** standar desviation =1  **/


    image_fft(fuIm2, dim2, dim2, DIRECT);

    /**** Centre "fimTil" ****/

    for (i = 0; i < dim; i ++)
        for (j = 0; j < dim; j ++)
            fuIm1[i+dim/2][j+dim/2] = fimTil[i][j];

    /***** Save "fuIm1" in "lectura"  ****/
    for (i = 0; i < dim2; i ++)
        for (j = 0; j < dim2; j ++)
            lectura[i][j] = fuIm1[i][j];


    image_fft(fuIm1, dim2, dim2, DIRECT);

    /******* Cross corelation ****************/

    for (i = 0; i < dim2; i++)
        for (j = 0; j <= dim; j++)
        {
            real1 = fuIm1[i][2*j];
            imag1 = fuIm1[i][2*j+1];
            real2 = fuIm2[i][2*j];
            imag2 = fuIm2[i][2*j+1];
            fuIm1[i][2*j]   = real1 * real2 + imag1 * imag2;
            fuIm1[i][2*j+1] = real2 * imag1 - real1 * imag2;
        }

    /******* Inverse F.F.T.   ********/

    image_fft(fuIm1, dim2, dim2, INVERSE);

    /******* Normalize the result *******/

    cornor(lectura, dim, dim, fuIm1, dim2, dim2, sum_patr);

    /******* Looking for the maximum ( 4 times ) ********/

    busca_maximos(fuIm1, dim + 1, dim + 1, maxi, imax, jmax, 4);

    /******* Calculate the gravity centre of the corelation        ****/
    /******* in a neighbourhood such as  maximum/sqrt(2) > value   ****/

    n_max = -1;
    para_gravedad = FALSE;
    while (!para_gravedad)
    {
        n_max ++;
        for (i = (-n_max); i <= n_max; i++)
            for (j = (-n_max); j <= n_max; j++)
            {
                i_actual = i + imax[0];
                j_actual = j + jmax[0];
                if (i_actual >= 0 && j_actual >= 0 &&
                    i_actual < dim + 1 && j_actual < dim + 1)
                    if (maxi[0] / 1.414 > fuIm1[i_actual][j_actual]) /*** Menor ***/
                        para_gravedad = TRUE;
            }
    }

    /*** We have the neighbourhood => looking for the gravity centre ***/

    xx_max = yy_max = sum_corr = 0;
    for (i = (-n_max); i <= n_max; i++)
        for (j = (-n_max); j <= n_max; j++)
        {
            i_actual = i + imax[0];
            j_actual = j + jmax[0];
            if (i_actual >= 0 && j_actual >= 0 &&
                i_actual < dim + 1 && j_actual < dim + 1)
            {
                yy_max += i_actual * fuIm1[i_actual][j_actual];
                xx_max += j_actual * fuIm1[i_actual][j_actual];
                sum_corr += fuIm1[i_actual][j_actual];
            }
        }
    xx_max /= sum_corr; /***This is the gravity centre ***/
    yy_max /= sum_corr;

    /******* Displace the image  with xx_max, yy_max ******/
    xx_max -= dim / 2;          /*** Centre the origin  ****/
    yy_max -= dim / 2;

    temporal[0][0] = 1;
    temporal[0][1] = 0;
    temporal[0][2] = 0;
    temporal[1][0] = 0;
    temporal[1][1] = 1;
    temporal[1][2] = 0;
    temporal[2][0] = xx_max;
    temporal[2][1] = yy_max;
    temporal[2][2] = 1;


    transforma(fimTil, lectura, dim, dim, temporal);

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            fimTil[i][j] = lectura[i][j];

}
/* **************************************************************************/
/* **************************************************************************/

