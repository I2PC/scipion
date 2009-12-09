/***************************************************************************
 *
 * Author:     Monica Chagoyen
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

      This file contains all the routines to call the corresponding
      fft routines, depending on the dimensions of the data (power of 2
      or not).

 ************************************************************************/

#include <cstdio>
#include <cstdlib>

#include "fftn.h"
#include "groe.h"

/*
# define fftn fftnf
*/
/* define some convenient macros for printing out values */
float *datamon = NULL;
/* element of a complex array stored [real0,imag0,real1,imag1, ...] */
# define Real_El(i,j) datamon [2*((j)*row + (i))]
# define Imag_El(i,j) datamon [2*((j)*row + (i)) + 1]

# define Real(i,j)       x[i][2*(j)]
# define Imag(i,j)       x[i][2*(j)+1]

/* define macros for volume */
# define Real_Vol(i,j,k)  datamon [2*((i)*row*col+(j)*col+(k))]
# define Imag_Vol(i,j,k)  datamon [2*((i)*row*col+(j)*col+(k))+1]

# define Realv(i,j,k)     vol[2*(i)][j][k]
# define Imagv(i,j,k)     vol[2*(i)+1][j][k]

/**************************************************************************
        This routine performs the fft (direct or inverse) of an image
**************************************************************************/

int image_fft(float **x, int row, int col, int kind)
{
    int ret, power;
    int dims [2];  /* pass fft dimensions */
    int m, n, col2;

    if (no_correct(row, &power) || no_correct(col, &power))
    {

        datamon = (float *) malloc(2 * row * col * sizeof(float));

        if (datamon == NULL)
        {
            fprintf(stderr, "Unable to allocate memory.\n");
            return 1;
        }
        dims [0] = row;
        dims [1] = col;

        if (kind == DIRECT)
        {
            /* 2D forward fft */
            for (m = 0; m < row; m++)
                for (n = 0; n < col; n++)
                {
                    Real_El(m, n) = x[m][n];
                    Imag_El(m, n) = .0;
                }

            ret = FFTN(2, dims,
                       &Real_El(0, 0), &Imag_El(0, 0), -2, 0.0);
            if (ret)
            {
                fprintf(stderr, "\n FFT coefficients error");
                return(ERROR);
            }

            col2 = (int)col / 2;
            for (m = 0; m < row; m++)
                for (n = 0; n <= col2; n++)
                {
                    Real(m, n) = Real_El(m, n);
                    Imag(m, n) = Imag_El(m, n);
                }
            fft_free();
            free(datamon);
        }

        else
        {
            /* -------------now do the inverse----------------- */

            col2 = (int)col / 2;
            for (m = 0; m < row; m++)
                for (n = 0; n < col; n++)
                {
                    if (n <= col2)
                    {
                        Real_El(m, n) = Real(m, n);
                        Imag_El(m, n) = Imag(m, n);
                    }
                    else if (m == 0)
                    {
                        Real_El(m, n) = Real(m, col - n);
                        Imag_El(m, n) = -Imag(m, col - n);
                    }
                    else
                    {
                        Real_El(m, n) = Real(row - m, col - n);
                        Imag_El(m, n) = -Imag(row - m, col - n);
                    }

                }


            ret = FFTN(2, dims,
                       &Real_El(0, 0), &Imag_El(0, 0),
                       2,
                       -1.0);
            if (ret) return ERROR;

            col2 = col / 2;
            for (m = 0; m < row; m++)
                for (n = 0; n < col; n++)
                {
                    if (n < col2)
                        x[m][n] = Real_El(m, n);
                    else
                        x[m][n] = Real_El(m, n) / (row * col);
                }

            fft_free();
            free(datamon);
        }
        return(OK);

    }
    else
    {
        ret = fft_pow2(x, row, col, kind);
        return (ret);
    }
}

/************************************************************************
        This routine perform the fft (direct or inverse of a volume
************************************************************************/

int volume_fft(float ***vol, int sli, int row, int col, int kind)
{
    /* vol: Fourier transform array       */
    /* sli, row, col: volume dimensions   */
    /* kind: of transform to be performed */
    int poweri, powerj, powerk;                  /* 2**power = nrow or ncol */
    int l, m, n, ret;
    int dims[3];

    /******************* Check image dimensions *********************************/

    if (no_correct(sli, &poweri) ||
        no_correct(row, &powerj) || no_correct(col, &powerk))
    {

        datamon = (float *)malloc(2 * row * col * sli * sizeof(float));
        if (datamon == NULL)
        {
            fprintf(stderr, "\nCan't allocate memory\n");
            return (ERROR);
        }

        dims[0] = sli;
        dims[1] = row;
        dims[2] = col;

        if (kind == DIRECT)
        {
            for (l = 0; l < sli; l++)
                for (m = 0; m < row; m++)
                    for (n = 0; n < col; n++)
                    {
                        Real_Vol(l, m, n) = vol[l][m][n];
                        Imag_Vol(l, m, n) = .0;
                    }
            /* Computes fft with fftn */
            ret = FFTN(3, dims, &Real_Vol(0, 0, 0), &Imag_Vol(0, 0, 0), -2, 0.0);
            if (ret)
            {
                printf("Error computing fft coefficients\n");
                return (ERROR);
            }
            for (l = 0; l <= sli / 2; l++)
                for (m = 0; m < row; m++)
                    for (n = 0; n < col; n++)
                    {
                        Realv(l, m, n) = Real_Vol(l, m, n);
                        Imagv(l, m, n) = Imag_Vol(l, m, n);
                    }
            fft_free();
            free(datamon);
            return (OK);
        }
        else
        {
            for (l = 0; l < sli; l++)
                for (m = 0; m < row; m++)
                    for (n = 0; n < col; n++)
                    {
                        if (l <= sli / 2)
                        {
                            Real_Vol(l, m, n) = Realv(l, m, n);
                            Imag_Vol(l, m, n) = Imagv(l, m, n);
                        }
                        else if (m == 0 && n == 0)
                        {
                            Real_Vol(l, m, n) = Realv(sli - l, m, n);
                            Imag_Vol(l, m, n) = -Imagv(sli - l, m, n);
                        }
                        else if (m == 0)
                        {
                            Real_Vol(l, m, n) = Realv(sli - l, m, col - n);
                            Imag_Vol(l, m, n) = -Imagv(sli - l, m, col - n);
                        }
                        else if (n == 0)
                        {
                            Real_Vol(l, m, n) = Realv(sli - l, row - m, n);
                            Imag_Vol(l, m, n) = -Imagv(sli - l, row - m, n);
                        }
                        else
                        {
                            Real_Vol(l, m, n) = Realv(sli - l, row - m, col - n);
                            Imag_Vol(l, m, n) = -Imagv(sli - l, row - m, col - n);
                        }
                    }
            ret = FFTN(3, dims, &Real_Vol(0, 0, 0), &Imag_Vol(0, 0, 0), 2, -1.0);
            if (ret)
            {
                printf("Error computing fft coefficients\n");
                return (ERROR);
            }
            for (l = 0; l < sli; l++)
                for (m = 0; m < row; m++)
                    for (n = 0; n < col; n++)
                    {
                        if (l < sli / 2)
                            vol[l][m][n] = Real_Vol(l, m, n);
                        else
                            vol[l][m][n] = Real_Vol(l, m, n) / (sli * row * col);
                    }
            fft_free();
            free(datamon);
            return(OK);
        }
    }
    else
    {
        ret = vfft_pow2(vol, sli, row, col, kind);
        return(ret);
    }
}
/* ---------------------- end-of-file (c source) ---------------------- */
