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
#include <cmath>
#include <values.h>

#include "groe.h"
#include "spider.h"

#ifdef __STDC__
int Resta_Plano(unsigned char ** , int , int, unsigned int);
#else
int Resta_Plano();
#endif


/*****************************************************************/
int Resta_Plano(unsigned char **Image , int rows, int cols, unsigned int io_Size)
{

    double x2, xy, x, y2, y, xz, yz, z, n, det, deta, detb, detc, i_rxi_r, i_r, j_c;
    double rctr, cctr;
    int i, j;
    double ch;
    double  f , ImagTemp, max, min ;
    int tipo = NATURAL;
    int temp;

    double a, b , c;


    x2 = 0.0;
    xy = 0.0;
    x = 0.0;
    y2 = 0.0;
    y = 0.0;
    xz = 0.0;
    yz = 0.0;
    z  = 0.0;

    rctr = (rows - 1) / 2.0;
    cctr = (cols - 1) / 2.0;
    if (rctr == 0 || cctr == 0)
    {
        fprintf(stderr, "Cannot find tilt on one pixel wide image!\n");
        return(0);
    }

    n = rows * cols;

    switch (tipo)
    {
    case NATURAL :
        for (i = 0; i < rows; i++)
        {
            i_r     = i - rctr;
            i_rxi_r = i_r * i_r;
            for (j = 0; j < cols; j++)
            {
                temp = (int)((j + i * cols + 1) / io_Size) ;
                ch = (double) Image[temp][j+i*cols+1 -temp*io_Size-1];
                j_c = j - cctr;
                x  += j_c;
                x2 += j_c * j_c;
                xy += j_c * i_r;
                y  += i_r;
                y2 += i_rxi_r;
                xz += j_c * ch;
                yz += i_r * ch;
                z  += ch;
            }
        }
        break;


    case FLOATFMT :
        for (i = 0; i < rows; i++)
        {
            i_r     = i - rctr;
            i_rxi_r = i_r * i_r;
            for (j = 0; j < cols; j++)
            {

                /* poner las lineas para el caso FLOAT  */

                j_c = j - cctr;
                x  += j_c;
                x2 += j_c * j_c;
                xy += j_c * i_r;
                y  += i_r;
                y2 += i_rxi_r;
                xz += j_c * f ;
                yz += i_r * f ;
                z  += f;
            }
        }
        break;


    default:
        fprintf(stderr, "Unkown data storage type\n");
        return(0);
        break;
    }
    /* Find the matrix determinant */
    det = x2 * (y2 * n - y * y) - xy * (xy * n - y * x) + x * (xy * y - y2 * x);

    /* Find the Cramer's determinants */
    deta = xz * (y2 * n - y * y) - xy * (yz * n - y * z) + x * (yz * y - y2 * z);
    detb = x2 * (yz * n - y * z) - xz * (xy * n - x * y) + x * (xy * z - yz * x);
    detc = x2 * (y2 * z - yz * y) - xy * (xy * z - yz * x) + xz * (xy * y - y2 * x);

    /* Compute the plane fit coefficients */
    a = deta / det;
    b = detb / det;
    c = detc / det;

    /* Resta a la imagen original el plano de coeficientes a, b, c.   */
    max = 0.;
    min = MAXDOUBLE;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            temp = (int)((j + i * cols + 1) / io_Size) ;
            ImagTemp = (double)(Image[temp][j+i*cols+1 -temp*io_Size-1]) - (double)(1 * (a * j + b * i + c));
            min = (min >= ImagTemp) * ImagTemp + (min < ImagTemp) * min;
            max = (max <= ImagTemp) * ImagTemp + (max > ImagTemp) * max;
        }
    }
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            temp = (int)((j + i * cols + 1) / io_Size) ;
            Image[temp][j+i*cols+1 -temp*io_Size-1]   = (unsigned char)(((double)Image[temp][j+i*cols+1-temp*io_Size-1]                 - (double)(1 * (a * j + b * i + c)) - min) / (max - min) * 255.);

        }
    }
    return(1);  /* success */

}
