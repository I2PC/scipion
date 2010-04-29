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
#include <cstring>
#include <cstdlib>

#ifndef __APPLE__
#include <malloc.h>
#endif

#include "groe.h"

/***************************************************************
  This routine does a Low Pass Filtering of an input image
  "im_ent". The resulting image is loaded in "im_sal".
 ***************************************************************/
void filtra_pb(float **im_ent, float **im_sal, int fil, int col)

{
    int i, j;

    /***** Copia la orla a la salida *****/

    for (i = 0; i < fil; i++)
    {
        im_sal[i][0] = im_ent[i][0];
        im_sal[i][col-1] = im_ent[i][col-1];
    }
    for (j = 0; j < col; j++)
    {
        im_sal[0][j] = im_ent[0][j];
        im_sal[fil-1][j] = im_ent[fil-1][j];
    }
    /***** Promedia salvajemente *********/

    for (i = 1; i < fil - 1; i++)
        for (j = 1; j < col - 1; j++)
            im_sal[i][j] = (im_ent[i-1][j-1] + im_ent[i-1][j] + im_ent[i-1][j+1] +
                            im_ent[i][j-1]   + im_ent[i][j]   + im_ent[i][j+1]   +
                            im_ent[i+1][j-1] + im_ent[i+1][j] + im_ent[i+1][j+1]) / 9.;

    /****** Copiamos al original *********/

    for (i = 0; i < fil; i++)
        for (j = 0; j < col; j++)
            im_ent[i][j] = im_sal[i][j];
}

/***************************************************************/
/***************************************************************/

#define Re(i,j) x[i][2*(j)]
#define Im(i,j) x[i][2*(j)+1]

void cos_alz(float **x, int dim, int r1, int r2, int r3, int r4)
{
    int i, j, k;
    long r12, r22, r32, r42;
    long radio;
    float aux, aux1, aux2;   /*** pi/(r2-r1), pi/(r4-r3) ***/
    int dim2;
    int ii;
    long i2;

    r12 = r1 * (long)r1;
    r22 = r2 * (long)r2;
    r32 = r3 * (long)r3;
    r42 = r4 * (long)r4;

    if (r1 != r2)
        aux1 = PI / (r2 - r1);
    else
        aux1 = 0;
    if (r3 != r4)
        aux2 = PI / (r4 - r3);
    else
        aux2 = 0;

    dim2 = dim / 2;

    if (r1 < 0 || r2 < 0 || r1 > r2 || r1 > dim || r2 > dim)
    {
        puts("Radios de la ventana err¢neos (0 <= r1 <= r2 < dimensi¢n) (mascara)");
        exit(1);
    }

    for (i = 0; i < dim; i++)
    {
        ii = i > dim2 ? i - dim : i;
        i2 = ii * (long)ii;
        for (j = 0; j <= dim2; j++)
        {
            radio = i2 + j * (long)j;
            if (radio < r12)
                Re(i, j) = Im(i, j) = 0;
            else if (radio < r22)   /*** Aplicamos la ventana ***/
            {
                aux = (1. + cos((r2 - sqrt((double)radio)) * aux1)) * 0.5;
                Re(i, j) *= aux;
                Im(i, j) *= aux;
            }
            else if (radio < r32)
                ;           /**** No hacemos nada, claro *****/
            else if (radio < r42)   /*** Aplicamos la ventana ***/
            {
                aux = (1. + cos((sqrt((double)radio) - r3) * aux2)) * 0.5;
                Re(i, j) *= aux;
                Im(i, j) *= aux;
            }
            else
                Re(i, j) = Im(i, j) = 0;
        }
    }
}
