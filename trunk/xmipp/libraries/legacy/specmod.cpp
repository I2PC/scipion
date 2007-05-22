/**************************************************************************
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

/**********************************************************************/
/* This program accepts a file containing a description of the rot-   */
/* ational spectrum of an image and generates a file with a clearer   */
/* description of the spectrum.    (J.P. Secilla - IBM/MSC)           */
/* See recons do a description of the spectrum file format.           */
/**********************************************************************/
/* Original in PL/I by A. Santisteban                                 */
/**********************************************************************/
/**********************************************************************/

#include <cstdio>
#include <cmath>
#include <cstring>

#ifndef __APPLE__
#include <malloc.h>
#else
#include <cstdlib>
#endif

#include "groe.h"

#define MAX_HARMONIC 51   /* Max. no of harmonics accepted */
#define eprintf(x) printf(x);fflush(stdout);
#define eprintf2(x,y) printf(x,y);fflush(stdout);
#define eprintf3(x,y,z) printf(x,y,z);fflush(stdout);
#define eprintf4(x,y,z,u) printf(x,y,z,u);fflush(stdout);


void spectro_m(char *nom, float xr1, float xr2, float xdr, float xr,
               float *result)

/** rad_min = xr1 **/
/** rad_max = xr2 **/
/** inc_int = xdr **/
/** long_int = xr **/

{
    float *e[MAX_HARMONIC], *rv, *st, *ep[MAX_HARMONIC], *erp [MAX_HARMONIC],
    *rp1, *rp2, *sp;
    FILE *in, *prt; /* input & output file */
    short ir, numin, numax;
    float xx0, yy0, rl, rh, dr, *c, *s;
    char nomin[128], nomprt[128];  /* names of files */
    int i1, n, m, i, j1, k, j, iflag, ir1, ir2, ndr, nr, ncol, nvez,
    irk, k1, k2;

    /******************** request some bytes of storage *******************/

    if ((c = (float *) calloc(5191, sizeof(float))) == NULL ||
        (s = (float *) calloc(5191, sizeof(float))) == NULL)
    {
        printf("\nERROR: no memory (spectro 1)");
        return;
    }

    /******************** Open input files *******************************/

    strcpy(nomin, nom);
    strcat(nomin, ".spc");
    if ((in = fopen(nomin, "rb")) == NULL)
    {
        printf("Error: Couldn't read the file %s\n", nomin);
        free((char *) c);
        free((char *) s);
        return;
    }

    /******************Read spectrum file into memory *********************/

    leespe(in, &ir, &numin, &numax, &xx0, &yy0, &rl, &rh, &dr, &c[1], &s[1]);
    fclose(in);

    /***************** Open output file **********************************/

    strcpy(nomprt, nom);
    strcat(nomprt, ".spt");
    if ((prt = fopen(nomprt, "w")) == NULL)
    {
        printf("Error: Couldn't create the file %s\n", nomprt);
        free((char *) c);
        free((char *) s);
        return;
    }

    fprintf(prt, "\n#Frecuency  minimum = %d, maximum = %d", numin, numax);
    fprintf(prt, "\n#Radius: low  = %f, high = %f, increment = %f", rl, rh, dr);

    /******************* Request more bytes ******************************/

    n = numax - numin + 1;
    m = (int)((rh - rl) / dr + 1);
    for (i = 1; i <= n; i++)
        if ((e[i] = (float *) calloc(m + 1, sizeof(float))) == NULL)
        {
            printf("\nERROR: no memory ");
            free((char *) c);
            free((char *) s);
            return;
        }
    if ((rv = (float *) calloc(m + 1, sizeof(float))) == NULL ||
        (st = (float *) calloc(m + 1, sizeof(float))) == NULL)
    {
        printf("\nERROR: no memory");
        for (i = 1; i <= n; i++)
            free((char *) ep[i]);
        free((char *) c);
        free((char *) s);
        fclose(prt);
        return;
    }

    /******************** The curre begins now (don't ask anything) ******/

    if (numin == 0)
        j1 = 2;
    else
        j1 = 1;
    k = 0;
    for (i = 1; i <= n; i++)
        for (j = 1; j <= m; j++)
        {
            k++;
            e[i][j] = c[k] * c[k] + s[k] * s[k];
        }

    for (i = 1; i <= m; i++)
    {
        rv[i] = rl + dr * (i - 1);
        st[i] = 0;
        for (j = j1; j <= n; j++)
            st[i] += e[j][i];
    }
    iflag = 0;

    ir1 = (int)((xr1 - rl) / dr + 1);
    if (ir1 < 1)
        ir1 = 1;
    ir2 = (int)((min(xr2, rh) - rl) / dr + 1);
    if (ir2 < ir1)
        ir2 = ir1;
    ndr = (int)(xdr / dr);
    if (ndr < 1)
        ndr = 1;
    if (xr < 0)
        nr = m;
    else
        nr = (int)(xr / dr + 1);
    ir2 = min(ir2, m);
    ncol = ir2 - nr + 1 - ir1;
    if (ncol < 0)
        ncol = 0;
    else
        ncol = 1 + ncol / ndr;
    if (ncol == 0)
    {
        printf("\n*** Incorrect data ***");
        for (i = 1; i <= n; i++)
            free((char *) ep[i]);
        free((char *) c);
        free((char *) s);
        free((char *) rv);
        free((char *) st);
        fclose(prt);
        return;
    }
    nvez = (ncol - 1) / 13 + 1;
    for (i = 1; i <= n; i++)
        if ((ep[i] = (float *) calloc(ncol + 1, sizeof(float))) == NULL ||
            (erp[i] = (float *) calloc(ncol + 1, sizeof(float))) == NULL)
        {
            printf("\nERROR: no memory ");
            for (i = 1; i <= n; i++)
                free((char *) ep[i]);
            free((char *) c);
            free((char *) s);
            free((char *) rv);
            free((char *) st);
            fclose(prt);
            return;
        }
    if ((rp1 = (float *) calloc(ncol + 1, sizeof(float))) == NULL ||
        (rp2 = (float *) calloc(ncol + 1, sizeof(float))) == NULL ||
        (sp  = (float *) calloc(ncol + 1, sizeof(float))) == NULL)
    {
        printf("\nERROR: no memory");
        for (i = 1; i <= n; i++)
        {
            free((char *) ep[i]);
            free((char *) erp [i]);
        }
        free((char *) c);
        free((char *) s);
        free((char *) rv);
        free((char *) st);
        fclose(prt);
        return;
    }
    for (k = 1; k <= ncol; k++)
    {
        irk = ir1 + (k - 1) * ndr - 1;
        rp1[k] = 10 * rv[irk+1];
        rp2[k] = 10 * rv[irk+nr];
        sp[k] = 0;
        for (k1 = 1; k1 <= nr; k1++)
            sp[k] += st[irk+k1];
        for (i = 1; i <= n; i++)
        {
            ep[i][k] = 0;
            for (k1 = 1; k1 <= nr; k1 ++)
                ep [i][k] += e[i][irk+k1];
            erp[i][k] = 1000000. * ep[i][k] / sp[k];
        }
    }
    for (k = 1; k <= nvez; k++)
    {
        k1 = 13 * (k - 1) + 1;
        k2 = min(ncol, 13 * k);
        if (iflag == 1)
        {
            fprintf(prt, "\n#Rotational Energy");
            fprintf(prt, "\n#   radius x 10\n#   ");
            for (j = k1; j <= k2; j++)
                fprintf(prt, "%4.0f-%4.0f ", rp1[j], rp2[j]);
            fprintf(prt, "\n#   ");
            for (j = k1; j <= k2; j++)
                fprintf(prt, "==========");
            fprintf(prt, "\n#nu");
            for (i = 1; i <= n; i++)
            {
                i1 = numin + i - 1;
                fprintf(prt, "\n%2d ", i1);
                for (j = k1; j <= k2; j++)
                    fprintf(prt, "%9.2f ", ep[i][j]);
            }
        }
        fprintf(prt, "\n\n#Fraction of the energy x 1.000.000");
        fprintf(prt, "\n#   radius x 10\n#   ");
        for (j = k1; j <= k2; j++)
            fprintf(prt, "%4.0f-%4.0f ", rp1[j], rp2[j]);
        fprintf(prt, "\n#   ");
        for (j = k1; j <= k2; j++)
            fprintf(prt, "==========");
        fprintf(prt, "\n#nu");
        for (i = j1; i <= n; i++)
        {
            i1 = numin + i - 1;
            fprintf(prt, "\n%2d ", i1);
            for (j = k1; j <= k2; j++)
            {
                fprintf(prt, "%9.0f ", erp[i][j]);
                if (j == k1)                    // para el .SIM
                    result[i-j1] = erp[i][j] / 10000;
                // para que salga bien con formato 5.2f
            }
        }
        fprintf(prt, "\n\n#   Whole energy of the non-null components\n#   ");
        for (j = k1; j <= k2; j++)
            fprintf(prt, "%9.2f ", sp [j]);
    }
    for (i = 1; i <= n; i++)
    {
        free((char *) ep[i]);
        free((char *) erp [i]);
    }
    free((char *) c);
    free((char *) s);
    free((char *) rv);
    free((char *) st);
    free((char *) rp1);
    free((char *) rp2);
    free((char *) sp);
    fclose(prt);
}

/**********************************************************************/

void leespe(FILE *in, short *ir, short *numin, short *numax,
            float *x0, float *y0, float *rl, float *rh, float *dr, float *ac, float *as)

{
    fread((char *)ir,    sizeof(*ir),    1,   in);
    fread((char *)numin, sizeof(*numin), 1,   in);
    fread((char *)numax, sizeof(*numax), 1,   in);
    fread((char *)x0,    sizeof(*x0),    1,   in);
    fread((char *)y0,    sizeof(*y0),    1,   in);
    fread((char *)rl,    sizeof(*rl),    1,   in);
    fread((char *)rh,    sizeof(*rh),    1,   in);
    fread((char *)dr,    sizeof(*dr),    1,   in);
    fread((char *)ac,   sizeof(*ac), (*numax - *numin + 1)*(*ir), in);
    fread((char *)as,   sizeof(*as), (*numax - *numin + 1)*(*ir), in);
}
/************************** The } ***********************************/
