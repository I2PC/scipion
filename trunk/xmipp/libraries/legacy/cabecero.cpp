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

This file contains all the routines to read, write and modify
the image or volume header called 'cabecero'.

 ************************************************************************/

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <sys/types.h>
#include <sys/stat.h>

#include "spider.h"
#include "groe.h"

CABECERO cabecero;

/***************************************************************
 This function open the file "fname" and check if it is a
 SPIDER file. It returns '1' and the image dimensions if success,
 '0' if failure, and '-1' if there's  a file open error .
 ***************************************************************/
int SPIDheader(char *fname, int *row, int *col)
{
    unsigned long tam;
    FILE *fp;
    struct stat info;

    fp = fopen(fname, "rb");
    if (!fp) return(-1);   /*unable to open file*/

    if (fread(&cabecero, SPID_HDR_SIZE , 1, fp) != 1)
    {
        fclose(fp);
        return(0) ;       /*not a SPIDER file*/
    }

    if (!fstat(fileno(fp), &info))
    {
        fclose(fp);
        tam = (unsigned long)(cabecero.fNsam * cabecero.fLabrec * 4) +
              (unsigned long)(cabecero.fNsam * cabecero.fNrow * 4);
        if (tam == info.st_size)
        {
            *row = (int)cabecero.fNrow;
            *col = (int)cabecero.fNsam;
            return 1;
        }
    }
    else
        fclose(fp);
    return 0;
}

/***************************************************************
 This function open the file "fname" and check if it is a BYTE
 file, that is, the file "fname" must be "filxcol" bytes long.
 It returns '1' if success, '0' if failure, and '-1' if there's
 a file open error .
 ***************************************************************/
int BYTEheader(char *fname, int fil, int col)
{
    unsigned long tam;
    FILE *fp;
    struct stat info;

    fp = fopen(fname, "rb");
    if (!fp) return(-1);

    if (!fstat(fileno(fp), &info))
    {
        fclose(fp);
        tam = fil * col * sizeof(char);
        if (tam == info.st_size)   return 1;
        else return 0;
    }

    fclose(fp);
    return 0;
}


/***************************************************************
 This routine read/write the geometric information from/to the
 header.
 ***************************************************************/
void header_geo(double matriz[3][3], float *angle, int rdwr)
{
    int i;

    switch (rdwr)
    {
    case WRITING :
        for (i = 0; i < 3; i++)
        {
            cabecero.fGeo_matrix[i][0] = matriz[i][0];
            cabecero.fGeo_matrix[i][1] = matriz[i][1];
            cabecero.fGeo_matrix[i][2] = matriz[i][2];
        }
        cabecero.fAngle1 = *angle;
        break;

    case READING :
        for (i = 0; i < 3; i++)
        {
            matriz[i][0] = cabecero.fGeo_matrix[i][0];
            matriz[i][1] = cabecero.fGeo_matrix[i][1];
            matriz[i][2] = cabecero.fGeo_matrix[i][2];
        }
        *angle = cabecero.fAngle1;
        break;
    }

}

/***************************************************************
 This function check if there's geometric information in the
 matrix ( non zero values ); returns '1' if this information
 exists and '0' in other case.
 ***************************************************************/
int IsInfogeo(double matriz[3][3])
{
    int i, j;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            if (matriz[i][j] != 0.)  return 1;

    return 0;
}


/***************************************************************
  This simple routine initializes the header 'cabecero' with
  zeros.
 ***************************************************************/
void ceroCabecero(void)
{
    int i;
    BYTE *punt;

    punt = (BYTE *) & cabecero;
    for (i = 0; i < 1024; i++) *punt++ = 0;
    /** *(punt+i)='\0'; Esto es mas lento **/
}


/***************************************************************
  And this one sets a 3x3 matrix to the Identity.
 ***************************************************************/

void identMatrix(double matrix[3][3])
{
    int i;

    for (i = 0; i < 9; i++)
    {
        if (i % 3 == i / 3)
            matrix[i/3][i%3] = 1;
        else
            matrix[i/3][i%3] = 0;
    }
}

/***************************************************************
 This kind routine sets the geometric information in the header
 to default values,that is, sets the angle=0 and the matrix with
 geometric information equal to the the Identity.
 ***************************************************************/
void defaultHeader(void)
{

    double matriz_geo[3][3];
    float angulo = 0.;

    identMatrix(matriz_geo);
    header_geo(matriz_geo, &angulo, WRITING);

}

/***************************************************************
 This routine normalize an image and modify the header setting
 up the fields with average, sigma, maximum and minimum values.
 If the I/O operation is 'read' check if the image has alredy
 been normalized , and in that case exits.
 ***************************************************************/
void normalize_io(float **image, int rdwr, char *name)
{
    float fmax, fmin;

    switch (rdwr)
    {
    case WRITING :
        normImg(image, (int)cabecero.fNrow, (int)cabecero.fNsam,
                &fmax, &fmin);
        cabecero.fAv = 0.;
        cabecero.fSig = 1.;
        cabecero.fImami = 1.;
        cabecero.fFmax = fmax;
        cabecero.fFmin = fmin;
        break;

    case READING :
        if (cabecero.fAv != 0. || cabecero.fSig != 1.)
        {
            normImg(image, (int)cabecero.fNrow, (int)cabecero.fNsam,
                    &fmax, &fmin);
            cabecero.fAv = 0.;
            cabecero.fSig = 1.;
            cabecero.fImami = 1.;
            cabecero.fFmax = fmax;
            cabecero.fFmin = fmin;
        }
        break;
    }

    if (cabecero.fFmax == cabecero.fFmin)
    {
        fprintf(stdout, "\n WARNING:  file %s : MAX and MIN values are identical\n", name);
        fflush(stdout);
    }
}

/***************************************************************
 This is the routine that does the work of normalize the image.
 The image is in a "iNrow*iNsam" matrix of floats called
 "fimagen".
 ***************************************************************/
void normImg(float **fImagen, int iNrow, int iNsam, float *ppmax, float *ppmin)
{

    int i, j;
    float fmax = -1e36 , fmin = 1e36;
    double aux;
    double fmedia, fsigma;
    double suma1 = 0;
    double suma2 = 0;
    double numpuntos;
    double theroot;

    for (i = 0; i < iNrow ; i++)
        for (j = 0; j < iNsam ; j++)
        {
            aux = fImagen[i][j];
            suma1 += aux;
            suma2 += aux * aux;
        }

    numpuntos = (double) iNrow * (double) iNsam;
    /****
      fmedia =  (suma1/(double)iNrow)/(double)iNsam ;
      fsigma =  sqrt ( ((suma2/(double)iNrow)/(double)iNsam)-(fmedia*fmedia));
    ****/

    fmedia = (double)(suma1 / numpuntos);
    theroot = (double)(suma2 / numpuntos) - (fmedia * fmedia) ;
    fsigma =  sqrt(fabs(theroot));


    /***** setting average=0 and sigma=1 *****/

    for (i = 0; i < iNrow ; i++)
        for (j = 0; j < iNsam ; j++)
        {
            fImagen[i][j] = (fImagen[i][j] - fmedia) / fsigma;
            if (fImagen[i][j] < fmin) fmin = fImagen[i][j];
            if (fImagen[i][j] > fmax) fmax = fImagen[i][j];
        }
    *ppmin = fmin;
    *ppmax = fmax;
}

/**************************************************************************/
/******************************* VOLUMES **********************************/

/***************************************************************
 This function open the file "fname" and check if it is a
 SPIDER volume. It returns '1' and the voume dimensions if
 success, '0' if failure, and '-1' if there's  a file open error .
 ***************************************************************/
int SPIDvolum(char *fname, int *slice, int *row, int *col)
{
    unsigned long tam;
    FILE *fp;
    struct stat info;

    fp = fopen(fname, "rb");
    if (!fp) return(-1);   /*unable to open file*/

    if (fread(&cabecero, SPID_HDR_SIZE , 1, fp) != 1)
    {
        fclose(fp);
        return(0) ;       /*not a SPIDER file*/
    }

    if ((fabs(cabecero.fIform)) != 3.)
        return (0);

    if (!fstat(fileno(fp), &info))
    {
        fclose(fp);
        tam = (unsigned long)(cabecero.fNsam * cabecero.fLabrec * 4) +
              (unsigned long)(abs((int)cabecero.fNslice) * cabecero.fNsam * cabecero.fNrow * 4);
        if (tam == info.st_size)
        {
            *slice = (int)abs((int)cabecero.fNslice) ;
            *row = (int)cabecero.fNrow;
            *col = (int)cabecero.fNsam;
            return 1;
        }
    }
    else
        fclose(fp);
    return 0;
}

/***************************************************************
 This function open the file "fname" and check if it is a BYTE
 voume, that is, the file "fname" must be "slice x row x col"
 bytes long. It returns '1' if success, '0' if failure, and '-1'
 if there's a file open error .
 ***************************************************************/
int BYTEvolum(char *fname, int slice, int row, int col)
{
    unsigned long tam;
    FILE *fp;
    struct stat info;

    fp = fopen(fname, "rb");
    if (!fp) return(-1);

    if (!fstat(fileno(fp), &info))
    {
        fclose(fp);
        tam = slice * row * col * sizeof(char);
        if (tam == info.st_size)   return 1;
        else return 0;
    }

    fclose(fp);
    return 0;
}

/***************************************************************
 This routine normalize a voulme and modify the header setting
 up the fields with average, sigma, maximum and minimum values.
 If the I/O operation is 'read' check if the volume has alredy
 been normalized , and in that case exits.
 ***************************************************************/
void norm_of_volume(float ***volumen, int rdwr, char *name)
{

    float fmax = -1e36, fmin = 1e36;

    switch (rdwr)
    {
    case WRITING :
        normVol(volumen, (int)cabecero.fNslice,
                (int)cabecero.fNrow , (int)cabecero.fNsam,
                &fmax , &fmin);
        cabecero.fAv = 0.;
        cabecero.fSig = 1.;
        cabecero.fImami = 1.;
        cabecero.fFmax = fmax;
        cabecero.fFmin = fmin;
        break;

    case READING :
        if (cabecero.fAv != 0. || cabecero.fSig != 1.)
        {
            normVol(volumen, (int)cabecero.fNslice,
                    (int)cabecero.fNrow , (int)cabecero.fNsam,
                    &fmax , &fmin);
            cabecero.fAv = 0.;
            cabecero.fSig = 1.;
            cabecero.fImami = 1.;
            cabecero.fFmax = fmax;
            cabecero.fFmin = fmin;
        }
        break;
    }

    if (cabecero.fFmax == cabecero.fFmin)
    {
        fprintf(stdout, "\n WARNING:  file %s : MAX and MIN values are identical\n"                      , name);
        fflush(stdout);
    }
}

/***************************************************************
 This is the routine that does the work of normalize the volume.
 The volume is in a "capas x fil x col " array of floats called
 "volu".
 ***************************************************************/
void normVol(float ***volu, int capas, int fil, int col,
             float *fvmax, float *fvmin)
{
    float themax = -1e36, themin = 1e36;
    double aux;
    double fmedia, fsigma;
    double suma1 = 0;
    double suma2 = 0;
    double numpuntos;
    double theroot;

    int i, j, k;

    for (i = 0; i < capas; i++)
        for (j = 0; j < fil ; j++)
            for (k = 0; k < col; k++)
            {
                aux = volu[i][j][k];
                suma1 += aux;
                suma2 += aux * aux;
            }

    numpuntos = (double)capas * (double) fil * (double) col;

    fmedia = (double)(suma1 / numpuntos);
    theroot = (double)(suma2 / numpuntos) - (fmedia * fmedia);
    fsigma =  sqrt(fabs(theroot));

    /***** Setting average=0 and sigma=1 *****/

    for (i = 0; i < capas; i++)
        for (j = 0; j < fil ; j++)
            for (k = 0; k < col; k++)
            {
                volu[i][j][k] = (float)((volu[i][j][k] - fmedia) / fsigma);
                if (volu[i][j][k] > themax) themax = volu[i][j][k];
                if (volu[i][j][k] < themin) themin = volu[i][j][k];

            }
    *fvmax = themax;
    *fvmin = themin;

}
/**********************************************************/
