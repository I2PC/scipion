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

/***********************************************************************/
/* Useful definitions for the programs that compose the GROE system    */
/*      Juan P. Secilla    IBM/MSC     Nov/86 (modyfied Jan/89)        */
/***********************************************************************/

/******************** Basic type definitions ***************************/
#ifndef H_GROE
#define H_GROE

#include "spider.h"

#include <cstdio>

typedef unsigned char  BYTE;       /*** Only this and float are used ***/
typedef unsigned short UWORD;
typedef unsigned long  ULONG;

/******************* Constant definitions ******************************/

#define ERROR        -1                /* to show an error condition   */
#define OK            0                /* to indicate a routine is OK  */
#define AYUDA         2                /* To indicate that help is req */

#define NATURAL       1                /* natural, 1 byte/pixel format */
#define INTFMT        2                /* integer, 2 byte/pixel format */
#define LONGFMT       3                /* long, 4 bytes/pixel format   */
#define FLOATFMT      4                /* float, 4 bytes/pixel format  */
#define FOURIER       5                /* Fourier transform format     */
#define SPIDER        6                /* Spider (header) format       */

#define READING       5                /* File read operation          */
#define WRITING       6                /* File write operation         */

#define TRUE          1                /* Guess what?                  */
#define FALSE         0

#define DIRECT        0                /* Direct FFT is to be performed*/
#define INVERSE       1                /* Inverse FFT "    "    "      */

#define NEIGHBOUR 0                    /* Kinds of interpolation       */
#define BILINEAR 1
#define ON 1                           /* Cursor status                */
#define OFF 0
#define CROSS 0                        /* Cursor shape                 */
#define SQUARE 1

#define LINEAR      0                  /* Dynamic range modifications  */
#define LOGARITHMIC 1
#define SQRT        2
#define POWER       3

#define ALLOC_SIZE  65000              /* Allocation unit size         */

#define ANCHO    1024                  /* Screen dimension, in pixels  */
#define ALTO      768

#define PI 3.141592654f                 /* Guess what?                  */

/******************* Macro definitions *********************************/

#define max(x, y)     ((x)>(y)?(x):(y))
#define min(x, y)     ((x)<(y)?(x):(y))
/*#define abs(x)        ((x)<0?(-(x)):(x))*//*esta ya definido*/

/******************* Type declarations *********************************/

typedef BYTE           ** byte_image;
typedef UWORD          ** int_image;   /* Several image formats        */
typedef ULONG          ** long_image;
typedef float          ** float_image;


/******************* Constantes para el fichero de selecci¢n ***********/

#define DESCARTADA                 -1
#define CORTADA                     0
#define PRECENTRADA                 1
#define CENTRADA(veces)     (veces*2)      /*** Centrada "veces" veces ***/
#define GIRADA(veces)     (veces*2+1)      /*** Girada "veces" veces   ***/


#define MAX_GREY 247.

/*** Funtions in cabecero.c **/
#define SPID_HDR_SIZE   sizeof(CABECERO)
/** this is because there is an old funtion with this name **/
#define Formato( a, b, c) SPIDheader( a,&b,&c )

extern CABECERO cabecero;

/**
 Funtions in cabecero.c :
**/

int SPIDheader(char *, int *, int *);
int BYTEheader(char *, int , int);
int SPIDvolum(char *, int *, int *, int *);
int BYTEvolum(char *, int, int , int);

void header_geo(double matriz[3][3], float *, int);
int IsInfogeo(double matriz[3][3]);
void ceroCabecero(void);
void identMatrix(double matrix[3][3]);
void defaultHeader(void);

void normalize_io(float **, int , char *);
void norm_of_volume(float ***, int , char *);
void normImg(float **, int, int, float *, float *);
void normVol(float ***, int, int, int, float *, float *);
/**
 Funtions in varios.c :
**/
int image_io(char **, char **, int, int, char *, int, int, int);
void **imalloc(int, int, int);
void imfree(char **, int, int, int);
int exists(char *);
void asigna_extension(char *, int);
void busca_maximos(float **, int, int, float *, int *, int *, int);
void Tiempo(void);
void Cabecera(void);
float ScanfMejorFloat(void);
/**
 Funtions en varios2.c :
**/
void cornor(float **, int, int, float **, int, int, double);
void transforma(float **, float **, int, int, double matriz[3][3]);
void compon(double m1[3][3], double m2[3][3], double m3[3][3]);
/**
 Funtions en varios3.c :
**/
int capa_io(char **, int, int, int, int, int);
void ***trialloc(int , int , int , int);
void trifree(void ***, int, int, int, int);
/*int volume_io (void ***,int,int,int,char *,int,int,int);*/
int volume_io(char ***, int, int, int, char *, int, int, int);
void vox_4_a_1(float ***, BYTE ***, int);
void vox_1_a_4(float ***, BYTE ***, int);
/**
 Funtions in dimension.c :
**/
int Count_element(FILE *, char *, int, int);
int  First_name(FILE *, char *, int, int);
int Get_dimension(char *, int *, int *);
int Get_Vol_dimension(char *, int *, int *, int *);
/**
 Funtions in floatnat.c :
**/
void oraCharoraFloat(BYTE **, float **, int, int, int);
void volCharvolFloat(BYTE ***, float ***, int, int, int, int);
void InterpBilinear(BYTE **, BYTE **, int, int, int, int);
void matriz_tira(unsigned char **, unsigned char *, int, int,
                 unsigned long *, unsigned long *, int, int);
void mFloat_tira(float **, unsigned char *, int, int, int,
                 unsigned long *, unsigned long *);
void tira_matriz(unsigned char *, unsigned char **, int, int);
/**
 Funtions in iocabecero.c :
**/
void rdwr_geo(double matriz[3][3], float *, GEO_INFO *, int);
void infog_cabecero(GEO_INFO *, int);
/**
 Funtions in usages.c
**/
int check_comm(char **);
void getTempoName(char *);
void getDate(char *, char *);
void getNameAbs(char *, char *);
/**
 Funtion in fft.c
**/
void fftcol(float **, int, float **, int, int, int, int);
int image_fft(float **, int , int , int);
int no_correct(int , int *);
void rfft(float *, int, int, int);
void cfft(float *, float *, int, int, int);
int volume_fft(float ***, int , int, int , int);
int fft_pow2(float **x, int row, int col, int kind);
int vfft_pow2(float ***vol, int sli, int row, int col, int kind);
/**
 Funtions in fftmod.c
**/
float  fmodulo(float, float);
int fast_fourier_trans_img(float **, int, int);
/**
 Funtions in filtros.c
**/
void filtra_pb(float **, float **, int, int);
void cos_alz(float **, int, int, int, int, int);

/**
 Funtions in specmod.c
**/
void spectro_m(char *nom, float xr1, float xr2, float xdr, float xr,
               float *result);
void leespe(FILE *in, short *ir, short *numin, short *numax,
            float *x0, float *y0, float *rl, float *rh, float *dr, float *ac, float *as);
/** End file **/
#endif
