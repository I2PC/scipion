/***************************************************************************
 *
 * Author:     Monica Chagoyen          monica@b12sg1.cnb.uam.es
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

/****************************************************************************/
/* Program for finding the center of an image. The original in Fortran was  */
/* written by A. Santisteban. For any information about the program, contact*/
/* him. Translated into C by Juan P. Secilla (MSC)  Jun/86		    */
/****************************************************************************/
/****************************************************************************/
/* Created a wrapper to fit NewXmipp 					    */
/* Alberto Pascual, October 2001 					    */
/* pascual@cnb.uam.es 					    		    */
/****************************************************************************/

/********************** Include's and macro definitions *********************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <OldXmipp/spider.h>
#include <OldXmipp/groe.h>
#include <XmippData/xmippArgs.hh>

//CABECERO cabecero;

#define DEF_IT     5  /*Default iterations number */
#define DEF_DEL    2
#define DEF_IN     2

/************************* Global variables *********************************/

float coseno [1281];
BYTE **imagen;
float r1, r2, r3, cxp1, cxp2, cxp3, cxm1, cxm2, cxm3, rh = 1000., xc0, yc0, del,
      rbajo, ralto, zzz0;
int largo, lancho, indmul, in, ni, mu, ncic;
int ir, m, mu1, mu4, ntot,  mt, idz, ncic2, ntot4;
/*
float conv1x (double, double);
void busca (), suaviza (), ergrot (double, double, float *), Usage(char *);
*/

#if defined (__STDC__) || defined (_CCAIX3)
  extern void busca();
  extern void suaviza();
  extern void ergrot (double, double, float *);
  extern float conv1x (double, double);
  static void Usage( char *);
#else
  extern void busca();
  extern void suaviza();
  extern void ergrot();
  extern float conv1x ();
  static void Usage();
#endif

/****************************************************************************/

main (int argc, char **argv)
{
    short i,j;
    FILE *inp_image;
    static char nom_imagen[30];
    char *ent, *tmp;
    short ancho, alto;
    int tipo;
    float **imgtmp, **imgtmp2;
    BYTE **imgbyte;

/*if (argc < 11 || argc > 14)
 {
   Usage(argv[0]);
   exit(1);
 }
*/

if (check_param(argc, argv, "-img"))
   ent = (char*) get_param(argc, argv, "-img", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
strcpy (nom_imagen, ent);


//xc0=atof(argv[2]);    /* Initial coordinates of the centre */
//yc0=atof(argv[3]);

if (check_param(argc, argv, "-x0"))
  tmp = (char*) get_param(argc, argv, "-x0", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
xc0 = atof(tmp);

if (check_param(argc, argv, "-y0"))
  tmp = (char*) get_param(argc, argv, "-y0", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
yc0 = atof(tmp);

xc0++;
yc0++;

if (check_param(argc, argv, "-low"))
  tmp = (char*) get_param(argc, argv, "-low", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
rbajo = atof(tmp);
if (check_param(argc, argv, "-high"))
  tmp = (char*) get_param(argc, argv, "-high", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
ralto = atof(tmp);

//rbajo=atof(argv[4]);  /* Smothing radius (low, high) */
//ralto=atof(argv[5]);

if (check_param(argc, argv, "-r1"))
  tmp = (char*) get_param(argc, argv, "-r1", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
r1 = atof(tmp);

if (check_param(argc, argv, "-r2"))
  tmp = (char*) get_param(argc, argv, "-r2", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
r2 = atof(tmp);

tmp = (char*) get_param(argc, argv, "-rinc", "1");
r3 = atof(tmp);


//r1=atof(argv[6]);     /* Integrating radius (low, high, increment */
//r2=atof(argv[7]);
//r3=atof(argv[8]);

tmp = (char*) get_param(argc, argv, "-harm", "1");
ncic = atoi(tmp);
tmp = (char*) get_param(argc, argv, "-opt", "-1");
indmul = atoi(tmp);

//ncic=atoi(argv[9]);    /* Harmonic to optimize */
//indmul=atoi(argv[10]); /* Maximize(+1) or minimize (-1)*/


if (indmul!=1 && indmul!=-1)
{
  printf("Error %s > Optimization not valid\n", argv[0]);
  Usage (argv[0]);
  exit(1);
}
del=DEF_DEL;
in=DEF_IN;
ni=DEF_IT;
/*for (i=11; i<argc;)
 {
  if (!strncmp(argv[i],"-i", 2))
   {
     ni=atoi(argv[i]+2);
     i++;
   }
  else if (!strncmp(argv[i],"-g",2))
   {
     del=atoi(argv[i]+2);
     i++;
   }
  else if (!strncmp(argv[i],"-p",2))
   {
     in=atoi(argv[i]+2);
     i++;
   }
  else
   {
     printf("\n %s: bad option \n",argv[i]);
     fflush(stdout);
     Usage(argv[0]);
     exit(1);
   }
 }*/
mu= (int) (PI/2*r2/ncic);

if (mu<3)
{
  printf("\nA HIGHER INTEGRATION RADIUS IS NEEDED (R2 > 6*ORDEN/PI)\n");
  exit(1);
}
//strcpy(nom_imagen, argv[1]);
if (!exists (nom_imagen))
{
  printf(" %s : File not found ( %s ).\n", nom_imagen, argv[0]);
  exit(1);
}

ceroCabecero();
if ((tipo=SPIDheader(nom_imagen, &largo, &lancho))<0)
  {

  printf ("There is a problem opening the file %s.\n", nom_imagen);
      exit(1);
    }
    if (!tipo) /** not SPIDER file **/
    {
       printf("\n %s > ERROR: %s is not in SPIDER format !!\n", argv[0],nom_imagen);
       exit(1);
    }
/************************** ALLOC MEMORY ********************************/

if ((imgtmp=(float **)imalloc (largo, lancho, FLOATFMT))==NULL)
 {
  printf("\n%s > Sorry, no memory", argv[0]);
  exit(1);
 }

if ((imgtmp2=(float **)imalloc (largo, lancho, FLOATFMT))==NULL)
 {
  printf("\n%s > Sorry, no memory", argv[0]);
  exit(1);
 }

if ((imgbyte=(unsigned char **)imalloc (largo, lancho, NATURAL))==NULL)
 {
  printf("\n%s > Sorry, no memory", argv[0]);
  exit(1);
 }
if ((imagen=(unsigned char **)imalloc (largo+1, lancho+1, NATURAL))==NULL)
 {
  printf("\n%s > Sorry, no memory", argv[0]);
  exit(1);
 }

/* This temporal image is needed due to the different indexing rules in
   Fortran  */
/************************* REDING THE FILE ****************************/

if ((image_io((char **)imgtmp2,(char **)imgtmp,largo,lancho,
          nom_imagen,READING,FLOATFMT,TRUE))==-1)
 {
  printf("%s > Error: there's a problem reading the file %s", argv[0],nom_imagen);
  imfree((char **)imgtmp,largo,lancho,FLOATFMT);
  imfree((char **)imgtmp2,largo,lancho,FLOATFMT);
  imfree((char **)imagen,largo+1, lancho+1, NATURAL);
  imfree((char **)imgbyte,largo,lancho,NATURAL);
  exit(1);
 }

oraCharoraFloat(imgbyte,imgtmp,FLOATFMT,largo,lancho);

for (i=0; i<largo; i++)
 for (j=0; j<lancho; j++)
   imagen[i+1][j+1]=imgbyte[i][j];

imfree((char **)imgtmp,largo,lancho,FLOATFMT);
imfree((char **)imgtmp2,largo,lancho,FLOATFMT);
imfree((char **)imgbyte,largo,lancho,NATURAL);

    ancho = lancho;
    alto = largo;

    busca ();
imfree((char **)imagen,largo+1, lancho+1, NATURAL);

}

/****************************************************************************/

void ergrot(double xc0, double yc0, float* zzz) {
//double xc0, yc0;
//float *zzz;	/* It hides global variable, handle with care */

    static float a[266], b[257];
    double za,zb;
    double xp1,xp2,xp3,xm1,xm2,xm3;
    double axp1,axp2,axp3,axm1,axm2,axm3;
    short i7, iz2, iz1, kv, l1, l2, i1, i, j;
    float r, x, y, zz, ai, bi;
    float bj, aj, z;

    r=r1-r3;
    axp1=0.;
    axm1=0.;
    axp2=0.;
    axm2=0.;
    axp3=0.;
    axm3=0.;
    for (i7 = 1; i7 <= ir; i7++) /* do 31 */
    {	r+=r3;
	iz2=ntot-m+mu4;
	iz1=iz2-1;
	for (i = 1; i <= mu; i++) /* do 13 */
	{   iz1+=idz;
	    iz2+=idz;
	    za=0.;
	    zb=0.;
	    for (kv = 1; kv <= ncic; kv++) /* do 11 */
	    {	iz1+=m;
		x=xc0+r*coseno[iz1];
		y=yc0+r*coseno[iz1-mu4];
		z = conv1x(y,x);
		iz1+=m;
		x=xc0+r*coseno[iz1];
		y=yc0+r*coseno[iz1-mu4];
		zz = conv1x(y,x);
		zb+=(z-zz);
		iz2+=m;
		x=xc0+r*coseno[iz2];
		y=yc0+r*coseno[iz2-mu4];
		z = conv1x(y,x);
		iz2+=m;
		x=xc0+r*coseno[iz2];
		y=yc0+r*coseno[iz2-mu4];
		zz = conv1x(y,x);
		za+=(z-zz);
	    }
	    b[i]=zb;
	    a[i]=za;
	}
	xp1=a[mu];
	xp2=b[mu];
	xm2=xp2*xp2;
	xp3=xp1*xp2*coseno[mu4-ncic];
	xm3=xp3;
	xp2=xm2*coseno[mu4-ncic2];
	xp1=xp1*xp1;
	xm1=xp1;
	for (i = 1; i<= mu1; i++) /* do 14 */
	{   l1=4*i*ncic+mu4;
	    l2=mu4;
	    ai=a[i];
	    bi=b[i];
	    xp1+=(ai*ai*coseno[l1]);
	    xp2+=(bi*bi*coseno[l1-ncic2]);
	    xp3+=(ai*bi*coseno[l1-ncic]);
	    xm1+=(ai*ai);
	    xm2+=(bi*bi);
	    xm3+=(ai*bi*coseno[mu4+ncic]);
	    i1=i+1;
	    for (j = i1; j <= mu; j++) /* do 15 */
	    {	l1+=ncic2;
		l2+=ncic2;
		ai=a[i];
		aj=a[j];
		bi=b[i];
		bj=b[j];
		xp1+=(2.*aj*ai*coseno[l1]);
		xm1+=(2.*aj*ai*coseno[l2]);
		xp2+=(2.*bj*bi*coseno[l1-ncic2]);
		xm2+=(2.*bj*bi*coseno[l2]);
		xp3+=((ai*bj+aj*bi)*coseno[l1-ncic]);
		xm3+=(ai*bj*coseno[l2-ncic]+aj*bi*coseno[l2+ncic]);
	    }
	}
	axp1+=xp1*r;
	axm1+=xm1*r;
	axp2+=xp2*r;
	axm2+=xm2*r;
	axp3+=xp3*r;
	axm3+=xm3*r;
    }
    (*zzz)=axp1*cxp1+axp2*cxp2+axp3*cxp3+axm1*cxm1+axm2*cxm2+axm3*cxm3;
    (*zzz)/=(zzz0*ir);
}

/****************************************************************************/


void suaviza()
{
    unsigned char pix;
    short i, j, k, n;
    long isuma;
    float racua, rbcua, dr, r;


    isuma=0;
    n=0;
    i=0;
    racua=ralto*ralto;
    rbcua=rbajo*rbajo;
    dr=3.141592653/(ralto-rbajo);
    for (k = 1; k <= lancho; k++)
    {	if (k%16 == 0)
	 /*   printf ("\n%d", k);*/
	for (j = 1; j <= largo; j++)
	{   r=(xc0-k)*(xc0-k)+(yc0-j)*(yc0-j);
	    if (r > racua)
	    {	isuma += imagen [k][j];
		n++;
	    }
	}
    }
    if (n != 0)
	i= (short int) (((float) isuma)/n + .5);
    m=i;
    pix=i;
    for (k = 1; k <= lancho; k++)
    {	if (k%16 == 0)
	  /*  printf ("\n%d", k);*/
	for (j = 1; j <= largo; j++)
	{   r=(xc0-k)*(xc0-k)+(yc0-j)*(yc0-j);
	    if (r < rbcua)
		continue;
	    if (r > racua)
	    {	imagen[k][j]=pix;
		continue;
	    }
	    r=sqrt(r);
	    r-=rbajo;
	    imagen[k][j]=(unsigned char) ((0.5+0.5*cos(dr*r))*(imagen[k][j]-m)+m);
	}
    }
    printf ("\nImage smoothed. Outer mean = %d\n", m);
}

/***************************************************************************/

float conv1x (double y, double x)

//double y, x;

/**************************************************************************/
/* Returns the value of the image at the point (y, x) using linear interp-*/
/* olation. For higher accuracy, use "conv3x" (available in 370 assembly  */
/* languaje.								  */
/**************************************************************************/
{   short i, j; 	 /* Row and column */
    float intfila1, intfila2;	/* Partial row interpolations */
    float escala;   /* Scale factor */

    j = (short int) y;     /* Trunc the x, y coordinates to short */
    i = (short int) x;

    escala =  y - j;
    /* 1st row interpolation */
    intfila1 = imagen[i][j] + escala*((short)imagen[i][j+1] - (short)imagen[i][j]);
    /* 2nd row interpolation */
    intfila2 = imagen[i+1][j] + escala*((short)imagen[i+1][j+1] - (short)imagen [i+1][j]);
    /* Column interpolation */
    return intfila1 + (x - i)*(intfila2 - intfila1);
}

void Usage (char *name)
{
   printf("\nUsage: \n");
   puts ("xmipp_findcenter  [parameters]");
   puts ("Purpose: Finds the center of symmetry of a given image");
   puts("Parameters:");
   puts ("       -img         : Input file name");
   puts ("       -x0  -y0     : center coordinates");
   puts ("       -low  -high  : smoothing radii");
   puts ("       -r1 -r2      : Integration radius (low high)");
   puts ("       [-rInc]      : Integration increment (1 by default)");
   puts ("       [-harm]      : Harmonic to optimize (1 by default)");
   puts ("       [-opt]       : Type of optimization");
   puts ("                      1 maximize ");
   puts ("                      -1 minimize (default)");
   puts ("");
   puts ("Example:");
   puts ("xmipp_findcenter -img file.ave -x0 25 -y0 25 -low 23 -high 25 -r1 10 -r2 20");
   puts ("");

/* printf("\nUsage: \n  %s",name);
 puts (" file_name x0 y0 low high R1 R2 RIn harm +1/-1 ");
 puts (" [-gGrid] [-pPoints] [-iIter]");
 puts (" where:");
 puts ("         - file_name      : input file (SPIDER)");
 puts ("         - x0 - y0        : starting centre coordinates");
 puts ("         - low - high     : smoothing radius");
 puts ("         - R1 R2 RIn      : integration radius (low high increment)");
 puts ("         - harm           : harmonic to optimize");
 puts ("         - +1/-1          : maximize (+1) or minimize (-1)");
 puts ("         - [-gGrid]       : interval grid (default: 2)");
 puts ("         - [-pPoints]     : grid points (default: 2)");
 puts ("         - [-iIter]       : number of iterations (default: 5)");*/
}

