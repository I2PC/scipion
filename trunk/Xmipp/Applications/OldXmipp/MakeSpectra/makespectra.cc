/***************************************************************************
 *
 * Authors:     Irene Martinez          irene@b12sg1.cnb.uam.es
 *              Roberto Marabini        roberto@b12sg1.cnb.uam.es
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
/* Program for rotational filtering  of images. Original in Fortran   */
/* by A. Santisteban. Translated into C by J. P. Secilla and Pp Acu¤a */
/* Santisteban is the one to ask any question to.                     */
/**********************************************************************/
/* Input image format: natural, largo x lancho BYTEs                 */
/* Spectrum files: 2 bytes : ir, numin, numax                         */
/*                 4 bytes : x0, y00, r1, r2, r3                      */
/*                 4 bytes : ampcos (numax-numin+1)*ir elements       */
/*                 4 bytes : ampcen (numax-numin+1)*ir elements       */
/* (ir elemets represent a certain rot. frequency)                    */
/*                                                                    */
/*                            Juan P. Secilla Sept. 86                */
/**********************************************************************/
/* Mas modificaciones para el sistema GroE. Mery M. Nov. 89           */
/* S¢lo realiza la primera parte del recons y ademas saca el fichero  */
/* en ascii para poder leerlo.                                        */
/* NOTA: S¢lo esta implementada la primera opcion del recons.         */
/**********************************************************************/

/**********************************************************************/
/* Modificaciones por Alberto Pascual (pascual@cnb.uam.es)            */   
/* para soportar cualquier nombre de sel y simplificar la salida.     */
/* 1 de Octubre de 2001                                               */   
/**********************************************************************/


/********************* macroes & constants definitions ****************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <OldXmipp/spider.h>
#include <OldXmipp/groe.h>
#include <XmippData/xmippArgs.hh>

#define eprintf(x) printf(x);fflush(stdout)
#define eprintf2(x,y) printf(x,y);fflush(stdout)
#define eprintf3(x,y,z) printf(x,y,z);fflush(stdout)

#define MAX_SIM 51      /*** concuerda con el n£mero m ximo de arm¢nicos ***/
                        /*** definido en el specmod (mod. del spectro    ***/


/***************** Global variables (originally a common) ********/

short numin,numax,largo,lancho;
int kind;
float r1,r2,r3,x0,y00;
BYTE **imagen;   /* Pointer to the rows of the image **/
FILE  *new_esp_file;


static void recons(void);
static float conv1x (double y,double x);
static void Usage( char *name);

/********************* Main routine ***********************************/

main (int argc, char **argv)

{
int i, j, total;
int nnumin,nnumax,llargo,llancho;
int dim, codigo;    
int imagen_terminada;
int sim;
int rf1, rf2, rf3, rf4;
int filtra = FALSE;          /**** Filtro paso bajo           ****/
int fil_cos = FALSE; 
int chek_path = FALSE;

float **im_float_1;                   /**** imagen en float (no sirve) ****/
float **im_float_2;                   /**** imagen en float (no sirve) ****/
float max_img, min_img;               /**** para escalar la img.       ****/
float aux;
float rad_min, rad_max, inc_int, long_int; /*** para pasar a ascii (spectro)**/
float res_sim[MAX_SIM];

char *plin, *pnumero;
char nom_entrada[200], linea[200], nombre[100];
char *ent, *sal, *tmp;
char geo_name[100];
char nom_path[200], tempo_nombre[200];
char nom_info[100];   /**** Nombre fichero informaci¢n ****/
char lin[100], numero[100];
char actDay[48], actHour[32];

FILE *info;
FILE *file_sel;


kind = 1;    /**** viene del recons ****/

if (check_param(argc, argv, "-sel"))
   ent = (char*) get_param(argc, argv, "-sel", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }

strcpy (nom_entrada, ent);

if (check_param(argc, argv, "-out"))
   sal = (char*) get_param(argc, argv, "-out", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
strcpy (nom_info, sal);
   
if (check_param(argc, argv, "-x0"))
  tmp = (char*) get_param(argc, argv, "-x0", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
y00 = atof(tmp);

if (check_param(argc, argv, "-y0"))
  tmp = (char*) get_param(argc, argv, "-y0", "");
else    
 {  printf ("Error in the arguments");
    Usage(argv[0]);
    exit(1);
 }
  
x0 = atof(tmp);


x0 ++;              /*** problemas del fortran ***/
y00 ++;

tmp = (char*) get_param(argc, argv, "-low", "1");
nnumin = atoi(tmp);
tmp = (char*) get_param(argc, argv, "-high", "15");
nnumax = atoi(tmp);

numin = nnumin; /* No se puede hacer scanf de short en CMS (JPS) */
numax = nnumax;

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

rad_min = r1;   /*** parametros para la rutina spectro ****/
rad_max = r2;
inc_int = r3;

//tmp = (char*) get_param(argc, argv, "-int", "1");
//long_int = atof(tmp); /*** paso a fich. ASCII ****/

long_int = r2-r1;

/* 
cout << "rad_min="  << rad_min  << endl
     << "rad_max="  << rad_max  << endl
     << "inc_int="  << inc_int  << endl
     << "long_int=" << long_int << endl
     << "nnumin="   << nnumin   << endl
     << "nnumax="   << nnumax   << endl;
*/
     
/* ****************  OPENING THE FILES  *********************/


if (!exists (nom_entrada))
{   printf (" %s : File not found.\n", nom_entrada);
    exit (1);
}

if ((file_sel = fopen (nom_entrada, "r")) == NULL)
{  printf ("There is a problem opening the file %s.\n", nom_entrada);
   exit (1);
}

if ((info = fopen (nom_info, "w")) == NULL) /*** Fich. .SIM ***/
{   printf ("There is a problem opening the file %s.\n", nom_info);
    exit (1);
}

sim = nnumax - nnumin + 1;
fprintf (info, "                                   \n");

/*getDate(actDay,actHour);
fprintf (info, "#\n");
fprintf (info, "# DATE: %s  TIME: %s\n",actDay,actHour);
fprintf (info, "SYM %d\n",sim);
fprintf (info, "********");
for (i=0; i < sim; i++)
fprintf (info, "********");
fprintf (info, "\n Image  ");
for (i=0; i < sim; i++)
    fprintf (info, " |  %2d  ", nnumin+i);
fprintf (info, "\n********");
for (i=0; i < sim; i++)
fprintf (info, "********");
putc ('\n', info);*/
fflush (info);



/* ****************   IMAGE DIMENSIONS   ********************/

        /* Open one of the files to get the dimensions */

if(! First_name( file_sel,nom_entrada,DESCARTADA,2))
 {
    fprintf(stderr,"\n Couldn't get the dimensions \n");
    exit(1);
 }
if (chek_path)
 {  strcpy (tempo_nombre, nom_path);
    strcat (tempo_nombre, nom_entrada);
    strcpy (nom_entrada, tempo_nombre);
 }
ceroCabecero();
if(SPIDheader(nom_entrada,&dim,&dim)<=0 )
 {
   fprintf(stdout,"%s not a SPIDER file\n", nom_entrada);
   exit( 1);
 }
llargo = dim; 
llancho = dim;
largo = llargo;
lancho = llancho;
if (r2 > x0-4 || r2 > y00-4 || r2 > largo-x0-3 || r2 > lancho-y00-3)
{   printf ("\n%f> Error: radius r2 is very long",r2);
 exit (1);
}

/* ****************   ALLOC MEMORY   ************************/


if ((im_float_1 = (float **) imalloc (dim, dim, FOURIER)) == NULL)
 {
    fprintf(stdout,"\n%s > Sorry, no memory ",argv[0]);
    fflush(stdout);
    exit(1);
 }

if ((im_float_2 = (float **) imalloc (dim, dim, FOURIER)) == NULL)
 {
    fprintf(stdout,"\n%s > Sorry, no memory ",argv[0]);
    fflush(stdout);
    exit(1);
 }

if ((imagen = (BYTE **) imalloc (dim+1, dim+1, NATURAL)) == NULL)
 {
    fprintf(stdout,"\n%s > Sorry, no memory ",argv[0]);
    fflush(stdout);
    exit(1);
 }

/*** se alloca dim +1 por problemas fortranescos...... ***/


/****************************************************************************/
/*** Empieza el trabajito en serio                                        ***/
/****************************************************************************/

total = 0;
while (fgets (linea, 200, file_sel) != NULL)
{  imagen_terminada = ' ';
   if( (i=sscanf (linea,"%s %d %c",nombre, &codigo, &imagen_terminada))
             >2 || i==0 )
       continue;
   if (codigo <=DESCARTADA) 
       continue;           

   strcpy (nom_entrada, nombre);
   if (chek_path)
    {  strcpy (tempo_nombre, nom_path);
       strcat (tempo_nombre, nom_entrada);
       strcpy (nom_entrada, tempo_nombre);
    }
    
   if (image_io ((char **)im_float_1,(char **)im_float_2, dim, dim, 
          nom_entrada, READING, FLOATFMT,TRUE) == -1)
    {   printf ("%s > Error: image %s not found \n",argv[0], nom_entrada);
        continue;
    }

	/***** filtra paso bajo la imagen transformada ********/

   if (filtra)
    {   cout << "Filtrado 1\n";
        filtra_pb (im_float_2, im_float_1, dim, dim);
        for (i = 0; i < dim; i++)    
            for (j = 0; j < dim; j++)
                im_float_2[i][j] = im_float_1[i][j];
    }

	/***** filtra en coseno alzado                  *******/
   if (fil_cos)
    {
	cout << "Filtrado 1\n";
	/*********** Transformada de Fourier directa **************/

        image_fft (im_float_2, dim, dim, DIRECT);

	/*********** Filtro en coseno alzado **********************/

        cos_alz (im_float_2, dim, rf1, rf2, rf3, rf4);

	/*********** Transformada de Fourier inversa **************/

        image_fft (im_float_2, dim, dim, INVERSE);
    }

         /****** Paso a BYTE  ********/

    max_img = -1e36;
    min_img = 1e36;

    for (i = 0; i < dim ; i++)
       for (j = 0; j < dim ; j++)
       {   if (im_float_2[i][j] > max_img)
              max_img = im_float_2[i][j];
           if (im_float_2[i][j] < min_img)
              min_img = im_float_2[i][j];
       }

    aux = 255/(max_img-min_img);

    for (i = 0; i < dim ; i++)
       for (j = 0; j < dim ; j++)
           imagen[i+1][j+1] =  (unsigned char)
              ((im_float_2[i][j] - min_img) * aux + .5);
                          /*** los indices empiezan en 1 (fortran) ***/

	/**************************/

    strcpy (nom_entrada, nombre);
    strcat (nom_entrada, ".spc");
    if ((new_esp_file = fopen (nom_entrada, "wb")) == NULL)
    {   printf ("\n%s > Couldn't write the file %s ",argv[0],nom_entrada);
        continue; 
    }

	/************************************/

    recons (); /* Do the curre now */
    fclose (new_esp_file);

	/******************************/

    for (i=0; i<sim; i++)  /**** inicializo a 0 ****/
        res_sim[i]= 0;

    /*
    cout << "rad_min=" << rad_min << endl;
    cout << "rad_max=" << rad_max << endl;
    cout << "inc_int=" << inc_int << endl;
    cout << "long_int=" << long_int << endl;
    */

    spectro_m (nombre, rad_min, rad_max, inc_int, long_int, res_sim);
    strcpy (nom_entrada, nombre);
    strcat (nom_entrada, ".spc");
    unlink (nom_entrada);

	/*******************************/
/*    fprintf (info, "%8s ", nombre);
    for (i=0; i < sim; i++)
         fprintf (info, "| %5.2f ", res_sim[i]);
    putc ('\n', info);*/
    
    for (i=0; i < sim; i++)
         fprintf (info, "%5.2f ", res_sim[i]);
    fprintf (info, "%8s ", nombre);
    putc ('\n', info);
    
    fflush (info);
    total++;

} /*** end while ***/

rewind(info);
fprintf (info, "%d %d", sim, total);

fclose (info);
fclose (file_sel);

}/** end MAIN **/

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
float coseno[1281],ampcos[5191],ampsen[5191],peso[51];

static void recons(void)
	/********* Does everything **************/
{
static double ac,as,bc,bs;
static float b1,coefca,coefcb,coefsa,coefsb,add;
static float ck,cx,cy,ccc1,ccc2,css1,css2,cc,cs;
static float d1,dr1,dr2,e1,f,fi,g1,h,hdpi;
static float rh,r11,rhh,r,rl,th, ys,x,y,zs,z,znul,zn, ys2, zs2;
static short ir,k0[51], vlt[2049];
static unsigned char imag[513];
static int nu[51];
static int ind,i1c,i1s,i2c,i2s,irep,irhhhh,indnfr,ick,irl;
static int ixx,iiii,iijjkk,jr,k,kk,larmax,l,ntot,ndim,numfrq;
static short iyy,larxxx;
static int my,my2,my3,my4,my5;
register int i,k5,j,iii,i7; /* Tocar en el AT (2 register) */
static BYTE imgbin[513],bin;
static short i0 = 0, i1 = 1, i2 = 2048, i3 = 3, i11 = 11, imgval = 0;
static int numin1, numax1,lanmax;

/********************** Initialise matrices **************************/

for (i = 1; i <= 50; i++)
    peso[i] = 1.;              /* Esto va por algunos datas originales */
for (i = 1; i <= 512; i ++)
    imgbin [i] = 0;

/****************** Begin to work (Ununderstandable) ******************/

rh=MIN(r2,x0-4.);
rh=MIN(rh,y00-4.);
rh=MIN(rh,largo-x0-3.);
rh=MIN(rh,lancho-y00-3.);

ir=(short int)((rh-r1)/r3+1);
ind=0;
numin1=numin+1;
numax1=numax+1;

for(kk=numin1;kk<=numax1;kk++)
{    k=kk-1;
     if (k != 0)
     {    my=(short int)(1+PI*rh/2./k);
          my2=2*my;
          my4=my*k;
          my5=my4-1;
          ntot=4*my4;
          h=2.*PI/ntot;
          hdpi=h/PI;
          th=k*h;
          ys=sin (th);
          zs=cos (th);
          ys2 = sin (2.*th);
          b1=2./(th*th)*(1.+ zs*zs - ys2/th);
          g1=4./(th*th)*(ys/th-zs);
          d1=2.*th/45.;
          e1=d1*ys*2.;
          d1*=ys2;
          coefca=(b1+e1)*hdpi;
          coefcb=(g1-d1)*hdpi;
          coefsa=(b1-e1)*hdpi;
          coefsb=(g1+d1)*hdpi;
     }
     else
     {
          my=(short int)(1+PI*rh/2.);
          my2=2*my;
          my4=my;
          my5=my4-1;
          ntot=4*my4;
          h=2.*PI/ntot;
          coefca=h/PI/2.;
     }
     for (i=1;i<=my5;i++)
     {   fi=i*h;
         coseno[i] = sin(fi);
     }
     coseno[my4] = 1.;
     my3=2*my4;
     for (i=1;i<=my5;i++)
     {
         coseno[my3-i] = coseno[i];
         coseno[my3+i] = -coseno[i];
         coseno[ntot-i] = -coseno[i];
         coseno[ntot+i] = coseno[i];
     }
     coseno[my3] =0.;
     coseno[my3+my4] = -1.;
     coseno[ntot] =0.;
     coseno[ntot+my4] = 1.;
     r11=r1-r3;
     for (jr=1;jr<=ir;jr++)
     {
         ind++;
         r=r11+r3*jr;
         ac=0.;
         i1c=my4;
         i1s=0;
         if (k != 0)
         {
          as=bc=bs=0.;
          for (i=1;i<=k;i++)
            {
                i2c=my4;
                i2s=0;
                for (j=1;j<=my2;j++)
                {
                    i1c++;
                    i1s++;
                    i2c+=k;
                    i2s+=k;
                    x=x0+r*coseno[i1c];
                    y=y00+r*coseno[i1s];
                    z = conv1x(y,x);
                    bc+=z*coseno[i2c];
                    bs+=z*coseno[i2s];
                    i1c++;
                    i1s++;
                    i2c+=k;
                    i2s+=k;
                    x=x0+r*coseno[i1c];
                    y=y00+r*coseno[i1s];
                    z = conv1x(y,x);
                    ac+=z*coseno[i2c];
                    as+=z*coseno[i2s];
                }
            }
            ampcos[ind] = coefca*ac+coefcb*bc;
            ampsen[ind] = coefsa*as+coefsb*bs;
        }
        else
            for (j=1;j<=ntot;j++)
            {
                i1c++;
                i1s++;
                x=x0+r*coseno[i1c];
                y=y00+r*coseno[i1s];
                z=conv1x(y,x);
                ac+=z;
                ampcos[ind] =coefca*ac;
                ampsen[ind] =0.;
            }
    }
    if (k == 0)
    {
        ac=0.;
        for (j=1;j<=ir;j++)
        {
            r=r11+j*r3;
            ac+=ampcos[j]*2.*PI*r;
        }
        ac/=(PI*(r*r-r1*r1));
        for (j=1;j<=ir;j++)
            ampcos[j] -= ac;
    }
}

/***************** Write new spectrum file *************************/

fwrite ((char *)&ir, sizeof (ir), 1, new_esp_file);
fwrite ((char *)&numin, sizeof (numin), 1, new_esp_file);
fwrite ((char *)&numax, sizeof (numax), 1, new_esp_file);
fwrite ((char *)&x0, sizeof (x0), 1, new_esp_file);
fwrite ((char *)&y00, sizeof (y00), 1, new_esp_file);
fwrite ((char *)&r1, sizeof (r1), 1, new_esp_file);
fwrite ((char *)&r2, sizeof (r2), 1, new_esp_file);
fwrite ((char *)&r3, sizeof (r3), 1, new_esp_file);
fwrite ((char *)&ampcos[1], sizeof (ampcos[1]), (numax-numin+1)*ir,
                new_esp_file);
fwrite ((char *)&ampsen[1], sizeof (ampsen[1]), (numax-numin+1)*ir,
                new_esp_file);

/*
cout << "ir=" << ir << endl
     << "numin=" << numin << endl
     << "numax=" << numax << endl
     << "x0=" << x0 << endl
     << "y0=" << y00 << endl
     << "r1=" << r1 << endl
     << "r2=" << r2 << endl
     << "r3=" << r3 << endl;
for (int i=0; i<(numax-numin+1)*ir; i++) {
//   cout << "ampcos sin[" << i << "]"
   cout << ampcos[i+1] << " " << ampsen[i+1]<< endl;
}*/
}

/***********************************************************************/

static float conv1x (double y,double x)

/***********************************************************************/
/* Returns the value of the image at the posit.  (y,x) using linear    */
/* interpolation. For higher accuracy, use "conv3x" (available in 370  */
/* assembly language).                                                 */
/***********************************************************************/
/** J.P.Secilla **/

{   int i, j;      /* Row and column */
    float intfila1, intfila2;   /* Partial row interpolations */
    float escala;   /* Scale factor */

    j = (int)y;     /* Trunc the x, y coordinates to int */
    i = (int)x;

    escala =  y-j;
    /* 1st row interpolation */
    intfila1 = (short) imagen[i][j] +
    escala*((short) imagen[i][j+1] - (short) imagen[i][j]);
    /* 2nd row interpolation */
    intfila2 = (short) imagen[i+1][j] +
    escala*((short) imagen[i+1][j+1] - (short) imagen [i+1][j]);
    /* Column interpolation */
    return intfila1 + (x - i)*(intfila2 - intfila1);
}
/**********************************************************************/

static void Usage( char *name)
{
   printf("\nUsage: \n");
   puts ("xmipp_makespectra  [parameters]");
   puts ("Purpose: creates the rotational spectra of a set of images");
   puts("Parameters:");
   puts ("       -sel         : Input sel file name");
   puts ("       -out         : Output file name");
   puts ("       -x0  -y0     : Simetry center ");
   puts ("       -low  -high  : minimum and maximum harmonics to be calculated");
   puts ("                      by default low = 1 and high = 15");
   puts ("       -r1 -r2      : Integration radius (low high)");
   puts ("       [-rInc]      : Integration increment (1 by default)");
//   puts ("       [-int]       : Lenght of integration interval (1 by default)");
   puts ("");
   puts ("Example:");
   puts ("xmipp_makespectra -sel file.sel -out myfile.sim -x0 25 -y0 25 -r1 10 -r2 20");
   puts ("");
   puts ("WARNING: THIS PROGRAM DOES NOT APPLY THE ALIGNMENT PARAMETERS IN THE IMAGE HEADERS! ");
}

/*****************************************************************************/
/********************************* The end ***********************************/
/*****************************************************************************/

