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

/**************************************************************************/
/* Input/output from a file to memory. Image is a pointer to pointers to  */
/* the rows of the image. (row, col) are the image dimensions, Name is    */
/* the name of the file to read/write from/to and rdwr = 0 ==> reading,   */

/* rdwr == 1 ==> writing. In this latter case, if the file exists, it is  */
/* erased before writing.                                                 */
/* Uses low-level I/O routines to speed up the operations.                */
/* format indicates the kind of file to read/write (see groe.h)           */
/* It uses a temporal buffer to speed up I/O. We allocate a block of      */
/* 32000 bytes in memory. Reading/writing is performed in this block,     */
/* which implies less calls to I/O routines. This block is copied from    */
/* or to the image using the routine "memcpy", which allows for fast data */
/* transfers.                                                             */
/* Returns: OK if everything O.K., ERROR if error detected                */
/**************************************************************************/
/*     Version 1.0: supports all image formats, image allocator added,    */
/*                  fast input/output using a temporal memory buffer.     */
/*                  Juan P. Secilla    IBM/MSC   Nov/86                   */
/**************************************************************************/

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fcntl.h>
#include <unistd.h>
#include <memory.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifndef __APPLE__
   #include <malloc.h>
#else
   #include <cstdlib>
#endif

#include "spider.h"
#include "groe.h"

#ifdef _PARAMID
void vuelta_cabecero( )
{
	int i, k;
	unsigned int new, *pnew, *ip;
	double ndou;
	unsigned char *p1, *p2;

	ip = (unsigned int *)&cabecero.fNslice;
	pnew = &new;
	for( i=0; i<36; i++ ) {
        	*pnew = (((*ip & 0xff000000)>>24) |
       		       ((*ip & 0x00ff0000)>>8) |
      		       ((*ip & 0x0000ff00)<<8) |
       		       ((*ip & 0x000000ff)<<24) );
		*ip++ = *pnew;
	}
	p2 = (unsigned char *)& ndou;
	p1 = (unsigned char *)& cabecero.fGeo_matrix[0][0];
	for( k=0; k<9; k++ ){
		for( i=0; i<8; i++) p2[i] = p1[ 7-i ];
		for( i=0; i<8; i++) p1[i] = p2[ i ];
	/**	cabecero.fGeo_matrix[i/3][i%3] = ndou;  **/	
		p1 += 8;
	}

	ip = (unsigned int *)&cabecero.fAngle1;
        *pnew = (((*ip & 0xff000000)>>24) |
       	       ((*ip & 0x00ff0000)>>8) |
      	       ((*ip & 0x0000ff00)<<8) |
       	       ((*ip & 0x000000ff)<<24) );
	*ip = *pnew;
}


void vuelta_imagen( unsigned int ** image, int rows, int cols )
{
	int i, j;
	unsigned int new, *pnew, *pn;
	
	pnew = &new;	
	for ( i=0 ; i<rows;i++) {
		pn = (unsigned int *) image[i];
		for ( j=0; j<cols; j++) {
        		*pnew =  (((*pn & 0xff000000)>>24) |
				((*pn & 0x00ff0000)>>8) |
				((*pn & 0x0000ff00)<<8) |
				((*pn & 0x000000ff)<<24) );
			*pn++ = *pnew;
		}
	}
}

/**
void vuelta2_imagen( char ** image, int rows, int cols, int format )
{
	int i, j;
	char **imagen_alreves;

	if ( (imagen_alreves = imalloc(rows, cols, format) ) == NULL) {
        	puts("\n Error : no puedo almacenar memoria ");
        	puts(" para hacer la conversion al i860 ");
               exit (1);
	}

	for ( i=0 ; i<rows;i++) {
		for ( j=0; j<cols*sizeof(float) ; j+=sizeof(float)) {
			imagen_alreves[i][j]   = image [i][j+3];
			imagen_alreves[i][j+1] = image [i][j+2];
			imagen_alreves[i][j+2] = image [i][j+1];
			imagen_alreves[i][j+3] = image [i][j];
		}
	}
	for ( i=0 ; i<rows;i++) {
		 for ( j=0; j<cols*sizeof(float) ; j++)
			image[i][j] = imagen_alreves[i][j];
        }

	imfree( imagen_alreves, rows, cols,format);
}
**/
#endif
/***************************************************************
   This function read/write the image "image" from/to the file
 called "name", depending on the value of the flag "rdwr".
 The image's dimensions are "row x col" and its kind is indicated
 in "format", that is, the image could be a matrix of chars
 ("format" values NATURAL) or a matrix of floats ("format" values
 FLOATFMT or SPIDER).
   If you want translate and rotate the image using its header
 information, you need another matrix (with the same format and
 dimensions) to load this "changes", and it'll be the parameter
 "image_geo".
   But if you don't want aplicate this geometric information to
 the image, or you have a NATURAL image,  you must call this
 function passing the NULL value to this parameter, as:

       image_io ( myimage, NULL, row,col, filename,.......)

   The float images ( format values FLOATFMT) could be normalized
 after reading or before writing. If you want use such flavour
 you must call this funtion setting '1' or 'TRUE' in the parameter
 "norma", and '0' or 'FALSE' in other case.
   We recommend use normalized images, that is, with average=0 and
 standard desviation = 1.
 ***************************************************************/

int image_io (char **image,char **image_geo,int row,int col,char *name,
              int rdwr,int format,int norma)

                    /* If it is not a natural format image, a cast has      */
                    /* been previously done by the calling routine          */
    /*  name       :   File name                                            */
    /* rdwr, format:   Read/write, kind of image                            */
    /* norma       :   Flag for normalize                                   */

{   int i, j, k;                    /* Counters                             */
    int fichero;                    /* File handle                          */
    int oflag, pmode = 0;           /* Open options (see library reference) */
    int size;                       /* Size of element to read/write        */
    int n_block;                    /* No. of rows that can fit in a block  */
    int n_rest;                     /* No. of rows in the last short block  */
    int no_blocks;                  /* No. of blocks to read/write          */
    int rest;                       /* Rest of rows to read                 */
    int headrec;                    /* more spider header fun               */
    int GEO=0;                      /* Flag for geometric information       */
    int NORMwr=0;                   /* Flags for normalized images          */
    int NORMrd=0;                   /* (one for writing , one for reading)  */
    char * temp_buffer;             /* Temporal storage buffer to speed up  */
    char * aux_ptr;                 /* Auxiliary pointer                    */
    unsigned long io_size;          /* Size of the block to read/write      */
    long tot_size;                  /* Total no. of bytes to read/write     */
    float angle;
    double matriz_geo[3][3];

/********************* Check that the user is honest ************************/

if (row <= 0 || col <= 0)
    return ERROR;   /* Allowed range: 1..? */

/****************** Assign appropriate values to I/O parameters *************/

switch (rdwr) {
    case READING:
	oflag = O_RDONLY;   /* read */
        break;
    case WRITING:
	oflag = O_WRONLY | O_TRUNC | O_CREAT;  /* write      */
#ifndef _PARAMID
	pmode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;/*Read-write permission*/
#endif
        break;
    default:
        return ERROR;                        /* Operation not allowed        */
}

/***************** Assign value to the size flag ***************************/

switch (format) {
    case NATURAL:
        size = sizeof (BYTE);
        break;
    case INTFMT:
        size = sizeof (UWORD);
        break;
    case LONGFMT:
        size = sizeof (ULONG);
        break;
    case FLOATFMT:
        size = sizeof (float);
        if(image_geo!=NULL  &&  rdwr==READING) GEO=1;
        if(norma==TRUE)
           if(rdwr==READING) NORMrd=1;
           else NORMwr=1;
        break;
    case FOURIER:
        size = sizeof (float);
        col += 2;               /* Column index ranges from (0..N/2) * 2   */
        break;
    default:
        return ERROR;            /* Unrecognized format                    */
}

size *= col;                     /* No. of bytes to read in each row       */

/*************** Allocate temp. buffer & adjust reading parameters *********/

n_block = 32000/size;  /* no. of rows that fit in 32 K (roughly 32000)     */
if (n_block == 0)
    return ERROR;    /* The user's pulling the routine's leg, too much     */
io_size = n_block * size;                    /* This is the block i/o size */

while ((temp_buffer = (char *) malloc (io_size)) == NULL)
{   io_size -=size;              /* Not enough memory, reduce requeriments */
    if (io_size <= 0)
        return ERROR;                         /* No memory at all, goodbye */
}

n_block = io_size / size;                   /* No. of rows per block       */
tot_size = (long) size * (long) row;        /* Total no. of bytes to read  */
no_blocks = tot_size / io_size;             /* no. of 32K blocks to read   */
rest = tot_size % io_size;                  /* no. of bytes that are left  */
n_rest = rest / size;                       /* no. of rows that are left   */

/*********************** Open file, checking errors ************************/

if ((fichero = open (name, oflag , pmode)) == -1)
{   free (temp_buffer);
    return ERROR;           /* File not found error */
}

/* THE HEADER*/
if ((format == FLOATFMT) || (format == FOURIER) )
  {
      i= sizeof(CABECERO);
      cabecero.fNlabel = (float)((int)(256/col+1));
      cabecero.fLabrec= (float)ceil((float)256/col);

      if(rdwr == WRITING)
         {
              cabecero.fNslice =1.0;
              cabecero.fNsam = (float) col;
              cabecero.fNrow = (float) row;
              headrec = (int) 1024 / ((int)cabecero.fNsam * 4);
              if (format == FOURIER)
                  cabecero.fIform = (-1);
              else
                  cabecero.fIform = 1;
              if(  (1024%(int)cabecero.fNsam !=0))
                    {
                    cabecero.fNrec= cabecero.fNrow+1;
                    headrec = headrec + 1;
                    }
              else
                    cabecero.fNrec=cabecero.fNrow;

              /*cabecero.fLabbyt = (float) headrec * cabecero.fNsam * 4 ; */
              cabecero.fLabbyt = cabecero.fNsam*cabecero.fLabrec*4;
              cabecero.fLenbyt = (float) cabecero.fNsam * 4;
              cabecero.fPhi = -cabecero.fAngle1;
              Tiempo(); /*calculate time & date and store them in the header*/
              i=(int)cabecero.fNsam*(int)cabecero.fLabrec*4;
         }

	if( NORMwr )  normalize_io((float **)image,rdwr,name);
#ifdef _PARAMID
	if(rdwr == WRITING) vuelta_cabecero( );
#endif
         switch (rdwr) {
             case READING:
	         if ( read(fichero, (char *)&cabecero,i) != i) {
		         close (fichero);           /* Perform I/O                 */
		         return ERROR;             /* EOF prematurelly reached    */
	         }
                 break;
             case WRITING:
	         if ( write(fichero, (char *)&cabecero,i) != i) {
		         close (fichero);           /* Perform I/O                 */
		         return ERROR;             /* EOF prematurelly reached    */
	         }
                 break;
         }

	if(rdwr == READING) {
#ifdef _PARAMID
		vuelta_cabecero( );
#endif
		lseek(fichero,col*(int)cabecero.fLabrec*4,SEEK_SET);
	/******************************************************************/
	/*     cambiar numero de planos de negativo a positivo            */
	/******************************************************************/
		cabecero.fNslice = fabs(cabecero.fNslice);
	}
  }

/********************** Read/write file/memory *****************************/
#ifdef _PARAMID
if ( rdwr == WRITING ) vuelta_imagen( (unsigned int **)image, row, col );
#endif
i = 0;
for (j = 0; j < no_blocks; j++)          /* Read all the blocks of 32000  */
{
    aux_ptr = temp_buffer;               /* Reset aux. pointer            */
    if (rdwr == WRITING)                 /* Copy rows to temp. buffer     */
        for (k = 0; k < n_block; k++, i++, aux_ptr += size)
            memcpy (aux_ptr, image [i], size);

   switch (rdwr) {
       case READING:
          if (read(fichero, temp_buffer, io_size) != io_size)
          {   close (fichero);
              free (temp_buffer);
	      return ERROR;
          }
           break;
       case WRITING:
          if (write(fichero, temp_buffer, io_size) != io_size)
          {   close (fichero);
              free (temp_buffer);
	      return ERROR;
          }
           break;
   }


    if (rdwr == READING)                 /* Copy temp. buffer to rows     */
        for (k = 0; k < n_block; k++, i++, aux_ptr += size)
            memcpy (image [i], aux_ptr, size);
}
aux_ptr = temp_buffer;                   /* Reset aux. pointer            */
if (rdwr == WRITING)
    for (k = 0; k < n_rest; k++, i++, aux_ptr += size)
        memcpy (aux_ptr, image[i], size);/* Copy rows to aux. memory      */

switch (rdwr) {
    case READING:
         if (read(fichero, temp_buffer, rest) != rest)
         {   close (fichero);
             free (temp_buffer);
             return ERROR;
         }
        break;
    case WRITING:
         if (write(fichero, temp_buffer, rest) != rest)
         {   close (fichero);
             free (temp_buffer);
             return ERROR;
         }
        break;
}

if (rdwr == READING)                     /* Copy temp. buffer to rows     */
    for (k = 0; k < n_rest; k++, i++, aux_ptr += size)
        memcpy (image [i], aux_ptr, size);
#ifdef _PARAMID
if ( rdwr == READING ) vuelta_imagen( (unsigned int **)image, row, col );
#endif
if( NORMrd ) normalize_io((float **)image,rdwr,name);

if( GEO) {
    header_geo(matriz_geo,&angle,READING);
    if( !IsInfogeo(matriz_geo) )
     {
        identMatrix(matriz_geo);
        angle=0.0;
        header_geo(matriz_geo,&angle,WRITING);
     }

    transforma ((float **)image,(float **)image_geo,row,col,matriz_geo);
}

/************* Close file, free temp. buffer & return OK ******************/
free (temp_buffer);
close (fichero);
return OK;
}

/**************************************************************************/
/* Image allocation routine. Returns a pointer to pointers to individually*/
/* allocated rows of the image. (row, col) are the image dimensions.      */
/* Format indicates the format of the resultant image.                    */
/* Returns: the pointer to the pointers if everything OK, NULL if there   */
/* is not enough memory. In this latter case, it leaves everything as was */
/* before the call.                                                       */
/* Version 2.0: allocates blocks of rows to improve speed                 */
/* Version 3.0: modified for OS/2 1.0 ==> halloc used instead of malloc   */
/*    just because malloc does not give more then 4Mb, and halloc does.   */
/*    should work in any environment just changing halloc for malloc      */
/**************************************************************************/

void **imalloc (int row, int col, int format)


{   int i, j, k;                    /* Counters                             */
    unsigned element_size;          /* Size of element to allocate          */
    unsigned pointer_size;          /* Id. of pointers                      */
    char **temp;                    /* Temporal value to work with          */
    unsigned no_blocks;             /* No. of ALLOC_SIZE blocks to allocate */
    long tot_size;                  /* Total allocation size                */
    unsigned row_alloc;             /* No. of rows to alloc at the same time*/
    unsigned all_size;              /* No. of bits to alloc at one time     */
    unsigned rest_size;             /* Rest of bytes to alloc               */
    unsigned rest_rows;             /* Rest of rows to alloc.               */
    char *aux;                      /* Aux. pointer                         */


/******************* Assign appropriate value to size flag ******************/

if (format == FOURIER)
    col += 2;      /* Special treatment for FFT format (see foutrans.c) */

switch (format) {
    case NATURAL:
        element_size = sizeof (BYTE);
        pointer_size = sizeof (BYTE *);
        break;
    case INTFMT:
        element_size = sizeof (UWORD);
        pointer_size = sizeof (UWORD *);
        break;
    case LONGFMT:
        element_size = sizeof (ULONG);
        pointer_size = sizeof (ULONG *);
        break;
    case FLOATFMT:
    case FOURIER:
        element_size = sizeof (float);
        pointer_size = sizeof (float *);
        break;
    default:
        return NULL;
}

row_alloc = ALLOC_SIZE/(col*element_size);  /* No. of rows to alloc */
all_size = row_alloc*col*element_size;
tot_size = ((long) element_size) * ((long) row) * ((long) col);
no_blocks = tot_size/all_size;
rest_size = tot_size - ((long) no_blocks) * ((long) all_size);
rest_rows = rest_size / (col*element_size);

/********************* Allocate base pointer ********************************/

if ((temp = (char **) malloc (row*pointer_size)) == NULL)
    return NULL;               /* Not even this little bit of memory */

/*********************** Allocate most blocks *******************************/

j = 0;
for (i = 0; i < no_blocks; i++)
{   if ((aux = (char *)malloc ((long)all_size)) == NULL)
    {   for (j = 0; j < i; j++)
	    free (temp[j*row_alloc]);
        free ((char *) temp);
        return NULL;
    }
    for (k = 0; k < row_alloc; k++, j++)
        temp [j] = aux + k*col*element_size;
}

/*************************** Alloc the last block ***************************/

if (rest_size != 0)
{   if ((aux = (char *)malloc ((long)rest_size)) == NULL)
    {   for (j = 0; j < no_blocks; j++)
	    free (temp[j*row_alloc]);
        free (temp);
        return NULL;
    }
    for (k = 0; k < rest_rows; k++, j++)
        temp [j] = aux + k*col*element_size;
}

/************************* return OK pointer value  **********************/

return (void **)temp;
}

/**************************************************************************/

/**************************************************************************/
/* This function frees an image previously allocated with image_alloc.    */
/* hfree used instead of free (see imalloc header). Change hfree to free  */
/* for portability                                                        */
/**************************************************************************/

void imfree (char **image, int row,int  col, int format)

{   int i;                          /* Counters                             */
    unsigned element_size;          /* Size of element to allocate          */
    unsigned pointer_size;          /* Id. of pointers                      */
    unsigned no_blocks;             /* No. of ALLOC_SIZE blocks to allocate */
    long tot_size;                  /* Total allocation size                */
    unsigned row_alloc;             /* No. of rows to alloc at the same time*/
    unsigned all_size;              /* No. of bits to alloc at one time     */

if (image == NULL)                  /* No allocation at the moment */
    return;

if (format == FOURIER)
    col += 2;      /* Special treatment for FFT format (see foutrans.c) */

switch (format) {
    case NATURAL:
        element_size = sizeof (BYTE);
        pointer_size = sizeof (BYTE *);
        break;
    case INTFMT:
        element_size = sizeof (UWORD);
        pointer_size = sizeof (UWORD *);
        break;
    case LONGFMT:
        element_size = sizeof (ULONG);
        pointer_size = sizeof (ULONG *);
        break;
    case FLOATFMT:
    case FOURIER:
        element_size = sizeof (float);
        pointer_size = sizeof (float *);
        break;
    default:
        return;
}

row_alloc = ALLOC_SIZE/(col*element_size);  /* No. of rows to free  */
all_size = row_alloc*col*element_size;
tot_size = ((long) element_size) * ((long) row) * ((long) col);
no_blocks = tot_size/all_size;

if (image == NULL)  /* No allocation at the moment */
    return;

/*************************** Free most blocks *******************************/

for (i = 0; i < no_blocks; i++)
    if (image [i*row_alloc] != NULL)
	free (image [i*row_alloc]);

/*************************** Free the last block ****************************/

if (image [i*row_alloc] != NULL)
    free (image [i*row_alloc]);

free (image);

}

/****************************************************************************/
/* Indicates wether a particular FILE exists or not                         */
/****************************************************************************/

int exists (char *filenam)

{   FILE *aux;

if ((aux = fopen (filenam, "r")) == NULL)
    return FALSE;
fclose (aux);
return TRUE;
}

/****************************************************************************/
/*                                                                          */
/****************************************************************************/

void asigna_extension (char *extension, int codigo)

{
switch (codigo) {
    case DESCARTADA:
        strcpy (extension, ".");
        return;
    case CORTADA:
        strcpy (extension, ".img");
        return;
    case PRECENTRADA:
        strcpy (extension, ".pre");
        return;
}
if (codigo < DESCARTADA)
{   strcpy (extension, ".");
    return;
}

if (codigo%2 == 0)              /**** Es una centrada ****/
    sprintf (extension, ".c%02d", codigo/2);
else                            /**** Es una girada ******/
    sprintf (extension, ".g%02d", codigo/2);
}


/***************************************************************
 This function finds the four local maximums of an image.
 ***************************************************************/
void busca_maximos (float **imagen, int fil, int col,
               float *maxi, int *imax, int *jmax, int nmax)

{   int i,j,k,l,m,im,jm;
    float maximo, aux;
    int continua;
    int ini_f, fin_f, ini_c, fin_c;

maximo = -1e38;
for (i=0; i < fil; i++)
    for (j=0; j < col; j++)
        if (imagen[i][j] > maximo)
        {   maximo = imagen[i][j];
            im = i;
            jm = j;
        }
maxi[0] = maximo;
imax[0] = im;
jmax[0] = jm;

for (k=1; k < nmax; k++)
{   maximo = -1e38;
    for (i=0; i < fil; i++)
        for (j=0; j < col; j++)
        if ((aux = imagen[i][j]) > maximo && aux < maxi[k-1])
            /******* No a los m ximos anteriores ********/
        {   continua = TRUE;
            if (i == 0)         /**** Para no salirnos de madre ****/
                ini_f = 0;
            else
                ini_f = -1;
            if (i == fil -1)
                fin_f =  0;
            else
                fin_f = 1;
            if (j == 0)
                ini_c = 0;
            else
                ini_c = -1;
            if (j == col -1)
                fin_c = 0;
            else
                fin_c = 1;
            for (l=ini_f; l <= fin_f; l++)
                for (m= ini_c; m <= fin_c; m++)
                    if (l != 0 || m != 0) /*** No mires el central ***/
                        if (aux <= imagen[i+l][j+m])
                        {   continua = FALSE;
                            break;
                        }
            if (continua)
            {   maximo = aux;
                im = i;
                jm = j;
            }
        }
    imax[k] = im;
    jmax[k] = jm;
    maxi[k] = maximo;
}
}

/*************************************************************************/
/*    Puts the data and time in the image header                         */
/*************************************************************************/

void Tiempo(void)
{
     time_t lTiempo;
     struct tm *tmTiempoGreng;


     time(&lTiempo);
     tmTiempoGreng=localtime(&lTiempo);
     tmTiempoGreng->tm_mon++;

sprintf(cabecero.szITim,"%d%s%d",tmTiempoGreng->tm_hour,":",tmTiempoGreng->tm_min);
sprintf(cabecero.szIDat,"%d%s%d%s%d",tmTiempoGreng->tm_mday,"-",tmTiempoGreng->tm_mon,"-",tmTiempoGreng->tm_year);
}

/**************************************************************************/
/* Fills the image header "cabecero"                                      */
/**************************************************************************/
void Cabecera(void)
{
float fSalida;

puts("CALCULO DE LA CABECERA");

if((fSalida=ScanfMejorFloat())!=-9999)
   cabecero.fNslice=fabs(fSalida);

printf("\n\n\nNUMBER OF SLICES IN VOLUMEN:[%f]",fabs(cabecero.fNslice));
if((fSalida=ScanfMejorFloat())!=-9999)
   cabecero.fNslice=fabs(fSalida);
printf ("\nDate=%f",cabecero.fNslice);

printf("\n\n\nNUMBER OF ROWS PER SLICE (X-AXIS):[%f]",(cabecero.fNrow));
if((fSalida=ScanfMejorFloat())!=-9999)
   cabecero.fNrow=fSalida;
printf ("\nDate=%f",cabecero.fNrow);

printf("\n\n\nNUMBER OF PIXELS PER LINE (Y-AXIS):[%f]",(cabecero.fNsam));
if((fSalida=ScanfMejorFloat())!=-9999)
   cabecero.fNsam=fSalida;
printf ("\nDate=%f",cabecero.fNsam);

puts("\n\nFILE TYPE.- \n\t\t\t +3 FOR A 3-D FILE");
puts("\t\t\t +1 FOR A 2-D IMAGE");
puts("\t\t\t -1 FOR A 2-D FOURIER TRANSFORM");
puts("\t\t\t -3 FOR A 3-D FOURIER TRANSFORM");
puts("\t\t\t -5 FOR A NEW 2-D FOURIER TRANSFORM");
puts("\t\t\t -7 FOR A NEW 3-D FOURIER TRANSFORM");
puts("\t\t\t 66 'B' FOR A BIESPECTRUN PLANE");
/*puts("\t\t\+11 FOR A 2-D EIGHT BIT COLOR IMAGE FILE");*/
printf("\n\n\nFILE TYPE:[%f]",(cabecero.fIform));

if((fSalida=ScanfMejorFloat())!=-9999)
   cabecero.fIform=fSalida;
printf ("\nDate=%f",cabecero.fIform);

printf("\n\n\nTITLE:[%s]",cabecero.szITit);
gets(cabecero.szITit);
/*****Now some calculations to keep Spider Compability***/

cabecero.fNlabel = (float)((int)(256/cabecero.fNsam+1));

cabecero.fLabrec=ceil(256/cabecero.fNsam);
	
cabecero.fSig = -1;
cabecero.fImami =0;
}

/**************************************************************************/
/*               Scanf mejorado, admite entrada nula y float              */
/**************************************************************************/
float ScanfMejorFloat (void)
{
char szNumero[16];
gets(szNumero);
if (szNumero[0]=='\0')
   return (-9999);
else
   return(atof(szNumero));
}



