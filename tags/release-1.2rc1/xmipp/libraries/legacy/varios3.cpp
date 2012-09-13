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
/* Input/output from a layer to memory. Image is a pointer to pointers to */
/* the rows of the image. (row, col) are the image dimensions, Name is    */
/* the name of the file to read/write from/to and rdwr = 0 ==> reading,   */
/* rdwr == 1 ==> writing. In this latter case, if the file exists, it is  */
/* erased before writing.                                                 */
/* Uses low-level I/O routines to speed up the operations.                */
/* format indicates the kind of file to read/write (see virus.h)          */
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

#include <fcntl.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifndef __APPLE__
#include <malloc.h>
#else
#include <cstdlib>
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "spider.h"
#include "groe.h"

/***************************************************************
   This function read/write a slice of a volume  from/to the
 file handled by "handle", depending on the value of the flag
 "rdwr".
   The slice's dimensions are "row x col" and its kind is indicated
 in "format", that is, the slice could be a matrix of chars
 ("format" values NATURAL) or a matrix of floats ("format" values
 FLOATFMT or SPIDER).
   It's calling from "volume_io" (see below) and it's similar to
 the rutine "image_io" in file "varios.c".
 ***************************************************************/

int capa_io(char **image, int row, int col, int handle, int rdwr, int format)

{
    int i, j, k;                    /* Counters                             */
    int size;                       /* Size of element to read/write        */
    int n_block;                    /* No. of rows that can fit in a block  */
    int n_rest;                     /* No. of rows in the last short block  */
    char *temp_buffer;              /* Temporal storage buffer to speed up  */
    int io_size;                    /* Size of the block to read/write      */
    int no_blocks;                  /* No. of blocks to read/write          */
    int rest;                       /* Rest of rows to read                 */
    long tot_size;                  /* Total no. of bytes to read/write     */
    char *aux_ptr;                  /* Auxiliary pointer                    */


    /********************* Check that the user is honest ************************/

    if (row <= 0 || col <= 0)
        return ERROR;   /* Allowed range: 1..? */

    /***************** Assign value to the size flag ***************************/

    switch (format)
    {
    case NATURAL:
        size = sizeof(BYTE);
        break;
    case INTFMT:
        size = sizeof(UWORD);
        break;
    case LONGFMT:
        size = sizeof(ULONG);
        break;
    case FLOATFMT:
        size = sizeof(float);
        break;
    case FOURIER:
        size = sizeof(float);
        col += 2;               /* Column index ranges from (0..N/2) * 2   */
        break;
    default:
        return ERROR;            /* Unrecognized format                    */
    }

    size *= col;                     /* No. of bytes to read in each row       */

    /*************** Allocate temp. buffer & adjust reading parameters *********/

    n_block = 32000 / size;  /* no. of rows that fit in 32 K (roughly 32000)     */
    if (n_block == 0)
        return ERROR;    /* The user's pulling the routine's leg, too much     */
    io_size = n_block * size;                    /* This is the block i/o size */

    while ((temp_buffer = (char *) malloc(io_size)) == NULL)
    {
        io_size -= size;              /* Not enough memory, reduce requeriments */
        if (io_size <= 0)
            return ERROR;                         /* No memory at all, goodbye */
    }

    n_block = io_size / size;                   /* No. of rows per block       */
    tot_size = (long) size * (long) row;        /* Total no. of bytes to read  */
    no_blocks = tot_size / io_size;             /* no. of 32K blocks to read   */
    rest = tot_size % io_size;                  /* no. of bytes that are left  */
    n_rest = rest / size;                       /* no. of rows that are left   */

    /********************** Read/write file/memory *****************************/

    i = 0;
    for (j = 0; j < no_blocks; j++)          /* Read all the blocks of 32000  */
    {
        aux_ptr = temp_buffer;               /* Reset aux. pointer            */
        if (rdwr == WRITING)                 /* Copy rows to temp. buffer     */
            for (k = 0; k < n_block; k++, i++, aux_ptr += size)
                memcpy(aux_ptr, image [i], size);

        switch (rdwr)
        {
        case READING:
            if (read(handle, temp_buffer, io_size) != io_size)
            {
                free(temp_buffer);
                return ERROR;                    /* EOF prematurelly reached      */
            }
            break;
        case WRITING:
            if (write(handle, temp_buffer, io_size) != io_size)
            {
                free(temp_buffer);
                return ERROR;                    /* EOF prematurelly reached      */
            }
            break;
        }
        if (rdwr == READING)                 /* Copy temp. buffer to rows     */
            for (k = 0; k < n_block; k++, i++, aux_ptr += size)
                memcpy(image [i], aux_ptr, size);
    }
    aux_ptr = temp_buffer;                   /* Reset aux. pointer            */
    if (rdwr == WRITING)
        for (k = 0; k < n_rest; k++, i++, aux_ptr += size)
            memcpy(aux_ptr, image[i], size);/* Copy rows to aux. memory      */
    switch (rdwr)
    {
    case READING:
        if (read(handle, temp_buffer, rest) != rest)
        {
            free(temp_buffer);
            return ERROR;                        /* EOF prematurelly reached      */
        }
        break;
    case WRITING:
        if (write(handle, temp_buffer, rest) != rest)
        {
            free(temp_buffer);
            return ERROR;                        /* EOF prematurelly reached      */
        }
        break;
    }
    if (rdwr == READING)                     /* Copy temp. buffer to rows     */
        for (k = 0; k < n_rest; k++, i++, aux_ptr += size)
            memcpy(image [i], aux_ptr, size);

    /************* Close file, free temp. buffer & return OK ******************/

    free(temp_buffer);
    return OK;
}


/***************************************************************
 This function allocs memory for a volume, using the function
 "imalloc" (in file "varios.c") which allocs memory for an image,
 and in this case, allocs memory for each slice of the volume.
 ***************************************************************/
void ***trialloc(int capas, int row, int col, int format)

{
    int i;                    /* Counters                             */
    void ***temp;

    if (format == FOURIER)         /**** 2 slices more ****/
        capas += 2;

    temp = (void ***) malloc(capas * sizeof(void **));
    for (i = 0; i < capas; i++)
        if ((temp[i] = imalloc(row, col, format)) == NULL)
            return NULL;

    return temp;

}

/***************************************************************
 This function frees the memory allocated by "trialloc"
 (see the previous function)
 ***************************************************************/
void trifree(void ***volumen, int capas, int row, int col, int format)

{
    int i;                    /* Counters                             */

    if (format == FOURIER)         /**** 2 slices more ****/
        capas += 2;

    for (i = 0; i < capas; i++)
        imfree((char **)volumen[i], row, col, format);

    free(volumen);

}

/***************************************************************
   This function read/write the volume "volumen" from/to the file
 called "name", depending on the value of the flag "rdwr".
 The volume's dimensions are "capas x fil x col" and its kind is
 indicated in "format", that is, the volume could be a
 three-dimensional array of chars ("format" values NATURAL)
 or of floats ("format" values FLOATFMT or SPIDER).
   For reading and writing this function calls to "capa_io" which
 is similar to "image_io" of the file "varios.c".
   The float volumes ( format values FLOATFMT) could be normalized
 after reading or before writing. If you want use such flavour
 you must call this funtion setting '1' or 'TRUE' in the parameter
 "norma", and '0' or 'FALSE' in other case.
   We recommend use normalized volumes, that is, with average=0 and
 standard desviation = 1.
 ***************************************************************/

int volume_io(char ***volumen, int capas, int fil, int col, char *name,
              int rdwr, int format, int norma)

{
    int handle, oflag, pmode, i;
    int  headrec;    /* more spider header fun               */
    int NORMwr = 0;                   /* Flags for normalized volume          */
    int NORMrd = 0;                   /* (one for writing , one for reading)  */


    if (capas <= 0 || fil <= 0 || col <= 0)
        return ERROR;   /* Allowed range: 1..? */

    if (format == FOURIER)         /**** 2 slices more ****/
        capas += 2;

    /****************** Assign appropriate values to I/O parameters *************/

    switch (rdwr)
    {
    case READING:
        oflag = O_RDONLY ;                  /* read                         */
        break;
    case WRITING:
        oflag = O_WRONLY | O_TRUNC | O_CREAT;             /* write          */
        pmode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
        /* Read-write permission */
        break;
    default:
        return ERROR;                       /* Operation not allowed        */
    }

    if ((format == FLOATFMT) && (norma == TRUE))
    {
        if (rdwr == READING) NORMrd = 1;
        else NORMwr = 1;
    }


    /*********************** Open file, checking errors ***************************/

    if ((handle = open(name, oflag, pmode)) == -1)
        return ERROR;                                     /* File not found error */

    /********       THE HEADER      *************/

    if ((format == FLOATFMT) || (format == FOURIER))
    {
        i = sizeof(CABECERO);
        cabecero.fNlabel = (float)((int)(256 / col + 1));
        cabecero.fLabrec = (float) ceil((float)256 / col);

        if (rdwr == WRITING)
        {
            cabecero.fNslice = (float) capas;
            cabecero.fNsam = (float) col;
            cabecero.fNrow = (float) fil;
            headrec = (int) 1024 / ((int)cabecero.fNsam * 4);
            if (format == FOURIER)
                cabecero.fIform = (-3);
            else
                cabecero.fIform = 3;
            if ((1024 % (int)cabecero.fNsam != 0))
            {
                cabecero.fNrec = cabecero.fNrow + 1;
                headrec = headrec + 1;
            }
            else
                cabecero.fNrec = cabecero.fNrow;

            cabecero.fLabbyt = (float) headrec * cabecero.fNsam * 4 ;
            cabecero.fLenbyt = (float) cabecero.fNsam * 4;
            i = (int)cabecero.fNsam * (int)cabecero.fLabrec * 4;
            Tiempo(); /*calculate time & date and store them in the header*/
        }

        if (NORMwr) norm_of_volume((float ***)volumen, rdwr, name);

        switch (rdwr)
        {
        case READING:
            if (read(handle, (char*)&cabecero, i) != i)
            {
                close(handle);           /* Perform I/O                 */
                return ERROR;             /* EOF prematurelly reached    */
            }
            break;
        case WRITING:
            if (write(handle, (char*)&cabecero, i) != i)
            {
                close(handle);           /* Perform I/O                 */
                return ERROR;             /* EOF prematurelly reached    */
            }
            break;
        }

        if (rdwr == READING)
        {
            lseek(handle, col*(int)cabecero.fLabrec*4, SEEK_SET);
        }
        /****************************************************************************/
        /*         cambiar numero de planos de negativo a positivo                  */
        /***************************************************************************/
        cabecero.fNslice = fabs(cabecero.fNslice);
    }


    for (i = 0; i < capas; i++)
        if (capa_io(volumen[i], fil, col, handle, rdwr, format) == ERROR)
        {
            close(handle);
            return ERROR;
        }

    if (NORMrd) norm_of_volume((float ***)volumen, rdwr, name);

    close(handle);
    return OK;
}

/*************************************************************************/
/* This routine converts a float volume to a char volume.                */
/*************************************************************************/

void vox_4_a_1(float ***volumen, BYTE ***volbyte, int dim)

{
    int i, j, k;
    float maxi = -1e38, mini = 1e38, aux;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            for (k = 0; k < dim; k++)
            {
                aux = volumen[i][j][k];
                if (aux < mini)
                    mini = aux;
                if (aux > maxi)
                    maxi = aux;
            }

    if (maxi == mini)
        aux = 0;
    else
        aux = 255. / (maxi - mini);

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            for (k = 0; k < dim; k++)
                volbyte[i][j][k] = (unsigned char)((volumen[i][j][k] - mini) * aux + 0.5);

}

/*************************************************************************/
/* This routine converts a char volume to a float volume.                */
/*************************************************************************/

void vox_1_a_4(float ***volumen, BYTE ***volbyte, int dim)

{
    int i, j, k;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            for (k = 0; k < dim; k++)
                volumen[i][j][k] = volbyte[i][j][k];
}

