/*******************************************************************************/
/* tiff2raw.c																											  08/08/2000 */
/*******************************************************************************/
/* TIFF to raw converter.																											 */
/* The TIFF file should be in the tile oriented format                         */
/* It will make two files: the raw file and another one (filename.inf) with    */
/* the TIFF image size, the bits per sample and the samples per pixel          */
/*******************************************************************************/
/* Centro Nacional de Biotecnología                                            */
/* Carlos Manzanares Sancho                                                    */
/*******************************************************************************/
/*
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "tiffio.h"

/*
// initFile
//
// Creates a file with the specified size
// (to speed up the process it makes a 100000 bytes buffer)
//
*/
int initFile( char *filename, size_t len ) {
	unsigned char* 	buffer = NULL;
	size_t         	i;
	FILE *					fd;

	buffer = (unsigned char*)malloc( sizeof(unsigned char) * 100000 );
        for (int i=0; i<100000; i++) buffer[i]=0;
	if ( buffer == NULL ) return( -1 );

	fd = fopen( filename, "w" );

  for (i=0; i< len/100000; i++)
      fwrite( buffer, sizeof(unsigned char), 100000, fd );
	fwrite( buffer, sizeof(unsigned char), len%100000 , fd );

	fclose( fd );

	return( 0 );
}

/*
//
// rawWriteTile
//
// write the content of a tile from a TIFF file to an image array
//
*/
void rawWriteTile( unsigned char* raw_buf, unsigned char* tif_buf, 
    		  unsigned int x, unsigned int y,
		  unsigned int imageWidth, unsigned int imageLength, 
		  unsigned int tileWidth, unsigned int tileLength,
                  unsigned int *minval, unsigned int *maxval) {
	unsigned int 	i, j;
	unsigned int 	x_max = x + tileWidth,
			y_max = y + tileLength;

	if ( x_max > imageWidth )  x_max = imageWidth;
	if ( y_max > imageLength ) y_max = imageLength;

	for( j = y; j < y_max; j++ )
		for( i = x; i < x_max; i++) {
                        unsigned char val=
                           ((unsigned char*)tif_buf)[(j-y)*tileWidth+(i-x)];
                       /* printf("%d %d %d\n",val,minval,maxval);*/
                        if (val<*minval) *minval=val;
                        if (val>*maxval) *maxval=val;
			raw_buf[j*imageWidth + i] = val;
                }
}

/*
//
// main
//
// usage: tiff2raw inputfile outputfile
// inputfile is the TIFF filename
// outputfile is the RAW filename
// (outputfile.inf will be the properties file)
//
*/
int main(int argc, char *argv[])
{
  TIFF*          	tif     = NULL;
	unsigned char*	tif_buf = NULL;
	int							raw_fd;
	unsigned char*  raw_buf = NULL;
	FILE*						inf     = NULL;

  unsigned short	bitsPerSample;
  unsigned short	samplesPerPixel;
  unsigned int		imageWidth;
  unsigned int		imageLength;
	unsigned int		tileWidth;
	unsigned int		tileLength;
	size_t          len;

	unsigned int		x, y;
        unsigned int    minval, maxval;

	if ( argc < 3 ) {
		fprintf( stderr, "Usage: tiff2raw inputfile outputfile\n" );
		exit(-1);
  }

         minval=255;
         maxval=0;

	/* Open TIFF image */
	tif = TIFFOpen( argv[1], "r" );
	if ( tif == NULL ) exit( -1 );

	/* Get TIFF image properties */
	TIFFGetField( tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample );
	TIFFGetField( tif, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel );
	TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &imageWidth );
	TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &imageLength );
	TIFFGetField( tif, TIFFTAG_TILEWIDTH, &tileWidth );
	TIFFGetField( tif, TIFFTAG_TILELENGTH, &tileLength );
        
	/* Calculates the TIFF image size */
	len = (imageLength * imageWidth * bitsPerSample * samplesPerPixel)/8;
        initFile( argv[2], len );

	/* Mapping from file to memory (it will save primary memory) */
	raw_fd  = open( argv[2], O_RDWR );
	raw_buf = (unsigned char*)mmap( NULL, len, PROT_READ|PROT_WRITE, MAP_SHARED, 
																	raw_fd, 0 );

	/* Writes the TIFF image properties to the .inf file */
	inf     = fopen( strcat( argv[2], ".inf" ), "w" );
	fprintf( inf, "# Bits per sample\n" );
	fprintf( inf, "bitspersample= %d\n", bitsPerSample );
	fprintf( inf, "# Samples per pixel\n" );
	fprintf( inf, "samplesperpixel= %d\n", samplesPerPixel );
	fprintf( inf, "# Image width\n" );
	fprintf( inf, "Xdim= %d\n", imageWidth );
	fprintf( inf, "# Image length\n" );
	fprintf( inf, "Ydim= %d\n", imageLength );

	/* Start to convert the TIFF image to the RAW file */
	tif_buf = (unsigned char*)_TIFFmalloc( TIFFTileSize(tif) );
        
	for( y = 0; y < imageLength; y +=tileLength )
		for( x = 0; x < imageWidth; x += tileWidth ) {
			TIFFReadTile( tif, tif_buf, x, y, 0, 0 );
			rawWriteTile( raw_buf, tif_buf, x, y, 
			   imageWidth, imageLength,
			   tileWidth, tileLength, &minval, &maxval );
		}

	_TIFFfree( tif_buf );
	munmap( (char *) raw_buf, len );

	TIFFClose( tif );
	close( raw_fd );

	fprintf( inf, "# Minimum value=%d\n",minval );
	fprintf( inf, "# Maximum value=%d\n",maxval );
	fclose( inf );
}

/* Colimate menu =========================================================== */
/*Colimate:
   PROGRAM Tiff2Raw {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Tiff2Raw/Help/tiff2raw.html";
      help="Convert TIFF micrographs into raw images suitable for marking";
      OPEN MENU Tiff2Raw;
      COMMAND LINES {
         + usual: xmipp_tiff2raw $FILEIN $FILEOUT
      }
      PARAMETER DEFINITIONS {
         $FILEIN {
            label="Input TIFF micrograph";
            type=file existing;
         }
         $FILEOUT {
            label="Output RAW micrograph";
            type=file;
         }
      }
   }
   MENU Tiff2Raw {
      "I/O parameters"
      $FILEIN
      $FILEOUT
   }
*/
