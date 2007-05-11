/***************************************************************************
 *
 * Authors:     DEbora gil
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
/*****************************************************************************/
/* APH Files: MRC                                                            */
/*****************************************************************************/

#ifndef _XVMRC_H
   #define _XVMRC_H

#include <data/funcs.h>
#include <data/matrix2d.h>
#include <data/matrix3d.h>
#include <data/geometry.h>
#include <data/image.h>
#include <data/volume.h>

#include <sys/stat.h>

/**@name cp4 Files */
//@{
/** CCP$ Files.
    This is a class to read/write ccp4 files as produced by MRC rograms
    */

#define MRC_LABEL_SIZE      80
#define MRC_USER            28
#define MRC_NUM_LABELS      10

/*
 *  The different modes supported by the MRC format.
 */
#define MODE_BYTE           0
#define MODE_SHORT          1
#define MODE_FLOAT          2
#define MODE_SHORT_COMPLEX  3
#define MODE_FLOAT_COMPLEX  4

/*
 * Axis.
 */
#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3

/*
 *	This is the official size of the mrc header.
 *	C requires the '\0' terminator for strings which would
 *	offset the size by MRC_NUM_LABELS bytes.
 */
#define SIZEOF_MRC_HEADER 1024
/*******************************************************************************
 *
 *  STRUCTURE   :   struct MRCheader / MRC_HEADER
 *
 *  PURPOSE     :   This is the MRC header as defined in file mrc.doc
 */
typedef struct MRCheader
{
    int         nx;                 /* Number of columns -      */
                                    /* Fastest changing in map. */
    int         ny;

    int         nz;                 /* Number of sections -     */
                                    /* slowest changing in map. */

    int         mode;               /* See modes above. */

    int         nxstart;            /* No. of first column in map   default 0.*/
    int         nystart;            /* No. of first row in map      default 0.*/
    int         nzstart;            /* No. of first section in map  default 0.*/

    int         mx;                 /* Number of intervals along X. */
    int         my;                 /* Number of intervals along Y. */
    int         mz;                 /* Number of intervals along Z. */

    float       xlen;               /* Cell dimensions (Angstroms). */
    float       ylen;               /* Cell dimensions (Angstroms). */
    float       zlen;               /* Cell dimensions (Angstroms). */

    float       alpha;              /* Cell angles (Degrees). */
    float       beta;               /* Cell angles (Degrees). */
    float       gamma;              /* Cell angles (Degrees). */

    int         mapc;               /* Which axis corresponds to Columns.  */
                                    /* (1,2,3 for X,Y,Z.                   */
    int         mapr;               /* Which axis corresponds to Rows.     */
                                    /* (1,2,3 for X,Y,Z.                   */
    int         maps;               /* Which axis corresponds to Sections. */
                                    /* (1,2,3 for X,Y,Z.                   */

    float       amin;               /* Minimum density value. */
    float       amax;               /* Maximum density value. */
    float       amean;              /* Mean density value.    */

    int         ispg;               /* Space group number (0 for images). */

    int         nsymbt;             /* Number of bytes used for storing   */
                                    /* symmetry operators. */
    unsigned long  user[MRC_USER];  /* For user - all set to zero by default. */
    char        map[4];              /* Char string 'MAP ' to identify file type */
    char        machst[4];           /* Machine stamp indicating the machine type
			               which wrote file
                                       The machine stamp is a 32-bit
                                       quantity containing a set of four
                                       `nibbles' (half-bytes)---only half
                                       the space is used. Each nibble is
                                       a number specifying the representation
                                       of (in C terms) double (d) , float (f),
                                       int (i) and unsigned char (c) types.
                                        Thus each stamp is of the form
                                        0xdfic0000. The values for the
                                        floating point nibbles may be taken
                                        from the list (following HDF):
                                       1 	Big-endian ieee
                                       2 	VAX
                                       3 	Cray
                                       4 	Little-endian ieee
                                       5 	Convex native
                                       6 	Fijitsu VP
                                       */


    float       xorigin;            /* X origin. */
    float       yorigin;            /* Y origin. */

    int         nlabl;              /* Number of labels being used. */

                                    /* 10 text labels of 80 characters each. */
    char        labels[MRC_NUM_LABELS][MRC_LABEL_SIZE+1];

} MRC_HEADER;

class CCP4 {
public:
   /** Estructure with ccp4 header
        copied from ccp4 code */
   MRC_HEADER my_mrc_header;


public:
/** write a ccp4 2D file*/
   void write(const FileName &fn_out, const ImageXmipp &I, bool reversed=false,
                                 double x_length=0., double y_length=0., double z_length=0.);

/** write a ccp4 3D file*/
   void write(const FileName &fn_out, const VolumeXmipp &I, bool reversed=false,
                                 double x_length=0., double y_length=0., double z_length=0.);

/** read  a ccp4 2D file*/
   void read (const FileName &fn_out,  ImageXmipp &I, bool reversed=false);

/** read  a ccp4 3D file*/
   void read (const FileName &fn_out,  VolumeXmipp &V, bool reversed=false);

/** Empties actual my_mrc_header structure. */
   void clear();

/** Fill mrc header from 2D xmipp image. */
   void fill_header_from_xmippimage(ImageXmipp I, bool reversed=false,
                                 double x_length=0., double y_length=0., double z_length=0.);

/** Fill mrc header from 3D xmipp image. */
   void fill_header_from_xmippvolume(VolumeXmipp V, bool reversed=false,
                                 double x_length=0., double y_length=0., double z_length=0.);

/** Fill mrc header from mrc file.
    Returns tre is either reversed is true or the native endianess is
    different from the one used to store the data.

    */
   bool read_header_from_file(const FileName &fn_in, bool reversed=false);
};



//@}


#endif
