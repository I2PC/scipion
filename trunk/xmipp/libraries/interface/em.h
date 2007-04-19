/***************************************************************************
 *
 * Authors:     Sjors Scheres
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

#include <data/funcs.h>
#include <data/matrix3d.h>
#include <data/geometry.h>
#include <data/volume.h>

#include <bilib/configs.h>

#include <sys/stat.h>

/**@name EM Files */
//@{
/** EM Files.
    This is a class to read/write EM-format volume files as produced
    by the EM-toolbox
*/

#define EM_COMMENT_SIZE    80
#define EM_USER_INTS       40
#define EM_USER_DATA       256

/*
 *  The different modes supported by the EM format.
 */
#define MODE_BYTE           1
#define MODE_SHORT          2
#define MODE_LONG_INT       4
#define MODE_FLOAT          5
#define MODE_COMPLEX        8
#define MODE_DOUBLE         9

/*
 * Axis.
 */
#define X_AXIS 1 /* fastest */
#define Y_AXIS 2
#define Z_AXIS 3

/*
 *	This is the official size of the EM header.
 */
#define SIZEOF_EM_HEADER 512
/*******************************************************************************
 *
 *  STRUCTURE   :   struct EM header / EM_HEADER
 *
 */
typedef struct EMheader
{
  unsigned char   machine;                  /* Machine coding  */
  unsigned char   general_use;              /* General purpose */
  unsigned char   not_used;                 /* Not used */
  unsigned char   mode;                     /* See modes above. */

  unsigned long   nx;                       /* Number of columns X (fastest) */
  unsigned long   ny;                       /* Number of columns Y  */
  unsigned long   nz;                       /* Number of columns Z (slowest) */

  unsigned char   comment[EM_COMMENT_SIZE]; /* User comment */
  unsigned long   user[EM_USER_INTS];       /* For user - all set to zero by default. */
  unsigned char   user_data[EM_USER_DATA];  /* User data  */


} EM_HEADER;

class EM {
public:
   /** Estructure with EM EM_header */
   EM_HEADER my_em_header;

public:

/** write a EM 3D file*/
   void write(const FileName &fn_out, const VolumeXmipp &I, bool reversed=false);

/** read a EM 3D file*/
   void read (const FileName &fn_out,  VolumeXmipp &V, bool reversed=false);

/** Empties actual my_em_header structure. */
   void clear();

/** Fill EM header from 3D xmipp image. */
   void fill_header_from_xmippvolume(VolumeXmipp V, bool reversed=false);

/** Fill EM header from EM file.
    Returns true if either reversed is true or the native endianess is
    different from the one used to store the data.
    */
   bool read_header_from_file(const FileName &fn_in, bool reversed=false);
};

//@}

