/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#ifndef RWEM_H_
#define RWEM_H_

///defgroup EM EM File format
///@ingroup ImageFormats

/** Information obtained from:
 *  http://www.biochem.mpg.de/doc_tom/TOM_Release_2008/IOfun/tom_isemfile.html
 */

typedef struct
{
    char machine;
    /**  Machine:    Value:
          OS-9         0
          VAX          1
          Convex       2
          SGI          3
          Mac          5
          PC           6 */
    char gPurpose;  //> General purpose. On OS-9 system: 0 old version 1 is new version
    char unused;   //> Not used in standard EM-format, if this byte is 1 the header is abandoned.
    char datatype;
    /** Data TypeCoding:
     *  Image Type:  No. of Bytes:   Value:
          byte            1           1
          short           2           2
          long int        4           4
          float           4           5
          complex         8           8
          double          8           9 */
    //>  Three long integers (3x4 bytes) are image size in x, y, z Dimension
    int xdim;
    int ydim;
    int zdim;
    char comment[80]; //> 80 Characters as comment
    long int params[40]; //> 40 long integers (4 x 40 bytes) are user defined parameters
    /**  The parameters are coded as follwing:
    No.  |  Name  |  Value  |  Factor  |  Comment
    1       U        Volt      1000       accelerating voltage
    2       COE      ???m        1000       Cs of objective lense
    3       APE      mrad      1000       aperture
    4       VE       x         1          end magnification
    5       VN       -         1000       postmagnification of CCD
    6       ET       s         1000       exposure time in seconds
    7       XG       -         1          pixelsize in object-plane
    8       DG       nm        1000       EM-Code:
    EM420=1;
    CM12=2;
    CM200=3;
    CM120/Biofilter=4;
    CM300=5;
    CM300/Tecnai=6;
    extern=0;
    9       APD      nm        1000       photometer aperture
    10      L        nm        1000       phys_pixel_size * nr_of_pixels
    11      DF       Angstr.   1          defocus, underfocus is neg.
    12      FA       Angstr.   1          astigmatism
    13      PHI      deg/1000  1000       angle of astigmatism
    14      DR       Angstr.   1          drift in Angstr.
    15      DELT     deg/1000  1000       direction of drift
    16      DDF      Angstr.   1          focusincr. for focus-series
    17      X0       -         1          obsolete
    18      Y0       -         1          obsolete
    19      KW       deg/1000  1000       tiltangle
    20      KR       deg/1000  1000       axis perpend. to tiltaxis
    21      -        Angstr.   1
    22      SC       ASCII     1
    23      -        -         -
    24      -        pixel     1          internal: subframe X0
    25      -        pixel     1          internal: subframe Y0
    26      -        Angstr.   1000       internal: resolution
    27      -        -         -          internal: density
    28      -        -         -          internal: contrast
    29      -        -         -          internal: unknown
    30      SP       -         1000       mass centre X
    31      SP       -         1000       mass centre Y
    32      SP       -         1000       mass centre Z
    33      H        -         1000       height
    34      -        -         1000       internal: unknown
    35      D1       -         1000       width 'Dreistrahlbereich'
    36      D2       -         1000       width 'Achrom. Ring'
    37      -        -         1          internal: lambda
    38      -        -         1          internal: delta theta
    39      -        -         1          internal: unknown
    40      -        -         1          internal: unknown   */
    char userdata[256]; //>  256 Byte with userdata, i.e. the username
}
EMHead;

/** Read image from EM file format
 *
 * @param select_img Index number of selected image
 * @return
 */
int readEM(size_t select_img);

/** Write image to EM file format
 *
 * @param select_img Number of selected image to write
 * @param isStack  Force to write as stack if possible
 * @param mode Write file type: overwrite, replace, append....
 * @return
 */
int writeEM(size_t select_img = ALL_IMAGES, bool isStack=false, int mode=WRITE_OVERWRITE);


#endif /* RWEM_H_ */
