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

#ifndef RWPIF_H_
#define RWPIF_H_

///defgroup PIF Portable Image Format for EM Data File format
///@ingroup ImageFormats

/** Information obtained from:
 *  http://cryoem.ucsd.edu/programDocs/pif.pdf
 */


typedef struct
{
    // One time file header - 60 bytes -->> 512 total
    long int file_id;       //> BYTE*8 Magic bytes - First 4 bytes ignored.
    char realScaleFactor[16];  //> char[16] How to convert floatInt*4 into real*4
    int numImages;          //> int*4 Number of images in this file
    int endianNess;
    /**  int*4 What endian type of machine created this file
    0 = VAX/DECstation
     (Little)
    1 = SGI etc... (Big) */
    char genProgram[32];   //> Program used to generate this file with version number
    int htype;
    /**  int*4- Possible values:
      1 - all projections have the same number of pixels and the same depth (number of bits per pixel)
      0 - otherwise */
    int nx;      //> int*4 Number of columns
    int ny;      //> int*4 Number of rows
    int nz;      //> int*4 Number of sections
    int mode;
    /**  int*4 EM data type e.g.
         0 = byte*1 - Zeiss scans
         1 = int*2
         2 = floatInt*4
         3 = complex int*2
         4 = complex floatInt*4
         5 = Structure Factors
         6 = Boxed data unfloated with background value placed past radius of box.
         7 = floatInt*2
         8 = complex floatInt*2
         9 = float*4
         10 = complex float*4
         20 = MAP floatInt*2
         21 = MAP floatInt*4
         22 = MAP floatInt*4 PFTS rot*4 dimension
         31 = Structure Factors Amp/Phase floatInt*4
         32 = Structure Factors Apart/Bpart floatInt*4
         88 = Accumulated TIF's in int*2 (st2pif)
         97 = DEPTHCUED etc... */
    int futureUse[107]; //> 428 Bytes Save some space for use later.
}
PIFMainHeader;

typedef struct
{
    int nx;      //> int*4 Number of columns
    int ny;      //> int*4 Number of rows
    int nz;      //> int*4 Number of sections
    int mode;
    /**  int*4 EM data type e.g.
         0 = byte*1 - Zeiss scans
         1 = int*2
         2 = floatInt*4
         3 = complex int*2
         4 = complex floatInt*4
         5 = Structure Factors
         6 = Boxed data unfloated with background value placed past radius of box.
         7 = floatInt*2
    8 = complex floatInt*2
         9 = float*4
         10 = complex float*4
         20 = MAP floatInt*2
         21 = MAP floatInt*4
         22 = MAP floatInt*4 PFTS rot*4 dimension
         31 = Structure Factors Amp/Phase floatInt*4
         32 = Structure Factors Apart/Bpart floatInt*4
         88 = Accumulated TIF's in int*2 (st2pif)
         97 = DEPTHCUED etc... */
    int bkgnd; //> int*4 Background value
    int packRadius; //> int*4 Radius of boxed image
    int nxstart; //> int*4 Number of first col in map
    int nystart; //>  int*4 Number of first row in map
    int nzstart; //>  int*4 Number of first section in map
    int mx; //> int*4 Number of intervals along x
    int my; //> int*4 Number of intervals along y
    int mz; //> int*4 Number of intervals along z
    float xlength; //> floatInt*4 Cell Dimensions (Angstroms)
    float ylength; //> floatInt*4           "
    float zlength; //>  floatInt*4              "
    float alpha; //>  floatInt*4 Cell Angles (degrees)
    float beta; //> floatInt*4          "
    float gamma; //> floatInt*4         "
    int mapc; //> int*4 Which axis is col(1,2,3 =x,y,z)
    int mapr; //> int*4 Which axis is row(1,2,3 =x,y,z)
    int maps; //> int*4 Which axis is sections(1,2,3 = x,y,z)
    float min; //> floatInt*4 Min density value
    float max; //> floatInt*4 Max density value
    float mean; //> floatInt*4 Mean density value
    float stdDev; //> floatInt*4 StdDev of GrayLevels
    int ispg; //> int*4 Space group number
    int nsymbt; //> int*4 Number of bytes for symetry ops
    float xorigin; //> floatInt*4 x origin
    float yorigin; //> floatInt*4 y origin
    char titleDescription[80]; //> User defined description
    char timeStamp[32]; //> Date/time data last modified
    char microGraphDesignation[16]; //> Unique Micrograph Number
    char scanNumber[8]; //> Scan Number of Micrograph
    float aoverb; //> floatInt*4    AOVERB
    float map_abang; //> floatInt*4    MAP_ABANG
    float dela; //>              floatInt*4   DELA
    float delb; //>                floatInt*4    DELB
    float delc; //>               floatInt*4    DELC
    int t_matrix[6]; //>  t_matrix (array of 6 ints)
    float dthe; //> floatInt*4 dthe
    float dphi_90; //> floatInt*4 dphi_90
    float symmetry; //> floatInt*4 symmetry
    int binFactor; //> int*4        Image compression factor
    float a_star; //> floatInt*4 emsf3dbt/emmap3dt stuff
    float b_star; //> floatInt*4
    float c_star; //> floatInt*4
    float alp_star; //> floatInt*4
    float bet_star; //> floatInt*4
    float gam_star; //> floatInt*4
    float pixelSize; //> floatInt*4 From em3dr
    int futureUse[40]; //> int*4 (160 bytes) Save some space for use later
}
PIFDataHeader;

/** Read image from PIF file format
 *
 * @param select_img Index number of selected image
 * @return
 */
int readPIF(size_t select_img);

/** Write image to PIF file format
 *
 * @param select_img Number of selected image to write
 * @param isStack  Force to write as stack if possible
 * @param mode Write file type: overwrite, replace, append....
 * @return
 */
int writePIF(size_t select_img = ALL_IMAGES, bool isStack=false, int mode=WRITE_OVERWRITE);


#endif /* RWPIF_H_ */
