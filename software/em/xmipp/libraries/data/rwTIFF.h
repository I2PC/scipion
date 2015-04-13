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

#ifndef RWTIFF_H_
#define RWTIFF_H_

///@defgroup TIFF TIFF File format
///@ingroup ImageFormats
//@{

/** TIFF Data Header
*/
struct TIFFDirHead
{                                   // Header for each Directory in TIFF
    unsigned short  bitsPerSample;
    unsigned short  samplesPerPixel;
    unsigned int   imageWidth;
    unsigned int   imageLength;
    uint16           imageSampleFormat;
    unsigned short  resUnit;
    float            xTiffRes,yTiffRes;
    unsigned int subFileType;
    uint16 pNumber, pTotal; // pagenumber and total number of pages of current directory
    TIFFDirHead()
    {
        	bitsPerSample=samplesPerPixel=0;
        	imageWidth=imageLength=subFileType=0;
            imageSampleFormat=0;
            xTiffRes=yTiffRes=0;
    }
};

/** castTiffTile2T
  * write the content of a tile from a TIFF file to an image array
  */
void castTiffTile2T(
    size_t offset ,
    char* tif_buf,
    unsigned int x, unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned int tileWidth, unsigned int tileLength,
    unsigned short samplesPerPixel,
    DataType datatype);

/** castTiffLine2T
  * write the content of a line from a TIFF file to an image array
  */
void castTiffLine2T(
    size_t offset,
    char* tif_buf,
    unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned short samplesPerPixel,
    DataType datatype);

/** Determine datatype of the TIFF format file.
  * @ingroup TIFF
  */
DataType datatypeTIFF(TIFFDirHead dHead);

/** Read TIFF format files.
  */
int readTIFF(size_t select_img, bool isStack=false);

/** Write TIFF format files.
  */
int writeTIFF(size_t select_img, bool isStack=false, int mode=WRITE_OVERWRITE, String bitDepth="", CastWriteMode castMode = CW_CAST);
//@}
#endif /* RWTIFF_H_ */
