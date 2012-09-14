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

#ifndef IMAGE_MACROS_H_
#define IMAGE_MACROS_H_

//@addtogroup Images

//@{
// Macros used to select slice in readPreview
#define CENTRAL_SLICE -1
#define ALL_SLICES 0
#define SLICE_INDEX(select_slice) ((select_slice == ALL_SLICES) ? 0 : select_slice - 1)
/** Macro for represent ALL IMAGES on stack index */
#define ALL_IMAGES 0
#define FIRST_IMAGE 1
#define IMG_INDEX(select_img) ((select_img == ALL_IMAGES) ? 0 : select_img - 1)
/** Macro for appending an image when doing mapFile2Write */
#define APPEND_IMAGE 0
//@}

#endif /* IMAGE_MACROS_H_ */
