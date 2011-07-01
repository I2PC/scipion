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

#ifndef RWRAW_H_
#define RWRAW_H_

///@defgroup RAW1 RAW Data type
///@ingroup ImageFormats

// I/O prototypes
/** RAW Reader
  * @ingroup RAW1
*/

DataType datatypeRAW(String strDT);

///@defgroup RAW RAW File format
///@ingroup ImageFormats

// I/O prototypes
/** RAW Reader
  * @ingroup RAW1
*/

int readRAW(size_t select_img,bool isStack=false);

#endif /* RWINF_H_ */
