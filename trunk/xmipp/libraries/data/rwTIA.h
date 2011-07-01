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

#ifndef RWTIA_H_
#define RWTIA_H_


///@defgroup TIA TIA File format
///@ingroup ImageFormats

// I/O prototypes
/** TIA Reader
  * @ingroup TIA
*/
int readTIA(int img_select,bool isStack=false);

/** TIA Writer
  * @ingroup TIA
*/
int writeTIA(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(ERR_IMG_NOWRITE, "ERROR: writeTIA is not implemented.");
    return(-1);
}

#endif /* RWTIA_H_ */
