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

#ifndef RWSPE_H_
#define RWSPE_H_

///@defgroup SPE Princeton Instruments File Format
///@ingroup ImageFormats

// I/O prototypes
/** SPE Reader
  * @ingroup SPE
*/
int readSPE(size_t select_img,bool isStack=false);

/** SPE Writer
  * @ingroup SPE
*/
int writeSPE(size_t select_img, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(ERR_IMG_NOWRITE, "writeSPE is not implemented.");
    return(-1);
}
#endif /* RWSPE_H_ */
