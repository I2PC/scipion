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

#ifndef RWDM3_H_
#define RWDM3_H_

///@defgroup DM3 DM3 File format
///@ingroup ImageFormats

/** DM3 Reader
  * @ingroup DM3
*/
int readDM3(size_t img_select,bool isStack=false);

/** DM3 Writer
  * @ingroup DM3
*/
int writeDM3(size_t img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(ERR_IO_NOWRITE, "ERROR: writeDM3 is not implemented.");
    return(-1);
}

#endif /* RWDM3_H_ */
