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

#ifndef XMIPP_HDF5_H_
#define XMIPP_HDF5_H_

#include "../../external/hdf5-1.8.10/src/hdf5.h"
#include "../../external/hdf5-1.8.10/c++/src/H5Cpp.h"
#include <iostream>


/** @defgroup Tools for General Purpose handling of hdf5 files
 *  @ingroup DataLibrary
 *
 *  @{
 */


class XmippH5File: public H5::H5File
{


public:

    void showTree(std::ostream &out = std::cout);


};

herr_t showObjectInfo(hid_t group, const char *name, void *op_data);


/** @}
 */
#endif /* XMIPP_HDF5_H_ */
