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
#include "matrix1d.h"
#include <iostream>



/** @defgroup Tools for General Purpose handling of hdf5 files
 *  @ingroup DataLibrary
 *
 *  @{
 */


class XmippH5File: public H5::H5File
{

public:

	/**
	 * Show the groups and dataset and print them in the output stream out.
	 * @param out Output stream
	 */
	void showTree(std::ostream &out = std::cout);

    /**
     * Open HDF5 file
     * @param name	File name
     * @param flags	tandard hdf5 flags
     * @param access_plist Standard hdf5 plist
     */
    void openFile(const H5std_string& name, unsigned int flags,
                  const H5::FileAccPropList& access_plist = H5::FileAccPropList::DEFAULT);

    /** Return the values in the dataset dsname and return them in a Matrix1D double data
     *
     * @param dsname Dataset name
     * @param data	 Vector data
     * @param reportError If true throw an exception in case of failure, otherwise it returns
     * a negative number
     */
    int getDataset(const char* dsname, Matrix1D<double> &data, bool reportError = true) const;
};

herr_t showObjectInfo(hid_t group, const char *name, void *op_data);


/** @}
 */
#endif /* XMIPP_HDF5_H_ */
