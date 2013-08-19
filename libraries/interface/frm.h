/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
/*****************************************************************************/
/* INTERACTION WITH Fast Rotational Matching                                 */
/*****************************************************************************/

#ifndef _XMIPP_FRM_HH
#define _XMIPP_FRM_HH

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <data/multidim_array.h>

/**@defgroup FRMInterface Fast Rotational Matching
   @ingroup InterfaceLibrary */
//@{
/// Initialize Python to be used from C++
void initializeXmippPython(String &xmippPython);

/// Convert from MultidimArray to numpy array
PyObject* convertToNumpy(const MultidimArray<double> &I);

/** Get pointer to FRM function.
 * This is done to avoid losing time in importing the module each time.
 * Remind to free it when you do not need it any longer with
 * Py_DECREF(pFunc);
 */
PyObject * getPointerToPythonFRMFunction();

/** Align two volumes using FRM.
 * The first argument is the pointer to the FRM python function. You may obtain it with
 * getPointerToPythonFRMFunction()
 */
void alignVolumesFRM(PyObject *pFunc, const MultidimArray<double> &Iref, const MultidimArray<double> &I,
		double &rot, double &tilt, double &psi, double &x, double &y, double &z, double &score);

//@}
#endif
