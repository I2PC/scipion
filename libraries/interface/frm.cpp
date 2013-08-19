/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "frm.h"

void initializeXmippPython(String &xmippPython)
{
	String xmippHome=getenv("XMIPP_HOME");
	if (xmippHome=="")
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot access to $XMIPP_HOME");
	String xmippPythonDir=xmippHome+"/external/python/Python-2.7.2";
	xmippPython=xmippPythonDir+"/python";
	Py_SetProgramName((char *)xmippPython.c_str());
	Py_Initialize();
	import_array(); // For working with numpy
}

PyObject* convertToNumpy(const MultidimArray<double> &I)
{
	npy_intp dim[3];
	dim[0]=ZSIZE(I);
	dim[1]=YSIZE(I);
	dim[2]=XSIZE(I);
	PyObject* pyI=PyArray_SimpleNewFromData(3, dim, NPY_DOUBLE, (void *)MULTIDIM_ARRAY(I));
	return pyI;
}

PyObject * getPointerToPythonFRMFunction()
{
	PyObject * pName = PyString_FromString("sh_alignment.frm"); // Import sh_alignment.frm
	PyObject * pModule = PyImport_Import(pName);
	Py_DECREF(pName);
	PyObject * pFunc = PyObject_GetAttrString(pModule, "frm_align");
	return pFunc;
}

void alignVolumesFRM(PyObject *pFunc, const MultidimArray<double> &I1, const MultidimArray<double> &I2,
		double &rot, double &tilt, double &psi, double &x, double &y, double &z, double &score)
{
	PyObject *pyI1=convertToNumpy(I1);
	PyObject *pyI2=convertToNumpy(I2);

	// Call frm
	int bandWidthSphericalHarmonics0=4;
	int bandWidthSphericalHarmonicsF=64;
	int frequency=20; // Pixels
	int maxshift=10;
	PyObject *arglistfrm = Py_BuildValue("OOOO(ii)i", pyI1, Py_None, pyI2, Py_None,
			bandWidthSphericalHarmonics0, bandWidthSphericalHarmonicsF, frequency, maxshift);
	PyObject *resultfrm = PyObject_CallObject(pFunc, arglistfrm);
	Py_DECREF(arglistfrm);
	if (resultfrm!=NULL)
	{
		PyObject *shift=PyTuple_GetItem(resultfrm,0);
		PyObject *euler=PyTuple_GetItem(resultfrm,1);
		score=PyFloat_AsDouble(PyTuple_GetItem(resultfrm,2));
		x=PyFloat_AsDouble(PyList_GetItem(shift,0));
		y=PyFloat_AsDouble(PyList_GetItem(shift,1));
		z=PyFloat_AsDouble(PyList_GetItem(shift,2));
		rot=PyFloat_AsDouble(PyList_GetItem(euler,0));
		tilt=PyFloat_AsDouble(PyList_GetItem(euler,1));
		psi=PyFloat_AsDouble(PyList_GetItem(euler,2));
		Py_DECREF(resultfrm);
	}
}
