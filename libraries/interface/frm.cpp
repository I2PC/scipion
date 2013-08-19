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
#include <data/geometry.h>

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

void alignVolumesFRM(PyObject *pFunc, const MultidimArray<double> &Iref, const MultidimArray<double> &I,
		double &rot, double &tilt, double &psi, double &x, double &y, double &z, double &score)
{
	PyObject *pyIref=convertToNumpy(Iref);
	PyObject *pyI=convertToNumpy(I);

	// Call frm
	int bandWidthSphericalHarmonics0=4;
	int bandWidthSphericalHarmonicsF=64;
	int frequency=20; // Pixels
	int maxshift=10;
	PyObject *arglistfrm = Py_BuildValue("OOOO(ii)i", pyI, Py_None, pyIref, Py_None,
			bandWidthSphericalHarmonics0, bandWidthSphericalHarmonicsF, frequency, maxshift);
	PyObject *resultfrm = PyObject_CallObject(pFunc, arglistfrm);
	Py_DECREF(arglistfrm);
	Py_DECREF(pyIref);
	Py_DECREF(pyI);
	if (resultfrm!=NULL)
	{
		PyObject *shift=PyTuple_GetItem(resultfrm,0);
		PyObject *euler=PyTuple_GetItem(resultfrm,1);
		score=PyFloat_AsDouble(PyTuple_GetItem(resultfrm,2));
		x=PyFloat_AsDouble(PyList_GetItem(shift,0))-XSIZE(Iref)/2;
		y=PyFloat_AsDouble(PyList_GetItem(shift,1))-YSIZE(Iref)/2;
		z=PyFloat_AsDouble(PyList_GetItem(shift,2))-ZSIZE(Iref)/2;
		double angz1=PyFloat_AsDouble(PyList_GetItem(euler,0));
		double angz2=PyFloat_AsDouble(PyList_GetItem(euler,1));
		double angx=PyFloat_AsDouble(PyList_GetItem(euler,2));
		Matrix2D<double> Efrm, E(3,3);
		Euler_anglesZXZ2matrix(-angz1, -angx, -angz2, Efrm); // -angles because FRM rotation definition
		                                                  // is the opposite of Xmipp

		// Reorganize the matrix because the Z and X axes in the coordinate system are reversed with
		// respect to Xmipp
		// Exmipp=[0 0 1; 0 1 0; 1 0 0]*Efrm*[0 0 1; 0 1 0; 1 0 0]
		MAT_ELEM(E,0,0)=MAT_ELEM(Efrm, 2, 2); MAT_ELEM(E,0,1)=MAT_ELEM(Efrm, 2, 1); MAT_ELEM(E,0,2)=MAT_ELEM(Efrm, 2, 0);
		MAT_ELEM(E,1,0)=MAT_ELEM(Efrm, 1, 2); MAT_ELEM(E,1,1)=MAT_ELEM(Efrm, 1, 1); MAT_ELEM(E,1,2)=MAT_ELEM(Efrm, 1, 0);
		MAT_ELEM(E,2,0)=MAT_ELEM(Efrm, 0, 2); MAT_ELEM(E,2,1)=MAT_ELEM(Efrm, 0, 1); MAT_ELEM(E,2,2)=MAT_ELEM(Efrm, 0, 0);
		E=E.inv();
		Euler_matrix2angles(E,rot,tilt,psi);
		Py_DECREF(resultfrm);

		Euler_angles2matrix(0,20,20, E);
	}
}
