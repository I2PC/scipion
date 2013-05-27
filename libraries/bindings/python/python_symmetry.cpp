/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include "xmippmodule.h"

/***************************************************************/
/*                            SymList                          */
/***************************************************************/

/* Constructor */
PyObject *
SymList_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    SymListObject *self = (SymListObject*) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->symlist = new SymList();
        //self->symlist->readSymmetryFile("i3");
    }
    return (PyObject *) self;
}

/* Destructor */
 void SymList_dealloc(SymListObject* self)
{
    delete self->symlist;
    self->ob_type->tp_free((PyObject*) self);
}

/* readSymmetryFile */
PyObject *
SymList_readSymmetryFile(PyObject * obj, PyObject *args, PyObject *kwargs)
{
    char * str = NULL;

    if (PyArg_ParseTuple(args, "s", &str))
    {
        try
        {
            SymListObject *self = (SymListObject*) obj;
            self->symlist->readSymmetryFile(str);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/* computeDistance */
PyObject *
SymList_computeDistance(PyObject * obj, PyObject *args, PyObject *kwargs)
{

    PyObject *pyMd = NULL;
    PyObject *pyProjdirMode = Py_False;
    PyObject *pyCheckMirrors = Py_False;
    PyObject *pyObjectRotation = Py_False;

    if (PyArg_ParseTuple(args, "O|OOO", &pyMd,
                         &pyProjdirMode,
                         &pyCheckMirrors,
                         &pyObjectRotation))
    {
        if (!MetaData_Check(pyMd))
            PyErr_SetString(PyExc_TypeError,
                            "Expected MetaData as first argument");
        try
        {
            bool projdir_mode    = false;
            bool check_mirrors   = false;
            bool object_rotation = false;
            if (PyBool_Check(pyProjdirMode))
                projdir_mode = (pyProjdirMode == Py_True);
            if (PyBool_Check(pyCheckMirrors))
                check_mirrors = (pyCheckMirrors == Py_True);
            if (PyBool_Check(pyObjectRotation))
                object_rotation = (pyObjectRotation == Py_True);
            SymListObject *self = (SymListObject*) obj;
            self->symlist->computeDistance(MetaData_Value(pyMd),projdir_mode,
                                                                check_mirrors,
                                                                object_rotation);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* computeDistance */
PyObject *
SymList_computeDistanceAngles(PyObject * obj, PyObject *args, PyObject *kwargs)
{
	double rot1, tilt1, psi1, rot2, tilt2, psi2;
    PyObject *pyProjdirMode = Py_False;
    PyObject *pyCheckMirrors = Py_False;
    PyObject *pyObjectRotation = Py_False;

    if (PyArg_ParseTuple(args, "dddddd|OOO", &rot1, &tilt1, &psi1, &rot2, &tilt2, &psi2,
                         &pyProjdirMode,
                         &pyCheckMirrors,
                         &pyObjectRotation))
    {
        try
        {
            bool projdir_mode    = false;
            bool check_mirrors   = false;
            bool object_rotation = false;
            if (PyBool_Check(pyProjdirMode))
                projdir_mode = (pyProjdirMode == Py_True);
            if (PyBool_Check(pyCheckMirrors))
                check_mirrors = (pyCheckMirrors == Py_True);
            if (PyBool_Check(pyObjectRotation))
                object_rotation = (pyObjectRotation == Py_True);
            SymListObject *self = (SymListObject*) obj;
            double dist=self->symlist->computeDistance(rot1,tilt1,psi1,rot2,tilt2,psi2,
            		projdir_mode,check_mirrors,object_rotation);
            return PyFloat_FromDouble(dist);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* SymList methods */
PyMethodDef SymList_methods[] =
{
   { "readSymmetryFile", (PyCFunction) SymList_readSymmetryFile,
     METH_VARARGS, "read symmetry file" },
   { "computeDistance", (PyCFunction) SymList_computeDistance,
       METH_VARARGS, "compute angular distance in a metadata" },
   { "computeDistanceAngles", (PyCFunction) SymList_computeDistanceAngles,
	   METH_VARARGS, "compute angular distance between two sets of angles" },
   { NULL } /* Sentinel */
};//SymList_methods

/*SymList Type */
PyTypeObject SymListType =
{
    PyObject_HEAD_INIT(NULL)
    0, /*ob_size*/
    "xmipp.SymList", /*tp_name*/
    sizeof(SymListObject), /*tp_basicsize*/
    0, /*tp_itemsize*/
    (destructor)SymList_dealloc, /*tp_dealloc*/
    0, /*tp_print*/
    0, /*tp_getattr*/
    0, /*tp_setattr*/
    0, /*tp_compare*/
    0, /*tp_repr*/
    0, /*tp_as_number*/
    0, /*tp_as_sequence*/
    0, /*tp_as_mapping*/
    0, /*tp_hash */
    0, /*tp_call*/
    0, /*tp_str*/
    0, /*tp_getattro*/
    0, /*tp_setattro*/
    0, /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT, /*tp_flags*/
    "Python wrapper to Xmipp SymList class",/* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    SymList_methods, /* tp_methods */
    0, /* tp_members */
    0, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    0, /* tp_init */
    0, /* tp_alloc */
    SymList_new, /* tp_new */
};//SymListType

