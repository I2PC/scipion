/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Roberto Marabini       (roberto@cnb.csic.es)
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

#ifndef _PYTHON_PROGRAM_H
#define _PYTHON_PROGRAM_H

#include "Python.h"

/***************************************************************/
/*                            Program                         */
/**************************************************************/

#define Program_Check(v) (((v)->ob_type == &ProgramType))
#define Program_Value(v) ((*((ProgramObject*)(v))->program))

/*Program Object*/
typedef struct
{
    PyObject_HEAD
    XmippProgramGeneric * program;
}
ProgramObject;

#define ProgramObject_New() (ProgramObject*)malloc(sizeof(ProgramObject))

/* Destructor */
void Program_dealloc(ProgramObject* self);

/* Constructor */
PyObject *
Program_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

/* addUsageLine */
PyObject *
Program_addUsageLine(PyObject *obj, PyObject *args, PyObject *kwargs);

/* addExampleLine */
PyObject *
Program_addExampleLine(PyObject *obj, PyObject *args, PyObject *kwargs);

/* addParamsLine */
PyObject *
Program_addParamsLine(PyObject *obj, PyObject *args, PyObject *kwargs);

/* usage */
PyObject *
Program_usage(PyObject *obj, PyObject *args, PyObject *kwargs);

/* endDefinition */
PyObject *
Program_endDefinition(PyObject *obj, PyObject *args, PyObject *kwargs);

/* read */
PyObject *
Program_read(PyObject *obj, PyObject *args, PyObject *kwargs);

/* checkParam */
PyObject *
Program_checkParam(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getParam */
PyObject *
Program_getParam(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getListParam */
PyObject *
Program_getListParam(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Program methods */
static PyMethodDef Program_methods[] =
{
    { "addUsageLine", (PyCFunction) Program_addUsageLine, METH_VARARGS,
      "Add a line to program usage" },
    { "addExampleLine", (PyCFunction) Program_addExampleLine, METH_VARARGS,
      "Add a line to program examples" },
    { "addParamsLine", (PyCFunction) Program_addParamsLine, METH_VARARGS,
      "Add a line to program params definition" },
    { "usage", (PyCFunction) Program_usage, METH_VARARGS,
      "Print program usage" },
    { "endDefinition", (PyCFunction) Program_endDefinition, METH_VARARGS,
      "End definition of params" },
    { "read", (PyCFunction) Program_read, METH_VARARGS,
      "read arguments" },
    { "checkParam", (PyCFunction) Program_checkParam, METH_VARARGS,
      "Check if a param was passed in arguments" },
    { "getParam", (PyCFunction) Program_getParam, METH_VARARGS,
      "Get the value passed of this param" },
    { "getListParam", (PyCFunction) Program_getListParam, METH_VARARGS,
      "Get the list of all values passed of this param" },
    { NULL } /* Sentinel */
};

/*Program Type */
static PyTypeObject ProgramType =
{
    PyObject_HEAD_INIT(NULL)
    0, /*ob_size*/
    "xmipp.Program", /*tp_name*/
    sizeof(ProgramObject), /*tp_basicsize*/
    0, /*tp_itemsize*/
    (destructor)Program_dealloc, /*tp_dealloc*/
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
    "Python wrapper to Xmipp Program class",/* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    Program_methods, /* tp_methods */
    0, /* tp_members */
    0, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    0, /* tp_init */
    0, /* tp_alloc */
    Program_new, /* tp_new */
};


#endif
