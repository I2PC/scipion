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
/*                            Program                         */
/**************************************************************/

/* Program methods */
PyMethodDef Program_methods[] =
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
PyTypeObject ProgramType =
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

/* Destructor */
void Program_dealloc(ProgramObject* self)
{
    delete self->program;
    self->ob_type->tp_free((PyObject*) self);
}

/* Constructor */
PyObject *
Program_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        self->program = new XmippProgramGeneric();
        PyObject * runWithoutArgs = Py_False;
        if (PyArg_ParseTuple(args, "|O", &runWithoutArgs))
        {
            try
            {
                if (PyBool_Check(runWithoutArgs))
                    self->program->runWithoutArgs = (runWithoutArgs == Py_True);
                else
                    PyErr_SetString(PyExc_TypeError, "MetaData::new: Expecting boolean value");
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return (PyObject *) self;
}

/* addUsageLine */
PyObject *
Program_addUsageLine(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    PyObject *verbatim = Py_False;
    if (self != NULL)
    {
        char * line = NULL;
        if (PyArg_ParseTuple(args, "s|O", &line, &verbatim))
        {
            try
            {
                if (PyBool_Check(verbatim))
                {
                    self->program->addUsageLine(line, verbatim == Py_True);
                    Py_RETURN_NONE;
                }
                else
                    PyErr_SetString(PyExc_TypeError, "Program::addUsageLine: Expecting boolean value");
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}

/* addExampleLine */
PyObject *
Program_addExampleLine(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    PyObject *verbatim = Py_True;
    if (self != NULL)
    {
        char * line = NULL;
        if (PyArg_ParseTuple(args, "s|O", &line, &verbatim))
        {
            try
            {
                if (PyBool_Check(verbatim))
                {
                    self->program->addExampleLine(line, verbatim == Py_True);
                    Py_RETURN_NONE;
                }
                else
                    PyErr_SetString(PyExc_TypeError, "Program::addExampleLine: Expecting boolean value");
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}

/* addParamsLine */
PyObject *
Program_addParamsLine(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        char * line = NULL;
        if (PyArg_ParseTuple(args, "s", &line))
        {
            try
            {
                self->program->addParamsLine(line);
                Py_RETURN_NONE;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}

/* usage */
PyObject *
Program_usage(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        int verbose = 0;
        if (PyArg_ParseTuple(args, "|i", &verbose))
        {
            try
            {
                self->program->usage(verbose);
                Py_RETURN_NONE;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}

/* endDefinition */
PyObject *
Program_endDefinition(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        try
        {
            self->program->endDefinition();
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* read */
PyObject *
Program_read(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        PyObject *list = NULL;
        if (PyArg_ParseTuple(args, "O", &list))
        {
            if (PyList_Check(list))
            {
                size_t size = PyList_Size(list);
                PyObject * item = NULL;
                char ** argv = new char*[size];
                std::vector<double> vValue(size);
                for (size_t i = 0; i < size; ++i)
                {
                    item = PyList_GetItem(list, i);
                    if (!PyString_Check(item))
                    {
                        PyErr_SetString(PyExc_TypeError,
                                        "Program arguments should be of type string");
                        return NULL;
                    }

                    argv[i] = PyString_AsString(item);

                    if (i == 0)
                    {
                        FileName temp(argv[i]);
                        temp = temp.removeDirectories();
                        argv[i] = strdup(temp.c_str());
                    }
                }
                self->program->read((int)size, (const char **)argv);
                if (self->program->doRun)
                    Py_RETURN_TRUE;
                else
                    Py_RETURN_FALSE;
            }
        }
    }
    return NULL;
}

/* checkParam */
PyObject *
Program_checkParam(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        char * param = NULL;
        if (PyArg_ParseTuple(args, "s", &param))
        {
            try
            {
                if (self->program->checkParam(param))
                    Py_RETURN_TRUE;
                else
                    Py_RETURN_FALSE;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}

/* getParam */
PyObject *
Program_getParam(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        char * param = NULL;
        int arg = 0;
        if (PyArg_ParseTuple(args, "s|i", &param, &arg))
        {
            try
            {
                const char * value = self->program->getParam(param, arg);
                return PyString_FromString(value);
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}

/* getListParam */
PyObject *
Program_getListParam(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        char * param = NULL;
        if (PyArg_ParseTuple(args, "s", &param))
        {
            try
            {
                StringVector list;
                self->program->getListParam(param, list);
                size_t size = list.size();
                PyObject * pylist = PyList_New(size);

                for (size_t i = 0; i < size; ++i)
                    PyList_SetItem(pylist, i, PyString_FromString(list[i].c_str()));
                return pylist;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
            }
        }
    }
    return NULL;
}
