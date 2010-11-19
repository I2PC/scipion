#include <Python.h>

#include "../filename.h"

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
} xmipp_NoddyObject;

static PyTypeObject xmipp_NoddyType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "xmipp.Noddy",             /*tp_name*/
    sizeof(xmipp_NoddyObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Noddy objects",           /* tp_doc */
};

static PyMethodDef xmipp_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

/* Trying to create the FileName type */
typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    FileName * filename;
} xmipp_FileNameObject;

static PyTypeObject xmipp_FileNameType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "xmipp.FileName",          /*tp_name*/
    sizeof(xmipp_FileNameObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Wrapper of Xmipp FileName class",/* tp_doc */
};

PyMODINIT_FUNC
initxmipp(void)
{
    PyObject* m;

    xmipp_NoddyType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&xmipp_NoddyType) < 0)
        return;
    xmipp_FileNameType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&xmipp_FileNameType) < 0)
        return;

    m = Py_InitModule3("xmipp", xmipp_methods,
                       "Example module that creates an extension type.");

    Py_INCREF(&xmipp_NoddyType);
    PyModule_AddObject(m, "Noddy", (PyObject *)&xmipp_NoddyType);

    Py_INCREF(&xmipp_FileNameType);
    PyModule_AddObject(m, "FileName", (PyObject *)&xmipp_FileNameType);
}
