#include <Python.h>

#include <data/filename.h>

static PyMethodDef xmipp_methods[] = {
    {NULL}  /* Sentinel */
};

/*************************** FileName ***********************/
/*FileName Object*/
typedef struct {
    PyObject_HEAD
    FileName * filename;
}FileNameObject;

/*FileName functions*/
static void
FileName_dealloc(FileNameObject* self)
{
    delete self->filename;
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
FileName_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    FileNameObject *self = (FileNameObject*)type->tp_alloc(type, 0);

    if (self != NULL)
    {
      self->filename = new FileName();
    }

    return (PyObject *)self;
}

static int
FileName_init(FileNameObject *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

/*FileName Type */
static PyTypeObject FileNameType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "xmipp.FileName",          /*tp_name*/
    sizeof(FileNameObject),   /*tp_basicsize*/
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

    FileNameType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&FileNameType) < 0)
        return;

    m = Py_InitModule3("xmipp", xmipp_methods,
                       "Example module that creates an extension type.");

    Py_INCREF(&FileNameType);
    PyModule_AddObject(m, "FileName", (PyObject *)&FileNameType);
}
