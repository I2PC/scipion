#include <Python.h>

#include <data/filename.h>

/***************************************************************/
/*                            FileName                         */
/**************************************************************/

/*FileName Object*/
typedef struct
{
    PyObject_HEAD
    FileName * filename;
}
FileNameObject;

/* Destructor */
static void
FileName_dealloc(FileNameObject* self)
{
    delete self->filename;
    self->ob_type->tp_free((PyObject*)self);
}

/* Constructor */
static PyObject *
FileName_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*)type->tp_alloc(type, 0);

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = "", *ext="";
        int number = -1;
        if (PyArg_ParseTuple(args, "|Ois", &input, &number, &ext)
            || PyArg_ParseTuple(args, "|Os", &input, &ext))
        {
            pyStr = PyObject_Str(input);
            if (pyStr != NULL)
                str = PyString_AsString(pyStr);
        }

        self->filename = new FileName(str, number, ext);
    }
    return (PyObject *)self;
}

/* String representation */
static PyObject *
FileName_repr(PyObject * obj)
{
    FileNameObject *self = (FileNameObject*)obj;
    return PyString_FromString(self->filename->c_str());
}

/* Initialization */
static int
FileName_init(FileNameObject *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

/* compose */
static PyObject *
FileName_compose(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*)obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = "", *ext="";
        int number = -1;
        if (PyArg_ParseTuple(args, "Ois", &input, &number, &ext))
        {
            pyStr = PyObject_Str(input);
            if (pyStr != NULL)
                str = PyString_AsString(pyStr);
            self->filename->compose(str, number, ext);
        }
        else if (PyArg_ParseTuple(args, "iO", &number, &input))
        {
          pyStr = PyObject_Str(input);
          if (pyStr != NULL)
              str = PyString_AsString(pyStr);
            self->filename->compose(number, str);
        }
        else
          return NULL;
    }
    return Py_BuildValue("");//Return None(similar to void in C)
}

/* isInStack */
static PyObject *
FileName_isInStack(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  FileNameObject *self = (FileNameObject*)obj;

  if (self->filename->isInStack())
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

/* isInStack */
static PyObject *
FileName_isMetaData(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  FileNameObject *self = (FileNameObject*)obj;

  if (self->filename->isMetaData(false))
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

static PyObject *
FileName_getExtension(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  FileNameObject *self = (FileNameObject*)obj;

  return PyString_FromString(self->filename->getExtension().c_str());
}

static PyObject *
FileName_getNumber(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  FileNameObject *self = (FileNameObject*)obj;

  return PyInt_FromLong(self->filename->getNumber());
}

static PyObject *
FileName_getBaseName(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  FileNameObject *self = (FileNameObject*)obj;

  return PyString_FromString(self->filename->getBaseName().c_str());
}

static PyObject *
FileName_decompose(PyObject *obj, PyObject *args, PyObject *kwargs)
{
  FileNameObject *self = (FileNameObject*)obj;
  int no;
  String str;
  self->filename->decompose(no, str);
  return Py_BuildValue("is", no, str.c_str());
}

/* FileName methods */
static PyMethodDef FileName_methods[] = {
    {"compose", (PyCFunction)FileName_compose, METH_VARARGS,
     "Compose from root, number and extension OR prefix with number @"
    },
    {"isInStack", (PyCFunction)FileName_isInStack, METH_NOARGS,
     "True if filename has stack format"
    },
    {"isMetaData", (PyCFunction)FileName_isMetaData, METH_NOARGS,
     "True if is a MetaData"
    },
    {"getExtension", (PyCFunction)FileName_getExtension, METH_NOARGS,
     "Get the last extension from a FileName"
    },
    {"getNumber", (PyCFunction)FileName_getNumber, METH_NOARGS,
     "Get the number from a FileName"
    },
    {"getBaseName", (PyCFunction)FileName_getBaseName, METH_NOARGS,
     "Get the base name from a FileName"
    },
    {"decompose", (PyCFunction)FileName_decompose, METH_NOARGS,
     "Decompose filenames with @. Mainly from selfiles"
    },
    {NULL}  /* Sentinel */
};

/*FileName Type */
static PyTypeObject FileNameType = {
   PyObject_HEAD_INIT(NULL)
   0,                         /*ob_size*/
   "xmipp.FileName",          /*tp_name*/
   sizeof(FileNameObject),   /*tp_basicsize*/
   0,                         /*tp_itemsize*/
   (destructor)FileName_dealloc, /*tp_dealloc*/
   0,                         /*tp_print*/
   0,                         /*tp_getattr*/
   0,                         /*tp_setattr*/
   0,                         /*tp_compare*/
   FileName_repr,             /*tp_repr*/
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
   "Python wrapper to Xmipp FileName class",/* tp_doc */
   0,                     /* tp_traverse */
   0,                     /* tp_clear */
   0,                     /* tp_richcompare */
   0,                     /* tp_weaklistoffset */
   0,                     /* tp_iter */
   0,                     /* tp_iternext */
   FileName_methods,  /* tp_methods */
   0,                      /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)FileName_init,      /* tp_init */
   0,                         /* tp_alloc */
   FileName_new,                 /* tp_new */
};


static PyMethodDef xmipp_methods[] =
{
   {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
initxmipp(void)
{
    PyObject* m;

    m = Py_InitModule3("xmipp", xmipp_methods,
                       "Example module that creates an extension type.");

    //FileNameType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&FileNameType) < 0)
        return;

    Py_INCREF(&FileNameType);
    PyModule_AddObject(m, "FileName", (PyObject *)&FileNameType);
}
