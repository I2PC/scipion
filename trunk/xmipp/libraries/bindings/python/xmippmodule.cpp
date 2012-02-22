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

#include <external/python/Python-2.7.2/Include/Python.h>

#include <data/xmipp_program.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image_generic.h>
#include <data/xmipp_image_extension.h>
#include <data/xmipp_color.h>
#include <reconstruction/ctf_estimate_from_micrograph.h>

static PyObject * PyXmippError;

#define FileName_Check(v)  (((v)->ob_type == &FileNameType))
#define FileName_Value(v)  ((*((FileNameObject*)(v))->filename))

#define MetaData_Check(v)  (((v)->ob_type == &MetaDataType))
#define MetaData_Value(v)  ((*((MetaDataObject*)(v))->metadata))

#define MDQuery_Check(v) (((v)->ob_type == &MDQueryType))
#define MDQuery_Value(v)  ((*((MDQueryObject*)(v))->query))

#define Image_Check(v) (((v)->ob_type == &ImageType))
#define Image_Value(v) ((*((ImageObject*)(v))->image))

#define Program_Check(v) (((v)->ob_type == &ProgramType))
#define Program_Value(v) ((*((ProgramObject*)(v))->program))

#define RETURN_MDOBJECT(value) return new MDObject((MDLabel)label, value)
/*Helper function to create an MDObject from a PyObject */
static MDObject *
createMDObject(int label, PyObject *pyValue);
static PyObject * getMDObjectValue(MDObject * obj);
static void setMDObjectValue(MDObject * obj, PyObject *pyobj);

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
static void FileName_dealloc(FileNameObject* self)
{
    delete self->filename;
    self->ob_type->tp_free((PyObject*) self);
}

/* Constructor */
static PyObject *
FileName_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) type->tp_alloc(type, 0);

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = "", *ext = "";
        int number = ALL_IMAGES;
        if (PyArg_ParseTuple(args, "|Ois", &input, &number, &ext))
            //|| PyArg_ParseTuple(args, "|Os", &input, &ext)) FIXME
        {
            pyStr = PyObject_Str(input);
            if (pyStr != NULL)
                str = PyString_AsString(pyStr);
        }
        if (number != ALL_IMAGES)
            self->filename = new FileName(str, number, ext);
        else
            self->filename = new FileName(str);

    }
    return (PyObject *) self;
}

/* String representation */
static PyObject *
FileName_repr(PyObject * obj)
{
    FileNameObject *self = (FileNameObject*) obj;
    return PyString_FromString(self->filename->c_str());
}

/* compose */
static PyObject *
FileName_compose(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = "", *ext = "";
        int number = -1;
        size_t n = PyTuple_Size(args);
        if (n == 3 && PyArg_ParseTuple(args, "Ois", &input, &number, &ext))
        {
            pyStr = PyObject_Str(input);
            if (pyStr != NULL)
                str = PyString_AsString(pyStr);
            self->filename->compose(str, number, ext);
        }
        else if (n == 2 && PyArg_ParseTuple(args, "iO", &number, &input))
        {
            pyStr = PyObject_Str(input);
            if (pyStr != NULL)
                str = PyString_AsString(pyStr);
            self->filename->compose(number, str);
        }
        else
            return NULL;
    }
    Py_RETURN_NONE;//Return None(similar to void in C)
}

/* composeBlock */
static PyObject *
FileName_composeBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *root = "", *ext = "", *block ="";
        int number = 1;
        PyArg_ParseTuple(args, "sis|s", &block, &number, &root, &ext);
        self->filename->composeBlock(block, number, root, ext);
    }
    Py_RETURN_NONE;//Return None(similar to void in C)
}

/* exists */
static PyObject *
FileName_exists(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (self->filename->exists())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isInStack */
static PyObject *
FileName_isInStack(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (self->filename->isInStack())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isMetadata */
static PyObject *
FileName_isMetaData(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (self->filename->isMetaData(false))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isImage */
static PyObject *
FileName_isImage(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (isImage(*(self->filename)))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isStar1 */
static PyObject *
FileName_isStar1(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    if (self->filename->isStar1(false))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
FileName_getExtension(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    return PyString_FromString(self->filename->getExtension().c_str());
}

static PyObject *
FileName_getNumber(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    return PyInt_FromLong(self->filename->getNumber());
}

static PyObject *
FileName_getBaseName(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;

    return PyString_FromString(self->filename->getBaseName().c_str());
}

static PyObject *
FileName_decompose(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;
    size_t no;
    String str;
    self->filename->decompose(no, str);
    return Py_BuildValue("is", no, str.c_str());
}

static PyObject *
FileName_withoutExtension(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;
    return PyString_FromString(self->filename->withoutExtension().c_str());
}

static PyObject *
FileName_removeBlockName(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self = (FileNameObject*) obj;
    return PyString_FromString(self->filename->removeBlockName().c_str());
}

/* FileName methods */
static PyMethodDef FileName_methods[] =
    {
        { "compose", (PyCFunction) FileName_compose, METH_VARARGS,
          "Compose from root, number and extension OR prefix with number @" },
        { "composeBlock", (PyCFunction) FileName_composeBlock, METH_VARARGS,
          "Compose from blockname, number, root and extension" },
        { "decompose", (PyCFunction) FileName_decompose, METH_NOARGS,
          "Decompose filenames with @. Mainly from selfiles" },
        { "getBaseName", (PyCFunction) FileName_getBaseName, METH_NOARGS,
          "Get the base name from a FileName" },
        { "getExtension", (PyCFunction) FileName_getExtension, METH_NOARGS,
          "Get the last extension from a FileName" },
        { "getNumber", (PyCFunction) FileName_getNumber, METH_NOARGS,
          "Get the number from a FileName" },
        { "isInStack", (PyCFunction) FileName_isInStack, METH_NOARGS,
          "True if filename has stack format" },
        { "exists", (PyCFunction) FileName_exists, METH_NOARGS,
          "True if FileName exists" },
        { "isMetaData", (PyCFunction) FileName_isMetaData, METH_NOARGS,
          "True if is a MetaData" },
        { "isImage", (PyCFunction) FileName_isImage, METH_NOARGS,
          "True if is an image" },
        { "isStar1", (PyCFunction) FileName_isStar1, METH_NOARGS,
          "True if is a Star1" },
        { "withoutExtension", (PyCFunction) FileName_withoutExtension, METH_NOARGS,
          "return filename without extension" },
        { "removeBlockName", (PyCFunction) FileName_removeBlockName, METH_NOARGS,
          "return filename without block" },
        { NULL } /* Sentinel */
    };

/*FileName Type */
static PyTypeObject FileNameType = {
                                       PyObject_HEAD_INIT(NULL)
                                       0, /*ob_size*/
                                       "xmipp.FileName", /*tp_name*/
                                       sizeof(FileNameObject), /*tp_basicsize*/
                                       0, /*tp_itemsize*/
                                       (destructor)FileName_dealloc, /*tp_dealloc*/
                                       0, /*tp_print*/
                                       0, /*tp_getattr*/
                                       0, /*tp_setattr*/
                                       0, /*tp_compare*/
                                       FileName_repr, /*tp_repr*/
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
                                       "Python wrapper to Xmipp FileName class",/* tp_doc */
                                       0, /* tp_traverse */
                                       0, /* tp_clear */
                                       0, /* tp_richcompare */
                                       0, /* tp_weaklistoffset */
                                       0, /* tp_iter */
                                       0, /* tp_iternext */
                                       FileName_methods, /* tp_methods */
                                       0, /* tp_members */
                                       0, /* tp_getset */
                                       0, /* tp_base */
                                       0, /* tp_dict */
                                       0, /* tp_descr_get */
                                       0, /* tp_descr_set */
                                       0, /* tp_dictoffset */
                                       0, /* tp_init */
                                       0, /* tp_alloc */
                                       FileName_new, /* tp_new */
                                   };

/***************************************************************/
/*                            Image                         */
/**************************************************************/

/*Image Object*/
typedef struct
{
    PyObject_HEAD
    ImageGeneric * image;
}
ImageObject;

#define ImageObject_New() (ImageObject*)malloc(sizeof(ImageObject))

/* Destructor */
static void Image_dealloc(ImageObject* self)
{
    delete self->image;
    self->ob_type->tp_free((PyObject*) self);
}

/* Constructor */
static PyObject *
Image_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) type->tp_alloc(type, 0);
    if (self != NULL)
    {
        PyObject *input = NULL;

        if (PyArg_ParseTuple(args, "|O", &input))
        {
            if (input != NULL)
            {
                try
                {
                    if (PyString_Check(input))
                        self->image = new ImageGeneric(PyString_AsString(input));
                    else if (FileName_Check(input))
                        self->image = new ImageGeneric(FileName_Value(input));
                    //todo: add copy constructor
                    else
                    {
                        PyErr_SetString(PyExc_TypeError,
                                        "Image_new: Expected string or FileName as first argument");
                        return NULL;
                    }
                }
                catch (XmippError &xe)
                {
                    PyErr_SetString(PyXmippError, xe.msg.c_str());
                    return NULL;
                }
            }
            else
                self->image = new ImageGeneric();
        }
        else
            return NULL;
    }
    return (PyObject *) self;
}

/* Image string representation */
static PyObject *
Image_repr(PyObject * obj)
{
    ImageObject *self = (ImageObject*) obj;
    String s;
    self->image->toString(s);
    return PyString_FromString(s.c_str());
}
/* Image compare function */
static int
Image_compare(PyObject * obj, PyObject * obj2)
{
    int result = -1;
    if (obj != NULL && obj2 != NULL)
    {
        try
        {
            if (Image_Value(obj) == Image_Value(obj2))
                result = 0;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return result;
}

/* write */
static PyObject *
Image_write(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        if (PyArg_ParseTuple(args, "O", &input))
        {
            try
            {
                if (PyString_Check(input))
                    self->image->write(PyString_AsString(input));
                else if (FileName_Check(input))
                    self->image->write(FileName_Value(input));
                else
                    return NULL;
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

/* read */
static PyObject *
Image_read(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;

    if (self != NULL)
    {
        int datamode = DATA;
        PyObject *input = NULL;
        if (PyArg_ParseTuple(args, "O|i", &input, &datamode))
        {

            try
            {
                if (PyString_Check(input))
                    self->image->read(PyString_AsString(input),(DataMode)datamode);
                else if (FileName_Check(input))
                    self->image->read(FileName_Value(input),(DataMode)datamode);
                else
                    return NULL;
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

/* read preview*/
static PyObject *
Image_readPreview(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL;
        int x;
        if (PyArg_ParseTuple(args, "Oi", &input, &x))
        {
            try
            {
                if (PyString_Check(input))
                    self->image->readPreview(PyString_AsString(input), x);
                else if (FileName_Check(input))
                    self->image->readPreview(FileName_Value(input), x);
                else
                    return NULL;
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


/* readApplyGeo */
//int ImageGeneric::readApplyGeo(const MetaData &md, size_t objId,
//bool only_apply_shifts, DataMode datamode, size_t select_img)
static PyObject *
Image_readApplyGeo(PyObject *obj, PyObject *args, PyObject *kwargs);

static PyObject *
Image_applyGeo(PyObject *obj, PyObject *args, PyObject *kwargs);


/* getPixel */
static PyObject *
Image_getPixel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    int i, j;

    if (self != NULL && PyArg_ParseTuple(args, "ii", &i, &j))
    {
        try
        {
            double value = self->image->getPixel(i, j);
            return PyFloat_FromDouble(value);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }

    return NULL;
}

/* setPixel */
static PyObject *
Image_setPixel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    int i, j;
    double value = -1;

    if (self != NULL && PyArg_ParseTuple(args, "iid", &i, &j, &value))
    {
        try
        {
            self->image->setPixel(i, j, value);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }

    return NULL;
}

/* initConstant */
static PyObject *
Image_initConstant(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    double value = -1;

    if (self != NULL && PyArg_ParseTuple(args, "d", &value))
    {
        try
        {
            self->image->initConstant(value);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }

    return NULL;
}

/* initRandom */
static PyObject *
Image_initRandom(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    double op1 = 0;
    double op2 = 1;
    RandomMode mode = RND_Uniform;
    int pyMode = -1;

    if (self != NULL && PyArg_ParseTuple(args, "|ddi", &op1, &op2, &pyMode))
    {
        try
        {
            if (pyMode != -1)
                mode = (RandomMode) pyMode;

            self->image->initRandom(op1, op2, mode);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }

    return NULL;

}

/* Resize Image */
static PyObject *
Image_resize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    double xDim = 0, yDim = 0, zDim = 1;
    size_t nDim = 1;

    if (self != NULL && PyArg_ParseTuple(args, "dd|dn", &xDim, &yDim, &zDim, &nDim))
    {
        try
        {
            self->image->resize(xDim, yDim, zDim, nDim, false); // TODO: Take care of copy mode if needed
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }

    return NULL;

}

/* Return image dimensions as a tuple */
static PyObject *
Image_getDimensions(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    if (self != NULL)
    {
        try
        {
            int xdim, ydim, zdim;
            size_t ndim;
            self->image->image->getDimensions(xdim, ydim, zdim, ndim);
            return Py_BuildValue("iiik", xdim, ydim, zdim, ndim);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* Return image dimensions as a tuple */
static PyObject *
Image_getEulerAngles(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    if (self != NULL)
    {
        try
        {
            double rot, tilt, psi;
            self->image->getEulerAngles(rot, tilt, psi);
            return Py_BuildValue("fff", rot, tilt, psi);

        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}//Image_getDimensions


/* Return image dimensions as a tuple */
static PyObject *
Image_computeStats(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;
    if (self != NULL)
    {
        try
        {
            double mean, dev, min, max;
            self->image->data->computeStats(mean, dev, min, max);


            return Py_BuildValue("ffff", mean, dev, min, max);

        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}//Image_getDimensions

static PyObject *
Image_add(PyObject *obj1, PyObject *obj2);
static PyObject *
Image_subtract(PyObject *obj1, PyObject *obj2);

static PyObject *
Image_iadd(PyObject *obj1, PyObject *obj2);
static PyObject *
Image_isubtract(PyObject *obj1, PyObject *obj2);


/* Image methods that behave like numbers */
static PyNumberMethods Image_NumberMethods = {
            Image_add, //binaryfunc  nb_add;
            Image_subtract, //binaryfunc  nb_subtract;
            0, //binaryfunc  nb_multiply;
            0, //binaryfunc  nb_divide;
            0, //binaryfunc  nb_remainder;
            0, //binaryfunc  nb_divmod;
            0, //ternaryfunc nb_power;
            0, //unaryfunc nb_negative;
            0, //unaryfunc nb_positive;
            0, //unaryfunc nb_absolute;
            0, //inquiry nb_nonzero;       /* Used by PyObject_IsTrue */
            0, //unaryfunc nb_invert;
            0, //binaryfunc  nb_lshift;
            0, //binaryfunc  nb_rshift;
            0, //binaryfunc  nb_and;
            0, //binaryfunc  nb_xor;
            0, //binaryfunc  nb_or;
            0, //coercion nb_coerce;       /* Used by the coerce() function */
            0, //unaryfunc nb_int;
            0, //unaryfunc nb_long;
            0, //unaryfunc nb_float;
            0, //unaryfunc nb_oct;
            0, //unaryfunc nb_hex;

            /* Added in release 2.0 */
            Image_iadd, //binaryfunc  nb_inplace_add;
            Image_isubtract, //binaryfunc  nb_inplace_subtract;
            0, //binaryfunc  nb_inplace_multiply;
            0, //binaryfunc  nb_inplace_divide;
            0, //binaryfunc  nb_inplace_remainder;
            0, //ternaryfunc nb_inplace_power;
            0, //binaryfunc  nb_inplace_lshift;
            0, //binaryfunc  nb_inplace_rshift;
            0, //binaryfunc  nb_inplace_and;
            0, //binaryfunc  nb_inplace_xor;
            0, //binaryfunc  nb_inplace_or;

            /* Added in release 2.2 */
            0, //binaryfunc  nb_floor_divide;
            0, //binaryfunc  nb_true_divide;
            0, //binaryfunc  nb_inplace_floor_divide;
            0, //binaryfunc  nb_inplace_true_divide;

            /* Added in release 2.5 */
            0 //unaryfunc nb_index;
        };

/* Image methods */
static PyMethodDef Image_methods[] =
    {
        { "applyGeo", (PyCFunction) Image_applyGeo, METH_VARARGS,
          "Apply geometry in refering metadata to an image" },
        { "read", (PyCFunction) Image_read, METH_VARARGS,
          "Read image from disk" },
        { "readPreview", (PyCFunction) Image_readPreview, METH_VARARGS,
          "Read image preview" },
        { "readApplyGeo", (PyCFunction) Image_readApplyGeo, METH_VARARGS,
          "Read image from disk applying geometry in refering metadata" },
        { "write", (PyCFunction) Image_write, METH_VARARGS,
          "Write image to disk" },
        { "getPixel", (PyCFunction) Image_getPixel, METH_VARARGS,
          "Return a pixel value" },
        { "initConstant", (PyCFunction) Image_initConstant, METH_VARARGS,
          "Initialize to value" },
        { "initRandom", (PyCFunction) Image_initRandom, METH_VARARGS,
          "Initialize to random value" },
        { "resize", (PyCFunction) Image_resize, METH_VARARGS,
          "resize the image dimensions" },
        { "setPixel", (PyCFunction) Image_setPixel, METH_VARARGS,
          "Set the value of some pixel" },
        { "getDimensions", (PyCFunction) Image_getDimensions, METH_VARARGS,
          "Return image dimensions as a tuple" },
        { "getEulerAngles", (PyCFunction) Image_getEulerAngles, METH_VARARGS,
          "Return euler angles as a tuple" },
        { "computeStats", (PyCFunction) Image_computeStats, METH_VARARGS,
          "Compute image statistics, return mean, dev, min and max" },
        { NULL } /* Sentinel */
    };

/*Image Type */
static PyTypeObject ImageType = {
                                    PyObject_HEAD_INIT(NULL)
                                    0, /*ob_size*/
                                    "xmipp.Image", /*tp_name*/
                                    sizeof(ImageObject), /*tp_basicsize*/
                                    0, /*tp_itemsize*/
                                    (destructor)Image_dealloc, /*tp_dealloc*/
                                    0, /*tp_print*/
                                    0, /*tp_getattr*/
                                    0, /*tp_setattr*/
                                    Image_compare, /*tp_compare*/
                                    Image_repr, /*tp_repr*/
                                    &Image_NumberMethods, /*tp_as_number*/
                                    0, /*tp_as_sequence*/
                                    0, /*tp_as_mapping*/
                                    0, /*tp_hash */
                                    0, /*tp_call*/
                                    0, /*tp_str*/
                                    0, /*tp_getattro*/
                                    0, /*tp_setattro*/
                                    0, /*tp_as_buffer*/
                                    Py_TPFLAGS_DEFAULT, /*tp_flags*/
                                    "Python wrapper to Xmipp Image class",/* tp_doc */
                                    0, /* tp_traverse */
                                    0, /* tp_clear */
                                    0, /* tp_richcompare */
                                    0, /* tp_weaklistoffset */
                                    0, /* tp_iter */
                                    0, /* tp_iternext */
                                    Image_methods, /* tp_methods */
                                    0, /* tp_members */
                                    0, /* tp_getset */
                                    0, /* tp_base */
                                    0, /* tp_dict */
                                    0, /* tp_descr_get */
                                    0, /* tp_descr_set */
                                    0, /* tp_dictoffset */
                                    0, /* tp_init */
                                    0, /* tp_alloc */
                                    Image_new, /* tp_new */
                                };

/* Add two images, operator + */
static PyObject *
Image_add(PyObject *obj1, PyObject *obj2)
{
    ImageObject * result = PyObject_New(ImageObject, &ImageType);
    if (result != NULL)
    {
        try
        {
            result->image = new ImageGeneric(Image_Value(obj1));
            Image_Value(result).add(Image_Value(obj2));
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return (PyObject *)result;
}//operator +

/** Image inplace add, equivalent to += operator */
static PyObject *
Image_iadd(PyObject *obj1, PyObject *obj2)
{
    ImageObject * result = NULL;
    try
    {
        Image_Value(obj1).add(Image_Value(obj2));
        if (result = PyObject_New(ImageObject, &ImageType))
            result->image = new ImageGeneric(Image_Value(obj1));
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return (PyObject *)result;
}//operator +=

/* Subtract two images, operator - */
static PyObject *
Image_subtract(PyObject *obj1, PyObject *obj2)
{
    ImageObject * result = PyObject_New(ImageObject, &ImageType);
    if (result != NULL)
    {
        try
        {
            result->image = new ImageGeneric(Image_Value(obj1));
            Image_Value(result).subtract(Image_Value(obj2));
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return (PyObject *)result;
}//operator +

/** Image inplace subtraction, equivalent to -= operator */
static PyObject *
Image_isubtract(PyObject *obj1, PyObject *obj2)
{
    ImageObject * result = NULL;
    try
    {
        Image_Value(obj1).subtract(Image_Value(obj2));
        if (result = PyObject_New(ImageObject, &ImageType))
            result->image = new ImageGeneric(Image_Value(obj1));
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return (PyObject *)result;
}//operator +=


/***************************************************************/
/*                            Program                         */
/**************************************************************/

/*Program Object*/
typedef struct
{
    PyObject_HEAD
    XmippProgramGeneric * program;
}
ProgramObject;

#define ProgramObject_New() (ProgramObject*)malloc(sizeof(ProgramObject))

/* Destructor */
static void Program_dealloc(ProgramObject* self)
{
    delete self->program;
    self->ob_type->tp_free((PyObject*) self);
}

/* Constructor */
static PyObject *
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
static PyObject *
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
static PyObject *
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
static PyObject *
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
static PyObject *
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
static PyObject *
Program_endDefinition(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ProgramObject *self = (ProgramObject*) obj;
    if (self != NULL)
    {
        int verbose = 1;
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
static PyObject *
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
                self->program->read((int)size, argv);
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
static PyObject *
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
static PyObject *
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
static PyObject *
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

                for (int i = 0; i < size; ++i)
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
static PyTypeObject ProgramType = {
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




/***************************************************************/
/*                            MDQuery                          */
/***************************************************************/

/*MDQuery Object*/
typedef struct
{
    PyObject_HEAD
    MDQuery * query;
}
MDQueryObject;

/* Destructor */
static void MDQuery_dealloc(MDQueryObject* self)
{
    delete self->query;
    self->ob_type->tp_free((PyObject*) self);
}

/* String representation */
static PyObject *
MDQuery_repr(PyObject * obj)
{
    MDQueryObject *self = (MDQueryObject*) obj;
    if (self->query)
    {
        String s = self->query->whereString() + self->query->limitString()
                   + self->query->orderByString();
        return PyString_FromString(s.c_str());
    }
    else
        return PyString_FromString("");
}

/* MDQuery methods */
static PyMethodDef MDQuery_methods[] = { { NULL } /* Sentinel */
                                       };

/*MDQuery Type */
static PyTypeObject MDQueryType = {
                                      PyObject_HEAD_INIT(NULL)
                                      0, /*ob_size*/
                                      "xmipp.MDQuery", /*tp_name*/
                                      sizeof(MDQueryObject), /*tp_basicsize*/
                                      0, /*tp_itemsize*/
                                      (destructor)MDQuery_dealloc, /*tp_dealloc*/
                                      0, /*tp_print*/
                                      0, /*tp_getattr*/
                                      0, /*tp_setattr*/
                                      0, /*tp_compare*/
                                      MDQuery_repr, /*tp_repr*/
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
                                      "Python wrapper to Xmipp MDQuery class",/* tp_doc */
                                      0, /* tp_traverse */
                                      0, /* tp_clear */
                                      0, /* tp_richcompare */
                                      0, /* tp_weaklistoffset */
                                      0, /* tp_iter */
                                      0, /* tp_iternext */
                                      MDQuery_methods, /* tp_methods */
                                      0, /* tp_members */
                                      0, /* tp_getset */
                                      0, /* tp_base */
                                      0, /* tp_dict */
                                      0, /* tp_descr_get */
                                      0, /* tp_descr_set */
                                      0, /* tp_dictoffset */
                                      0, /* tp_init */
                                      0, /* tp_alloc */
                                      0, /* tp_new */
                                  };

/***************************************************************/
/*                            MetaData                         */
/**************************************************************/

/*MetaData Object*/
typedef struct
{
    PyObject_HEAD
    MetaData * metadata;
    MDIterator * iter;
}
MetaDataObject;

/* Destructor */
static void MetaData_dealloc(MetaDataObject* self)
{
    delete self->metadata;
    delete self->iter;
    self->ob_type->tp_free((PyObject*) self);
}

/* Constructor */
static PyObject *
MetaData_aggregate(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_aggregateSingle(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_join(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_importObjects(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_intersection(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_merge(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_operate(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_replace(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_readPlain(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_setComment(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_getComment(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_unionAll(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_randomize(PyObject *obj, PyObject *args, PyObject *kwargs);
static PyObject *
MetaData_selectPart(PyObject *obj, PyObject *args, PyObject *kwargs);

static int MetaData_print(PyObject *obj, FILE *fp, int flags)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        std::stringstream ss;
        self->metadata->write(ss);
        fprintf(fp, "%s", ss.str().c_str());
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
        return -1;
    }
    return 0;
}

/* String representation */
static PyObject *
MetaData_repr(PyObject * obj)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyString_FromString(
               (self->metadata->getFilename() + "(MetaData)").c_str());
}

/* MetaData compare function */
static int
MetaData_compare(PyObject * obj, PyObject * obj2)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    MetaDataObject *md2 = (MetaDataObject*) obj2;
    int result = -1;

    if (self != NULL && md2 != NULL)
    {
        try
        {
            if (*(self->metadata) == *(md2->metadata))
                result = 0;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return result;
}

/* read */
static PyObject *
MetaData_read(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *list = NULL;
        PyObject *input = NULL, *pyStr = NULL;
        char *str = NULL;
        int number = -1;
        if (PyArg_ParseTuple(args, "O|O", &input,  &list))
        {
            try
            {
                if ((pyStr = PyObject_Str(input)) != NULL)
                {
                    str = PyString_AsString(pyStr);
                    if (list!=NULL)
                    {
                        if (PyList_Check(list))
                        {
                            size_t size = PyList_Size(list);
                            PyObject * item = NULL;
                            int iValue = 0;
                            std::vector<MDLabel> vValue(size);
                            for (size_t i = 0; i < size; ++i)
                            {
                                item = PyList_GetItem(list, i);
                                if (!PyInt_Check(item))
                                {
                                    PyErr_SetString(PyExc_TypeError,
                                                    "MDL labels must be integers (MDLABEL)");
                                    return NULL;
                                }
                                iValue = PyInt_AsLong(item);
                                vValue[i] = (MDLabel)iValue;
                            }
                            self->metadata->read(str,&vValue);
                        }
                    }
                    else
                        self->metadata->read(str);
                    Py_RETURN_NONE;
                }
                else
                    return NULL;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
    }
    return NULL;
}

/* read */
static PyObject *
MetaData_readPlain(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *input2 = NULL;
        PyObject *pyStr = NULL, *pyLabels = NULL, *pySep = NULL;
        char *str = NULL, *sep = NULL, *labels = NULL;
        int number = -1;

        if (PyArg_ParseTuple(args, "OO|O", &input, &input2, &pySep))
        {
            try
            {
                if ((pyStr = PyObject_Str(input)) != NULL && (pyLabels
                        = PyObject_Str(input2)))
                {
                    str = PyString_AsString(pyStr);
                    labels = PyString_AsString(pyLabels);
                    self->metadata->readPlain(str, labels);
                    Py_RETURN_NONE;
                }
                else
                    return NULL;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
    }
    return NULL;
}

/* read block */
static PyObject *
MetaData_readBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *blockName = NULL, *pyStr = NULL, *pyStrBlock =
                                                 NULL;
        char *str = NULL, *strBlock = NULL;
        int number = -1;
        if (PyArg_ParseTuple(args, "OO", &input, &blockName))
        {
            try
            {
                if ((pyStr = PyObject_Str(input)) != NULL && (pyStrBlock
                        = PyObject_Str(blockName)) != NULL)
                {
                    str = PyString_AsString(pyStr);
                    strBlock = PyString_AsString(pyStrBlock);
                    self->metadata->read((std::string) (strBlock) + "@" + str,
                                         NULL);
                    Py_RETURN_NONE;
                }
                else
                    return NULL;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
    }
    return NULL;
}

/* write */
static PyObject *
MetaData_write(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    WriteModeMetaData wmd;
    wmd = MD_OVERWRITE;
    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = NULL;
        int number = -1;
        if (PyArg_ParseTuple(args, "O|i", &input, &wmd))
        {
            try
            {
                if ((pyStr = PyObject_Str(input)) != NULL)
                {
                    str = PyString_AsString(pyStr);
                    self->metadata->write(str, (WriteModeMetaData) wmd);
                    Py_RETURN_NONE;
                }
                else
                    return NULL;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
    }
    return NULL;
}

/* write */
static PyObject *
MetaData_writeBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *blockName = NULL;
        int number = -1;
        if (PyArg_ParseTuple(args, "OO", &input, &blockName))
        {
            try
            {
                String fn, block;
                if (PyString_Check(input))
                    fn = PyString_AsString(input);
                else if (FileName_Check(input))
                    fn = FileName_Value(input);
                else
                    return NULL;
                if (PyString_Check(blockName))
                    block = PyString_AsString(blockName);
                else if (FileName_Check(blockName))
                    block = FileName_Value(blockName);
                else
                    return NULL;
                self->metadata->_write(fn, block, MD_APPEND);
                Py_RETURN_NONE;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
    }
    return NULL;
}

/* append */
static PyObject *
MetaData_append(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = NULL;
        int number = -1;
        if (PyArg_ParseTuple(args, "O", &input))
        {
            try
            {
                if (PyString_Check(input))
                    self->metadata->append(PyString_AsString(input));
                else if (FileName_Check(input))
                    self->metadata->append(FileName_Value(input));
                else
                    return NULL;
                Py_RETURN_NONE;
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
    }
    return NULL;
}
/* addObject */
static PyObject *
MetaData_addObject(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyLong_FromUnsignedLong(self->metadata->addObject());
}
/* firstObject */
static PyObject *
MetaData_firstObject(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyLong_FromUnsignedLong(self->metadata->firstObject());
}
/* lastObject */
static PyObject *
MetaData_lastObject(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyLong_FromUnsignedLong(self->metadata->lastObject());
}
/* size */
static PyObject *
MetaData_size(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        return PyLong_FromUnsignedLong(self->metadata->size());
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
        return NULL;
    }
}
/* isEmpty */
static PyObject *
MetaData_isEmpty(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        if (self->metadata->isEmpty())
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
        return NULL;
    }
}
/* getColumnFormat */
static PyObject *
MetaData_getColumnFormat(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        if (self->metadata->getColumnFormat())
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
        return NULL;
    }
}
/* setColumnFormat */
static PyObject *
MetaData_setColumnFormat(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *input = NULL;
    if (PyArg_ParseTuple(args, "O", &input))
    {
        try
        {
            if (PyBool_Check(input))
            {
                MetaDataObject *self = (MetaDataObject*) obj;
                self->metadata->setColumnFormat(input == Py_True);
                Py_RETURN_NONE;
            }
            else
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::setColumnFormat: Expecting boolean value");
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* setValue */
static PyObject *
MetaData_setValue(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    size_t objectId = BAD_OBJID;
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "iOk", &label, &pyValue, &objectId))
    {
        try
        {
            MDObject * object = createMDObject(label, pyValue);
            if (!object)
                return NULL;
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->setValue(*object, objectId);
            delete object;
            Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* setValueCol */
static PyObject *
MetaData_setValueCol(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "iO", &label, &pyValue))
    {
        try
        {
            MDObject * object = createMDObject(label, pyValue);
            if (!object)
                return NULL;
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->setValueCol(*object);
            delete object;
            Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* drop column */
static PyObject *
MetaData_removeLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            if (self->metadata->removeLabel((MDLabel) label))
                Py_RETURN_TRUE;
            else
                Py_RETURN_FALSE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* getValue */
static PyObject *
MetaData_getValue(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    size_t objectId = BAD_OBJID;
    PyObject *pyValue;

    if (PyArg_ParseTuple(args, "ik", &label, &objectId))
    {
        try
        {
            MDObject * object = new MDObject((MDLabel) label);
            MetaDataObject *self = (MetaDataObject*) obj;
            if (self->metadata->getValue(*object, objectId))
            {
                pyValue = getMDObjectValue(object);
                delete object;
                return pyValue;
            }
            else
            {
                delete object;
                Py_RETURN_NONE;
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* getValue */
static PyObject *
MetaData_getColumnValues(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    if (PyArg_ParseTuple(args, "i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;

            std::vector<MDObject> v;
            self->metadata->getColumnValues((MDLabel) label,v);

            size_t size=v.size();
            PyObject * list = PyList_New(size);

            for (int i = 0; i < size; ++i)
                PyList_SetItem(list, i, getMDObjectValue(&(v[i])));

            return list;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* setValue */
static PyObject *
MetaData_setColumnValues(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *list = NULL;
    if (PyArg_ParseTuple(args, "iO", &label, &list))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            size_t size=PyList_Size(list);
            bool addObjects=(self->metadata->size()==0);
            MDObject object((MDLabel) label);
            if (addObjects)
            {
                for (size_t i=0; i<size; ++i)
                {
                    size_t id=self->metadata->addObject();
                    setMDObjectValue(&object,PyList_GetItem(list,i));
                    self->metadata->setValue(object,id);
                }
            }
            else
            {
                if (self->metadata->size()!=size)
                    PyErr_SetString(PyXmippError, "Metadata size different from list size");
                size_t i=0;
                FOR_ALL_OBJECTS_IN_METADATA(*(self->metadata))
                {
                    setMDObjectValue(&object,PyList_GetItem(list,i++));
                    self->metadata->setValue(object,__iter.objId);
                }
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
            return NULL;
        }
    }
    Py_RETURN_NONE;
}

/* containsLabel */
static PyObject *
MetaData_getActiveLabels(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyValue;
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        std::vector<MDLabel>* labels = self->metadata->getActiveLabelsAddress();
        int size = labels->size();
        PyObject * list = PyList_New(size);

        for (int i = 0; i < size; ++i)
            PyList_SetItem(list, i, PyInt_FromLong(labels->at(i)));

        return list;

    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}
/* containsLabel */
static PyObject *
xmipp_getBlocksInMetaDataFile(PyObject *obj, PyObject *args)
{
    int label;
    PyObject *input;
    FileName fn;
    StringVector blocks;

    try
    {
        if (PyArg_ParseTuple(args, "O", &input))
        {

            if (PyString_Check(input))
                fn = PyString_AsString(input);
            else if (FileName_Check(input))
                fn = FileName_Value(input);
            else
                return NULL;
            getBlocksInMetaDataFile(fn, blocks);
            int size = blocks.size();
            PyObject * list = PyList_New(size);

            for (int i = 0; i < size; ++i)
            {
                PyList_SetItem(list, i, PyString_FromString(blocks[i].c_str()));
            }
            return list;
        }
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}

/* containsLabel */
static PyObject *
MetaData_getMaxStringLength(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyValue;
    if (PyArg_ParseTuple(args, "i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            int length = self->metadata->getMaxStringLength((MDLabel) label);

            return PyInt_FromLong(length);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* containsLabel */
static PyObject *
MetaData_containsLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyValue;
    if (PyArg_ParseTuple(args, "i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            if (self->metadata->containsLabel((MDLabel) label))
                Py_RETURN_TRUE;
            else
                Py_RETURN_FALSE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/* addLabel */
static PyObject *
MetaData_addLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label, pos = -1;
    PyObject *pyValue;
    if (PyArg_ParseTuple(args, "i|i", &label, &pos))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->addLabel((MDLabel) label);
            Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/* removeObjects */
static PyObject *
MetaData_removeObjects(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyQuery = NULL;

    if (PyArg_ParseTuple(args, "O", &pyQuery))
    {
        try
        {
            if (!MDQuery_Check(pyQuery))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::removeObjects: Expecting MDQuery as second arguments");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->removeObjects(MDQuery_Value(pyQuery));
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* removeObjects */
static PyObject *
MetaData_removeDisabled(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyQuery = NULL;

    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        self->metadata->removeDisabled();
        Py_RETURN_NONE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}

/* Make absolute path */
static PyObject *
MetaData_makeAbsPath(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label = (MDLabel) MDL_IMAGE;
    if (PyArg_ParseTuple(args, "|i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->makeAbsPath((MDLabel) label);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* clear */
static PyObject *
MetaData_clear(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        self->metadata->clear();
        Py_RETURN_NONE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
        return NULL;
    }
}

/** Iteration functions */
static PyObject *
MetaData_iter(PyObject *obj)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        self->iter = new MDIterator(*(self->metadata));
        Py_INCREF(self);
        return (PyObject *) self;
        //return Py_BuildValue("l", self->metadata->iteratorBegin());
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}
static PyObject *
MetaData_iternext(PyObject *obj)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        size_t objId = self->iter->objId;
        self->iter->moveNext();
        if (objId == BAD_OBJID)
            return NULL;
        //type format should be "n" instead of "i" but I put i since python 2.4 does not support n
        return Py_BuildValue("i", objId);
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}
/** Sort Metadata */
static PyObject *
MetaData_sort(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label = (MDLabel) MDL_IMAGE;
    if (PyArg_ParseTuple(args, "|i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            MetaData MDaux = *(self->metadata);
            self->metadata->clear();
            self->metadata->sort(MDaux, (MDLabel) label);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/** remove Duplicate rows Metadata */
static PyObject *
MetaData_removeDuplicates(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        MetaData MDaux = *(self->metadata);
        self->metadata->clear();
        self->metadata->removeDuplicates(MDaux);
        Py_RETURN_NONE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}
/* MetaData methods */
static PyMethodDef
MetaData_methods[] =
    {
        { "read", (PyCFunction) MetaData_read, METH_VARARGS,
          "Read data from file" },
        { "write", (PyCFunction) MetaData_write, METH_VARARGS,
          "Write MetaData content to disk" },
        { "readBlock", (PyCFunction) MetaData_readBlock,
          METH_VARARGS, "Read block from a metadata file" },
        { "writeBlock", (PyCFunction) MetaData_writeBlock,
          METH_VARARGS, "Append block to a metadata" },
        { "append", (PyCFunction) MetaData_append,
          METH_VARARGS, "Append MetaData content to disk" },
        { "addObject", (PyCFunction) MetaData_addObject,
          METH_NOARGS,
          "Add a new object and return its id" },
        { "firstObject", (PyCFunction) MetaData_firstObject,
          METH_NOARGS,
          "Goto first metadata object, return its object id" },
        { "lastObject", (PyCFunction) MetaData_lastObject,
          METH_NOARGS,
          "Goto last metadata object, return its object id" },
        { "size", (PyCFunction) MetaData_size, METH_NOARGS,
          "Return number of objects in MetaData" },
        { "isEmpty", (PyCFunction) MetaData_isEmpty,
          METH_NOARGS,
          "Check whether the MetaData is empty" },
        { "clear", (PyCFunction) MetaData_clear, METH_NOARGS,
          "Clear MetaData" },
        { "getColumnFormat",
          (PyCFunction) MetaData_getColumnFormat,
          METH_NOARGS, "Get column format info" },
        { "setColumnFormat",
          (PyCFunction) MetaData_setColumnFormat,
          METH_VARARGS, "Set column format info" },
        { "setValue", (PyCFunction) MetaData_setValue,
          METH_VARARGS,
          "Set the value for column(label) for a given object" },
        { "setValueCol", (PyCFunction) MetaData_setValueCol,
          METH_VARARGS,
          "Set the same value for column(label) for all objects" },
        { "removeLabel", (PyCFunction) MetaData_removeLabel,
          METH_VARARGS,
          "Remove a label if exists. The values are still in the table." },
        { "getValue", (PyCFunction) MetaData_getValue,
          METH_VARARGS, "Get the value for column(label)" },
        { "getColumnValues", (PyCFunction) MetaData_getColumnValues,
          METH_VARARGS, "Get all values value from column(label)" },
        { "setColumnValues", (PyCFunction) MetaData_setColumnValues,
          METH_VARARGS, "Set all values value from column(label)" },
        { "getActiveLabels",
          (PyCFunction) MetaData_getActiveLabels,
          METH_VARARGS,
          "Return a list with the labels of the Metadata" },
        { "getMaxStringLength",
          (PyCFunction) MetaData_getMaxStringLength,
          METH_VARARGS,
          "Return the maximun lenght of a value on this column(label)" },
        { "containsLabel",
          (PyCFunction) MetaData_containsLabel,
          METH_VARARGS,
          "True if this metadata contains this label" },
        { "addLabel", (PyCFunction) MetaData_addLabel,
          METH_VARARGS, "Add a new label to MetaData" },
        { "makeAbsPath", (PyCFunction) MetaData_makeAbsPath,
          METH_VARARGS,
          "Make filenames with absolute paths" },
        { "importObjects",
          (PyCFunction) MetaData_importObjects,
          METH_VARARGS,
          "Import objects from another metadata" },
        { "removeObjects",
          (PyCFunction) MetaData_removeObjects,
          METH_VARARGS, "Remove objects from metadata" },
        { "removeDisabled",
          (PyCFunction) MetaData_removeDisabled,
          METH_VARARGS, "Remove disabled objects from metadata" },
        { "aggregateSingle",
          (PyCFunction) MetaData_aggregateSingle,
          METH_VARARGS,
          "Aggregate operation in metadata (single value result)" },
        { "aggregate", (PyCFunction) MetaData_aggregate,
          METH_VARARGS,
          "Aggregate operation in metadata. The results is stored in self." },
        { "unionAll", (PyCFunction) MetaData_unionAll,
          METH_VARARGS,
          "Union of two metadatas. The results is stored in self." },
        { "merge", (PyCFunction) MetaData_merge, METH_VARARGS,
          "Merge columns of two metadatas. The results is stored in self." },
        {
            "join",
            (PyCFunction) MetaData_join,
            METH_VARARGS,
            "join between two metadatas, use MDL_UNDEFINED as label. The results is stored in self." },
        { "readPlain", (PyCFunction) MetaData_readPlain,
          METH_VARARGS,
          "Import metadata from a plain text file." },
        {"intersection",
         (PyCFunction) MetaData_intersection,
         METH_VARARGS,
         "Intersection of two metadatas using a common label. The results is stored in self." },
        { "setComment", (PyCFunction) MetaData_setComment,
          METH_VARARGS, "Set comment in Metadata." },
        { "getComment", (PyCFunction) MetaData_getComment,
          METH_VARARGS, "Get comment in Metadata." },
        { "operate", (PyCFunction) MetaData_operate,
          METH_VARARGS, "Replace strings values in some column." },
        { "replace", (PyCFunction) MetaData_replace,
          METH_VARARGS, "Basic operations on columns data." },
        {"randomize",
         (PyCFunction) MetaData_randomize,
         METH_VARARGS,
         "Randomize another metadata and keep in self." },
        {"selectPart",
         (PyCFunction) MetaData_selectPart,
         METH_VARARGS,
         "select a part of another metadata starting from start and with a number of objects" },
        { "removeDuplicates", (PyCFunction) MetaData_removeDuplicates,
          METH_VARARGS, "Remove duplicate rows" },
        {
            "sort", (PyCFunction) MetaData_sort,
            METH_VARARGS,
            "Sort metadata according to a label" },
        { NULL } /* Sentinel */
    };

/*MetaData Type */
static PyTypeObject MetaDataType = {
                                       PyObject_HEAD_INIT(NULL)
                                       0, /*ob_size*/
                                       "xmipp.MetaData", /*tp_name*/
                                       sizeof(MetaDataObject), /*tp_basicsize*/
                                       0, /*tp_itemsize*/
                                       (destructor)MetaData_dealloc, /*tp_dealloc*/
                                       MetaData_print, /*tp_print*/
                                       0, /*tp_getattr*/
                                       0, /*tp_setattr*/
                                       MetaData_compare, /*tp_compare*/
                                       MetaData_repr, /*tp_repr*/
                                       0, /*tp_as_number*/
                                       0, /*tp_as_sequence*/
                                       0, /*tp_as_mapping*/
                                       0, /*tp_hash */
                                       0, /*tp_call*/
                                       0, /*tp_str*/
                                       0, /*tp_getattro*/
                                       0, /*tp_setattro*/
                                       0, /*tp_as_buffer*/
                                       Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /*tp_flags*/
                                       "Python wrapper to Xmipp MetaData class",/* tp_doc */
                                       0, /* tp_traverse */
                                       0, /* tp_clear */
                                       0, /* tp_richcompare */
                                       0, /* tp_weaklistoffset */
                                       MetaData_iter, /* tp_iter */
                                       MetaData_iternext, /* tp_iternext */
                                       MetaData_methods, /* tp_methods */
                                       0, /* tp_members */
                                       0, /* tp_getset */
                                       0, /* tp_base */
                                       0, /* tp_dict */
                                       0, /* tp_descr_get */
                                       0, /* tp_descr_set */
                                       0, /* tp_dictoffset */
                                       0, /* tp_init */
                                       0, /* tp_alloc */
                                       MetaData_new, /* tp_new */
                                   };

PyObject *
MetaData_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) type->tp_alloc(type, 0);

    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        PyArg_ParseTuple(args, "|O", &input);
        if (input != NULL)
        {
            try
            {
                if (MetaData_Check(input))
                    self->metadata = new MetaData(MetaData_Value(input));
                else if ((pyStr = PyObject_Str(input)) != NULL)
                {
                    char * str = PyString_AsString(pyStr);
                    self->metadata = new MetaData(str);
                }
                else
                {
                    PyErr_SetString(PyExc_TypeError,
                                    "MetaData_new: Bad string value for reading metadata");
                    return NULL;
                }
            }
            catch (XmippError &xe)
            {
                PyErr_SetString(PyXmippError, xe.msg.c_str());
                return NULL;
            }
        }
        else
        {
            self->metadata = new MetaData();
        }
    }
    return (PyObject *) self;
}

/* importObjects */
static PyObject *
MetaData_importObjects(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyMd = NULL;
    PyObject *pyQuery = NULL;

    if (PyArg_ParseTuple(args, "OO", &pyMd, &pyQuery))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::importObjects: Expecting MetaData as first argument");
                return NULL;
            }
            if (!MDQuery_Check(pyQuery))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::importObjects: Expecting MDQuery as second argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->importObjects(MetaData_Value(pyMd),
                                          MDQuery_Value(pyQuery));
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* aggregateSingle */
static PyObject *
MetaData_aggregateSingle(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    AggregateOperation op;
    MDLabel label;
    PyObject *pyValue;

    if (PyArg_ParseTuple(args, "ii", &op, &label))
    {
        try
        {
            MDObject * object = new MDObject((MDLabel) label);
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->aggregateSingle(*object, (AggregateOperation) op,
                                            (MDLabel) label);
            pyValue = getMDObjectValue(object);
            delete object;
            return pyValue;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/*aggregate*/
static PyObject *
MetaData_aggregate(PyObject *obj, PyObject *args, PyObject *kwargs)
{

    AggregateOperation op;
    MDLabel aggregateLabel;
    MDLabel operateLabel;
    MDLabel resultLabel;
    PyObject *pyMd = NULL;

    if (PyArg_ParseTuple(args, "Oiiii", &pyMd, &op, &aggregateLabel,
                         &operateLabel, &resultLabel))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::aggregate: Expecting MetaData as first argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->aggregate(MetaData_Value(pyMd),
                                      (AggregateOperation) op, (MDLabel) aggregateLabel,
                                      (MDLabel) operateLabel, (MDLabel) resultLabel);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* UnionAll */
static PyObject *
MetaData_unionAll(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyMd = NULL;
    PyObject *pyQuery = NULL;

    if (PyArg_ParseTuple(args, "O", &pyMd))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::unionAll: Expecting MetaData as first argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->unionAll(MetaData_Value(pyMd));
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* merge */
static PyObject *
MetaData_merge(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyMd = NULL;
    PyObject *pyQuery = NULL;

    if (PyArg_ParseTuple(args, "O", &pyMd))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::merge: Expecting MetaData as first argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->merge(MetaData_Value(pyMd));
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* setComment */
static PyObject *
MetaData_setComment(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    char *str = "", *ext = "";
    if (PyArg_ParseTuple(args, "s", &str))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->setComment(str);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* getComment */
static PyObject *
MetaData_getComment(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyString_FromString(self->metadata->getComment().c_str());
}

/* join */
static PyObject *
MetaData_join(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int labelLeft;
    int labelRight;
    PyObject *pyMdLeft = NULL;
    PyObject *pyMdright = NULL;
    PyObject *pyQuery = NULL;
    JoinType jt;

    if (PyArg_ParseTuple(args, "OOiii", &pyMdLeft, &pyMdright, &labelLeft,&labelRight, &jt))
    {
        try
        {
            if (!MetaData_Check(pyMdLeft))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::join: Expecting MetaData as first argument");
                return NULL;
            }
            if (!MetaData_Check(pyMdright))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::join: Expecting MetaData as second argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->join(MetaData_Value(pyMdLeft),
                                 MetaData_Value(pyMdright), (MDLabel) labelLeft,
                                 (MDLabel) labelRight, (JoinType) jt);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* Intersection */
static PyObject *
MetaData_intersection(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyMd = NULL;

    if (PyArg_ParseTuple(args, "Oi", &pyMd, &label))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::intersection: Expecting MetaData as first argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->intersection(MetaData_Value(pyMd), (MDLabel) label);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* Basic operations on columns data. */
static PyObject *
MetaData_operate(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    char *str = "";
    if (PyArg_ParseTuple(args, "s", &str))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->operate(str);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/*Replace string values in some label. */
static PyObject *
MetaData_replace(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    char *oldStr = "";
    char *newStr = "";
    if (PyArg_ParseTuple(args, "iss", &label, &oldStr, &newStr))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->replace((MDLabel)label, oldStr, newStr);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* Randomize */
static PyObject *
MetaData_randomize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyMd = NULL;

    if (PyArg_ParseTuple(args, "O", &pyMd))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::randomize: Expecting MetaData as first argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->randomize(MetaData_Value(pyMd));
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* Select Part */
static PyObject *
MetaData_selectPart(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label=(int)MDL_OBJID;
    size_t start, numberOfObjects;
    PyObject *pyMd = NULL;

    if (PyArg_ParseTuple(args, "Okk|i", &pyMd, &start, &numberOfObjects, &label))
    {
        try
        {
            if (!MetaData_Check(pyMd))
            {
                PyErr_SetString(PyExc_TypeError,
                                "MetaData::selectPart: Expecting MetaData as first argument");
                return NULL;
            }
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->selectPart(MetaData_Value(pyMd), start, numberOfObjects, (MDLabel) label);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

MDObject *
createMDObject(int label, PyObject *pyValue)
{
    try
    {
        if (PyInt_Check(pyValue))
        {
            int iValue = PyInt_AS_LONG(pyValue);
            RETURN_MDOBJECT(iValue);
        }
        if (PyLong_Check(pyValue))
        {
            size_t value = PyLong_AsUnsignedLong(pyValue);
            RETURN_MDOBJECT(value);
        }
        if (PyString_Check(pyValue))
        {
            RETURN_MDOBJECT(std::string(PyString_AsString(pyValue)));
        }
        if (FileName_Check(pyValue))
        {
            RETURN_MDOBJECT(*((FileNameObject*)pyValue)->filename);
        }
        if (PyFloat_Check(pyValue))
        {
            double dValue = PyFloat_AS_DOUBLE(pyValue);
            RETURN_MDOBJECT(double(dValue));
        }
        if (PyBool_Check(pyValue))
        {
            bool bValue = (pyValue == Py_True);
            RETURN_MDOBJECT(bValue);
        }
        if (PyList_Check(pyValue))
        {
            size_t size = PyList_Size(pyValue);
            PyObject * item = NULL;
            double dValue = 0.;
            std::vector<double> vValue(size);
            for (size_t i = 0; i < size; ++i)
            {
                item = PyList_GetItem(pyValue, i);
                if (!PyFloat_Check(item))
                {
                    PyErr_SetString(PyExc_TypeError,
                                    "Vectors are only supported for double");
                    return NULL;
                }
                dValue = PyFloat_AS_DOUBLE(item);
                vValue[i] = dValue;
            }
            RETURN_MDOBJECT(vValue);
        }
        PyErr_SetString(PyExc_TypeError, "Unrecognized type to create MDObject");
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}

void setMDObjectValue(MDObject *obj, PyObject *pyValue)
{
    try
    {
        if (PyInt_Check(pyValue))
            obj->setValue((int)PyInt_AS_LONG(pyValue));
        else if (PyLong_Check(pyValue))
            obj->setValue((size_t)PyLong_AsUnsignedLong(pyValue));
        else if (PyString_Check(pyValue))
            obj->setValue(std::string(PyString_AsString(pyValue)));
        else if (FileName_Check(pyValue))
            obj->setValue((*((FileNameObject*)pyValue)->filename));
        else if (PyFloat_Check(pyValue))
            obj->setValue(PyFloat_AS_DOUBLE(pyValue));
        else if (PyBool_Check(pyValue))
            obj->setValue((pyValue == Py_True));
        else if (PyList_Check(pyValue))
        {
            size_t size = PyList_Size(pyValue);
            PyObject * item = NULL;
            double dValue = 0.;
            std::vector<double> vValue(size);
            for (size_t i = 0; i < size; ++i)
            {
                item = PyList_GetItem(pyValue, i);
                if (!PyFloat_Check(item))
                {
                    PyErr_SetString(PyExc_TypeError,
                                    "Vectors are only supported for double");
                }
                dValue = PyFloat_AS_DOUBLE(item);
                vValue[i] = dValue;
            }
            obj->setValue(vValue);
        }
        else
            PyErr_SetString(PyExc_TypeError, "Unrecognized type to create MDObject");
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
}

static PyObject *
getMDObjectValue(MDObject * obj)
{
    if (obj->label == MDL_UNDEFINED) //if undefine label, store as a literal string
        return NULL;
    switch (MDL::labelType(obj->label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
        if (obj->data.boolValue)
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    case LABEL_INT:
        return PyInt_FromLong(obj->data.intValue);
    case LABEL_LONG:
        return PyLong_FromLong(obj->data.longintValue);
    case LABEL_DOUBLE:
        return PyFloat_FromDouble(obj->data.doubleValue);
    case LABEL_STRING:
        return PyString_FromString(obj->data.stringValue->c_str());
    case LABEL_VECTOR:
        std::vector<double> & vector = *(obj->data.vectorValue);
        int size = vector.size();
        PyObject * list = PyList_New(size);
        for (int i = 0; i < size; ++i)
            PyList_SetItem(list, i, PyFloat_FromDouble(vector[i]));
        return list;
    }//close switch
    return NULL;
}

static PyObject *
Image_readApplyGeo(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;

    if (self != NULL)
    {
        PyObject *md = NULL;
        PyObject *only_apply_shifts = Py_False;
        PyObject *wrap ;
        if (WRAP)
            wrap = Py_True;
        else
            wrap = Py_False;
        size_t objectId = BAD_OBJID;
        bool boolOnly_apply_shifts = false;
        bool boolWrap = WRAP;
        int datamode = DATA;
        size_t select_img = ALL_IMAGES;

        if (PyArg_ParseTuple(args, "Ok|OikO", &md, &objectId, &only_apply_shifts, &datamode, &select_img, &wrap))
        {
            try
            {
                if (PyBool_Check(only_apply_shifts))
                    boolOnly_apply_shifts = (only_apply_shifts == Py_True);
                else
                    PyErr_SetString(PyExc_TypeError, "ImageGeneric::readApplyGeo: Expecting boolean value");
                if (PyBool_Check(wrap))
                    boolWrap = (wrap == Py_True);
                else
                    PyErr_SetString(PyExc_TypeError, "ImageGeneric::readApplyGeo: Expecting boolean value");
                self->image->readApplyGeo(MetaData_Value(md), objectId, boolOnly_apply_shifts,
                                          (DataMode)datamode, select_img,boolWrap);
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

static PyObject *
Image_applyGeo(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    ImageObject *self = (ImageObject*) obj;

    if (self != NULL)
    {
        PyObject *md = NULL;
        PyObject *only_apply_shifts = Py_False;
        PyObject *wrap ;
        if (WRAP)
            wrap = Py_True;
        else
            wrap = Py_False;
        size_t objectId = BAD_OBJID;
        bool boolOnly_apply_shifts = false;
        bool boolWrap = WRAP;

        if (PyArg_ParseTuple(args, "Ok|OO", &md, &objectId, &only_apply_shifts, &wrap))
        {
            try
            {
                if (PyBool_Check(only_apply_shifts))
                    boolOnly_apply_shifts = (only_apply_shifts == Py_True);
                else
                    PyErr_SetString(PyExc_TypeError, "ImageGeneric::applyGeo: Expecting boolean value");
                if (PyBool_Check(wrap))
                    boolWrap = (wrap == Py_True);
                else
                    PyErr_SetString(PyExc_TypeError, "ImageGeneric::applyGeo: Expecting boolean value");
                self->image->applyGeo(MetaData_Value(md), objectId, boolOnly_apply_shifts,boolWrap);
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
/***************************************************************/
/*                            Global methods                   */
/***************************************************************/
static PyObject *
xmipp_str2Label(PyObject *obj, PyObject *args)
{
    char * str;
    if (PyArg_ParseTuple(args, "s", &str))
        return Py_BuildValue("i", (int) MDL::str2Label(str));
    return NULL;
}

static PyObject *
xmipp_label2Str(PyObject *obj, PyObject *args)
{
    int label;
    if (PyArg_ParseTuple(args, "i", &label))
    {
        String labelStr = MDL::label2Str((MDLabel) label);
        return PyString_FromString(labelStr.c_str());
    }
    return NULL;
}

static PyObject *
xmipp_colorStr(PyObject *obj, PyObject *args)
{
    char *str;
    int color;
    if (PyArg_ParseTuple(args, "is", &color, &str))
    {
        String labelStr = colorString(str, color);
        return PyString_FromString(labelStr.c_str());
    }
    return NULL;
}

static PyObject *
xmipp_labelType(PyObject *obj, PyObject *args)
{
    PyObject * input;
    if (PyArg_ParseTuple(args, "O", &input))
    {
        if (PyString_Check(input))
            return Py_BuildValue("i",
                                 (int) MDL::labelType(PyString_AsString(input)));
        else if (PyInt_Check(input))
            return Py_BuildValue("i",
                                 (int) MDL::labelType((MDLabel) PyInt_AsLong(input)));
        else
            PyErr_SetString(PyExc_TypeError,
                            "labelType: Only int or string are allowed as input");
    }
    return NULL;
}

/* isInStack */
static PyObject *
xmipp_isValidLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    char *str;
    int label;
    if (PyArg_ParseTuple(args, "s", &str))
        label = MDL::str2Label(str);
    else if (PyArg_ParseTuple(args, "i", &label))
        ;
    else
        return NULL;
    if (MDL::isValidLabel((MDLabel) label))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* Methods for constructing concrete queries */

/* Helper function to create relational queries */
static PyObject *
createMDValueRelational(PyObject *args, int op)
{
    int label, limit = -1, offset = 0, orderLabel = (int) MDL_OBJID;
    PyObject *pyValue; //Only used to skip label and value

    if ((op == -1 && PyArg_ParseTuple(args, "iO|iiii", &label, &pyValue, &op,
                                      &limit, &offset, &orderLabel)) || PyArg_ParseTuple(args, "iO|iii",
                                              &label, &pyValue, &limit, &offset, &orderLabel))
    {
        MDObject * object = createMDObject(label, pyValue);
        if (!object)
            return NULL;
        MDQueryObject * pyQuery = PyObject_New(MDQueryObject, &MDQueryType);
        pyQuery->query = new MDValueRelational(*object, (RelationalOp) op,
                                               limit, offset, (MDLabel) orderLabel);
        delete object;
        return (PyObject *) pyQuery;
    }
    return NULL;
}
/* MDValue Relational */
static PyObject *
xmipp_MDValueRelational(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, -1);
}
/* MDValueEQ */
static PyObject *
xmipp_MDValueEQ(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, EQ);
}
/* MDValueEQ */
static PyObject *
xmipp_MDValueNE(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, NE);
}
/* MDValueLT */
static PyObject *
xmipp_MDValueLT(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, LT);
}
/* MDValueLE */
static PyObject *
xmipp_MDValueLE(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, LE);
}
/* MDValueLT */
static PyObject *
xmipp_MDValueGT(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, GT);
}
/* MDValueLE */
static PyObject *
xmipp_MDValueGE(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, GE);
}
/* MDValueRange */
static PyObject *
xmipp_MDValueRange(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label, limit = -1, offset = 0, orderLabel = (int) MDL_OBJID;
    PyObject *pyValue1, *pyValue2; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "iOO|iii", &label, &pyValue1, &pyValue2, &limit,
                         &offset, &orderLabel))
    {
        MDObject * object1 = createMDObject(label, pyValue1);
        MDObject * object2 = createMDObject(label, pyValue2);
        if (!object1 || !object2)
            return NULL;
        MDQueryObject * pyQuery = PyObject_New(MDQueryObject, &MDQueryType);
        pyQuery->query = new MDValueRange(*object1, *object2, limit, offset,
                                          (MDLabel) orderLabel);
        delete object1;
        delete object2;
        return (PyObject *) pyQuery;
    }
    return NULL;
}
/* createEmptyFile */
static PyObject *
xmipp_createEmptyFile(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyValue1, *pyValue2; //Only used to skip label and value

    int Xdim,Ydim,Zdim;
    size_t Ndim;
    Zdim=1;
    Ndim=1;
    PyObject * input;
    if (PyArg_ParseTuple(args, "Oii|ii", &input, &Xdim, &Ydim, &Zdim,
                         &Ndim))
    {
        createEmptyFile(PyString_AsString(input),Xdim,Ydim,Zdim,Ndim,true,WRITE_REPLACE);
        Py_RETURN_NONE;
    }
    return NULL;
}
/* SingleImgSize */
static PyObject *
xmipp_SingleImgSize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "O", &pyValue))
    {
        try
        {

            PyObject * pyStr = PyObject_Str(pyValue);
            char * str = PyString_AsString(pyStr);
            int xdim, ydim, zdim;
            size_t ndim;
            getImageSize(str, xdim, ydim, zdim, ndim);
            Py_DECREF(pyStr);
            return Py_BuildValue("iiik", xdim, ydim, zdim, ndim);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/* ImgSize (from metadata filename)*/
static PyObject *
xmipp_ImgSize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "O", &pyValue))
    {
        try
        {
            MetaData md;
            if (PyString_Check(pyValue))
            {
                char * str = PyString_AsString(pyValue);
                md.read(str);
            }
            else if (FileName_Check(pyValue))
            {
                md.read(FileName_Value(pyValue));
            }
            else if (MetaData_Check(pyValue))
            {
                md = MetaData_Value(pyValue);
            }
            else
            {
                PyErr_SetString(PyXmippError, "Invalid argument: expected String, FileName or MetaData");
                return NULL;
            }
            int xdim, ydim, zdim;
            size_t ndim;
            getImageSize(md, xdim, ydim, zdim, ndim);
            return Py_BuildValue("iiik", xdim, ydim, zdim, ndim);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}/* ImgSize (from metadata filename)*/

static PyObject *
xmipp_ImgCompare(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *filename1, *filename2;

    if (PyArg_ParseTuple(args, "OO", &filename1, &filename2))
    {
        try
        {
            PyObject * pyStr1 = PyObject_Str(filename1);
            PyObject * pyStr2 = PyObject_Str(filename2);
            char * str1 = PyString_AsString(pyStr1);
            char * str2 = PyString_AsString(pyStr2);
            bool result = compareImage(str1, str2);
            Py_DECREF(pyStr1);
            Py_DECREF(pyStr2);
            if (result)
                Py_RETURN_TRUE;
            else
                Py_RETURN_FALSE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

static PyObject *
xmipp_compareTwoFiles(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *filename1, *filename2;
    size_t offset=0;
    if (PyArg_ParseTuple(args, "OO|i", &filename1, &filename2, &offset))
    {
        try
        {
            PyObject * pyStr1 = PyObject_Str(filename1);
            PyObject * pyStr2 = PyObject_Str(filename2);
            char * str1 = PyString_AsString(pyStr1);
            char * str2 = PyString_AsString(pyStr2);
            bool result = compareTwoFiles(str1, str2, offset);
            Py_DECREF(pyStr1);
            Py_DECREF(pyStr2);
            if (result)
                Py_RETURN_TRUE;
            else
                Py_RETURN_FALSE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/***************************************************************/
/*                   Some specific utility functions           */
/***************************************************************/
/* readMetaDataWithTwoPossibleImages */
static PyObject *
xmipp_readMetaDataWithTwoPossibleImages(PyObject *obj, PyObject *args,
                                        PyObject *kwargs)
{
    PyObject *pyStr, *pyMd; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "OO", &pyStr, &pyMd))
    {
        try
        {
            if (!MetaData_Check(pyMd))
                PyErr_SetString(PyExc_TypeError,
                                "Expected MetaData as second argument");
            else
            {
                if (PyString_Check(pyStr))
                    readMetaDataWithTwoPossibleImages(PyString_AsString(pyStr),
                                                      MetaData_Value(pyMd));
                else if (FileName_Check(pyStr))
                    readMetaDataWithTwoPossibleImages(FileName_Value(pyStr),
                                                      MetaData_Value(pyMd));
                else
                    PyErr_SetString(PyExc_TypeError,
                                    "Expected string or FileName as first argument");
                Py_RETURN_NONE;
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* substituteOriginalImages */
static PyObject *
xmipp_substituteOriginalImages(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyStrFnOrig, *pyStrFnOut;
    int label, skipFirstBlock;

    if (PyArg_ParseTuple(args, "OOOii", &pyStrFn, &pyStrFnOrig, &pyStrFnOut,
                         &label, &skipFirstBlock))
    {
        try
        {
            FileName fn, fnOrig, fnOut;
            if (PyString_Check(pyStrFn))
                fn = PyString_AsString(pyStrFn);
            else if (FileName_Check(pyStrFn))
                fn = FileName_Value(pyStrFn);
            else
                PyErr_SetString(PyExc_TypeError,
                                "Expected string or FileName as first argument");

            if (PyString_Check(pyStrFnOrig))
                fnOrig = PyString_AsString(pyStrFnOrig);
            else if (FileName_Check(pyStrFnOrig))
                fnOrig = FileName_Value(pyStrFnOrig);
            else
                PyErr_SetString(PyExc_TypeError,
                                "Expected string or FileName as second argument");

            if (PyString_Check(pyStrFnOut))
                fnOut = PyString_AsString(pyStrFnOut);
            else if (FileName_Check(pyStrFnOut))
                fnOut = FileName_Value(pyStrFnOut);
            else
                PyErr_SetString(PyExc_TypeError,
                                "Expected string or FileName as third argument");

            substituteOriginalImages(fn, fnOrig, fnOut, (MDLabel) label,
                                     (bool) skipFirstBlock);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

bool validateInputImageString(PyObject * pyImage, PyObject *pyStrFn, FileName &fn)
{
    if (!Image_Check(pyImage))
    {
        PyErr_SetString(PyExc_TypeError,
                        "bad argument: Expected Image as first argument");
        return false;
    }
    if (PyString_Check(pyStrFn))
        fn = PyString_AsString(pyStrFn);
    else if (FileName_Check(pyStrFn))
        fn = FileName_Value(pyStrFn);
    else
    {
        PyErr_SetString(PyExc_TypeError,
                        "bad argument:Expected string or FileName as second argument");
        return false;
    }
    return true;
}

/* calculate enhanced psd and return preview*/
static PyObject *
xmipp_fastEstimateEnhancedPSD(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyImage;
    //ImageObject *pyImage;
    double downsampling;
    int dim, Nthreads;
    FileName fn;

    if (PyArg_ParseTuple(args, "OOdii", &pyImage, &pyStrFn, &downsampling, &dim, &Nthreads))
    {
        try
        {
            if (validateInputImageString(pyImage, pyStrFn, fn))
            {
                MultidimArray<double> data;
                fastEstimateEnhancedPSD(fn, downsampling, data, Nthreads);
                selfScaleToSize(LINEAR, data, dim, dim);
                Image_Value(pyImage).setDatatype(Double);
                Image_Value(pyImage).data->setImage(data);
                Py_RETURN_NONE;
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/** Some helper macros repeated in filter functions*/
#define FILTER_TRY()\
try {\
if (validateInputImageString(pyImage, pyStrFn, fn)) {\
Image<double> img;\
img.read(fn);\
MultidimArray<double> &data = MULTIDIM_ARRAY(img);\
ArrayDim idim;\
data.getDimensions(idim);

#define FILTER_CATCH()\
int w = dim, h = dim, &x = idim.xdim, &y = idim.ydim;\
double ddim = dim;\
if (x > y) h = y * (dim/x);\
else if (y > x)\
  w = x * (dim/y);\
selfScaleToSize(LINEAR, data, w, h);\
Image_Value(pyImage).setDatatype(Double);\
MULTIDIM_ARRAY_GENERIC(Image_Value(pyImage)).setImage(data);\
Py_RETURN_NONE;\
}} catch (XmippError &xe)\
{ PyErr_SetString(PyXmippError, xe.msg.c_str());}\


/* calculate enhanced psd and return preview
* used for protocol preprocess_particles*/
static PyObject *
xmipp_bandPassFilter(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyImage;
    double w1, w2, raised_w;
    int dim;
    FileName fn;

    if (PyArg_ParseTuple(args, "OOdddi", &pyImage, &pyStrFn, &w1, &w2, &raised_w, &dim))
    {
        FILTER_TRY()
        bandpassFilter(data, w1, w2, raised_w);
        FILTER_CATCH()
    }
    return NULL;
}

/* calculate enhanced psd and return preview
 * used for protocol preprocess_particles*/
static PyObject *
xmipp_gaussianFilter(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyImage;
    double freqSigma;
    int dim;
    FileName fn;

    if (PyArg_ParseTuple(args, "OOdi", &pyImage, &pyStrFn, &freqSigma, &dim))
    {
        FILTER_TRY()
        gaussianFilter(data, freqSigma);
        FILTER_CATCH()
    }
    return NULL;
}

/* calculate enhanced psd and return preview
 * used for protocol preprocess_particles*/
static PyObject *
xmipp_badPixelFilter(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyImage;
    double factor;
    int dim;
    FileName fn;

    if (PyArg_ParseTuple(args, "OOdi", &pyImage, &pyStrFn, &factor, &dim))
    {
        FILTER_TRY()
        BadPixelFilter filter;
        filter.type = BadPixelFilter::OUTLIER;
        filter.factor = factor;
        filter.apply(data);
        FILTER_CATCH()
    }
    return NULL;
}

static PyMethodDef
xmipp_methods[] =
    {

        { "getBlocksInMetaDataFile",
          xmipp_getBlocksInMetaDataFile, METH_VARARGS,
          "return list with metadata blocks in a file" },
        { "label2Str", xmipp_label2Str, METH_VARARGS,
          "Convert MDLabel to string" },
        { "colorStr", xmipp_colorStr, METH_VARARGS,
          "Create a string with color characters sequence for print in console" },
        { "labelType", xmipp_labelType, METH_VARARGS,
          "Return the type of a label" },
        { "str2Label", xmipp_str2Label, METH_VARARGS,
          "Convert an string to MDLabel" },
        { "isValidLabel", (PyCFunction) xmipp_isValidLabel,
          METH_VARARGS,
          "Check if the label is a valid one" },
        { "MDValueRelational",
          (PyCFunction) xmipp_MDValueRelational,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueEQ", (PyCFunction) xmipp_MDValueEQ,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueNE", (PyCFunction) xmipp_MDValueNE,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueLT", (PyCFunction) xmipp_MDValueLT,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueLE", (PyCFunction) xmipp_MDValueLE,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueGT", (PyCFunction) xmipp_MDValueGT,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueGE", (PyCFunction) xmipp_MDValueGE,
          METH_VARARGS, "Construct a relational query" },
        { "MDValueRange", (PyCFunction) xmipp_MDValueRange,
          METH_VARARGS, "Construct a range query" },
        { "createEmptyFile", (PyCFunction) xmipp_createEmptyFile,
          METH_VARARGS, "create empty stack (speed up things)" },
        { "SingleImgSize", (PyCFunction) xmipp_SingleImgSize,
          METH_VARARGS, "Get image dimensions" },
        { "ImgSize", (PyCFunction) xmipp_ImgSize, METH_VARARGS,
          "Get image dimensions of first metadata entry" },
        { "ImgCompare", (PyCFunction) xmipp_ImgCompare,  METH_VARARGS,
          "return true if both files are identical" },
        { "compareTwoFiles", (PyCFunction) xmipp_compareTwoFiles, METH_VARARGS,
          "return true if both files are identical" },
        { "readMetaDataWithTwoPossibleImages", (PyCFunction) xmipp_readMetaDataWithTwoPossibleImages, METH_VARARGS,
          "Read a 1 or two column list of micrographs" },
        { "substituteOriginalImages", (PyCFunction) xmipp_substituteOriginalImages, METH_VARARGS,
          "Substitute the original images into a given column of a metadata" },
        { "fastEstimateEnhancedPSD", (PyCFunction) xmipp_fastEstimateEnhancedPSD, METH_VARARGS,
          "Utility function to calculate PSD preview" },
        { "bandPassFilter", (PyCFunction) xmipp_bandPassFilter, METH_VARARGS,
          "Utility function to apply bandpass filter" },
        { "gaussianFilter", (PyCFunction) xmipp_gaussianFilter, METH_VARARGS,
          "Utility function to apply gaussian filter in Fourier space" },
        { "badPixelFilter", (PyCFunction) xmipp_badPixelFilter, METH_VARARGS,
          "Bad pixel filter" },

        { NULL } /* Sentinel */
    };

void addIntConstant(PyObject * dict, const char * name, const long &value)
{
    PyObject * pyValue = PyInt_FromLong(value);
    PyDict_SetItemString(dict, name, pyValue);
    Py_DECREF(pyValue);
}

PyMODINIT_FUNC initxmipp(void)
{
    //Initialize module variable
    PyObject* module;
    module = Py_InitModule3("xmipp", xmipp_methods,
                            "Xmipp module as a Python extension.");

    // Add FileName type
    if (PyType_Ready(&FileNameType) < 0)
        return;
    Py_INCREF(&FileNameType);
    PyModule_AddObject(module, "FileName", (PyObject *) &FileNameType);

    // Add Image type
    if (PyType_Ready(&ImageType) < 0)
        return;
    Py_INCREF(&ImageType);
    PyModule_AddObject(module, "Image", (PyObject *) &ImageType);

    // Add Program type
    if (PyType_Ready(&ProgramType) < 0)
        return;
    Py_INCREF(&ProgramType);
    PyModule_AddObject(module, "Program", (PyObject *) &ProgramType);

    // Add MDQuery type, as no specific new is create, use the generic one
    MDQueryType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MDQueryType) < 0)
        return;
    Py_INCREF(&MDQueryType);
    PyModule_AddObject(module, "MDQuery", (PyObject *) &MDQueryType);

    //Add MetaData type
    if (PyType_Ready(&MetaDataType) < 0)
        return;
    Py_INCREF(&MetaDataType);
    PyModule_AddObject(module, "MetaData", (PyObject *) &MetaDataType);

    //Add PyXmippError
    PyXmippError = PyErr_NewException("xmipp.XmippError", NULL, NULL);
    Py_INCREF(PyXmippError);
    PyModule_AddObject(module, "XmippError", PyXmippError);

    PyObject * dict = PyModule_GetDict(module);

    //Add constants
    addIntConstant(dict, "AGGR_COUNT", (long) AGGR_COUNT);
    addIntConstant(dict, "AGGR_MAX", (long) AGGR_MAX);
    addIntConstant(dict, "AGGR_MIN", (long) AGGR_MIN);
    addIntConstant(dict, "AGGR_SUM", (long) AGGR_SUM);
    addIntConstant(dict, "AGGR_AVG", (long) AGGR_AVG);
    addIntConstant(dict, "UNION", (long) UNION);
    addIntConstant(dict, "UNION_DISTINCT", (long) UNION_DISTINCT);
    addIntConstant(dict, "INTERSECTION", (long) INTERSECTION);
    addIntConstant(dict, "SUBSTRACTION", (long) SUBSTRACTION);
    addIntConstant(dict, "INNER_JOIN", (long) INNER_JOIN);
    addIntConstant(dict, "LEFT_JOIN", (long) LEFT_JOIN);
    addIntConstant(dict, "NATURAL_JOIN", (long) NATURAL_JOIN);
    addIntConstant(dict, "OUTER_JOIN", (long) OUTER_JOIN);
    addIntConstant(dict, "INNER", (long) INNER);
    addIntConstant(dict, "LEFT", (long) LEFT);
    addIntConstant(dict, "OUTER", (long) OUTER);
    addIntConstant(dict, "NATURAL", (long) NATURAL);
    addIntConstant(dict, "EQ", (long) EQ);
    addIntConstant(dict, "NE", (long) NE);
    addIntConstant(dict, "GT", (long) GT);
    addIntConstant(dict, "LT", (long) LT);
    addIntConstant(dict, "GE", (long) GE);
    addIntConstant(dict, "LE", (long) LE);
    addIntConstant(dict, "MD_OVERWRITE", (long) MD_OVERWRITE);
    addIntConstant(dict, "MD_APPEND", (long) MD_APPEND);
    addIntConstant(dict, "MDL_UNDEFINED", (long) MDL_UNDEFINED);
    addIntConstant(dict, "MDL_FIRST_LABEL", (long) MDL_FIRST_LABEL);
    //Metadata Labels
    addIntConstant(dict, "MDL_OBJID", (size_t) MDL_OBJID);
    addIntConstant(dict, "MDL_ANGLE_COMPARISON", (long) MDL_ANGLE_COMPARISON);
    addIntConstant(dict, "MDL_ANGLEPSI2", (long) MDL_ANGLEPSI2);
    addIntConstant(dict, "MDL_ANGLEPSI", (long) MDL_ANGLEPSI);
    addIntConstant(dict, "MDL_ANGLEROT2", (long) MDL_ANGLEROT2);
    addIntConstant(dict, "MDL_ANGLEROT", (long) MDL_ANGLEROT);
    addIntConstant(dict, "MDL_ANGLETILT2", (long) MDL_ANGLETILT2);
    addIntConstant(dict, "MDL_ANGLETILT", (long) MDL_ANGLETILT);
    addIntConstant(dict, "MDL_ASSOCIATED_IMAGE1", (long) MDL_ASSOCIATED_IMAGE1);
    addIntConstant(dict, "MDL_ASSOCIATED_IMAGE2", (long) MDL_ASSOCIATED_IMAGE2);
    addIntConstant(dict, "MDL_ASSOCIATED_IMAGE3", (long) MDL_ASSOCIATED_IMAGE3);
    addIntConstant(dict, "MDL_ASSOCIATED_IMAGE4", (long) MDL_ASSOCIATED_IMAGE4);
    addIntConstant(dict, "MDL_ASSOCIATED_IMAGE5", (long) MDL_ASSOCIATED_IMAGE5);
    addIntConstant(dict, "MDL_AVG", (long) MDL_AVG);
    addIntConstant(dict, "MDL_AZIMUTALANGLE", (long) MDL_AZIMUTALANGLE);
    addIntConstant(dict, "MDL_BGMEAN", (long) MDL_BGMEAN);
    addIntConstant(dict, "MDL_BLOCK", (long) MDL_BLOCK);
    addIntConstant(dict, "MDL_CELLX", (long) MDL_CELLX);
    addIntConstant(dict, "MDL_CELLY", (long) MDL_CELLY);
    addIntConstant(dict, "MDL_COMMENT", (long) MDL_COMMENT);
    addIntConstant(dict, "MDL_COST", (long) MDL_COST);
    addIntConstant(dict, "MDL_COUNT", (long) MDL_COUNT);
    addIntConstant(dict, "MDL_CTFINPUTPARAMS", (long) MDL_CTFINPUTPARAMS);
    addIntConstant(dict, "MDL_CTFMODEL", (long) MDL_CTFMODEL);
    addIntConstant(dict, "MDL_CTFMODEL2", (long) MDL_CTFMODEL2);
    addIntConstant(dict, "MDL_CTF_SAMPLING_RATE", (long) MDL_CTF_SAMPLING_RATE);
    addIntConstant(dict, "MDL_CTF_VOLTAGE", (long) MDL_CTF_VOLTAGE);
    addIntConstant(dict, "MDL_CTF_DEFOCUSU", (long) MDL_CTF_DEFOCUSU);
    addIntConstant(dict, "MDL_CTF_DEFOCUSV", (long) MDL_CTF_DEFOCUSV);
    addIntConstant(dict, "MDL_CTF_DEFOCUS_ANGLE", (long) MDL_CTF_DEFOCUS_ANGLE);
    addIntConstant(dict, "MDL_CTF_CS", (long) MDL_CTF_CS);
    addIntConstant(dict, "MDL_CTF_CA", (long) MDL_CTF_CA);
    addIntConstant(dict, "MDL_CTF_GROUP", (long) MDL_CTF_GROUP);
    addIntConstant(dict, "MDL_CTF_ENERGY_LOSS", (long) MDL_CTF_ENERGY_LOSS);
    addIntConstant(dict, "MDL_CTF_LENS_STABILITY",
                   (long) MDL_CTF_LENS_STABILITY);
    addIntConstant(dict, "MDL_CTF_CONVERGENCE_CONE",
                   (long) MDL_CTF_CONVERGENCE_CONE);
    addIntConstant(dict, "MDL_CTF_LONGITUDINAL_DISPLACEMENT",
                   (long) MDL_CTF_LONGITUDINAL_DISPLACEMENT);
    addIntConstant(dict, "MDL_CTF_TRANSVERSAL_DISPLACEMENT",
                   (long) MDL_CTF_TRANSVERSAL_DISPLACEMENT);
    addIntConstant(dict, "MDL_CTF_Q0", (long) MDL_CTF_Q0);
    addIntConstant(dict, "MDL_CTF_K", (long) MDL_CTF_K);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN_K", (long) MDL_CTFBG_GAUSSIAN_K);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN_SIGMAU",
                   (long) MDL_CTFBG_GAUSSIAN_SIGMAU);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN_SIGMAV",
                   (long) MDL_CTFBG_GAUSSIAN_SIGMAV);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN_CU", (long) MDL_CTFBG_GAUSSIAN_CU);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN_CV", (long) MDL_CTFBG_GAUSSIAN_CV);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN_ANGLE",
                   (long) MDL_CTFBG_GAUSSIAN_ANGLE);
    addIntConstant(dict, "MDL_CTFBG_SQRT_K", (long) MDL_CTFBG_SQRT_K);
    addIntConstant(dict, "MDL_CTFBG_SQRT_U", (long) MDL_CTFBG_SQRT_U);
    addIntConstant(dict, "MDL_CTFBG_SQRT_V", (long) MDL_CTFBG_SQRT_V);
    addIntConstant(dict, "MDL_CTFBG_SQRT_ANGLE", (long) MDL_CTFBG_SQRT_ANGLE);
    addIntConstant(dict, "MDL_CTFBG_BASELINE", (long) MDL_CTFBG_BASELINE);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN2_K", (long) MDL_CTFBG_GAUSSIAN2_K);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN2_SIGMAU",
                   (long) MDL_CTFBG_GAUSSIAN2_SIGMAU);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN2_SIGMAV",
                   (long) MDL_CTFBG_GAUSSIAN2_SIGMAV);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN2_CU",
                   (long) MDL_CTFBG_GAUSSIAN2_CU);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN2_CV",
                   (long) MDL_CTFBG_GAUSSIAN2_CV);
    addIntConstant(dict, "MDL_CTFBG_GAUSSIAN2_ANGLE",
                   (long) MDL_CTFBG_GAUSSIAN2_ANGLE);
    addIntConstant(dict, "MDL_CTF_CRITERION_PSDCORRELATION90",
                   (long) MDL_CTF_CRITERION_PSDCORRELATION90);
    addIntConstant(dict, "MDL_CTF_CRITERION_FIRSTZERORATIO",
                   (long) MDL_CTF_CRITERION_FIRSTZERORATIO);
    addIntConstant(dict, "MDL_CTF_CRITERION_FIRSTZEROAVG",
                   (long) MDL_CTF_CRITERION_FIRSTZEROAVG);
    addIntConstant(dict, "MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT",
                   (long) MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT);
    addIntConstant(dict, "MDL_CTF_CRITERION_DAMPING",
                   (long) MDL_CTF_CRITERION_DAMPING);
    addIntConstant(dict, "MDL_CTF_CRITERION_PSDRADIALINTEGRAL",
                   (long) MDL_CTF_CRITERION_PSDRADIALINTEGRAL);
    addIntConstant(dict, "MDL_CTF_CRITERION_FITTINGSCORE",
                   (long) MDL_CTF_CRITERION_FITTINGSCORE);
    addIntConstant(dict, "MDL_CTF_CRITERION_FITTINGCORR13",
                   (long) MDL_CTF_CRITERION_FITTINGCORR13);
    addIntConstant(dict, "MDL_CTF_CRITERION_PSDVARIANCE",
                   (long) MDL_CTF_CRITERION_PSDVARIANCE);
    addIntConstant(dict, "MDL_CTF_CRITERION_PSDPCA1VARIANCE",
                   (long) MDL_CTF_CRITERION_PSDPCA1VARIANCE);
    addIntConstant(dict, "MDL_CTF_CRITERION_PSDPCARUNSTEST",
                   (long) MDL_CTF_CRITERION_PSDPCARUNSTEST);
    addIntConstant(dict, "MDL_DATATYPE", (long) MDL_DATATYPE);
    addIntConstant(dict, "MDL_DEFGROUP", (long) MDL_DEFGROUP);
    addIntConstant(dict, "MDL_DM3_IDTAG", (long) MDL_DM3_IDTAG);
    addIntConstant(dict, "MDL_DM3_NODEID", (long) MDL_DM3_NODEID);
    addIntConstant(dict, "MDL_DM3_NUMBER_TYPE", (long) MDL_DM3_NUMBER_TYPE);
    addIntConstant(dict, "MDL_DM3_PARENTID", (long) MDL_DM3_PARENTID);
    addIntConstant(dict, "MDL_DM3_TAGCLASS", (long) MDL_DM3_TAGCLASS);
    addIntConstant(dict, "MDL_DM3_TAGNAME", (long) MDL_DM3_TAGNAME);
    addIntConstant(dict, "MDL_DM3_SIZE", (long) MDL_DM3_SIZE);
    addIntConstant(dict, "MDL_DM3_VALUE", (long) MDL_DM3_VALUE);
    addIntConstant(dict, "MDL_ENABLED", (long) MDL_ENABLED);
    addIntConstant(dict, "MDL_FLIP", (long) MDL_FLIP);
    addIntConstant(dict, "MDL_IMAGE_CLASS_COUNT", (long) MDL_IMAGE_CLASS_COUNT);
    addIntConstant(dict, "MDL_IMAGE_CLASS_GROUP", (long) MDL_IMAGE_CLASS_GROUP);
    addIntConstant(dict, "MDL_IMAGE_CLASS", (long) MDL_IMAGE_CLASS);
    addIntConstant(dict, "MDL_IMAGE", (long) MDL_IMAGE);
    addIntConstant(dict, "MDL_IMAGE_ORIGINAL", (long) MDL_IMAGE_ORIGINAL);
    addIntConstant(dict, "MDL_IMAGE_TILTED", (long) MDL_IMAGE_TILTED);
    addIntConstant(dict, "MDL_IMGMD", (long) MDL_IMGMD);
    addIntConstant(dict, "MDL_INTSCALE", (long) MDL_INTSCALE);
    addIntConstant(dict, "MDL_ITER", (long) MDL_ITER);
    addIntConstant(dict, "MDL_K", (long) MDL_K);
    addIntConstant(dict, "MDL_KSTEST", (long) MDL_KSTEST);
    addIntConstant(dict, "MDL_LL", (long) MDL_LL);
    addIntConstant(dict, "MDL_MAGNIFICATION", (long) MDL_MAGNIFICATION);
    addIntConstant(dict, "MDL_MASK", (long) MDL_MASK);
    addIntConstant(dict, "MDL_MAXCC", (long) MDL_MAXCC);
    addIntConstant(dict, "MDL_MAX", (long) MDL_MAX);
    addIntConstant(dict, "MDL_MICROGRAPH", (long) MDL_MICROGRAPH);
    addIntConstant(dict, "MDL_MICROGRAPH_TILTED", (long) MDL_MICROGRAPH_TILTED);
    addIntConstant(dict, "MDL_MIN", (long) MDL_MIN);
    addIntConstant(dict, "MDL_MIRRORFRAC", (long) MDL_MIRRORFRAC);
    addIntConstant(dict, "MDL_MISSINGREGION_NR", (long) MDL_MISSINGREGION_NR);
    addIntConstant(dict, "MDL_MISSINGREGION_TYPE",(long) MDL_MISSINGREGION_TYPE);
    addIntConstant(dict, "MDL_MISSINGREGION_THY0",(long) MDL_MISSINGREGION_THY0);
    addIntConstant(dict, "MDL_MISSINGREGION_THYF",(long) MDL_MISSINGREGION_THYF);
    addIntConstant(dict, "MDL_MISSINGREGION_THX0",(long) MDL_MISSINGREGION_THX0);
    addIntConstant(dict, "MDL_MISSINGREGION_THXF",(long) MDL_MISSINGREGION_THXF);
    addIntConstant(dict, "MDL_MODELFRAC", (long) MDL_MODELFRAC);
    addIntConstant(dict, "MDL_NEIGHBORS", (long)LABEL_VECTOR_LONG);
    addIntConstant(dict, "MDL_NEIGHBORHOOD_RADIUS", (long)LABEL_DOUBLE);
    addIntConstant(dict, "MDL_NMA", (long) MDL_NMA);
    addIntConstant(dict, "MDL_ORDER", (long) MDL_ORDER);
    addIntConstant(dict, "MDL_ORIGINX", (long) MDL_ORIGINX);
    addIntConstant(dict, "MDL_ORIGINY", (long) MDL_ORIGINY);
    addIntConstant(dict, "MDL_ORIGINZ", (long) MDL_ORIGINZ);
    addIntConstant(dict, "MDL_PICKING_FAMILY", (long) MDL_PICKING_FAMILY);
    addIntConstant(dict, "MDL_PICKING_COLOR", (long) MDL_PICKING_COLOR);
    addIntConstant(dict, "MDL_PICKING_PARTICLE_SIZE", (long) MDL_PICKING_PARTICLE_SIZE);
    addIntConstant(dict, "MDL_PICKING_FAMILY_STATE", (long) MDL_PICKING_FAMILY_STATE);
    addIntConstant(dict, "MDL_PICKING_MICROGRAPH_FAMILY_STATE", (long) MDL_PICKING_MICROGRAPH_FAMILY_STATE);
    addIntConstant(dict, "MDL_PMAX", (long) MDL_PMAX);
    addIntConstant(dict, "MDL_PSD", (long) MDL_PSD);
    addIntConstant(dict, "MDL_PSD_ENHANCED", (long) MDL_PSD_ENHANCED);
    addIntConstant(dict, "MDL_RANDOMSEED", (long) MDL_RANDOMSEED);
    addIntConstant(dict, "MDL_REF3D", (long) MDL_REF3D);
    addIntConstant(dict, "MDL_REF", (long) MDL_REF);
    addIntConstant(dict, "MDL_REFMD", (long) MDL_REFMD);
    addIntConstant(dict, "MDL_RESOLUTION_DPR", (long) MDL_RESOLUTION_DPR);
    addIntConstant(dict, "MDL_RESOLUTION_ERRORL2",(long) MDL_RESOLUTION_ERRORL2);
    addIntConstant(dict, "MDL_RESOLUTION_FRC", (long) MDL_RESOLUTION_FRC);
    addIntConstant(dict, "MDL_RESOLUTION_FRCRANDOMNOISE",(long) MDL_RESOLUTION_FRCRANDOMNOISE);
    addIntConstant(dict, "MDL_RESOLUTION_FREQ", (long) MDL_RESOLUTION_FREQ);
    addIntConstant(dict, "MDL_RESOLUTION_FREQREAL",(long) MDL_RESOLUTION_FREQREAL);
    addIntConstant(dict, "MDL_SAMPLINGRATE", (long) MDL_SAMPLINGRATE);
    addIntConstant(dict, "MDL_SAMPLINGRATEX", (long) MDL_SAMPLINGRATEX);
    addIntConstant(dict, "MDL_SAMPLINGRATEY", (long) MDL_SAMPLINGRATEY);
    addIntConstant(dict, "MDL_SAMPLINGRATEZ", (long) MDL_SAMPLINGRATEZ);
    addIntConstant(dict, "MDL_SCALE", (long) MDL_SCALE);
    addIntConstant(dict, "MDL_SELFILE", (long) MDL_SELFILE);
    addIntConstant(dict, "MDL_SERIE", (long) MDL_SERIE);
    addIntConstant(dict, "MDL_SHIFTX", (long) MDL_SHIFTX);
    addIntConstant(dict, "MDL_SHIFTY", (long) MDL_SHIFTY);
    addIntConstant(dict, "MDL_SHIFTZ", (long) MDL_SHIFTZ);
    addIntConstant(dict, "MDL_SHIFT_CRYSTALX", (long) MDL_SHIFT_CRYSTALX);
    addIntConstant(dict, "MDL_SHIFT_CRYSTALY", (long) MDL_SHIFT_CRYSTALY);
    addIntConstant(dict, "MDL_SHIFT_CRYSTALZ", (long) MDL_SHIFT_CRYSTALZ);
    addIntConstant(dict, "MDL_SIGMANOISE", (long) MDL_SIGMANOISE);
    addIntConstant(dict, "MDL_SIGMAOFFSET", (long) MDL_SIGMAOFFSET);
    addIntConstant(dict, "MDL_SIGNALCHANGE", (long) MDL_SIGNALCHANGE);
    addIntConstant(dict, "MDL_SPHERICALABERRATION",(long) MDL_SPHERICALABERRATION);
    addIntConstant(dict, "MDL_STDDEV", (long) MDL_STDDEV);
    addIntConstant(dict, "MDL_SUM", (long) MDL_SUM);
    addIntConstant(dict, "MDL_SUMWEIGHT", (long) MDL_SUMWEIGHT);
    addIntConstant(dict, "MDL_SYMNO", (long) MDL_SYMNO);
    addIntConstant(dict, "MDL_TRANSFORMATIONMTRIX",(long) MDL_TRANSFORMATIONMTRIX);
    addIntConstant(dict, "MDL_VOLTAGE", (long) MDL_VOLTAGE);
    addIntConstant(dict, "MDL_WEIGHT", (long) MDL_WEIGHT);
    addIntConstant(dict, "MDL_WROBUST", (long) MDL_WROBUST);
    addIntConstant(dict, "MDL_XINT", (long) MDL_XINT);
    addIntConstant(dict, "MDL_XINTTILT", (long) MDL_XINTTILT);
    addIntConstant(dict, "MDL_X", (long) MDL_X);
    addIntConstant(dict, "MDL_YINT", (long) MDL_YINT);
    addIntConstant(dict, "MDL_YINTTILT", (long) MDL_YINTTILT);
    addIntConstant(dict, "MDL_Y", (long) MDL_Y);
    addIntConstant(dict, "MDL_ZINT", (long) MDL_ZINT);
    addIntConstant(dict, "MDL_Z", (long) MDL_Z);
    addIntConstant(dict, "MDL_ZSCORE", (long) MDL_ZSCORE);
    addIntConstant(dict, "MDL_LAST_LABEL", (long) MDL_LAST_LABEL);
    addIntConstant(dict, "LABEL_NOTYPE", (long) LABEL_NOTYPE);
    addIntConstant(dict, "LABEL_INT", (long) LABEL_INT);
    addIntConstant(dict, "LABEL_BOOL", (long) LABEL_BOOL);
    addIntConstant(dict, "LABEL_DOUBLE", (long) LABEL_DOUBLE);
    addIntConstant(dict, "LABEL_FLOAT", (long) LABEL_FLOAT);
    addIntConstant(dict, "LABEL_STRING", (long) LABEL_STRING);
    addIntConstant(dict, "LABEL_VECTOR", (long) LABEL_VECTOR);
    addIntConstant(dict, "LABEL_LONG", (long) LABEL_LONG);
    addIntConstant(dict, "_NONE", (long) _NONE);
    addIntConstant(dict, "HEADER", (long) HEADER);
    addIntConstant(dict, "_HEADER_ALL", (long) _HEADER_ALL);
    addIntConstant(dict, "DATA", (long) DATA);
    addIntConstant(dict, "_DATA_ALL", (long) _DATA_ALL);
    addIntConstant(dict, "WRAP", (long) WRAP);
    addIntConstant(dict, "ALL_IMAGES", (long) ALL_IMAGES);
    addIntConstant(dict, "FILENAMENUMBERLENGTH", (long) FILENAMENUMBERLENGTH);
    addIntConstant(dict, "XMIPP_BLACK", (long) BLACK);
    addIntConstant(dict, "XMIPP_RED", (long) RED);
    addIntConstant(dict, "XMIPP_GREEN", (long) GREEN);
    addIntConstant(dict, "XMIPP_YELLOW", (long) YELLOW);
    addIntConstant(dict, "XMIPP_BLUE", (long) BLUE);
    addIntConstant(dict, "XMIPP_MAGENTA", (long) MAGENTA);
    addIntConstant(dict, "XMIPP_CYAN", (long) CYAN);
    addIntConstant(dict, "XMIPP_WHITE", (long) WHITE);
    addIntConstant(dict, "XMIPP_RND_UNIFORM", (long) RND_Uniform);
    addIntConstant(dict, "XMIPP_RND_GAUSSIAN", (long) RND_Gaussian);


}
