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
#include "data/ctf.h"
#include "data/xmipp_image_macros.h"
PyObject * PyXmippError;


/***************************************************************/
/*                            Global methods                   */
/***************************************************************/
PyObject *
xmipp_str2Label(PyObject *obj, PyObject *args)
{
    char * str;
    if (PyArg_ParseTuple(args, "s", &str))
        return Py_BuildValue("i", (int) MDL::str2Label(str));
    return NULL;
}

PyObject *
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

PyObject *
xmipp_colorStr(PyObject *obj, PyObject *args)
{
    char *str;
    int color;
    int attrib = BRIGHT;
    if (PyArg_ParseTuple(args, "is|i", &color, &str, &attrib))
    {
        String labelStr = colorString(str, color, attrib);
        return PyString_FromString(labelStr.c_str());
    }
    return NULL;
}

PyObject *
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

PyObject *
xmipp_labelHasTag(PyObject *obj, PyObject *args)
{
    PyObject * input;
    int tag;

    if (PyArg_ParseTuple(args, "Oi", &input, &tag))
    {
        MDLabel label = MDL_UNDEFINED;

        if (PyString_Check(input))
            label = MDL::str2Label(PyString_AsString(input));
        else if (PyInt_Check(input))
            label = (MDLabel) PyInt_AsLong(input);

        if (label != MDL_UNDEFINED)
        {
            if (MDL::hasTag(label, tag))
                Py_RETURN_TRUE;
            else
                Py_RETURN_FALSE;
        }

        PyErr_SetString(PyExc_TypeError,
                        "labelHasTag: Input label should be int or string");
    }
    return NULL;
}

PyObject *
xmipp_labelIsImage(PyObject *obj, PyObject *args)
{
    PyObject * input;
    int tag = TAGLABEL_IMAGE;

    if (PyArg_ParseTuple(args, "O", &input))
    {
        MDLabel label = MDL_UNDEFINED;

        if (PyString_Check(input))
            label = MDL::str2Label(PyString_AsString(input));
        else if (PyInt_Check(input))
            label = (MDLabel) PyInt_AsLong(input);

        if (label != MDL_UNDEFINED)
        {
            if (MDL::hasTag(label, tag))
                Py_RETURN_TRUE;
            else
                Py_RETURN_FALSE;
        }

        PyErr_SetString(PyExc_TypeError,
                        "labelIsImage: Input label should be int or string");
    }
    return NULL;
}

/* isInStack */
PyObject *
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

/* createEmptyFile */
//ROB: argument of this python function do not much arguments of C funcion Ndim does not exists in C
PyObject *
xmipp_createEmptyFile(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int Xdim,Ydim,Zdim;
    size_t Ndim;
    Zdim=1;
    Ndim=1;
    DataType dataType = DT_Float;

    PyObject * input;
    if (PyArg_ParseTuple(args, "Oii|iii", &input, &Xdim, &Ydim, &Zdim,
                         &Ndim, &dataType))
    {
    try
        {
        String inputStr = PyString_AsString(input);
        inputStr += "%";
        inputStr += datatype2Str(dataType);
        createEmptyFile(inputStr, Xdim, Ydim, Zdim, Ndim, true, WRITE_REPLACE);
//        createEmptyFile(PyString_AsString(input),Xdim,Ydim,Zdim,APPEND_IMAGE,true,WRITE_REPLACE);
        Py_RETURN_NONE;
        }
     catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/* getImageSize */
PyObject *
xmipp_getImageSize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "O", &pyValue))
    {
        try
        {

            PyObject * pyStr = PyObject_Str(pyValue);
            char * str = PyString_AsString(pyStr);
            size_t xdim, ydim, zdim, ndim;
            getImageSize(str, xdim, ydim, zdim, ndim);
            Py_DECREF(pyStr);
            return Py_BuildValue("kkkk", xdim, ydim, zdim, ndim);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

PyObject * xmipp_MetaDataInfo(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyValue; //Only used to skip label and value

    if (PyArg_ParseTuple(args, "O", &pyValue))
    {
        try
        {
            MetaData *md = NULL;
            size_t size; //number of elements in the metadata
            bool destroyMd = true;

            if (PyString_Check(pyValue))
            {
                char * str = PyString_AsString(pyValue);
                md = new MetaData();
                md->setMaxRows(1);
                md->read(str);
                size = md->getParsedLines();
            }
            else if (FileName_Check(pyValue))
            {
                md = new MetaData();
                md->setMaxRows(1);
                md->read(FileName_Value(pyValue));
                size = md->getParsedLines();
            }
            else if (MetaData_Check(pyValue))
            {
                md = ((MetaDataObject*)pyValue)->metadata;
                destroyMd = false;
                size = md->size();
            }
            else
            {
                PyErr_SetString(PyXmippError, "Invalid argument: expected String, FileName or MetaData");
                return NULL;
            }
            size_t xdim, ydim, zdim, ndim;
            getImageSize(*md, xdim, ydim, zdim, ndim);

            if (destroyMd)
                delete md;
            return Py_BuildValue("iiikk", xdim, ydim, zdim, ndim, size);
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}/* Metadata info (from metadata filename)*/
/* check if block exists in file*/

PyObject *
xmipp_existsBlockInMetaDataFile(PyObject *obj, PyObject *args, PyObject *kwargs)
{

    PyObject *input = NULL, *pyStr = NULL;
    char *str = NULL, *strBlock = NULL;
    if (PyArg_ParseTuple(args, "O", &input))
    {
        try
        {
            if ((pyStr = PyObject_Str(input)) != NULL )
            {
                str = PyString_AsString(pyStr);
                if (existsBlockInMetaDataFile( (std::string) str))
                    Py_RETURN_TRUE;
                else
                    Py_RETURN_FALSE;
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

    return NULL;
}

PyObject *
xmipp_CheckImageFileSize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *filename;

    if (PyArg_ParseTuple(args, "O", &filename))
    {
        try
        {
            PyObject * pyStr = PyObject_Str(filename);
            char * str = PyString_AsString(pyStr);
            bool result = checkImageFileSize(str);
            Py_DECREF(pyStr);
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

PyObject *
xmipp_CheckImageCorners(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *filename;

    if (PyArg_ParseTuple(args, "O", &filename))
    {
        try
        {
            PyObject * pyStr = PyObject_Str(filename);
            char * str = PyString_AsString(pyStr);
            bool result = checkImageCorners(str);
            Py_DECREF(pyStr);
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

PyObject *
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

PyObject *
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


PyObject *
xmipp_bsoftRemoveLoopBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
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
            bsoftRemoveLoopBlock(str1, str2);
            Py_DECREF(pyStr1);
            Py_DECREF(pyStr2);
			Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}


PyObject *
xmipp_bsoftRestoreLoopBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
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
            bsoftRestoreLoopBlock(str1, str2);
            Py_DECREF(pyStr1);
            Py_DECREF(pyStr2);
			Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

PyObject *
xmipp_compareTwoImageTolerance(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *input1, *input2;
    double tolerance=0.;
    int index1=0, index2=0;

    if (PyArg_ParseTuple(args, "OO|d", &input1, &input2, &tolerance))

    {
        try
        {
            size_t index1 = 0;
            size_t index2 = 0;
            char * fn1, *fn2;

            // If the inputs objects are tuples, consider it (index, filename)
            if (PyTuple_Check(input1))
            {
              // Get the index and filename from the Python tuple object
              index1 = PyInt_AsSsize_t(PyTuple_GetItem(input1, 0));
              fn1 = PyString_AsString(PyObject_Str(PyTuple_GetItem(input1, 1)));
            }
            else
              fn1 = PyString_AsString(PyObject_Str(input1));

            if (PyTuple_Check(input2))
            {
              // Get the index and filename from the Python tuple object
              index2 = PyInt_AsSsize_t(PyTuple_GetItem(input2, 0));
              fn2 = PyString_AsString(PyObject_Str(PyTuple_GetItem(input2, 1)));
            }
            else
              fn2 = PyString_AsString(PyObject_Str(input2));

            bool result = compareTwoImageTolerance(fn1, fn2, tolerance, index1, index2);

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
PyObject *
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
PyObject *
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

PyObject *
xmipp_compareTwoMetadataFiles(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn1, *pyStrFn2, *pyStrAux;

    if (PyArg_ParseTuple(args, "OO", &pyStrFn1, &pyStrFn2))
    {
        try
        {
            FileName fn1, fn2;

            pyStrAux = PyObject_Str(pyStrFn1);

            if (pyStrAux != NULL)
                fn1 = PyString_AsString(pyStrAux);
            else
                PyErr_SetString(PyExc_TypeError,
                                "Expected string or FileName as first argument");
            pyStrAux = PyObject_Str(pyStrFn2);
            if (pyStrAux != NULL)
                fn2 = PyString_AsString(pyStrAux);
            else
                PyErr_SetString(PyExc_TypeError,
                                "Expected string or FileName as first argument");

            if (compareTwoMetadataFiles(fn1, fn2))
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

/* calculate enhanced psd and return preview*/
PyObject *
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
                Image_Value(pyImage).setDatatype(DT_Double);
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
size_t w = dim, h = dim, &x = idim.xdim, &y = idim.ydim;\
if (x > y) h = y * (dim/x);\
else if (y > x)\
  w = x * (dim/y);\
selfScaleToSize(LINEAR, data, w, h);\
Image_Value(pyImage).setDatatype(DT_Double);\
data.resetOrigin();\
MULTIDIM_ARRAY_GENERIC(Image_Value(pyImage)).setImage(data);\
Py_RETURN_NONE;\
}} catch (XmippError &xe)\
{ PyErr_SetString(PyXmippError, xe.msg.c_str());}\


/* calculate enhanced psd and return preview
* used for protocol preprocess_particles*/
PyObject *
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
PyObject *
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
PyObject *
xmipp_realGaussianFilter(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyImage;
    double realSigma;
    int dim;
    FileName fn;

    if (PyArg_ParseTuple(args, "OOdi", &pyImage, &pyStrFn, &realSigma, &dim))
    {
        FILTER_TRY()
        realGaussianFilter(data, realSigma);
        FILTER_CATCH()
    }
    return NULL;
}

/* calculate enhanced psd and return preview
 * used for protocol preprocess_particles*/
PyObject *
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

/* dump metadatas to database*/
PyObject *
xmipp_dumpToFile(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyStrFn, *pyStrAux;
    FileName fn;

    if (PyArg_ParseTuple(args, "O", &pyStrFn))
    {
        pyStrAux = PyObject_Str(pyStrFn);
        if (pyStrAux != NULL)
        {
            fn = PyString_AsString(pyStrAux);
            MDSql::dumpToFile(fn);
            Py_RETURN_NONE;
        }
    }
    return NULL;
}
PyObject *
xmipp_Euler_angles2matrix(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    //PyObject *pyStrFn, *pyStrAux;
    //FileName fn;
    double rot, tilt, psi;
    if (PyArg_ParseTuple(args, "ddd", &rot,&tilt,&psi))
    {
        npy_intp dims[2];
        dims[0] = 3;
        dims[1] = 3;
        PyArrayObject * arr = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        void * data = PyArray_DATA(arr);
        Matrix2D<double> euler(3,3);
        Euler_angles2matrix(rot, tilt, psi,euler,false);
        memcpy(data, (euler.mdata), 9 * sizeof(double));
        return (PyObject*)arr;
    }
    return NULL;
}



PyObject *
xmipp_Euler_matrix2angles(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject * input;
    if (PyArg_ParseTuple(args, "O", &input))
    {
        PyArrayObject * arr = (PyArrayObject*) input;
        //this is 3*4 matrix so he need to delete last column
        //try first 3x3
        //IS DE DATA DOUBLE? CREATE NUMPY DOUBLE
        void * data = PyArray_DATA(arr);
        Matrix2D<double> euler(3,3);
        memcpy((euler.mdata),data, 9 * sizeof(double));
        double rot, tilt, psi;
        Euler_matrix2angles(euler,rot, tilt, psi);
        return Py_BuildValue("fff", rot, tilt, psi);//fff three real
    }
    return NULL;
}


PyObject *
xmipp_Euler_direction(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject * input;
    double rot, tilt, psi;
    if (PyArg_ParseTuple(args, "ddd", &rot,&tilt,&psi))
    {
        Matrix1D<double> direction(3);
        Euler_direction(rot, tilt, psi, direction);
        return Py_BuildValue("fff", VEC_ELEM(direction, 0), VEC_ELEM(direction, 1), VEC_ELEM(direction, 2));//fff three real
    }
    return NULL;
}
/* activateMathExtensions */
PyObject *
xmipp_activateMathExtensions(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        if (MDSql::activateMathExtensions())
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}

/* activateRegExtensions */
PyObject *
xmipp_activateRegExtensions(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        if (MDSql::activateRegExtensions())
            Py_RETURN_TRUE;
        else
            Py_RETURN_FALSE;
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return NULL;
}

PyObject *
xmipp_errorBetween2CTFs(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyMd1, *pyMd2;
    double minFreq=0.05,maxFreq=0.25;
    size_t dim=256;

    if (PyArg_ParseTuple(args, "OO|idd"
                         ,&pyMd1, &pyMd2
                         ,&dim,&minFreq,&maxFreq))
    {
        try
        {
            if (!MetaData_Check(pyMd1))
                PyErr_SetString(PyExc_TypeError,
                                "Expected MetaData as first argument");
            else if (!MetaData_Check(pyMd2))
                PyErr_SetString(PyExc_TypeError,
                                "Expected MetaData as second argument");
            else
            {
                double error = errorBetween2CTFs(MetaData_Value(pyMd1),
                                                 MetaData_Value(pyMd2),
                                                 dim,
                                                 minFreq,maxFreq);
                return Py_BuildValue("f", error);
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;

}

PyObject *
xmipp_errorMaxFreqCTFs(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyMd1;
    double phaseDiffRad;

    if (PyArg_ParseTuple(args, "Od", &pyMd1, &phaseDiffRad))
    {
        try
        {
            if (!MetaData_Check(pyMd1))
                PyErr_SetString(PyExc_TypeError,
                                "Expected MetaData as first argument");
            else
            {
                double resolutionA = errorMaxFreqCTFs(MetaData_Value(pyMd1),
                                                      phaseDiffRad
                                                     );
                return Py_BuildValue("f", resolutionA);
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;

}

PyObject *
xmipp_errorMaxFreqCTFs2D(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyMd1, *pyMd2;

    if (PyArg_ParseTuple(args, "OO", &pyMd1, &pyMd2))
    {
        try
        {
            if (!MetaData_Check(pyMd1))
                PyErr_SetString(PyExc_TypeError,
                                "Expected MetaData as first argument");
            else if (!MetaData_Check(pyMd2))
                PyErr_SetString(PyExc_TypeError,
                                "Expected MetaData as second argument");
            else
            {
                double resolutionA = errorMaxFreqCTFs2D(MetaData_Value(pyMd1),
                                                        MetaData_Value(pyMd2)
                                                       );
                return Py_BuildValue("f", resolutionA);
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
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
        { "labelHasTag", xmipp_labelHasTag, METH_VARARGS,
          "Return the if the label has a specific tag" },
        { "labelIsImage", xmipp_labelIsImage, METH_VARARGS,
          "Return if the label has the TAGLABEL_IMAGE tag" },
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
        { "addLabelAlias", (PyCFunction) xmipp_addLabelAlias,
          METH_VARARGS, "Add lable alias dinamically in run time. Use for reading non xmipp star files" },
        { "createEmptyFile", (PyCFunction) xmipp_createEmptyFile,
          METH_VARARGS, "create empty stack (speed up things)" },
        { "getImageSize", (PyCFunction) xmipp_getImageSize,
          METH_VARARGS, "Get image dimensions" },
        { "MetaDataInfo", (PyCFunction) xmipp_MetaDataInfo, METH_VARARGS,
          "Get image dimensions of first metadata entry and the number of entries" },
        { "existsBlockInMetaDataFile", (PyCFunction) xmipp_existsBlockInMetaDataFile, METH_VARARGS,
          "Does block exists in file" },
        { "ImgCompare", (PyCFunction) xmipp_ImgCompare,  METH_VARARGS,
          "return true if both files are identical" },
        { "checkImageFileSize", (PyCFunction) xmipp_CheckImageFileSize,  METH_VARARGS,
          "return true if the file has at least as many bytes as needed to read the image" },
        { "checkImageCorners", (PyCFunction) xmipp_CheckImageCorners,  METH_VARARGS,
          "return false if the image has repeated pixels at some corner" },
        { "compareTwoFiles", (PyCFunction) xmipp_compareTwoFiles, METH_VARARGS,
          "return true if both files are identical" },
        { "bsoftRemoveLoopBlock", (PyCFunction) xmipp_bsoftRemoveLoopBlock, METH_VARARGS,
          "convert bsoft star files to xmipp star files" },
        { "bsoftRestoreLoopBlock", (PyCFunction) xmipp_bsoftRestoreLoopBlock, METH_VARARGS,
          "convert xmipp star files to bsoft star files" },
        { "compareTwoImageTolerance", (PyCFunction) xmipp_compareTwoImageTolerance, METH_VARARGS,
          "return true if both images are very similar" },
        { "readMetaDataWithTwoPossibleImages", (PyCFunction) xmipp_readMetaDataWithTwoPossibleImages, METH_VARARGS,
          "Read a 1 or two column list of micrographs" },
        { "substituteOriginalImages", (PyCFunction) xmipp_substituteOriginalImages, METH_VARARGS,
          "Substitute the original images into a given column of a metadata" },
        { "fastEstimateEnhancedPSD", (PyCFunction) xmipp_fastEstimateEnhancedPSD, METH_VARARGS,
          "Utility function to calculate PSD preview" },
        { "compareTwoMetadataFiles", (PyCFunction) xmipp_compareTwoMetadataFiles, METH_VARARGS,
          "Compare two metadata files" },
        { "bandPassFilter", (PyCFunction) xmipp_bandPassFilter, METH_VARARGS,
          "Utility function to apply bandpass filter" },
        { "gaussianFilter", (PyCFunction) xmipp_gaussianFilter, METH_VARARGS,
          "Utility function to apply gaussian filter in Fourier space" },
        { "realGaussianFilter", (PyCFunction) xmipp_realGaussianFilter, METH_VARARGS,
          "Utility function to apply gaussian filter in Real space" },
        { "badPixelFilter", (PyCFunction) xmipp_badPixelFilter, METH_VARARGS,
          "Bad pixel filter" },
        { "dumpToFile", (PyCFunction) xmipp_dumpToFile, METH_VARARGS,
          "dump metadata to sqlite database" },
        { "Euler_angles2matrix", (PyCFunction) xmipp_Euler_angles2matrix, METH_VARARGS,
          "convert euler angles to transformation matrix" },
        { "Euler_matrix2angles", (PyCFunction) xmipp_Euler_matrix2angles, METH_VARARGS,
          "convert transformation matrix to euler angles" },
        { "Euler_direction", (PyCFunction) xmipp_Euler_direction, METH_VARARGS,
          "converts euler angles to direction" },
        { "activateMathExtensions", (PyCFunction) xmipp_activateMathExtensions,
          METH_VARARGS, "activate math function in metadatas" },
        { "activateRegExtensions", (PyCFunction) xmipp_activateRegExtensions,
          METH_VARARGS, "activate regular expressions in metadatas" },
        { "errorBetween2CTFs", (PyCFunction) xmipp_errorBetween2CTFs,
          METH_VARARGS, "difference between two metadatas" },
        { "errorMaxFreqCTFs", (PyCFunction) xmipp_errorMaxFreqCTFs,
          METH_VARARGS, "resolution at which CTFs phase differs more than 90º" },
        { "errorMaxFreqCTFs2D", (PyCFunction) xmipp_errorMaxFreqCTFs2D,
          METH_VARARGS, "resolution at which CTFs phase differs more than 90º, 2D case" },
        { NULL } /* Sentinel */
    };//xmipp_methods

#define INIT_TYPE(type) if (PyType_Ready(&type##Type) < 0) return; Py_INCREF(&type##Type);\
    PyModule_AddObject(module, #type, (PyObject *) &type##Type);

PyMODINIT_FUNC initxmipp(void)
{
    //Initialize module variable
    PyObject* module;
    module = Py_InitModule3("xmipp", xmipp_methods,
                            "Xmipp module as a Python extension.");
    import_array();

    //Check types and add to module
    INIT_TYPE(FileName);
    INIT_TYPE(Image);
    INIT_TYPE(MDQuery);
    INIT_TYPE(MetaData);
    INIT_TYPE(Program);
    INIT_TYPE(SymList);
    INIT_TYPE(FourierProjector);

    //Add PyXmippError
    char message[32]="xmipp.XmippError";
    PyXmippError = PyErr_NewException(message, NULL, NULL);
    Py_INCREF(PyXmippError);
    PyModule_AddObject(module, "XmippError", PyXmippError);

    //Add MDLabel constants
    PyObject * dict = PyModule_GetDict(module);
    addLabels(dict);
}
