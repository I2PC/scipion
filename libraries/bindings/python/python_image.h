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

#ifndef _PYTHON_IMAGE_H
#define _PYTHON_IMAGE_H

#include "Python.h"
#include "python_metadata.h"

/***************************************************************/
/*                            Image                         */
/**************************************************************/

#define Image_Check(v) (((v)->ob_type == &ImageType))
#define Image_Value(v) ((*((ImageObject*)(v))->image))

/*Image Object*/
typedef struct
{
    PyObject_HEAD
    ImageGeneric * image;
}
ImageObject;

#define ImageObject_New() (ImageObject*)malloc(sizeof(ImageObject))

/* Destructor */
void Image_dealloc(ImageObject* self);

/* Constructor */
PyObject *
Image_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

/* Image string representation */
PyObject *
Image_repr(PyObject * obj);

/* Image compare function */
int
Image_compare(PyObject * obj, PyObject * obj2);


/* Compare two images up to a precision */
PyObject *
Image_equal(PyObject *obj, PyObject *args, PyObject *kwargs);

/* write */
PyObject *
Image_write(PyObject *obj, PyObject *args, PyObject *kwargs);


/* read */
PyObject *
Image_read(PyObject *obj, PyObject *args, PyObject *kwargs);

/* read preview*/
PyObject *
Image_readPreview(PyObject *obj, PyObject *args, PyObject *kwargs);

/* convert to psd */
PyObject *
Image_convertPSD(PyObject *obj, PyObject *args, PyObject *kwargs);

/* readApplyGeo */
PyObject *
Image_readApplyGeo(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
Image_applyGeo(PyObject *obj, PyObject *args, PyObject *kwargs);

NPY_TYPES datatype2NpyType(DataType dt);

DataType npyType2Datatype(int npy);

/* getData */
PyObject *
Image_getData(PyObject *obj, PyObject *args, PyObject *kwargs);

/* projectVolumeDouble */
PyObject *
Image_projectVolumeDouble(PyObject *obj, PyObject *args, PyObject *kwargs);


/* setData */
PyObject *
Image_setData(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getPixel */
PyObject *
Image_getPixel(PyObject *obj, PyObject *args, PyObject *kwargs);

/* setPixel */
PyObject *
Image_setPixel(PyObject *obj, PyObject *args, PyObject *kwargs);

/* initConstant */
PyObject *
Image_initConstant(PyObject *obj, PyObject *args, PyObject *kwargs);

/* initRandom */
PyObject *
Image_initRandom(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Resize Image */
PyObject *
Image_resize(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Scale Image */
PyObject *
Image_scale(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Patch with other image */
PyObject *
Image_patch(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Set Data Type */
PyObject *
Image_setDataType(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Set Data Type */
PyObject *
Image_convert2DataType(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Return image dimensions as a tuple */
PyObject *
Image_getDimensions(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Return image dimensions as a tuple */
PyObject *
Image_getEulerAngles(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Return value from MainHeader*/
PyObject *
Image_getMainHeaderValue(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Set value to MainHeader*/
PyObject *
Image_setMainHeaderValue(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Return value from Header, now using only the first image*/
PyObject *
Image_getHeaderValue(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Set value to Header, now using only the first image*/
PyObject *
Image_setHeaderValue(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Return image dimensions as a tuple */
PyObject *
Image_computeStats(PyObject *obj, PyObject *args, PyObject *kwargs);

/* I1-adjusted(I2) */
PyObject *
Image_adjustAndSubtract(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
Image_add(PyObject *obj1, PyObject *obj2);
PyObject *
Image_subtract(PyObject *obj1, PyObject *obj2);

PyObject *
Image_iadd(PyObject *obj1, PyObject *obj2);
PyObject *
Image_isubtract(PyObject *obj1, PyObject *obj2);

extern PyNumberMethods Image_NumberMethods;
extern PyMethodDef Image_methods[];
extern PyTypeObject ImageType;


#endif
