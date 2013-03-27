/***************************************************************************
 *
 * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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
/*                            FourierProjector                         */
/**************************************************************/



/* Constructor */
PyObject *
FourierProjector_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{

    FourierProjectorObject *self = (FourierProjectorObject*) type->tp_alloc(type, 0);
    XMIPP_TRY

    if (self != NULL)
    {
        PyObject *image = NULL;
        if (PyArg_ParseTuple(args, "O", &image))
        {

          MultidimArray<double> * pVolume;
          Image_Value(image).data->getMultidimArrayPointer(pVolume);
          self->fourier_projector = new FourierProjector(*pVolume, 1, 0.5, NEAREST);

        }
    }
    XMIPP_CATCH

    return (PyObject *) self;
}

/* Destructor */
 void FourierProjector_dealloc(FourierProjectorObject* self)
{
    delete self->fourier_projector;
    self->ob_type->tp_free((PyObject*) self);
}

/* readSymmetryFile */
PyObject *
FourierProjector_projectVolume(PyObject * obj, PyObject *args, PyObject *kwargs)
{
      FourierProjectorObject *self = (FourierProjectorObject*) obj;
      double rot, tilt, psi;
      if (PyArg_ParseTuple(args, "O|ddd", &rot,&tilt,&psi))
      {
          try
          {
              Projection P;
              MultidimArray<double> * pVolume;
              //self->image->data->getMultidimArrayPointer(pVolume);
              ArrayDim aDim;
              pVolume->getDimensions(aDim);
              pVolume->setXmippOrigin();
              projectVolume(FourierProjector_Value(self), P, aDim.xdim, aDim.ydim,rot, tilt, psi);
              ImageObject * result = PyObject_New(ImageObject, &ImageType);
              Image <double> I;

              result->image = new ImageGeneric();
              result->image->setDatatype(DT_Double);
              result->image->data->setImage(MULTIDIM_ARRAY(P));
              return (PyObject *)result;
          }
          catch (XmippError &xe)
          {
              PyErr_SetString(PyXmippError, xe.msg.c_str());
          }
      }
      return NULL;
}


/* FourierProjector methods */
PyMethodDef FourierProjector_methods[] =
{
   { "projectVolume", (PyCFunction) FourierProjector_projectVolume,
     METH_VARARGS, "read symmetry file" },
   { NULL } /* Sentinel */
};//FourierProjector_methods

/*FourierProjector Type */
PyTypeObject FourierProjectorType =
{
>>>>>>> 5a12ca07fc5c116e804b6c117b15fe78c18ccfef
    PyObject_HEAD_INIT(NULL)
    0, /*ob_size*/
    "xmipp.FourierProjector", /*tp_name*/
    sizeof(FourierProjectorObject), /*tp_basicsize*/
    0, /*tp_itemsize*/
    (destructor)FourierProjector_dealloc, /*tp_dealloc*/
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
    "Python wrapper to Xmipp FourierProjector class",/* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    FourierProjector_methods, /* tp_methods */
    0, /* tp_members */
    0, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    0, /* tp_init */
    0, /* tp_alloc */
    FourierProjector_new, /* tp_new */
};//FourierProjectorType


