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




 /* FourierProjector methods */
 PyMethodDef FourierProjector_methods[] =
 {
    { "projectVolume", (PyCFunction) FourierProjector_projectVolume,
      METH_VARARGS, "projects Volume" },
    { NULL } /* Sentinel */
 };//FourierProjector_methods

 /*FourierProjector Type */
 PyTypeObject FourierProjectorType =
 {
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


 /* Constructor */
 PyObject *
 FourierProjector_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
 {

     FourierProjectorObject *self = (FourierProjectorObject*) type->tp_alloc(type, 0);
     XMIPP_TRY

     if (self != NULL)
     {
         PyObject *image = NULL;
         double padding_factor, max_freq, spline_degree;
         if (PyArg_ParseTuple(args, "O|ddd", &image, &padding_factor, &max_freq, &spline_degree))
         {
        	 ArrayDim dims;
        	 Image_Value(image).getDimensions(dims);
        	 self->dims = &dims;
        	 MultidimArray<double> *pdata;
        	 Image_Value(image).data->getMultidimArrayPointer(pdata);
        	 pdata->setXmippOrigin();
        	 self->fourier_projector = new FourierProjector(*(pdata), padding_factor, max_freq, spline_degree);

         }
     }
     XMIPP_CATCH

     return (PyObject *) self;
 }

 /* Destructor */
  void FourierProjector_dealloc(FourierProjectorObject* self)
 {
     delete self->fourier_projector;
     delete self->dims;
     self->ob_type->tp_free((PyObject*) self);
 }

/* projectVolume */

  PyObject * FourierProjector_projectVolume(PyObject * obj, PyObject *args, PyObject *kwargs)
{

      FourierProjectorObject *self = (FourierProjectorObject*) obj;
      double rot, tilt, psi;
      int i;
      PyObject *projection_image = NULL;
      if (self != NULL && PyArg_ParseTuple(args, "O|ddd", &projection_image, &rot, &tilt, &psi))
      {
          try
          {
        	  Projection P;
              projectVolume(FourierProjector_Value(self), P, self->dims->xdim, self->dims->ydim, rot, tilt, psi);
              Image_Value(projection_image).data->setImage(MULTIDIM_ARRAY(P));
          }
          catch (XmippError &xe)
          {
              PyErr_SetString(PyXmippError, xe.msg.c_str());
          }
          Py_RETURN_NONE;
      }


}






