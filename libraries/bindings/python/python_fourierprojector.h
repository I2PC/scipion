/***************************************************************************
 *
 * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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

#ifndef _PYTHON_FOURIERPROJECTOR_H
#define _PYTHON_FOURIERPROJECTOR_H

#include "Python.h"


/***************************************************************/
/*                            FourierProjector                         */
/**************************************************************/

#define FourierProjector_Check(v) (((v)->ob_type == &FourierProjectorType))
#define FourierProjector_Value(v) ((*((FourierProjectorObject*)(v))->fourierprojector))

/*FourierProjector Object*/
typedef struct
{
    PyObject_HEAD
    FourierProjector* fourierprojector;
}
FourierProjectorObject;

#define FourierProjectorObject_New() (FourierProjectorObject*)malloc(sizeof(FourierProjectorObject))

/* Destructor */
void FourierProjector_dealloc(FourierProjectorObject* self);

/* Constructor */
PyObject *
FourierProjector_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

/* FourierProjector string representation */
PyObject *
FourierProjector_repr(PyObject * obj);

/* projectVolumeDouble */
PyObject *
FourierProjector_projectVolumeDouble(PyObject *obj, PyObject *args, PyObject *kwargs);

extern PyMethodDef FourierProjector_methods[];
extern PyTypeObject FourierProjectorType;


#endif
