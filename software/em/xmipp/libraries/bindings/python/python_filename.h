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

#ifndef _PYTHON_FILENAME_H
#define _PYTHON_FILENAME_H

#include "Python.h"

/***************************************************************/
/*                            FileName                         */
/**************************************************************/

#define FileName_Check(v)  (((v)->ob_type == &FileNameType))
#define FileName_Value(v)  ((*((FileNameObject*)(v))->filename))

/*FileName Object*/
typedef struct
{
    PyObject_HEAD
    FileName * filename;
}
FileNameObject;

/* Destructor */
void FileName_dealloc(FileNameObject* self);

/* Constructor */
PyObject *
FileName_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

/* String representation */
PyObject *
FileName_repr(PyObject * obj);

/* compose */
PyObject *
FileName_compose(PyObject *obj, PyObject *args, PyObject *kwargs);

/* composeBlock */
PyObject *
FileName_composeBlock(PyObject *obj, PyObject *args, PyObject *kwargs);

/* exists */
PyObject *
FileName_exists(PyObject *obj, PyObject *args, PyObject *kwargs);

/* isInStack */
PyObject *
FileName_isInStack(PyObject *obj, PyObject *args, PyObject *kwargs);

/* isMetadata */
PyObject *
FileName_isMetaData(PyObject *obj, PyObject *args, PyObject *kwargs);

/* isImage */
PyObject *
FileName_isImage(PyObject *obj, PyObject *args, PyObject *kwargs);

/* isStar1 */
PyObject *
FileName_isStar1(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
FileName_getExtension(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
FileName_getNumber(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
FileName_getBaseName(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
FileName_decompose(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
FileName_withoutExtension(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
FileName_removeBlockName(PyObject *obj, PyObject *args, PyObject *kwargs);

extern PyMethodDef FileName_methods[];
extern PyTypeObject FileNameType;

#endif
