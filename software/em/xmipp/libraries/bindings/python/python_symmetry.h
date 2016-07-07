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

#ifndef _PYTHON_SYMMETRY_H
#define _PYTHON_SYMMETRY_H

#include "Python.h"

/***************************************************************/
/*                            SymList                          */
/***************************************************************/

/*SymList Object*/
typedef struct
{
    PyObject_HEAD
    SymList * symlist;
}
SymListObject;

PyObject *
SymList_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

/* Destructor */
void SymList_dealloc(SymListObject* self);

/* readSymmetryFile */
PyObject *
SymList_readSymmetryFile(PyObject * obj, PyObject *args, PyObject *kwargs);

/* computeDistance */
PyObject *
SymList_computeDistance(PyObject * obj, PyObject *args, PyObject *kwargs);

/* getSymmetryMatrices */
PyObject *
SymList_getSymmetryMatrices(PyObject * obj, PyObject *args, PyObject *kwargs);

/* SymList methods */
extern PyMethodDef SymList_methods[];
/*SymList Type */
extern PyTypeObject SymListType;

#endif
