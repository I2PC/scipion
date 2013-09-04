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

#ifndef _PYTHON_METADATA_H
#define _PYTHON_METADATA_H

#include "Python.h"

//#define MetaData_Check(v)  (((v)->ob_type == &MetaDataType))
#define MetaData_Check(v)  (((v)->ob_type == &MetaDataType))
#define MetaData_Value(v)  ((*((MetaDataObject*)(v))->metadata))

#define MDQuery_Check(v) (((v)->ob_type == &MDQueryType))
#define MDQuery_Value(v)  ((*((MDQueryObject*)(v))->query))

#define RETURN_MDOBJECT(value) return new MDObject((MDLabel)label, value)

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
void MDQuery_dealloc(MDQueryObject* self);

/* String representation */
PyObject *
MDQuery_repr(PyObject * obj);

/* Helper function to create relational queries */
 PyObject *
createMDValueRelational(PyObject *args, int op);

/* MDValue Relational */
 PyObject *
xmipp_MDValueRelational(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDValueEQ */
 PyObject *
xmipp_MDValueEQ(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDValueEQ */
 PyObject *
xmipp_MDValueNE(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDValueLT */
 PyObject *
xmipp_MDValueLT(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDValueLE */
 PyObject *
xmipp_MDValueLE(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDValueLT */
 PyObject *
xmipp_MDValueGT(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDValueLE */
 PyObject *
xmipp_MDValueGE(PyObject *obj, PyObject *args, PyObject *kwargs);

 /* MDValueRange */
  PyObject *
xmipp_MDValueRange(PyObject *obj, PyObject *args, PyObject *kwargs);

  /* addLabelAlias */
   PyObject *
xmipp_addLabelAlias(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MDQuery methods */
extern PyMethodDef MDQuery_methods[];

/*MDQuery Type */
extern PyTypeObject MDQueryType;

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
void MetaData_dealloc(MetaDataObject* self);
/* Constructor */
 PyObject *
MetaData_aggregate(PyObject *obj, PyObject *args, PyObject *kwargs);

 PyObject *
MetaData_aggregateMdGroupBy(PyObject *obj, PyObject *args, PyObject *kwargs);

 PyObject *
MetaData_aggregateSingle(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_aggregateSingleInt(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_addIndex(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_join(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_importObjects(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_intersection(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_merge(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
MetaData_new(PyTypeObject *type, PyObject *args, PyObject *kwargs);

 PyObject *
MetaData_operate(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_replace(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_readPlain(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_setComment(PyObject *obj, PyObject *args, PyObject *kwargs);
 PyObject *
MetaData_getComment(PyObject *obj, PyObject *args, PyObject *kwargs);
PyObject *
MetaData_unionAll(PyObject *obj, PyObject *args, PyObject *kwargs);
PyObject *
MetaData_randomize(PyObject *obj, PyObject *args, PyObject *kwargs);
PyObject *
MetaData_selectPart(PyObject *obj, PyObject *args, PyObject *kwargs);


int MetaData_print(PyObject *obj, FILE *fp, int flags);

/* String representation */
PyObject *
MetaData_repr(PyObject * obj);

/* MetaData compare function */
int
MetaData_compare(PyObject * obj, PyObject * obj2);

/* read */
PyObject *
MetaData_read(PyObject *obj, PyObject *args, PyObject *kwargs);

/* read */
PyObject *
MetaData_readPlain(PyObject *obj, PyObject *args, PyObject *kwargs);

/* read block */
PyObject *
MetaData_readBlock(PyObject *obj, PyObject *args, PyObject *kwargs);

/* write */
PyObject *
MetaData_write(PyObject *obj, PyObject *args, PyObject *kwargs);

/* append */
PyObject *
MetaData_append(PyObject *obj, PyObject *args, PyObject *kwargs);

/* addObject */
PyObject *
MetaData_addObject(PyObject *obj, PyObject *args, PyObject *kwargs);

/* firstObject */
PyObject *
MetaData_firstObject(PyObject *obj, PyObject *args, PyObject *kwargs);

/* lastObject */
PyObject *
MetaData_lastObject(PyObject *obj, PyObject *args, PyObject *kwargs);

/* size */
PyObject *
MetaData_size(PyObject *obj, PyObject *args, PyObject *kwargs);

/* size */
PyObject *
MetaData_getParsedLines(PyObject *obj, PyObject *args, PyObject *kwargs);

/* isEmpty */
PyObject *
MetaData_isEmpty(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getColumnFormat */
PyObject *
MetaData_getColumnFormat(PyObject *obj, PyObject *args, PyObject *kwargs);

/* setColumnFormat */
PyObject *
MetaData_setColumnFormat(PyObject *obj, PyObject *args, PyObject *kwargs);

/* setValue */
PyObject *
MetaData_setValue(PyObject *obj, PyObject *args, PyObject *kwargs);

/* setValueCol */
PyObject *
MetaData_setValueCol(PyObject *obj, PyObject *args, PyObject *kwargs);

/* drop column */
PyObject *
MetaData_removeLabel(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getValue */
PyObject *
MetaData_getValue(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getValue */
PyObject *
MetaData_getColumnValues(PyObject *obj, PyObject *args, PyObject *kwargs);

/* setValue */
PyObject *
MetaData_setColumnValues(PyObject *obj, PyObject *args, PyObject *kwargs);

/* containsLabel */
PyObject *
MetaData_getActiveLabels(PyObject *obj, PyObject *args, PyObject *kwargs);

/* containsLabel */
PyObject *
xmipp_getBlocksInMetaDataFile(PyObject *obj, PyObject *args);

/* containsLabel */
PyObject *
MetaData_getMaxStringLength(PyObject *obj, PyObject *args, PyObject *kwargs);

/* containsLabel */
PyObject *
MetaData_containsLabel(PyObject *obj, PyObject *args, PyObject *kwargs);

/* addLabel */
PyObject *
MetaData_addLabel(PyObject *obj, PyObject *args, PyObject *kwargs);


/* fillConstant */
PyObject *
MetaData_fillConstant(PyObject *obj, PyObject *args, PyObject *kwargs);

/* fillRandom */
PyObject *
MetaData_fillRandom(PyObject *obj, PyObject *args, PyObject *kwargs);

/* fillExpand */
PyObject *
MetaData_fillExpand(PyObject *obj, PyObject *args, PyObject *kwargs);

/* copyColumn */
PyObject *
MetaData_copyColumn(PyObject *obj, PyObject *args, PyObject *kwargs);

/* copyColumn */
PyObject *
MetaData_renameColumn(PyObject *obj, PyObject *args, PyObject *kwargs);

/* copyColumnTo */
PyObject *
MetaData_copyColumnTo(PyObject *obj, PyObject *args, PyObject *kwargs);

/* removeObjects */
PyObject *
MetaData_removeObjects(PyObject *obj, PyObject *args, PyObject *kwargs);

/* removeObjects */
PyObject *
MetaData_removeDisabled(PyObject *obj, PyObject *args, PyObject *kwargs);

/* Make absolute path */
PyObject *
MetaData_makeAbsPath(PyObject *obj, PyObject *args, PyObject *kwargs);

/* clear */
PyObject *
MetaData_clear(PyObject *obj, PyObject *args, PyObject *kwargs);


/** Iteration functions */
PyObject *
MetaData_iter(PyObject *obj);

PyObject *
MetaData_iternext(PyObject *obj);

/** Sort Metadata */
PyObject *
MetaData_sort(PyObject *obj, PyObject *args, PyObject *kwargs);

/** remove Duplicate rows Metadata */
PyObject *
MetaData_removeDuplicates(PyObject *obj, PyObject *args, PyObject *kwargs);

/* MetaData methods */
extern PyMethodDef MetaData_methods[];

/*MetaData Type */
extern PyTypeObject MetaDataType;

/*Helper function to create an MDObject from a PyObject */
MDObject * createMDObject(int label, PyObject *pyValue);

void setMDObjectValue(MDObject *obj, PyObject *pyValue);

PyObject * getMDObjectValue(MDObject * obj);

#endif
