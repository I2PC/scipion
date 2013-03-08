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

/***************************************************************/
/*                            MDQuery                          */
/***************************************************************/

PyTypeObject MDQueryType =
    {
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
    }; //MDQueryType

PyMethodDef MDQuery_methods[] = { { NULL } /* Sentinel */
                                };

/* Destructor */
void MDQuery_dealloc(MDQueryObject* self)
{
    delete self->query;
    self->ob_type->tp_free((PyObject*) self);
}

/* String representation */
PyObject *
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

/* Methods for constructing concrete queries */

/* Helper function to create relational queries */
PyObject *
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
PyObject *
xmipp_MDValueRelational(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, -1);
}
/* MDValueEQ */
PyObject *
xmipp_MDValueEQ(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, EQ);
}
/* MDValueEQ */
PyObject *
xmipp_MDValueNE(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, NE);
}
/* MDValueLT */
PyObject *
xmipp_MDValueLT(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, LT);
}
/* MDValueLE */
PyObject *
xmipp_MDValueLE(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, LE);
}
/* MDValueLT */
PyObject *
xmipp_MDValueGT(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, GT);
}
/* MDValueLE */
PyObject *
xmipp_MDValueGE(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    return createMDValueRelational(args, GE);
}
/* MDValueRange */
PyObject *
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


/***************************************************************/
/*                            MetaData                         */
/**************************************************************/
PyMethodDef MetaData_methods[] =
    {
        { "read", (PyCFunction) MetaData_read, METH_VARARGS,
          "Read data from file" },
        { "write", (PyCFunction) MetaData_write, METH_VARARGS,
          "Write MetaData content to disk" },
        { "readBlock", (PyCFunction) MetaData_readBlock,
          METH_VARARGS, "Read block from a metadata file" },
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
        { "fillConstant", (PyCFunction) MetaData_fillConstant,
          METH_VARARGS, "Fill a column with constant value" },
        { "fillRandom", (PyCFunction) MetaData_fillRandom,
          METH_VARARGS, "Fill a column with random value" },
        { "copyColumn", (PyCFunction) MetaData_copyColumn,
          METH_VARARGS, "Copy the values of one column to another" },
        { "copyColumnTo", (PyCFunction) MetaData_copyColumnTo,
          METH_VARARGS, "Copy the values of one column to another in other md" },
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
          "Aggregate operation in metadata (double single value result)" },
        { "aggregateSingleInt",
          (PyCFunction) MetaData_aggregateSingleInt,
          METH_VARARGS,
          "Aggregate operation in metadata (int single value result)" },
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
        {
            "addIndex",
            (PyCFunction) MetaData_addIndex,
            METH_VARARGS,
            "Create index so search is faster." },
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
          METH_VARARGS, "Replace values in some column." },
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
        { "renameColumn", (PyCFunction) MetaData_renameColumn,
          METH_VARARGS, "Rename one column" },
        {
            "sort", (PyCFunction) MetaData_sort,
            METH_VARARGS,
            "Sort metadata according to a label" },
        { NULL } /* Sentinel */
    };//MetaData_methods


PyTypeObject MetaDataType =
    {
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
    };//MetaDataType

/* Destructor */
void MetaData_dealloc(MetaDataObject* self)
{
    delete self->metadata;
    delete self->iter;
    self->ob_type->tp_free((PyObject*) self);
}

int MetaData_print(PyObject *obj, FILE *fp, int flags)
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
PyObject *
MetaData_repr(PyObject * obj)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyString_FromString(
               (self->metadata->getFilename() + "(MetaData)").c_str());
}

/* MetaData compare function */
int
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
PyObject *
MetaData_read(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *list = NULL;
        PyObject *input = NULL, *pyStr = NULL;
        char *str = NULL;
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
PyObject *
MetaData_readPlain(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *input2 = NULL;
        PyObject *pyStr = NULL, *pyLabels = NULL, *pySep = NULL;
        char *str = NULL, *labels = NULL;

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
PyObject *
MetaData_readBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL, *blockName = NULL, *pyStr = NULL, *pyStrBlock =
                                                 NULL;
        char *str = NULL, *strBlock = NULL;
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
PyObject *
MetaData_write(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    WriteModeMetaData wmd;
    wmd = MD_OVERWRITE;
    if (self != NULL)
    {
        PyObject *input = NULL, *pyStr = NULL;
        char *str = NULL;
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

/* append */
PyObject *
MetaData_append(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;

    if (self != NULL)
    {
        PyObject *input = NULL;
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
PyObject *
MetaData_addObject(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyLong_FromUnsignedLong(self->metadata->addObject());
}
/* firstObject */
PyObject *
MetaData_firstObject(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyLong_FromUnsignedLong(self->metadata->firstObject());
}
/* lastObject */
PyObject *
MetaData_lastObject(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyLong_FromUnsignedLong(self->metadata->lastObject());
}
/* size */
PyObject *
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
    }
    return NULL;
}
/* isEmpty */
PyObject *
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
    }
    return NULL;
}
/* getColumnFormat */
PyObject *
MetaData_getColumnFormat(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    try
    {
        MetaDataObject *self = (MetaDataObject*) obj;
        if (self->metadata->isColumnFormat())
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
/* setColumnFormat */
PyObject *
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
PyObject *
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
PyObject *
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
PyObject *
MetaData_removeLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
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
PyObject *
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
PyObject *
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

            for (size_t i = 0; i < size; ++i)
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
PyObject *
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
PyObject *
MetaData_getActiveLabels(PyObject *obj, PyObject *args, PyObject *kwargs)
{
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
PyObject *
xmipp_getBlocksInMetaDataFile(PyObject *obj, PyObject *args)
{
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
PyObject *
MetaData_getMaxStringLength(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
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
PyObject *
MetaData_containsLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
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
PyObject *
MetaData_addLabel(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label, pos = -1;
    if (PyArg_ParseTuple(args, "i|i", &label, &pos))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->addLabel((MDLabel) label, pos);
            Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* fillConstant */
PyObject *
MetaData_fillConstant(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    PyObject *pyValue = NULL, *pyStr = NULL;
    if (PyArg_ParseTuple(args, "i|O", &label, &pyValue))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            if ((pyStr = PyObject_Str(pyValue)) != NULL)
            {
                char * str = PyString_AsString(pyStr);
                if (str != NULL)
                {
                    self->metadata->fillConstant((MDLabel) label, str);
                    Py_RETURN_TRUE;
                }
            }
            PyErr_SetString(PyXmippError, "MetaData.fillConstant: couldn't convert second argument to string");
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* fillConstant */
PyObject *
MetaData_fillRandom(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    double op1, op2, op3 = 0.;
    PyObject *pyValue = NULL, *pyStr = NULL;

    if (PyArg_ParseTuple(args, "iOdd|d", &label, &pyValue, &op1, &op2, &op3))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            if ((pyStr = PyObject_Str(pyValue)) != NULL)
            {
                char * str = PyString_AsString(pyStr);
                if (str != NULL)
                {
                    self->metadata->fillRandom((MDLabel) label, str, op1, op2, op3);
                    Py_RETURN_TRUE;
                }
            }
            PyErr_SetString(PyXmippError, "MetaData.fillRandom: couldn't convert second argument to string");
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* copyColumn */
PyObject *
MetaData_copyColumn(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int labelDst, labelSrc;
    if (PyArg_ParseTuple(args, "ii", &labelDst, &labelSrc))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->copyColumn((MDLabel)labelDst, (MDLabel)labelSrc);
            Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* copyColumn */
PyObject *
MetaData_renameColumn(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject * oldLabel = NULL;
    PyObject * newLabel = NULL;
    if (PyArg_ParseTuple(args, "OO", &oldLabel, &newLabel))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            if(PyInt_Check ( oldLabel ) && PyInt_Check ( newLabel ))
            {
                self->metadata->renameColumn((MDLabel) PyInt_AsLong (oldLabel),
                                             (MDLabel) PyInt_AsLong (newLabel));
            }
            else if (PyList_Check(oldLabel)&& PyList_Check ( newLabel ))
            {
                size_t size = PyList_Size(oldLabel);
                PyObject * itemOld = NULL;
                PyObject * itemNew = NULL;
                int iOldValue = 0;
                int iNewValue = 0;
                std::vector<MDLabel> vOldValue(size);
                std::vector<MDLabel> vNewValue(size);
                for (size_t i = 0; i < size; ++i)
                {
                    itemOld = PyList_GetItem(oldLabel, i);
                    itemNew = PyList_GetItem(newLabel, i);
                    if (!PyInt_Check(itemOld) || !PyInt_Check(itemNew))
                    {
                        PyErr_SetString(PyExc_TypeError,
                                        "MDL labels must be integers (MDLABEL)");
                        return NULL;
                    }
                    iOldValue = PyInt_AsLong(itemOld);
                    iNewValue = PyInt_AsLong(itemNew);
                    vOldValue[i] = (MDLabel)iOldValue;
                    vNewValue[i] = (MDLabel)iNewValue;
                }
                //self->metadata->read(str,&vValue);
                self->metadata->renameColumn(vOldValue,vNewValue);
            }

            Py_RETURN_TRUE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}

/* copyColumnTo */
PyObject *
MetaData_copyColumnTo(PyObject *obj, PyObject *args, PyObject *kwargs);

/* removeObjects */
PyObject *
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
PyObject *
MetaData_removeDisabled(PyObject *obj, PyObject *args, PyObject *kwargs)
{
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
PyObject *
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
PyObject *
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
PyObject *
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
PyObject *
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
PyObject *
MetaData_sort(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label   =  MDL_IMAGE;
    PyObject *ascPy = Py_True;
    bool asc=true;
    int limit   = -1;
    int offset  =  0;
    if (PyArg_ParseTuple(args, "|iOii", &label,&ascPy,&limit,&offset))
    {
        try
        {
            if (PyBool_Check(ascPy))
                asc = (ascPy == Py_True);
            MetaDataObject *self = (MetaDataObject*) obj;
            MetaData MDaux = *(self->metadata);
            self->metadata->clear();
            self->metadata->sort(MDaux, (MDLabel) label,asc,limit,offset);
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
PyObject *
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
PyObject *
MetaData_importObjects(PyObject *obj, PyObject *args, PyObject *kwargs)
{
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
PyObject *
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

/* aggregateSingleInt */
PyObject *
MetaData_aggregateSingleInt(PyObject *obj, PyObject *args, PyObject *kwargs)
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
            self->metadata->aggregateSingleInt(*object, (AggregateOperation) op,
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
PyObject *
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
PyObject *
MetaData_unionAll(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyMd = NULL;

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
PyObject *
MetaData_merge(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    PyObject *pyMd = NULL;

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
PyObject *
MetaData_setComment(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    char * str = NULL;
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
PyObject *
MetaData_getComment(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    MetaDataObject *self = (MetaDataObject*) obj;
    return PyString_FromString(self->metadata->getComment().c_str());
}

/* addIndex MetaData_join*/
PyObject *
MetaData_addIndex(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    if (PyArg_ParseTuple(args, "i", &label))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->addIndex((MDLabel)label);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    return NULL;
}
/* join */
PyObject *
MetaData_join(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int labelLeft;
    int labelRight;
    PyObject *pyMdLeft = NULL;
    PyObject *pyMdright = NULL;
    JoinType jt=LEFT;

    if (PyArg_ParseTuple(args, "OOii|i", &pyMdLeft, &pyMdright, &labelLeft,&labelRight, &jt))
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
PyObject *
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
PyObject *
MetaData_operate(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    char * str = NULL;

    if (PyArg_ParseTuple(args, "s", &str))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->operate(str);
            //free(str);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
   // free(str);
    return NULL;
}

/*Replace string values in some label. */
PyObject *
MetaData_replace(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int label;
    char * oldStr = NULL;
    char * newStr = NULL;
    if (PyArg_ParseTuple(args, "iss", &label, &oldStr, &newStr))
    {
        try
        {
            MetaDataObject *self = (MetaDataObject*) obj;
            self->metadata->replace((MDLabel)label, oldStr, newStr);
            //free(oldStr);
            //free(newStr);
            Py_RETURN_NONE;
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
    }
    //free(oldStr);
    //free(newStr);
    return NULL;
}

/* Randomize */
PyObject *
MetaData_randomize(PyObject *obj, PyObject *args, PyObject *kwargs)
{
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
PyObject *
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
        if (PyBool_Check(pyValue))
        {
            bool bValue = (pyValue == Py_True);
            RETURN_MDOBJECT(bValue);
        }
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

/* copyColumnTo */
PyObject *
MetaData_copyColumnTo(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    int labelDst, labelSrc;
    PyObject * pyMd;
    if (PyArg_ParseTuple(args, "Oii", &pyMd, &labelDst, &labelSrc))
    {
        try
        {
            if (MetaData_Check(pyMd))
            {
                MetaData_Value(obj).copyColumnTo(MetaData_Value(pyMd), (MDLabel)labelDst, (MDLabel)labelSrc);
                Py_RETURN_TRUE;
            }
        }
        catch (XmippError &xe)
        {
            PyErr_SetString(PyXmippError, xe.msg.c_str());
        }
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

PyObject *
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
    case LABEL_SIZET:
        return PyLong_FromLong(obj->data.longintValue);
    case LABEL_DOUBLE:
        return PyFloat_FromDouble(obj->data.doubleValue);
    case LABEL_STRING:
        return PyString_FromString(obj->data.stringValue->c_str());
    case LABEL_VECTOR_DOUBLE:
        {
        std::vector<double> & vector = *(obj->data.vectorValue);
        int size = vector.size();
        PyObject * list = PyList_New(size);
        for (int i = 0; i < size; ++i)
            PyList_SetItem(list, i, PyFloat_FromDouble(vector[i]));
        return list;
        }
    default:
    	return NULL;
    }//close switch
    return NULL;
}

