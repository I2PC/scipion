# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from pyworkflow.em.convert import ImageHandler
from pyworkflow.mapper.sqlite import SqliteDb, SqliteFlatDb
from pyworkflow.mapper.sqlite_db import SqliteDb

"""
This modules implements a DataSet, a container of several Tables.
This will serve as an abstraction layer from where the data will be taken.
"""
from collections import OrderedDict, namedtuple

class DataSet(object):
    """ Holds several Tables
    All tables should have an unique tableName. 
    """
    
    def __init__(self, tables, tableName=None, volumeName=None, numberSlices=0):
        self._tables = list(tables)
        self._tableName = tableName
        #NAPA de LUXE: Hay que ver si el volumen name se usa en algun lado        
        self._volumeName = volumeName 
        self._numberSlices = numberSlices 
        self.projectPath = None
        
        
    def currentTable(self):
        """ Returns the name of the last selected table. """    
        return self._tableName
        
    def setVolumeName(self, volumeName):    
        self._volumeName = volumeName
        
    def getVolumeName(self):
        return self._volumeName
    
    def setNumberSlices(self, numberSlices):
        self._numberSlices = numberSlices
        
    def getNumberSlices(self):
        return self._numberSlices    
    
    def getNumberSlicesForTemplate(self):
        return range(self._numberSlices)
        
    def listTables(self):
        """ List the actual table names on the DataSet. """
        return self._tables
    
    def getTable(self, tableName=None):
        if tableName is None:
            tableName = self.listTables()[0]
        if not tableName in self._tables:
            raise Exception("DataSet: table '%s' not found.\n   Current tables: %s" % 
                            (tableName, self._tables))
            
        table = self._loadTable(tableName)
        self._tableName = tableName
        
        return table
    
    def getTypeOfColumn(self, label):
        """ this method should be implemented by subclasses. """
        pass
    
    def _loadTable(self, tableName):
        """ this method should be implemented by subclasses. """
        pass


class Table(object):
    """ Table to hold rows of data. 
    A table contains a list of columns. """
    def __init__(self, *columns):
        self._columns = OrderedDict()
        self._rowDict = OrderedDict()
        self._addColumn(Column('id', int))
        
        for col in columns:
            self._addColumn(col)
            
        colNames = [col.getName() for col in self.iterColumns()]
        # This imply that columns can only be added in __init__
        self.Row = namedtuple('Row', colNames, verbose=False)
    
    def setLabelToRender(self, labelToRender):
        self._labelToRender = labelToRender
        
    def _addColumn(self, col):
        self._columns[col.getName()] = col
        
    def getColumnValues(self, columnName):
        if (self.hasColumn(columnName)):
            return [getattr(row, columnName) for row in self.iterRows()] 
        else:
            return [None] * self.getSize()
        
    def iterColumns(self):
        return self._columns.itervalues()
    
    def hasColumn(self, columnName):
        """ Return true if column exists """
        return columnName in self._columns
    
    def getColumn(self, columnName):
        if not columnName in self._columns:
            raise Exception('Table: column "%s" not found.\Current columns: %s' % (columnName, '\n'.join(self._columns.keys())))
        return self._columns[columnName] 
    
    def hasEnabledColumn(self):
        """ Return true if enabled column exists """
        return self.hasColumn('enabled')
        
    def getColumns(self):
        """ Return all columns. """
        return self._columns.values()
    
    def getNumberOfColumns(self):
        return len(self._columns)
    
    def getSize(self):
        """ Return the number of rows. """
        return len(self._rowDict)
    
    def getRows(self):
        """ Return all rows. """
        return [row for row in self.iterRows()]
    
    def getRow(self, rowId):
        return self._rowDict[rowId]
    
    def _setRow(self, rowId, row):
        self._rowDict[rowId] = row
    
    def getDataToRenderAndExtra(self):
        return zip(self.getIdColumn(),
                   self.getColumnValues("enabled"),
                   self.getDataToRender(),
                   self.getTransformationMatrix())
    
    def getDataToRender(self):
        return self.getColumnValues(self._labelToRender)
    
    def getIdColumn(self):
        return self.getColumnValues("id")
        
    def getTransformationMatrix(self):    
        return self.getColumnValues(self._labelToRender+"_transformationMatrix")
    
    def _convertValues(self, values):
        """ Convert the input values to the actual
        expected type of each column. 
        """
        cValues = {}
        for k, v in values.iteritems():
            col = self.getColumn(k)
            cValues[k] = col.convert(v)
        
        return cValues
        
    def addRow(self, rowId, **values):
        """ With this implementation the rowId should be provided.
        We need to work around to also allow automatic generation of id's
        """
        values['id'] = rowId
        
        for col in self.iterColumns():
            if col.getName() not in values:
                if col.hasDefault():
                    values[col.getName()] = col.getDefault()
                else:
                    raise Exception('Table: value for column "%s" not provided.' % col.getName())
                
        row = self.Row(**self._convertValues(values))
        self._setRow(rowId, row)
        
    def updateRow(self, rowId, **values):
        """ Update a row given its rowId and some values to update. """
        row = self.getRow(rowId)
        self._setRow(rowId, row._replace(**self._convertValues(values)))
    
    def iterRows(self):
        """ Iterate over the rows. """
        return self._rowDict.itervalues()
    
    def getValueFromIndex(self, index, label):
        """ Return the value of the property 'label'
        in the element that has this 'index'.
        """
        value = self._rowDict.values()[index]._asdict()[label]
        return value
    
    def getIndexFromValue(self, value, label):
        """ Search the element that has property 'label'
        equals to value and returns its index.
        """
        for index, row in enumerate(self.iterRows()):
            if value == row._asdict()[label]:
                return index            
        return -1
    
    def __str__(self):
        return '\n'.join([str(row) for row in self.iterRows()])
    
        
COL_RENDER_NONE = 0
COL_RENDER_ID = 1
COL_RENDER_TEXT = 2
COL_RENDER_IMAGE = 3
COL_RENDER_CHECKBOX = 4
COL_RENDER_VOLUME = 5


class Column(object):
    def __init__(self, colName, colType=None, default=None, 
                 label=None, renderType=COL_RENDER_NONE):
        self._name = colName
        self._type = colType
        self._default = default
        self._label = label or colName
        self._renderType = renderType
        
    def getName(self):
        return self._name

    def getLabel(self):
        return self._label
        
    def getType(self):
        return self._type
    
    def convert(self, value):
        """ Try to convert the value to the column type. """
        return self._type(value)
    
    def hasDefault(self):
        return self._default is not None
    
    def getDefault(self):
        return self._default
    
    def getRenderType(self):
        return self._renderType
    
    def setRenderType(self, renderType):
        self._renderType = renderType
    

class SqliteDataSet(DataSet):
    """ Provide a DataSet implementation based on sqlite file.
    The tables of the dataset will be the object tables in database.
    Each block is a table on the dataset. 
    """
    
    def __init__(self, filename):
        self._dbName = filename
        db = SqliteDb()
        db._createConnection(filename, 1000)
        # Tables should be at pairs:
        # PREFIX_Classes
        # PREFIX_Objects  
        # where PREFIX can be empty
        self.tablePrefixes = OrderedDict()
        tables = db.getTables()
        for t in tables:
            if t.endswith('Classes'):
                prefix = t.replace('Classes', '')
                to = prefix + 'Objects'
                if to not in tables:
                    raise Exception('SqliteDataSet: table "%s" found, but not "%s"' % (t, to))
                flatDb = SqliteFlatDb(filename, tablePrefix=prefix)
                tableName = prefix + self._getPlural(flatDb.getSelfClassName())
                self.tablePrefixes[tableName] = prefix
                #tablePrefixes.append(prefix)
        DataSet.__init__(self, self.tablePrefixes.keys())
        db.close()
        
    def _getPlural(self, className):
        """ Get the plural of word for tables labels. """
        if className.startswith('Class'):
            return className.replace('Class', 'Classes')
        return className + 's'
        
    def _loadTable(self, tableName):
        """ Load information from tables PREFIX_Classes, PREFIX_Objects. """
        
        tableName = self.tablePrefixes[tableName]
        
        BASIC_COLUMNS = [Column('id', int, renderType=COL_RENDER_ID), 
                         Column('enabled', bool ,renderType=COL_RENDER_CHECKBOX),
                         Column('label', str), 
                         Column('comment', str),
                         Column('creation', str)]
        # Load columns from PREFIX_Classes table
        columns = list(BASIC_COLUMNS)
        db = SqliteDb()
        db._createConnection(self._dbName, 1000)
        db.executeCommand("SELECT * FROM %sClasses;" % tableName)
        # This will store the images columsn to join
        # the _index and the _filename
        imgCols = {}
        for row in db._iterResults():
            renderType = COL_RENDER_NONE
            colName = row['column_name']
            colLabel = row['label_property']
            
            if colLabel != 'self':
                # Keep track of _index and _filename pairs to mark as renderable images
                if colLabel.endswith('_index'):
                    imgCols[colLabel.replace('_index', '')] = colName
                
                elif colLabel.endswith('_filename'):
                    
                    # TODO: Maybe not all the labels endswith "_filename" 
                    # have to be rendered. 
                    # for example in the RotSpectra with '_representative._filename'

                    prefix = colLabel.replace('_filename', '')
                    if prefix in imgCols:
                        renderType = COL_RENDER_IMAGE
                        imgCols[colName] = imgCols[prefix]
                
                #CTF FIX
                elif (colLabel.endswith('_psdFile') or 
                      colLabel.endswith('_enhanced_psd') or 
                      colLabel.endswith('_ctfmodel_quadrant') or 
                      colLabel.endswith('_ctfmodel_halfplane')):
                    
                    renderType = COL_RENDER_IMAGE
                
                if row['class_name'] == 'Boolean':
                    renderType = COL_RENDER_CHECKBOX   
                columns.append(Column(colName, str, label=colLabel, renderType=renderType))
        table = Table(*columns)
        
        checkedImgCols = {} # Check if the image columns are volumes
        ih = ImageHandler() 
        
        # Populate the table in the DataSet
        db.executeCommand("SELECT * FROM %sObjects;" % tableName)
        for row in db._iterResults():
            rowDict = dict(row)
            for k, v in rowDict.iteritems():
                if v is None:
                    rowDict[k] = ''
                # Set the index@filename for images columns values
                if k in imgCols:
                    colName = imgCols[k]
                    index = rowDict[colName]
                    
                    filename = os.path.join(self.projectPath, rowDict[k])
                    filepath = filename.replace(":mrc", "")
                    if not checkedImgCols.get(colName, False):
                        if os.path.exists(filepath):
                            #print "Fn to get dims: %s@%s" % (index,filename)
                            x, y, z, n = ih.getDimensions((index, filename))
                            if z > 1:
                                table.getColumn(k).setRenderType(COL_RENDER_VOLUME)
                        checkedImgCols[colName] = True
                    if index:
                        rowDict[k] = '%06d@%s' % (index, filename)
            table.addRow(row['id'], **rowDict)
            
        return table
        
        
class SingleFileDataSet(DataSet):
    """ DataSet implementation for single files such as Images or Volumes. 
    """
    
    def __init__(self, filename):
        self._filename = filename
        self._tableName = ""
        DataSet.__init__(self, [self._tableName])
        self._table = self._createSingleTable()
        
    def _createSingleTable(self):
        table = Table(Column('filename', str, 
                            renderType=COL_RENDER_VOLUME)) #FIXME: for single images we need to read the dimensions
        table.addRow(1, filename=self._filename)
        
        return table
        
    def _loadTable(self, tableName):
        return self._table
