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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules implements a DataSet, a container of several Tables.
This will serve as an abstraction layer from where the data will be taken.
"""
from collections import OrderedDict, namedtuple

class DataSet(object):
    """ Holds several Tables
    All tables should have an unique tableName. 
    """
    
    def __init__(self, tables, tableName=None, volumeName=None, numberSlices=0, labelToRender=None):
        self._tables = list(tables)
        self._tableName = tableName
        #NAPA de LUXE: Hay que ver si el volumen name se usa en algun lado        
        self._volumeName = volumeName 
        self._numberSlices = numberSlices 
        self._labelToRender = labelToRender
        
    def setTableName(self, tableName):    
        self._tableName = tableName
        
    def setVolumeName(self, volumeName):    
        self._volumeName = volumeName
        
    def getVolumeName(self):
        return self._volumeName
    
    def setLabelToRender(self, labelToRender):
        self._labelToRender = labelToRender
        
    def getDataToRender(self):
        return self.getTable().getColumnValues(self._labelToRender)
    
    def getDataToRenderAndExtra(self):
        return zip(self.getIdColumn(),
                   self.getTable().getColumnValues("enabled"),
                   self.getDataToRender(),
                   self.getTransformationMatrix())
    
    def getIdColumn(self):
        return self.getTable().getColumnValues("id")
    
    def setNumberSlices(self, numberSlices):
        self._numberSlices = numberSlices
        
    def getNumberSlices(self):
        return self._numberSlices    
    
    def getNumberSlicesForTemplate(self):
        return range(self._numberSlices)
        
    def getTransformationMatrix(self):    
        return self.getTable().getColumnValues(self._labelToRender+"_transformationMatrix")
        
    def listTables(self):
        """ List the actual table names on the DataSet. """
        return self._tables
    
    def getTable(self, tableName=None):
        if tableName == None:
            tableName = self._tableName
            
        if not tableName in self._tables:
            raise Exception("DataSet: table '%s' not found.\n   Current tables: %s" % 
                            (tableName, self._tables))
        table = self._loadTable(tableName)
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
    
    def getElementById(self, index, label):
        return self._rowDict.values()[index]._asdict()[label]
    
    def __str__(self):
        return '\n'.join([str(row) for row in self.iterRows()])
    
    
        
class Column(object):
    def __init__(self, colName, colType=None, default=None):
        self._name = colName
        self._type = colType
        self._default = default
        
    def getName(self):
        return self._name
    
    def getType(self):
        return self._type
    
    def convert(self, value):
        """ Try to convert the value to the column type. """
        return self._type(value)
    
    def hasDefault(self):
        return self._default is not None
    
    def getDefault(self):
        return self._default
    
    
