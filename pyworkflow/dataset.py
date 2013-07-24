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
    
    def __init__(self):
        self._tables = []
    
    def listTables(self):
        """ List the actual tables in the DataSet. """
        return self._tables
    
    def getTable(self, tableName):
        if not tableName in self._tables:
            raise Exception("DataSet: table '%s' not found.\n   Current tables: %s" % (tableName, ))


class Table(object):
    """ Table to hold rows of data. 
    A table contains a list of columns. """
    def __init__(self, *columns):
        self._columns = []
        self._rowDict = OrderedDict()
        self._addColumn(Column('id', int))
        
        for col in columns:
            self._addColumn(col)
            
        # This imply that columns can only be added in __init__
        self.Row = namedtuple('Row', ['x', 'y'], verbose=True)
        
    def _addColumn(self, column):
        self._columns.append(column)
        
    def getColumns(self):
        """ Return all columns. """
        return self._columns
    
    def getRows(self):
        """ Return all rows. """
        return [row for row in self.iterRows()]
    
    def getRow(self, rowId):
        return self._rowDict[rowId]
    
    def addRow(self, *values):
        row = self.Row(values)
        self._rowDict[row.id] = row
    
    def iterRows(self):
        """ Iterate over the rows. """
        return self._rowDict.itervalues()
    
        
class Column(object):
    def __init__(self, colName, colType, default=None):
        self._name = colName
        self._type = colType
        self._default = default    
    
    
