# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *             Adrian Quintana (aquintana@cnb.csic.es)   
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
from collections import OrderedDict
import json


class ColumnsConfig():
    """ Store the configuration of the columsn for a given table in a dataset.
    The order of the columns will be stored and configuration for each columns.
    For each column, we store properties:
    - visible
    - allowSetVisible
    - renderable
    - allowSetRenderable
    - editable
    - allowSetEditable
    - renderFunc
    - renderFuncExtra
    """
    def __init__(self, ds, table, allowRender=True, defaultColumnsLayoutProperties=None):
        
        self._columnsDict = OrderedDict() 
         
        for col in table.iterColumns():
            self._columnsDict[col.getName()] = ColumnProperties(col, ds, allowRender, defaultColumnsLayoutProperties[col.getName()] if defaultColumnsLayoutProperties != None else {})
        
    def getRenderableColumns(self):
        """ Return a list with the name of renderable columns. """
        columns = [col.getLabel() for col in self._columnsDict.values() if col.isRenderable()]
        return columns
    
    def hasEnableColumn(self):
        for columnLayout in self._columnsDict.values():
            if "enable" == columnLayout.label:
                return True
        return False;
    
    def getColumnProperty(self, colName, propName):
        """ Get some property value of a given column. """
        col = self._columnsDict[colName]
        return getattr(col, propName)
    
    def configColumn(self, colName, **kwargs):
        """ Configure properties of a given column. """
        col = self._columnsDict[colName]
        for k, v in kwargs.iteritems():
            setattr(col, k, v)
            
    def printColumns(self):
        for col in self._columnsDict.values():
            print "column: ", col.getLabel()
            print "  values: ", col.getValues()
        
            
class ColumnProperties():
    def __init__(self, col, ds, allowRender, defaultColumnLayoutProperties):
        self._column = col        
        self.columnType = ds.getTypeOfColumn(col.getName())
        
        self.visible = not(self.columnType == 'id')
        self.allowSetVisible = True 
        
        self.editable = (self.columnType == 'text')
        self.allowSetEditable = self.editable
        
        self.renderable = False
        self.allowSetRenderable = (self.columnType == 'image' and allowRender)

        self.renderFunc = "get_image"
        self.extraRenderFunc = ""
        
    def getLabel(self):
        return self._column.getName()
    
    def getColumnType(self):
        return self.columnType
    
    def isRenderable(self):
        return self.renderable or self.allowSetRenderable
        
    def setValues(self, defaultColumnLayoutProperties):
        for key in defaultColumnLayoutProperties:
            setattr(self, key, defaultColumnLayoutProperties[key])
    
    def getValues(self):
        return {"visible":self.visible,
                "allowSetVisible":self.allowSetVisible,
                "editable":self.editable,
                "allowSetEditable":self.allowSetEditable,
                "renderable":self.renderable,
                "allowSetRenderable":self.allowSetRenderable,
                "renderFunc":self.renderFunc,
                "extraRenderFunc":self.extraRenderFunc,
                'columnType': self.columnType
                }
        

class ColumnPropertiesEncoder(json.JSONEncoder):
    def default(self, columnProperties):
        
        return {"typeOfColumn":columnProperties.getColumnType(),
                "columnLayoutProperties":columnProperties.getValues()
                }
        
    
