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

class TableLayoutConfiguration():
    def __init__(self, ds, tableDataset, allowRender=True, defaultColumnsLayoutProperties=None):
        
        self.columnsLayout = OrderedDict() 
         
        for col in tableDataset.iterColumns():
#            self.columnsLayout[col.getName()]=ColumnLayoutConfiguration(col, ds.getTypeOfColumn(col.getName()), allowRender, defaultColumnsLayoutProperties[col.getName()] if defaultColumnsLayoutProperties != None else {})
            self.columnsLayout[col.getName()]=ColumnLayoutConfiguration(col, ds, allowRender, defaultColumnsLayoutProperties[col.getName()] if defaultColumnsLayoutProperties != None else {})
            
        self.colsOrder = defineColsLayout(self.columnsLayout.keys())
        
    def getLabelsToRenderComboBoxValues(self):
        labelsToRender = [columnLayout.label for columnLayout in self.columnsLayout.values() if (columnLayout.columnLayoutProperties.renderable or columnLayout.columnLayoutProperties.allowSetRenderable)]
        return tuple(zip(labelsToRender,labelsToRender))
    
    def hasEnableColumn(self):
        for columnLayout in self.columnsLayout.values():
            if "enable" == columnLayout.label:
                return True
        return False;
        
            
class ColumnLayoutConfiguration():
    def __init__(self, col, ds, allowRender, defaultColumnLayoutProperties):
        self.columns = col
        
        self.label = col.getName()
        self.typeOfColumn = ds.getTypeOfColumn(col.getName())
        
        self.columnLayoutProperties = ColumnLayoutProperties(self.typeOfColumn, allowRender)
        self.columnLayoutProperties.setValues(defaultColumnLayoutProperties)
        
        print "columnLayoutProperties",self.columnLayoutProperties 
        
        
#NAPA DE LUXE: Y si lo pasamos a una namedtupple
class ColumnLayoutProperties():
    def __init__(self, typeOfColumn, allowRender=True):
        #self.visible = defaultColumnLayoutProperties["visible"] if "visible" in defaultColumnLayoutProperties else not(typeOfColumn == 'id')
        self.visible = not(typeOfColumn == 'id')
        self.allowSetVisible = True 
        
        self.editable = (typeOfColumn == 'text')
        self.allowSetEditable = self.editable
        
        self.renderable = False
        self.allowSetRenderable = (typeOfColumn == 'image' and allowRender)

        self.renderFunc = "get_image"
        self.extraRenderFunc = ""
        
    def setValues(self, defaultColumnLayoutProperties):
        print "defaultColumnLayoutProperties",defaultColumnLayoutProperties
        for key in defaultColumnLayoutProperties:
            setattr(self, key, defaultColumnLayoutProperties[key])
        
def defineColsLayout(labels):
    colsOrder = range(len(labels))
    if 'enabled' in labels:
        colsOrder.insert(0, colsOrder.pop(labels.index('enabled')))
    return colsOrder


class ColumnLayoutConfigurationEncoder(json.JSONEncoder):
    def default(self, columnLayoutConfiguration):
        columnLayoutConfigurationCoded={}
        columnLayoutConfigurationCoded={"typeOfColumn":columnLayoutConfiguration.typeOfColumn,
                                        "columnLayoutProperties":{"visible":columnLayoutConfiguration.columnLayoutProperties.visible,
                                                                  "allowSetVisible":columnLayoutConfiguration.columnLayoutProperties.allowSetVisible,
                                                                  "editable":columnLayoutConfiguration.columnLayoutProperties.editable,
                                                                  "allowSetEditable":columnLayoutConfiguration.columnLayoutProperties.allowSetEditable,
                                                                  "renderable":columnLayoutConfiguration.columnLayoutProperties.renderable,
                                                                  "allowSetRenderable":columnLayoutConfiguration.columnLayoutProperties.allowSetRenderable,
                                                                  "renderFunc":columnLayoutConfiguration.columnLayoutProperties.renderFunc,
                                                                  "extraRenderFunc":columnLayoutConfiguration.columnLayoutProperties.extraRenderFunc
                                                                  }
                                        }
        return columnLayoutConfigurationCoded  
