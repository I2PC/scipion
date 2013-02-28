# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini       (roberto@cnb.csic.es)
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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

VERSION=1.0
ROOTNAME = 'EMX'
HEADER = '''
  ##########################################################################
  #               EMX Exchange file 
  #               Produced by the prolib_emx module
  # 
  #  This is a EMX file.
  #
  #  Information on this file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  #  One of the best ways you can help us to improve this software
  #  is to let us know about any problems you find with it.
  #  Please report bugs to: emx@cnb.csic.es
  ##########################################################################
  '''

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from pyworkflow.mapper import Mapper
from pyworkflow.emx import EmxData

from os.path import exists

class XmlMapper(Mapper):
    '''Mapper for XML'''
    def __init__(self,emxData,rootName=ROOTNAME):
        Mapper.__init__(self)
        self.emxData = emxData
        self.foreignKeyAuxList=[]

    def emxDataToXML(self, header=HEADER, 
                           rootName=ROOTNAME, 
                           version = VERSION):
        self.root = ET.Element(rootName)
        comment = ET.Comment(header)
        self.root.append(comment)
        #Create empty tree
        self.root.set("version", str(version))
        for k,v in self.emxData.objLists.iteritems():
            for obj in v:
                self.insert(obj)
            
    def indent(self, elem, level=0):
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for _elem in elem:
                self.indent(_elem, level+1)
            if not _elem.tail or not _elem.tail.strip():
                _elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                if elem.tag == 'm11':
                        elem.tail = " "
                        return
                elem.tail = i
#        if elem.tag=='EMX':
#            elem.text='\n'
    def read(self, filename):
        self.filename = filename
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
    
    def convertToEmxData(self, emxData):
        '''Select object meetings some criterias'''
        for child in self.root:
            obj = self.buildObject(child.tag, id=child.attrib)
            #objList.append(obj)
            emxData.addObject(obj)
            self.fillObject(obj, child)
        #set properlly primary keys
        for obj in self.foreignKeyAuxList:
            for key, attr in obj.getAttributesToStore():
                if attr.isPointer():
                    attr.set(self.emxData.getObjectwithID(key, attr.id))

    def fillObject(self, obj, objElem):
        for child in objElem:
            childObj = getattr(obj, child.tag)
            if child.text and len(child.text.strip()):
                childObj.set(child.text)
            else:
                #does have attributes?
                attributes = child.attrib
                if (attributes and childObj.isPointer()):
                    childObj.id = attributes
                    #save obj in auxiliary file
                    self.foreignKeyAuxList.append(obj)
                self.fillObject(childObj, child)
    def write(self, filename):
        self.filename = filename
        self.commit()
        
    def commit(self):
        self.indent(self.root)
        tree = ET.ElementTree(self.root)
        tree.write(self.filename)

    def setObjectId(self, objElem, obj):
        # Set attributes of this object element
        # The id is assumed to be a dictionary
        if obj.id:
            for k, v in obj.id.iteritems():
                objElem.set(k, str(v))
        
    def insert(self, obj):
        '''Insert a new object into the system'''
        objElem = self.addSubElement(self.root, obj.getClassName(), obj.value) 
        self.setObjectId(objElem, obj)
        # Insert object childs
        self.insertObjectWithChilds(obj, objElem )
        
    def addSubElement(self, parentElem, name, value=None):
        childElem = ET.SubElement(parentElem, name)
        if value != None:
            childElem.text = str(value)
        return childElem
            
    def insertObjectWithChilds(self, obj, parentElem):
        for key, attr in obj.getAttributesToStore():
            if attr.hasValue():
                if attr.tag == 'attribute':
                    parentElem.set(key, str(attr))
                else:
                    if attr.isPointer():
                        childElem = self.addSubElement(parentElem, key)
                        self.setObjectId(childElem, attr.get())
                    else:
                        childElem = self.addSubElement(parentElem, key, attr.value)
                        self.insertObjectWithChilds(attr, childElem)
            
    def updateFrom(self, obj):
        '''Update object data with storage info'''
        pass
    
    def updateTo(self, obj):
        '''Update storage with object info'''
        pass
            
    def get(self, objId):
        '''Return the object which id is objId'''
        pass

            
    
