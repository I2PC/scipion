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

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from pyworkflow.mapper import Mapper, buildObject

class XmlMapper(Mapper):
    '''Mapper for XML'''
    def __init__(self, filename=None, rootName='ROOT', version=None, header=''):
        Mapper.__init__(self)
        self.filename = filename
        self.root = ET.Element(rootName)
        comment = ET.Comment(header)
        self.root.append(comment)
        self.root.set("version", str(version))
        
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
    
    def getNonNullObjectsToStoreXML(self):
        pass

            
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
    
    def select(self, **args):
        '''Select object meetings some criterias'''
        pass 
            
    
