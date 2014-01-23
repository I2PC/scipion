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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from mapper import Mapper


class XmlMapper(Mapper):
    """Mapper for XML"""
    def __init__(self, filename, dictClasses=None, rootName='ALL', **args):
        self.filename = filename        
        if exists(filename):
            self._read()           
        else:
            self._create(rootName, **args)
        Mapper.__init__(self, dictClasses)
        # Objects map (id should be defined)  
        self.objDict = {}      
        # Store some objects during parsing
        # that have some Pointer attributes and need to be fixed later 
        self.pendingPtrDict = {}
        # This dictionary serve to define how to write classes
        # for example if the pair 'Integer': 'attribute' is present
        # all Integer will be store as attributes in the xml
        # possible values are: 
        #   'attribute', 'class_only', 'class_name', 'name_class', 'name_only'
        self.classTags = {} 
        # Counter to provide default objects id's
        self.objCount = 0

    def setClassTag(self, classname, tag):
        self.classTags[classname] = tag
        
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
            # TODO: generalize the matching of consecutive items
            if level and (not elem.tail or not elem.tail.strip()):
                for ii in ['t11','t12','t13','t21','t22','t23','t31','t32','t33']:
                    if elem.tag == ii:
                        elem.tail = " "
                        return
                elem.tail = i
    
    def _create(self, rootName, **args):
            self.root = ET.Element(rootName)
            if 'header' in args:
                self.root.append(ET.Comment(args.get('header')))
            if 'version' in args:
                self.root.set("version", str(args.get('version')))
    
    def _read(self):
        if self.filename is None:
            raise Exception("XmlMapper.read: filename is None.")
        if not exists(self.filename):
            raise Exception("XmlMapper.read: filename '%s' does not exists.")
        self.tree = ET.parse(self.filename)
        self.root = self.tree.getroot()

    def strId(self, objId):
        """Create an unique id from the objId,
        specially for dictionaries, where all keys
        and values will be converted to string"""
        if type(objId) is dict:
            ordKeys = objId.keys()
            ordKeys.sort()
            sid = ''
            for k in ordKeys:
                sid += '%s:%s ' % (k, objId[k])
            return sid
        return str(objId)
        
    def selectById(self, objId):
        """Retrieve an object give the id"""
        return self.objDict.get(self.strId(objId), None)
    
    def _addObjectToDict(self, obj, objDict):
        """Add object if have id"""
        if obj.hasObjId():
            key = self.strId(obj.getObjId())
        else:
            self.objCount += 1
            key = '%s_%06d' % (obj.getClassName(), self.objCount)
        objDict[key] = obj
        
    def selectAll(self):
        """Select object from storage"""
        for child in self.root:
            obj = self._buildObject(child.tag)
            self.fillObject(obj, child)
            self._addObjectToDict(obj, self.objDict)
        #set properlly primary keys
        for obj in self.pendingPtrDict.values():
            for key, attr in obj.getAttributesToStore():
                if attr.isPointer(): 
                    ptr = self.selectById(attr._objId)
                    attr._objId = None
                    attr.set(ptr)
        return self.objDict.values()
                   
        
    def setChildObject(self, obj, childName, childClass=None):
        childObj = getattr(obj, childName, None)
        if childObj is None:
            if childClass is None:
                raise Exception('Attribute "%s" was not found in parent object and classname not stored' % childName)
            childObj = self._buildObject(childClass)
            setattr(obj, childName, childObj)
        return childObj
                  
    def fillObject(self, obj, objElem):
        # Set attributes first
        if not obj.isPointer():
            for k, v in objElem.attrib.iteritems():
                if k not in ['id', 'classname', 'attrname']:
                    childObj = self.setChildObject(obj, k)
                    childObj.set(v)
        # Set tree childs attributes
        for child in objElem:
            if 'classname' in child.attrib:
                childName = child.tag
                childClass = child.attrib['classname']
            elif 'attrname' in child.attrib:
                childName = child.attrib['attrname']
                childClass = child.tag
            else:
                tagKey = '%s.%s' % (obj.getClassName(), child.tag)
                
                tag = self.classTags.get(tagKey, None)
                if tag == 'name_only' or not child.tag in self.dictClasses:
                    childName = child.tag
                    childClass = None
                else:
                    childName = ''
                    childClass = child.tag
                
            childObj = self.setChildObject(obj, childName, childClass)
            
            if child.text and len(child.text.strip()):
                childObj.set(child.text.strip())
            else:
                #does have attributes?
                attributes = child.attrib
                if (attributes and childObj.isPointer()):
                    childObj._objId = attributes
                    #save obj in auxiliary file
                    self._addObjectToDict(obj, self.pendingPtrDict)
                self.fillObject(childObj, child)
        
    def commit(self):
        if self.filename is None:
            raise Exception("XmlMapper.commit: filename is None")
        self.indent(self.root)
        tree = ET.ElementTree(self.root)
        tree.write(self.filename,
                   xml_declaration=True,
                   encoding='utf-8')

    def setObjectId(self, objElem, objId):
        # Set attributes of this object element
        if objId:
            if type(objId) is dict:
                for k, v in objId.iteritems():
                    if v:
                        objElem.set(k, str(v))
            else:        
                objElem.set('id', str(objId))
        
    def insert(self, obj):
        """Insert a new object into the system"""
        objElem = self.addSubElement(self.root, obj.getClassName(), obj._objValue) 
        self.setObjectId(objElem, obj.getObjId())
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
                attrClass = attr.getClassName()
                objClass = obj.getClassName()
                allKey = '%s.ALL' % objClass # First try with .ALL
                tag = self.classTags.get(allKey, '')
                allKey = 'ALL.%s' % attrClass 
                tag = self.classTags.get(allKey, tag)
                allKey = 'ALL.%s' % key # also with ALL.attribute
                tag = self.classTags.get(allKey, tag)
                classKey = '%s.%s' % (objClass, attrClass) # Second, with .childClass
                tag = self.classTags.get(classKey, tag)
                attrKey = '%s.%s' % (objClass, key) # Finally with .attribute
                tag = self.classTags.get(attrKey, tag)
                
                if tag == 'attribute':
                    parentElem.set(key, str(attr))
                else:
                    if attr.isPointer():
                        childElem = self.addSubElement(parentElem, key)
                        self.setObjectId(childElem, attr.get().getObjId())
                    else:
                        # Select if we want to use the classname or elem name
                        # as element tagname
                        if tag.startswith('class'):
                            name = attrClass
                        else:
                            name = key
                        childElem = self.addSubElement(parentElem, name, attr._objValue)
                        if tag.endswith('class'):
                            childElem.set('classname', attrClass)
                        elif tag.endswith('name'):
                            childElem.set('attrname', key)
                        
                        
                        self.insertObjectWithChilds(attr, childElem)
            
    def updateFrom(self, obj):
        """Update object data with storage info"""
        pass
    
    def updateTo(self, obj):
        """Update storage with object info"""
        pass

            
    
