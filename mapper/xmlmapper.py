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

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from mapper import Mapper


class XmlMapper(Mapper):
    """Mapper for XML"""
    def __init__(self, dictClasses=None, rootName='ALL', **args):
        self.root = ET.Element(rootName)
        if 'header' in args:
            self.root.append(ET.Comment(args.get('header')))
        if 'version' in args:
            self.root.set("version", str(args.get('version')))
        Mapper.__init__(self, dictClasses)        
        self.pointerList=[]
        # This dictionary serve to define how to write classes
        # for example if the pair 'Integer': 'attribute' is present
        # all Integer will be store as attributes in the xml
        # possible values are: 
        #   'attribute', 'class_only', 'class_name', 'name_class', 'name_only'
        self.classTags = {} 

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

    def read(self, filename):
        self.filename = filename
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()

    def get(self, objId):
        for obj in self.objList:
            if obj.id == objId:
                return obj
        return None
        
    def getAll(self):
        """Select object from storage"""
        self.objList = []
        
        for child in self.root:
            obj = self.buildObject(child.tag, id=child.attrib)
            #objList.append(obj)
            self.objList.append(obj)
            self.fillObject(obj, child)
        #set properlly primary keys
        for obj in self.pointerList:
            for key, attr in obj.getAttributesToStore():
                if attr.isPointer():
                    attr.set(self.get(attr.id))
        return self.objList
                    
        
    def setChildObject(self, obj, childName, childClass=None):
        childObj = getattr(obj, childName)
        if childObj is None:
            if childClass is None:
                raise Exception('Attribute "%s" was not found in parent object and classname not stored' % childName)
            childObj = self.buildObject(childClass)
            setattr(obj, childName, childObj)
        return childObj
                  
    def fillObject(self, obj, objElem):
        # Set attributes first
        for k, v in objElem.attrib.iteritems():
            if k not in ['id', 'classname', 'name']:
                childObj = self.setChildObject(obj, k)
                childObj.set(v)
        # Set tree childs attributes
        for child in objElem:
            #print 'tag: ', child.tag, 'attrib: ', child.attrib
            if 'classname' in child.attrib:
                childName = child.tag
                childClass = child.attrib['classname']
            elif 'name' in child.attrib:
                childName = child.attrib['name']
                childClass = child.tag
            else:
                childName = child.tag
                childClass = None
                
            childObj = self.setChildObject(obj, childName, childClass)
            
            if child.text and len(child.text.strip()):
                childObj.set(child.text)
            else:
                #does have attributes?
                attributes = child.attrib
                if (attributes and childObj.isPointer()):
                    childObj.id = attributes
                    #save obj in auxiliary file
                    self.pointerList.append(obj)
                self.fillObject(childObj, child)
    
    def write(self, filename):
        self.filename = filename
        self.commit()
        
    def commit(self):
        self.indent(self.root)
        tree = ET.ElementTree(self.root)
        tree.write(self.filename,
                   xml_declaration=True,
                   encoding='utf-8')

    def setObjectId(self, objElem, obj):
        # Set attributes of this object element
        # The id is assumed to be a dictionary
        if type(obj.id) is dict:
            for k, v in obj.id.iteritems():
                objElem.set(k, str(v))
        else:
            objElem.set('id', str(obj.id))
        
    def insert(self, obj):
        """Insert a new object into the system"""
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
                attrClass = attr.getClassName()
                objClass = obj.getClassName()
                allKey = '%s.ALL' % objClass
                tag = self.classTags.get(allKey, '')
                classKey = '%s.%s' % (objClass, attrClass)
                tag = self.classTags.get(classKey, tag)
                if attr.tag == 'attribute' or tag == 'attribute':
                    parentElem.set(key, str(attr))
                else:
                    if attr.isPointer():
                        childElem = self.addSubElement(parentElem, key)
                        self.setObjectId(childElem, attr.get())
                    else:
                        # Select if we want to use the classname or elem name
                        # as element tagname
                        if tag.startswith('class'):
                            name = attrClass
                        else:
                            name = key
                        childElem = self.addSubElement(parentElem, name, attr.value)
                        if tag.endswith('class'):
                            childElem.set('classname', attrClass)
                        elif tag.endswith('name'):
                            childElem.set('name', key)
                        
                        
                        self.insertObjectWithChilds(attr, childElem)
            
    def updateFrom(self, obj):
        """Update object data with storage info"""
        pass
    
    def updateTo(self, obj):
        """Update storage with object info"""
        pass

            
    
