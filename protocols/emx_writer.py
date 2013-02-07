'''
/***************************************************************************
 * Authors:     RobertoMarabini (roberto@cnb.csic.es)
 *
 *
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
MODIFICATION ADVICE:

Please,  do not  generate or  distribute 
a modified version of this file under its original name. 
 '''
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from emx_data_model import emxDataTypes , EMX_SEP, FILENAME, INDEX
class EmxWriter:
    '''Base class for writting EmxData to disk'''


    def __init__(self):
        self.version = 1.0
        self.inFileName=None 

    def write(self, outFileName, emxData):
        self.outFileName = outFileName
        self.emxData = emxData
        self.writeData()
    
    def writeData(self):
        pass    

    def writeParticle(self, emxParticle):
        pass
    
    def writeMicrograph(self, emxMicrograph):
        pass
    
    def writeList(self, emxList, emxWriteFunc):
        if not emxList is None and len(emxList) > 0:
            for emxObject in emxList:
                emxWriteFunc(emxObject) 

class EmxXmlWriter (EmxWriter):
    '''This class writes emx-xml format
    '''
    def __init__(self,inFileName=None):
        EmxWriter.__init__(self)

        
    def writeData(self):
        self.writeHeader()
        if self.emxData.listMicrographs:
            self.writeList(self.emxData.listMicrographs, self.writeMicrograph)
        if self.emxData.listParticles:
            self.writeList(self.emxData.listParticles,   self.writeParticle)
        self.writeFooter()    

    def writeHeader(self):
        self.root = ET.Element('EMX')
        header=\
'''
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
        if self.inFileName != None:
            header +=\
'''#  Inputfile: ''' + self.inFileName +'''
##########################################################################
'''
        self.root.set("version", str(self.version))
        self.root.append(ET.Comment(header))

    def writeFooter(self):
        self.indent(self.root)
        tree = ET.ElementTree(self.root)
        tree.write(self.outFileName)
        

    def writeObjectPK(self, object, entity):
        ''' Write primary key. So far all entities has the
        same PK. We may need to specialize or use dictPrimaryKeys
        in the future
        '''
        entity.set(FILENAME, object.get(FILENAME))
        index = object.get(INDEX)
        if not index is None:
            entity.set(INDEX, str(index))
        
    def writeObject(self, object):
        #primary key
        entity = ET.SubElement(self.root, object.name)
        self.writeObjectPK(object, entity)
        #foreign key
        for key, value in object.dictForeignKeys.iteritems():
            if value is None:
                continue
            _entity = ET.SubElement(entity, key)
            self.writeObjectPK(value, _entity)
        #other attributes
        for key, value in object.iterAttributes():
            if value is None:
                continue
            #is this an special case, that is, does the label contain '_'?
            if EMX_SEP in key:
                keyParts = key.split(EMX_SEP)
                _entity = entity
                for keyPart in keyParts[:-1]:
                    parent = entity.find(keyPart)
                    if parent is None:
                        parent = ET.SubElement(_entity, keyPart)
                    _entity = parent
                _key = keyParts[-1] # Take last key                    
            else:
                parent = entity
                _key = key
            #create branch
            child = ET.SubElement(parent, _key)
            label = emxDataTypes[key]            
            if label.hasUnit():
                child.set('unit', label.getUnit())
            child.text = str(object.get(key))


    def writeMicrograph(self,micrograph):
        self.writeObject(micrograph)
            
    def writeParticle(self,particle):
        self.writeObject(particle)
    
    def indent(self,elem, level=0):
        '''Prints a tree with each node indented according to its depth. 
           This custom made function will help to pretty print the xml output:
           in particular we can format matrices here.
           INPUT: XML tree
        '''
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self.indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                for ii in ['t11','t12','t13','t21','t22','t23','t31','t32','t33']:
                    if elem.tag == ii:
                        elem.tail = " "
                        return
                elem.tail = i
