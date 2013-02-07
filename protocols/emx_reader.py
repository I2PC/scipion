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

from emx_data_model import *

class EmxReader:
    '''Base class for reading EmxData from disk'''

    def __init__(self):
        self.version = 1.0

    def read(self, inFileName, emxData):
        self.inFileName = inFileName
        self.emxData = emxData
        self.readData()
    
    def readData(self):
        pass
    
    def readParticle(self, emxParticle):
        pass
    
    def readMicrograph(self, emxMicrograph):
        pass
    
    def readList(self, emxList, emxReadFunc):
        if not emxList is None and len(emxList) > 0:
            for emxObject in emxList:
                emxReadFunc(emxObject)
    
class EmxXmlReader (EmxReader):
    '''This reads xml-emx format
       Other class may be implemented for other formats
    '''
    def __init__(self,inFileName=None):
        EmxReader.__init__(self)

    def readData(self):
        '''read all objects in file
        '''
        self.tree = ET.parse(self.inFileName)
        self.root = self.tree.getroot()
        #micrograph must be read first
        self.readList(self.root.findall("micrograph"), self.readMicrograph)
        self.readList(self.root.findall("particle"),   self.readParticle)

    def readObjectPK(self, branch):
        ''' read primary key. So far all entities has the
        same PK. We may need to specialize or use dictPrimaryKeys
        in the future
        '''
        fileName = branch.get(FILENAME)
        if fileName is None:
            raise Exception("readObjectPK: No fileName" )
        index = branch.get(INDEX)
        if not index is None:
            index    = emxDataTypes[INDEX].type (index)
        return {FILENAME:fileName, INDEX:index}
        
    def setObjectFK(self, branch, object):
        ''' read foreign key. So far this is only relevant
        for particles. We may need to specialize or use dictPrimaryKeys
        in the future
        '''
        for key in object.dictForeignKeys.keys():
            pk = branch.find(key)
            if not pk is None:
                dictFK = self.readObjectPK(pk)
                for k, v in self.emxData.dictLists.iteritems():
                    if key.startswith(k):
                        # Find object in current data model matching found foreign keys dictionary
                        fk = self.emxData.findObject(v, **dictFK)
                        if fk is None: #At this point foreign object should be found
                            raise Exception("setObjectFK: Can not find foreign key with filename=%(fileName)s and index=%(index)d" % dictFK )
                        else: # Update foreign key
                            object.dictForeignKeys[key] = fk
                        return 
                #TODO: Print the dictFK instead of filename and index
                #so far there are no other foreign keys

    def readObject(self,branch,object,tag=''):
        ''' read a single object
        '''
        children = list(branch)
        for child in children:
            #if foreign key skip it
            if child.tag in object.dictForeignKeys.keys():
                continue
            #are there grandsons? then follow them
            if tag != '':
                _tag = tag + EMX_SEP + child.tag
            else:
                _tag = child.tag
            grandsons = list(child)
            nGranSons=len(grandsons)
            if(len(grandsons)>0):
                self.readObject(child,object,_tag)
            else:
                object.set(_tag,emxDataTypes[_tag].type(child.text))

    def readMicrograph(self,micrographBranch):
        #primary key
        dict = self.readObjectPK(micrographBranch)
        emxMicrograph = EmxMicrograph(dict[FILENAME], dict[INDEX])
        #other attributes
        self.readObject(micrographBranch,emxMicrograph)
        self.emxData.appendMicrograph(emxMicrograph)
            
    def readParticle(self,particleBranch):
        #primary key
        dict = self.readObjectPK(particleBranch)
        emxParticle = EmxParticle(dict[FILENAME], dict[INDEX])
        #other attributes
        self.setObjectFK(particleBranch, emxParticle)
        self.readObject(particleBranch,emxParticle)
        self.emxData.appendParticle(emxParticle)

