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
# *  e-mail address 'xmipp@cnb.csic.es' Eso eso ...las quejas a XMIPP
# *
# **************************************************************************

from pyworkflow.object import *

class Mapper():
    """This class will serves as a Data Mapper pattern.
    It will store/retrieve objects from some storage enviroment.
    (like SQL, XML or others)
    The mapper should have access to class dictionary
    in order to build any give class by name"""
    def __init__(self, dictClasses=None):
        #print '=========> Mapper:', dictClasses
        if dictClasses:
            self.dictClasses = dictClasses 
        else:
            self.dictClasses = globals()
    
    def _buildObject(self, className, **args):
        """Build an instance of an object
        given the class name, it should be in 
        the classes dictionary"""
        if className in self.dictClasses:
            return self.dictClasses[className](**args)
        raise Exception('Mapper._buildObject: Unknown class: %s' % className)
    
    def commit(self):
        """Commit changes made to the storage"""
        pass
    
    def insert(self, obj):
        """Insert a new object into the system, the id will be set"""
        pass
    
    def update(self, obj, direction='to'):
        """Update an existing object, the id should not be None
        direction can be "to" or "from" which indicates the 
        priority of the update.
        If "to" is used, object data is put in storage.
        If "from", object data is retrieved from storage"""
        if direction == 'to': # Update db with object changes
            self.updateTo(obj)
        elif direction == 'from': # Update object with db info
            self.updateFrom(obj)
        else:
            raise Exception('Invalid option %s for Sqlite.updateObject' % direction)
    
    def updateFrom(self, obj):
        """Update object data with storage info"""
        pass
    
    def updateTo(self, obj):
        """Update storage with object info"""
        pass
    
    def store(self, obj):
        """Stores an object, it can be inserted or updated"""
        if obj._objId is None:
            self.insert(obj)
        else:
            self.updateDict = {}
            self.updateTo(obj)
            
    def selectById(self, objId):
        """Return the object which id is objId"""
        pass
    
    def selectAll(self):
        """Return all object from storage"""
        pass
    
    def selectFirst(self):
        """Return only the first element"""
        allObj = self.selectAll()
        if len(allObj):
            return allObj[0]
        return None
    
    def selectBy(self, **args):
        """Select object meetings some criterias"""
        pass 
    
    def selectByClass(self, className, includeSubclasses=True):
        """Select all object of a give class.
        By default object of subclasses will be retrieved also.
        This behaviour can be changed by passing includeSubclass=False"""
        return self.selectBy(classname=className)
    
    def addRelation(self, relName, creatorObj, parentObj, childObj):
        """ This function will add a new relation between two objects.
        Params:
            relName: the name of the relation to be added.
            creatorObj: this object will be the one who register the relation.
            parentObj: this is "parent" in the relation
            childObj: this is "child" in the relation
        """
        pass
    
    def getRelationChilds(self, relName, parentObj):
        """ Return all "child" objects for a given relation.
        Params:
            relName: the name of the relation.
            parentObj: this is "parent" in the relation
        Returns: 
            a list of "child" objects.
        """
        pass  
            
    def getRelationParents(self, relName, childObj):
        """ Return all "parent" objects for a given relation.
        Params:
            relName: the name of the relation.
            childObj: this is "child" in the relation
        Returns: 
            a list of "parent" objects.
        """
        pass   