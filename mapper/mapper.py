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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************


class Mapper():
    '''This class will serves as a Data Mapper pattern.
    It will store/retrieve objects from some storage enviroment.
    (like SQL, XML or others)'''
    def __init__(self):
        pass
    
    def commit(self):
        '''Commit changes made to the storage'''
        pass
    
    def insert(self, obj):
        '''Insert a new object into the system, the id will be set'''
        pass
    
    def update(self, obj, direction='to'):
        '''Update an existing object, the id should not be None
        direction can be "to" or "from" which indicates the 
        priority of the update.
        If "to" is used, object data is put in storage.
        If "from", object data is retrieved from storage'''
        if direction == 'to': # Update db with object changes
            self.updateTo(obj)
        elif direction == 'from': # Update object with db info
            self.updateFrom(obj)
        else:
            raise Exception('Invalid option %s for Sqlite.updateObject' % direction)
    
    def updateFrom(self, obj):
        '''Update object data with storage info'''
        pass
    
    def updateTo(self, obj):
        '''Update storage with object info'''
        pass
    
    def store(self, obj):
        '''Stores an object, it can be inserted or updated'''
        if obj.id is None:
            self.insert(obj)
        else:
            self.updateTo(obj)
            
    def get(self, objId):
        '''Return the object which id is objId'''
        pass
    
    def select(self, **args):
        '''Select object meetings some criterias'''
        pass 
            
    