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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.utils.path import replaceExt, joinExt
from mapper import Mapper

class SqliteMapper(Mapper):
    """Specific Mapper implementation using Sqlite database"""
    def __init__(self, dbName, dictClasses=None):
        Mapper.__init__(self, dictClasses)
        self.db = SqliteDb(dbName)
    
    def commit(self):
        self.db.commit()
        
    def __getObjectValue(self, obj):
        if obj.isPointer() and obj.hasValue():
            return obj.get().strId() # For pointers store the id of referenced object
        return obj.getInternalValue()
        
    def __insert(self, obj, namePrefix=None):
        obj._objId = self.db.insertObject(obj._objName, obj.getClassName(), 
                                          self.__getObjectValue(obj), obj._objParentId)
        sid = obj.strId()
        if namePrefix is None:
            namePrefix = sid
        else:
            namePrefix = joinExt(namePrefix, sid)
        self.insertChilds(obj, namePrefix)
        
    def insert(self, obj):
        """Insert a new object into the system, the id will be set"""
        self.__insert(obj)
        
    def insertChild(self, obj, key, attr, namePrefix=None):
        if namePrefix is None:
            namePrefix = self.__getNamePrefix(obj)
        attr._objName = joinExt(namePrefix, key)
        attr._objParentId = obj._objId
        self.__insert(attr, namePrefix)
        
    def insertChilds(self, obj, namePrefix=None):
        """ Insert childs of an object, if namePrefix is None,
        the it will be deduced from obj. """
        # This is also done in insertChild, but avoid 
        # doing the same for every child element
        if namePrefix is None:
            namePrefix = self.__getNamePrefix(obj)
        for key, attr in obj.getAttributesToStore():
            self.insertChild(obj, key, attr, namePrefix)
        
    def deleteChilds(self, obj):
        namePrefix = self.__getNamePrefix(obj)
        self.db.deleteChildObjects(namePrefix)
        
    def delete(self, obj):
        """Delete an object and all its childs"""
        self.deleteChilds(obj)
        self.db.deleteObject(obj.getId())
    
    def __getNamePrefix(self, obj):
        if len(obj._objName) > 0 and '.' in obj._objName:
            return replaceExt(obj._objName, obj.strId())
        return obj.strId()
    
    def updateTo(self, obj, level=1):
        self.db.updateObject(obj._objId, obj._objName, obj.getClassName(), 
                             self.__getObjectValue(obj), obj._objParentId)
        for key, attr in obj.getAttributesToStore():
            if attr._objId is None: # Insert new items from the previous state
                attr._objParentId = obj._objId
                #path = obj._objName[:obj._objName.rfind('.')] # remove from last .
                namePrefix = self.__getNamePrefix(obj)
                attr._objName = joinExt(namePrefix, key)
                self.__insert(attr, namePrefix)
            else:                
                self.updateTo(attr, level+2)
        
    def updateFrom(self, obj):
        objRow = self.db.selectObjectById(obj._objId)
        self.fillObject(obj, objRow)
            
    def selectById(self, objId):
        """Build the object which id is objId"""
        objRow = self.db.selectObjectById(objId)
        if objRow is None:
            obj = None
        else:
            obj = self._buildObject(objRow['classname'])
            self.fillObject(obj, objRow)
        return obj
    
    def fillObjectWithRow(self, obj, objRow):
        """Fill the object with row data"""
        obj._objId = objRow['id']
        name = objRow['name']
        if name is None or len(name) <= 0:
            obj._objName = ''#obj.strId()
        else:
            obj._objName = name
        objValue = objRow['value']
        obj._objParentId = objRow['parent_id']
        if obj.isPointer():
            objValue = self.selectById(objValue)
        obj.set(objValue)
        
    def fillObject(self, obj, objRow):
        self.fillObjectWithRow(obj, objRow)
        namePrefix = self.__getNamePrefix(obj)
        childs = self.db.selectObjectsByAncestor(namePrefix)
        childsDict = {obj._objId: obj}
        
        for childRow in childs:
            childParts = childRow['name'].split('.')
            childName = childParts[-1]
            parentId = int(childParts[-2])
            # Here we are assuming that always the parent have
            # been processed first, so it will be in the dictiorary
            parentObj = childsDict[parentId]
            
            childObj = getattr(parentObj, childName, None)
            if childObj is None:
                childObj = self._buildObject(childRow['classname'])
                setattr(parentObj, childName, childObj)
            self.fillObjectWithRow(childObj, childRow)  
            childsDict[childObj._objId] = childObj  
              
    def __objectsFromRows(self, objRows):
        """Create a set of object from a set of rows"""
        objs = []
        for objRow in objRows:
            obj = self._buildObject(objRow['classname'])
            self.fillObject(obj, objRow)
            objs.append(obj)
        return objs
        
    def selectBy(self, **args):
        """Select object meetings some criterias"""
        objRows = self.db.selectObjectsBy(**args)
        return self.__objectsFromRows(objRows)
    
    def selectByClass(self, className, includeSubclasses=True):
        if includeSubclasses:
            from pyworkflow.utils.reflection import getSubclasses
            whereStr = "classname='%s'" % className
            base = self.dictClasses.get(className)
            subDict = getSubclasses(base, self.dictClasses)
            for k, v in subDict.iteritems():
                if issubclass(v, base):
                    whereStr += " OR classname='%s'" % k
            objRows = self.db.selectObjectsWhere(whereStr)
            return self.__objectsFromRows(objRows)
        else:
            return self.selectBy(classname=className)
            
    
    def selectAll(self):
        objRows = self.db.selectObjectsByParent(parent_id=None)
        return self.__objectsFromRows(objRows)


class SqliteDb():
    """Class to handle a Sqlite database.
    It will create connection, execute queries and commands"""
    
    SELECT = "SELECT id, parent_id, name, classname, value FROM Objects WHERE "
    DELETE = "DELETE FROM Objects WHERE "
    
    def selectCmd(self, whereStr, orderByStr=' ORDER BY id'):
        return self.SELECT + whereStr + orderByStr
    
    def __init__(self, dbName, timeout=1000):
        self.__createConnection(dbName, timeout)
        self.__createTables()

    def __createConnection(self, dbName, timeout):
        """Establish db connection"""
        from sqlite3 import dbapi2 as sqlite
        self.connection = sqlite.Connection(dbName, timeout)
        self.connection.row_factory = sqlite.Row
        self.cursor = self.connection.cursor()
        # Define some shortcuts functions
        self.executeCommand = self.cursor.execute
        self.commit = self.connection.commit
        
#    def executeCommand(self, cmd, *args, **kargs):
#        print "executing: ", cmd
#        self.cursor.execute(cmd, *args, **kargs)
        
    def __createTables(self):
        """Create requiered tables if don't exists"""
        # Enable foreings keys
        self.executeCommand("PRAGMA foreign_keys=ON")
    
        self.executeCommand("""CREATE TABLE IF NOT EXISTS Objects
                     (id        INTEGER PRIMARY KEY AUTOINCREMENT,
                      parent_id INTEGER REFERENCES Objects(id),
                      name      TEXT,               -- object name 
                      classname TEXT,               -- object's class name
                      value     TEXT DEFAULT NULL   -- object value, used for Scalars
                      )""")
        self.commit()
        
    def insertObject(self, name, classname, value, parent_id):
        """Execute command to insert a new object. Return the inserted object id"""
        self.executeCommand("INSERT INTO Objects (parent_id, name, classname, value) VALUES (?, ?, ?, ?)", \
                            (parent_id, name, classname, value))
        return self.cursor.lastrowid
    
    def updateObject(self, objId, name, classname, value, parent_id):
        """Update object data """
        self.executeCommand("UPDATE Objects SET parent_id=?, name=?,classname=?, value=? WHERE id=?", \
                            (parent_id, name, classname, value, objId))
        
    def selectObjectById(self, objId):
        """Select an object give its id"""
        self.executeCommand(self.selectCmd("id=?"), (objId,))  
        return self.cursor.fetchone()
    
    def selectObjectsByParent(self, parent_id=None):
        """Select object with a given parent
        if the parent_id is None, all object with parent_id NULL
        will be returned"""
        if parent_id is None:
            self.executeCommand(self.selectCmd("parent_id is NULL"))
        else:
            self.executeCommand(self.selectCmd("parent_id=?"), (parent_id,))
        return self.cursor.fetchall()  
    
    def selectObjectsByAncestor(self, ancestor_namePrefix):
        """Select all objects in the hierachy of ancestor_id"""
        self.executeCommand(self.selectCmd("name LIKE '%s.%%'" % ancestor_namePrefix))
        return self.cursor.fetchall()          
    
    def selectObjectsBy(self, **args):     
        """More flexible select where the constrains can be passed
        as a dictionary, the concatenation is done by an AND"""
        whereList = ['%s=?' % k for k in args.keys()]
        whereStr = ' AND '.join(whereList)
        whereTuple = tuple(args.values())
        self.executeCommand(self.selectCmd(whereStr), whereTuple)
        return self.cursor.fetchall()
    
    def selectObjectsWhere(self, whereStr):
        self.executeCommand(self.selectCmd(whereStr))
        return self.cursor.fetchall()   
    
    def deleteObject(self, objId):
        """Delete an existing object"""
        self.executeCommand(self.DELETE + "id=?", (objId,))
        
    def deleteChildObjects(self, ancestor_namePrefix):
        """ Delete from db all objects that are childs 
        of an ancestor, now them will have the same starting prefix"""
        self.executeCommand(self.DELETE + "name LIKE '%s.%%'" % ancestor_namePrefix)
        
        
        
