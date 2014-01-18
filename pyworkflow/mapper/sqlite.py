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
        self.__initObjDict()
        self.__initUpdateDict()
        try:
            self.db = SqliteDb(dbName)
        except Exception, ex:
            raise Exception('Error creating SqliteMapper, dbName: %s\n error: %s' % (dbName, ex))
    
    def commit(self):
        self.db.commit()
        
    def __getObjectValue(self, obj):
        if obj.isPointer() and obj.hasValue():
            if obj.get().hasObjId(): # Check the object has been stored previously
                return obj.get().strId() # For pointers store the id of referenced object
            else:
                self.updatePendingPointers.append(obj)
                return "Pending update" 
        return obj.getObjValue()
        
    def __insert(self, obj, namePrefix=None):
        obj._objId = self.db.insertObject(obj._objName, obj.getClassName(),
                                          self.__getObjectValue(obj), obj._objParentId,
                                          obj._objLabel, obj._objComment)
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
        
    def deleteAll(self):
        """ Delete all objects stored """
        self.db.deleteAll()
                
    def delete(self, obj):
        """Delete an object and all its childs"""
        self.deleteChilds(obj)
        self.db.deleteObject(obj.getObjId())
    
    def __getNamePrefix(self, obj):
        if len(obj._objName) > 0 and '.' in obj._objName:
            return replaceExt(obj._objName, obj.strId())
        return obj.strId()
    
    def __printObj(self, obj):
        print "obj._objId", obj._objId
        print "obj._objParentId", obj._objParentId
        print "obj._objName", obj._objName
        print "obj.getObjValue()", obj.getObjValue()
    
    def updateTo(self, obj, level=1):
        self.__initUpdateDict()
        self.__updateTo(obj, level)
        # Update pending pointers to objects
        for ptr in self.updatePendingPointers:
            self.db.updateObject(ptr._objId, ptr._objName, ptr.getClassName(),
                             self.__getObjectValue(obj), ptr._objParentId, 
                             ptr._objLabel, ptr._objComment)
        
    def __updateTo(self, obj, level):
        self.db.updateObject(obj._objId, obj._objName, obj.getClassName(),
                             self.__getObjectValue(obj), obj._objParentId, 
                             obj._objLabel, obj._objComment)
        if obj.getObjId() in self.updateDict:
            for k, v in self.updateDict.iteritems():
                print "%d -> %s" % (k, v.getName())
            raise Exception('Circular reference, object: %s found twice' % obj.getName())
        
        self.updateDict[obj._objId] = obj
        for key, attr in obj.getAttributesToStore():
            if attr._objId is None: # Insert new items from the previous state
                attr._objParentId = obj._objId
                #path = obj._objName[:obj._objName.rfind('.')] # remove from last .
                namePrefix = self.__getNamePrefix(obj)
                attr._objName = joinExt(namePrefix, key)
                self.__insert(attr, namePrefix)
            else:  
                self.__updateTo(attr, level + 2)
        
    def updateFrom(self, obj):
        objRow = self.db.selectObjectById(obj._objId)
        self.fillObject(obj, objRow)
            
    def selectById(self, objId):
        """Build the object which id is objId"""
        if objId in self.objDict:
            obj = self.objDict[objId]
        else:
            objRow = self.db.selectObjectById(objId)
            if objRow is None:
                obj = None
            else:
                obj = self._buildObject(objRow['classname'])
                self.fillObject(obj, objRow)
        return obj
    
    def getParent(self, obj):
        """ Retrieve the parent object of another. """
        return self.selectById(obj._objParentId)
    
    def _getStrValue(self, value):
        """ Return empty string if value is None or empty. """
        if value is None or len(value) <= 0:
            return ''
        return value
        
    def fillObjectWithRow(self, obj, objRow):
        """Fill the object with row data"""
        obj._objId = objRow['id']
        self.objDict[obj._objId] = obj
        obj._objName = self._getStrValue(objRow['name'])
        obj._objLabel = self._getStrValue(objRow['label'])
        obj._objComment = self._getStrValue(objRow['comment'])
        objValue = objRow['value']
        obj._objParentId = objRow['parent_id']
        
        if obj.isPointer():
            if objValue is not None:
                objValue = self.selectById(int(objValue))
            else:
                objValue = None
        obj.set(objValue)
        
    def fillObject(self, obj, objRow):
        self.fillObjectWithRow(obj, objRow)
        namePrefix = self.__getNamePrefix(obj)
        childs = self.db.selectObjectsByAncestor(namePrefix)
        #childsDict = {obj._objId: obj}
        
        for childRow in childs:
            childParts = childRow['name'].split('.')
            childName = childParts[-1]
            parentId = int(childParts[-2])
            # Here we are assuming that always the parent have
            # been processed first, so it will be in the dictiorary
            parentObj = self.objDict.get(parentId, None)
            if parentObj is None: # Something went wrong
                raise Exception("Parent object (id=%d) was not found, object: %s" % (parentId, childRow['name']))
            
            childObj = getattr(parentObj, childName, None)
            if childObj is None:
                childObj = self._buildObject(childRow['classname'])
                setattr(parentObj, childName, childObj)
            self.fillObjectWithRow(childObj, childRow)  
            #childsDict[childObj._objId] = childObj  
              
    def __objFromRow(self, objRow):
        obj = self._buildObject(objRow['classname'])
        self.fillObject(obj, objRow)
        return obj
        
    def __iterObjectsFromRows(self, objRows, objectFilter=None):
        for objRow in objRows:
            obj = self.__objFromRow(objRow)
            if objectFilter is None or objectFilter(obj):
                yield obj
        
    def __objectsFromRows(self, objRows, iterate=False, objectFilter=None):
        """Create a set of object from a set of rows
        Params:
            objRows: rows result from a db select.
            iterate: if True, iterates over all elements, if False the whole list is returned
            objectFilter: function to filter some of the objects of the results. 
        """
        if not iterate:
            #return [self.__objFromRow(objRow) for objRow in objRows]
            return [obj for obj in self.__iterObjectsFromRows(objRows, objectFilter)]
        else:
            return self.__iterObjectsFromRows(objRows, objectFilter)
               
    def __initObjDict(self):
        """ Clear the objDict cache """        
        self.objDict = {}
        
    def __initUpdateDict(self):
        """ Clear the updateDict cache """        
        self.updateDict = {}
        # This is used to store pointers that pointed object are not stored yet
        self.updatePendingPointers = [] 
         
    def selectBy(self, iterate=False, objectFilter=None, **args):
        """Select object meetings some criterias"""
        self.__initObjDict()
        objRows = self.db.selectObjectsBy(**args)
        return self.__objectsFromRows(objRows, iterate, objectFilter)
    
    def selectByClass(self, className, includeSubclasses=True, iterate=False, objectFilter=None):
        self.__initObjDict()
        if includeSubclasses:
            from pyworkflow.utils.reflection import getSubclasses
            whereStr = "classname='%s'" % className
            base = self.dictClasses.get(className)
            subDict = getSubclasses(base, self.dictClasses)
            for k, v in subDict.iteritems():
                if issubclass(v, base):
                    whereStr += " OR classname='%s'" % k
            objRows = self.db.selectObjectsWhere(whereStr)
            return self.__objectsFromRows(objRows, iterate, objectFilter)
        else:
            return self.selectBy(iterate=iterate, classname=className)
            
    def selectAll(self, iterate=False, objectFilter=None):
        self.__initObjDict()
        objRows = self.db.selectObjectsByParent(parent_id=None)
        return self.__objectsFromRows(objRows, iterate, objectFilter)    
    
    def insertRelation(self, relName, creatorObj, parentObj, childObj):
        """ This function will add a new relation between two objects.
        Params:
            relName: the name of the relation to be added.
            creatorObj: this object will be the one who register the relation.
            parentObj: this is "parent" in the relation
            childObj: this is "child" in the relation
        """
        for o in [creatorObj, parentObj, childObj]:
            if not o.hasObjId():
                raise Exception("Before adding a relation, the object should be stored in mapper")
        self.db.insertRelation(relName, creatorObj.getObjId(), parentObj.getObjId(), childObj.getObjId())
    
    def __objectsFromIds(self, objIds):
        """Return a list of objects, given a list of id's
        """
        return [self.selectById(rowId['id']) for rowId in objIds]
        
    def getRelationChilds(self, relName, parentObj):
        """ Return all "child" objects for a given relation.
        Params:
            relName: the name of the relation.
            parentObj: this is "parent" in the relation
        Returns: 
            a list of "child" objects.
        """
        childIds = self.db.selectRelationChilds(relName, parentObj.getObjId())
        
        return self.__objectsFromIds(childIds)  
            
    def getRelationParents(self, relName, childObj):
        """ Return all "parent" objects for a given relation.
        Params:
            relName: the name of the relation.
            childObj: this is "child" in the relation
        Returns: 
            a list of "parent" objects.
        """
        parentIds = self.db.selectRelationParents(relName, childObj.getObjId())
        
        return self.__objectsFromIds(parentIds)  

    def getRelations(self, creatorObj):
        """ Return all relations created by creatorObj. """
        return self.db.selectRelationsByCreator(creatorObj.getObjId())
    
    def deleteRelations(self, creatorObj):
        """ Delete all relations created by object creatorObj """
        pass
    
    def insertRelationData(self, relName, creatorId, parentId, childId):
        self.db.insertRelation(relName, creatorId, parentId, childId)
    
    
class SqliteDb():
    """Class to handle a Sqlite database.
    It will create connection, execute queries and commands"""
    
    SELECT = "SELECT id, parent_id, name, classname, value, label, comment FROM Objects WHERE "
    DELETE = "DELETE FROM Objects WHERE "
    
    SELECT_RELATION = "SELECT object_%s_id AS id FROM Relations WHERE name=? AND object_%s_id=?"
    SELECT_RELATIONS = "SELECT * FROM Relations WHERE parent_id=?"
    
    def selectCmd(self, whereStr, orderByStr=' ORDER BY id'):
        return self.SELECT + whereStr + orderByStr
    
    def __init__(self, dbName, timeout=1000):
        self.__createConnection(dbName, timeout)
        self.__createTables()

    def __createConnection(self, dbName, timeout):
        """Establish db connection"""
        from sqlite3 import dbapi2 as sqlite
        self.connection = sqlite.Connection(dbName, timeout, check_same_thread = False)
        self.connection.row_factory = sqlite.Row
        self.cursor = self.connection.cursor()
        # Define some shortcuts functions
        self.executeCommand = self.cursor.execute
        #self.executeCommand = self.__debugExecute
        self.commit = self.connection.commit
        
    def __debugExecute(self, *args):
        print args
        self.cursor.execute(*args)
        
    def __createTables(self):
        """Create requiered tables if don't exists"""
        # Enable foreings keys
        self.executeCommand("PRAGMA foreign_keys=ON")
        # Create the Objects table
        self.executeCommand("""CREATE TABLE IF NOT EXISTS Objects
                     (id        INTEGER PRIMARY KEY AUTOINCREMENT,
                      parent_id INTEGER REFERENCES Objects(id),
                      name      TEXT,                -- object name 
                      classname TEXT,                -- object's class name
                      value     TEXT DEFAULT NULL,   -- object value, used for Scalars
                      label     TEXT DEFAULT NULL,   -- object label, text used for display
                      comment   TEXT DEFAULT NULL    -- object comment, text used for annotations
                      )""")
        # Create the Relations table
        self.executeCommand("""CREATE TABLE IF NOT EXISTS Relations
                     (id        INTEGER PRIMARY KEY AUTOINCREMENT,
                      parent_id INTEGER REFERENCES Objects(id), -- object that created the relation
                      name      TEXT,               -- relation name 
                      classname TEXT DEFAULT NULL,  -- relation's class name
                      value     TEXT DEFAULT NULL,  -- relation value
                      label     TEXT DEFAULT NULL,  -- relation label, text used for display
                      comment   TEXT DEFAULT NULL,  -- relation comment, text used for annotations
                      object_parent_id  INTEGER REFERENCES Objects(id) ON DELETE CASCADE,
                      object_child_id  INTEGER REFERENCES Objects(id) ON DELETE CASCADE 
                      )""")
        self.commit()
        
    def insertObject(self, name, classname, value, parent_id, label, comment):
        """Execute command to insert a new object. Return the inserted object id"""
        try:
            self.executeCommand("INSERT INTO Objects (parent_id, name, classname, value, label, comment) VALUES (?, ?, ?, ?, ?, ?)",
                                (parent_id, name, classname, value, label, comment))
            return self.cursor.lastrowid
        except Exception, ex:
            print "insertObject: ERROR "
            print "INSERT INTO Objects (parent_id, name, classname, value, label, comment) VALUES (?, ?, ?, ?, ?, ?)", (parent_id, name, classname, value, label, comment)
            raise ex
        
    def insertRelation(self, relName, parent_id, object_parent_id, object_child_id, **args):
        """Execute command to insert a new object. Return the inserted object id"""
        self.executeCommand("INSERT INTO Relations (parent_id, name, object_parent_id, object_child_id) VALUES (?, ?, ?, ?)",
                            (parent_id, relName, object_parent_id, object_child_id))
        return self.cursor.lastrowid
    
    def insertRelationRow(self, row):
        """Execute command to insert a new object. Return the inserted object id"""
        return self.insertRelation(row['name'], row['parent_id'], 
                                   row['object_parent_id'], row['object_child_id'])
    
    def updateObject(self, objId, name, classname, value, parent_id, label, comment):
        """Update object data """
        self.executeCommand("UPDATE Objects SET parent_id=?, name=?,classname=?, value=?, label=?, comment=? WHERE id=?",
                            (parent_id, name, classname, value, label, comment, objId))
        
    def selectObjectById(self, objId):
        """Select an object give its id"""
        self.executeCommand(self.selectCmd("id=?"), (objId,))  
        return self.cursor.fetchone()
    
    def _iterResults(self):
        row = self.cursor.fetchone()
        while row is not None:
            yield row
            row = self.cursor.fetchone()
        
    def _results(self, iterate=False):
        """ Return the results to which cursor, point to. 
        If iterates=True, iterate yielding each result independenly"""
        if not iterate:
            return self.cursor.fetchall()
        else:
            return self._iterResults()
        
    def selectObjectsByParent(self, parent_id=None, iterate=False):
        """Select object with a given parent
        if the parent_id is None, all object with parent_id NULL
        will be returned"""
        if parent_id is None:
            self.executeCommand(self.selectCmd("parent_id is NULL"))
        else:
            self.executeCommand(self.selectCmd("parent_id=?"), (parent_id,))
        return self._results(iterate)  
    
    def selectObjectsByAncestor(self, ancestor_namePrefix, iterate=False):
        """Select all objects in the hierachy of ancestor_id"""
        self.executeCommand(self.selectCmd("name LIKE '%s.%%'" % ancestor_namePrefix))
        return self._results(iterate)          
    
    def selectObjectsBy(self, iterate=False, **args):     
        """More flexible select where the constrains can be passed
        as a dictionary, the concatenation is done by an AND"""
        whereList = ['%s=?' % k for k in args.keys()]
        whereStr = ' AND '.join(whereList)
        whereTuple = tuple(args.values())
        self.executeCommand(self.selectCmd(whereStr), whereTuple)
        return self._results(iterate)
    
    def selectObjectsWhere(self, whereStr, iterate=False):
        self.executeCommand(self.selectCmd(whereStr))
        return self._results(iterate)   
    
    def deleteObject(self, objId):
        """Delete an existing object"""
        self.executeCommand(self.DELETE + "id=?", (objId,))
        
    def deleteChildObjects(self, ancestor_namePrefix):
        """ Delete from db all objects that are childs 
        of an ancestor, now them will have the same starting prefix"""
        self.executeCommand(self.DELETE + "name LIKE '%s.%%'" % ancestor_namePrefix)
        
    def deleteAll(self):
        """ Delete all objects from the db. """
        self.executeCommand(self.DELETE + "1")
        
    def selectRelationChilds(self, relName, object_parent_id):
        self.executeCommand(self.SELECT_RELATION % ('child', 'parent'), 
                            (relName, object_parent_id))
        return self._results()
        
    def selectRelationParents(self, relName, object_child_id):
        self.executeCommand(self.SELECT_RELATION % ('parent', 'child'), 
                            (relName, object_child_id))
        return self._results()
    
    def selectRelationsByCreator(self, parent_id):
        self.executeCommand(self.SELECT_RELATIONS, (parent_id,))
        return self._results()
    
    def deleteRelationsByCreator(self, parent_id):
        self.executeCommand("DELETE FROM Relations where parent_id=?", (parent_id,))



        
        
        
