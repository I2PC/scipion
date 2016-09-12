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
from sqlite_db import SqliteDb



class SqliteMapper(Mapper):
    """Specific Mapper implementation using Sqlite database"""
    def __init__(self, dbName, dictClasses=None):
        Mapper.__init__(self, dictClasses)
        self.__initObjDict()
        self.__initUpdateDict()
        try:
            self.db = SqliteObjectsDb(dbName)
        except Exception, ex:
            raise Exception('Error creating SqliteMapper, dbName: %s\n error: %s' % (dbName, ex))
    
    def close(self):
        self.db.close()
        
    def commit(self):
        self.db.commit()
        
    def __getObjectValue(self, obj):
        """ Get the value of the object to be stored.
        We need to handle the special case of pointer, where we should
        store as value the object id of the pointed object.
        """
        value = obj.getObjValue()

        if obj.isPointer() and obj.hasValue():
            if value.hasObjId(): # Check the object has been stored previously
                value = value.strId() # For pointers store the id of referenced object
            else:
                self.updatePendingPointers.append(obj)
                value = "Pending update" 
            
        return value
        
    def __insert(self, obj, namePrefix=None):
        if not hasattr(obj, '_objDoStore'):
            print "MAPPER: object '%s' doesn't seem to be an Object subclass," % obj
            print "       it does not have attribute '_objDoStore'. Insert skipped."
            return 
        obj._objId = self.db.insertObject(obj._objName, obj.getClassName(),
                                          self.__getObjectValue(obj), obj._objParentId,
                                          obj._objLabel, obj._objComment)
        self.updateDict[obj._objId] = obj
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
        if not hasattr(attr, '_objDoStore'):
            print "MAPPER: object '%s' doesn't seem to be an Object subclass," % attr
            print "       it does not have attribute '_objDoStore'. Insert skipped."
            return 

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

        # Delete any child objects that have not been found.
        # This could be the case if some elements (such as inside List)
        # were stored in the database and were removed from the object
        self.db.deleteMissingObjectsByAncestor(self.__getNamePrefix(obj),
                                               self.updateDict.keys())

    def __updateTo(self, obj, level):
        self.db.updateObject(obj._objId, obj._objName, obj.getClassName(),
                             self.__getObjectValue(obj), obj._objParentId, 
                             obj._objLabel, obj._objComment)

        if obj.getObjId() in self.updateDict:
            raise Exception('Circular reference, object: %s found twice'
                            % obj.getName())
        
        self.updateDict[obj._objId] = obj

        for key, attr in obj.getAttributesToStore():
            if attr._objId is None: # Insert new items from the previous state
                attr._objParentId = obj._objId
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
                if obj is not None:
                    self.fillObject(obj, objRow)
        return obj
    
    def getParent(self, obj):
        """ Retrieve the parent object of another. """
        return self.selectById(obj._objParentId)
        
    def fillObjectWithRow(self, obj, objRow):
        """Fill the object with row data"""
        obj._objId = objRow['id']
        self.objDict[obj._objId] = obj
        obj._objName = self._getStrValue(objRow['name'])
        obj._objLabel = self._getStrValue(objRow['label'])
        obj._objComment = self._getStrValue(objRow['comment'])
        obj._objCreation = self._getStrValue(objRow['creation'])
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
                #print "WARNING: Parent object (id=%d) was not found, object: %s. Ignored." % (parentId, childRow['name'])
                continue
            childObj = getattr(parentObj, childName, None)
            if childObj is None:
                childObj = self._buildObject(childRow['classname'])
                # If we have any problem building the object, just ignore it
                if childObj is None:
                    continue
                setattr(parentObj, childName, childObj)
            self.fillObjectWithRow(childObj, childRow)  
            #childsDict[childObj._objId] = childObj  
              
    def __objFromRow(self, objRow):
        objClassName = objRow['classname']
        
        obj = self._buildObject(objClassName)
        if obj is not None:
            self.fillObject(obj, objRow)
        
        return obj
        
    def __iterObjectsFromRows(self, objRows, objectFilter=None):
        for objRow in objRows:
            obj = self.__objFromRow(objRow)
            if (obj is not None and 
                objectFilter is None or objectFilter(obj)):
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
        """Select object meetings some criteria"""
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
    
    def insertRelation(self, relName, creatorObj, parentObj, childObj, 
                       parentExt=None, childExt=None):
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
        self.db.insertRelation(relName, creatorObj.getObjId(), parentObj.getObjId(), childObj.getObjId(), 
                               parentExt, childExt)
    
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

    def getRelationsByCreator(self, creatorObj):
        """ Return all relations created by creatorObj. """
        return self.db.selectRelationsByCreator(creatorObj.getObjId())
    
    def getRelationsByName(self, relationName):
        """ Return all relations stored of a given type. """
        return self.db.selectRelationsByName(relationName)

    def deleteRelations(self, creatorObj):
        """ Delete all relations created by object creatorObj """
        self.db.deleteRelationsByCreator(creatorObj.getObjId())
    
    def insertRelationData(self, relName, creatorId, parentId, childId,
                           parentExtended=None, childExtended=None):
        self.db.insertRelation(relName, creatorId, parentId, childId,
                               parentExtended, childExtended)
    
    
class SqliteObjectsDb(SqliteDb):
    """Class to handle a Sqlite database.
    It will create connection, execute queries and commands"""
    # Maintain the current version of the DB schema
    # useful for future updates and backward compatibility
    # version should be an integer number
    VERSION = 1
    
    SELECT = "SELECT id, parent_id, name, classname, value, label, comment, datetime(creation, 'localtime') as creation FROM Objects WHERE "
    DELETE = "DELETE FROM Objects WHERE "
    DELETE_SEQUENCE = "DELETE FROM SQLITE_SEQUENCE WHERE name='Objects'"
    
    SELECT_RELATION = "SELECT object_%s_id AS id FROM Relations WHERE name=? AND object_%s_id=?"
    SELECT_RELATIONS = "SELECT * FROM Relations WHERE "
    
    
    def selectCmd(self, whereStr, orderByStr=' ORDER BY id'):
        return self.SELECT + whereStr + orderByStr
    
    def __init__(self, dbName, timeout=1000):
        SqliteDb.__init__(self)
        self._createConnection(dbName, timeout)
        self._initialize()

    def _initialize(self):
        """ Create the required tables if needed. """
        tables = self.getTables()
        # Check if the tables have been created or not
        if not tables:
            self.__createTables()
        else:
            self.__updateTables()
        
    def __createTables(self):
        """Create required tables if don't exists"""
        # Enable foreings keys
        self.setVersion(self.VERSION)
        self.executeCommand("PRAGMA foreign_keys=ON")
        # Create the Objects table
        self.executeCommand("""CREATE TABLE IF NOT EXISTS Objects
                     (id        INTEGER PRIMARY KEY AUTOINCREMENT,
                      parent_id INTEGER REFERENCES Objects(id),
                      name      TEXT,                -- object name 
                      classname TEXT,                -- object's class name
                      value     TEXT DEFAULT NULL,   -- object value, used for Scalars
                      label     TEXT DEFAULT NULL,   -- object label, text used for display
                      comment   TEXT DEFAULT NULL,   -- object comment, text used for annotations
                      creation  DATE                 -- creation date and time of the object
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
                      object_child_id  INTEGER REFERENCES Objects(id) ON DELETE CASCADE,
                      creation  DATE,                 -- creation date and time of the object
                      object_parent_extended TEXT DEFAULT NULL, -- extended property to consider internal objects
                      object_child_extended TEXT DEFAULT NULL
                      )""")
        self.commit()
        
    def __updateTables(self):
        """ This method is intended to update the table schema
        in the case of dealing with old database version.
        """
        if self.getVersion() < self.VERSION: # This applies for version 1
            # Add the extra column for pointer extended attribute in Relations table
            # from version 1 on, there is not needed since the table will 
            # already contains this column
            columns = [c[1] for c in self.getTableColumns('Relations')]
            if not 'object_parent_extended' in columns:
                self.executeCommand("ALTER TABLE Relations "
                                    "ADD COLUMN object_parent_extended  TEXT DEFAULT NULL")
            if not 'object_child_extended' in columns:    
                self.executeCommand("ALTER TABLE Relations "
                                    "ADD COLUMN object_child_extended  TEXT DEFAULT NULL")
            self.setVersion(self.VERSION)
        
        
    def insertObject(self, name, classname, value, parent_id, label, comment):
        """Execute command to insert a new object. Return the inserted object id"""
        try:
            self.executeCommand("INSERT INTO Objects (parent_id, name, classname, value, label, comment, creation) VALUES (?, ?, ?, ?, ?, ?, datetime('now'))",
                                (parent_id, name, classname, value, label, comment))
            return self.cursor.lastrowid
        except Exception, ex:
            print "insertObject: ERROR "
            print "INSERT INTO Objects (parent_id, name, classname, value, label, comment, creation) VALUES (?, ?, ?, ?, ?, ?, datetime('now'))", (parent_id, name, classname, value, label, comment)
            raise ex
        
    def insertRelation(self, relName, parent_id, object_parent_id, object_child_id, 
                       object_parent_extended=None, object_child_extended=None, **kwargs):
        """Execute command to insert a new object. Return the inserted object id"""
        self.executeCommand("INSERT INTO Relations "
                            "(parent_id, name, object_parent_id, object_child_id, creation, "
                            "object_parent_extended, object_child_extended) "
                            " VALUES (?, ?, ?, ?, datetime('now'), ?, ?)",
                            (parent_id, relName, object_parent_id, object_child_id,
                             object_parent_extended, object_child_extended))
        return self.cursor.lastrowid
    
    def updateObject(self, objId, name, classname, value, parent_id, label, comment):
        """Update object data """
        self.executeCommand("UPDATE Objects SET parent_id=?, name=?, classname=?, value=?, label=?, comment=? WHERE id=?",
                            (parent_id, name, classname, value, label, comment, objId))
        
    def selectObjectById(self, objId):
        """Select an object give its id"""
        self.executeCommand(self.selectCmd("id=?"), (objId,))  
        return self.cursor.fetchone()
    
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
        """Select all objects in the hierarchy of ancestor_id"""
        self.executeCommand(self.selectCmd("name LIKE '%s.%%'"
                                           % ancestor_namePrefix))
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
        self.executeCommand(self.DELETE + "name LIKE '%s.%%'"
                            % ancestor_namePrefix)

    def selectMissingObjectsByAncestor(self, ancestor_namePrefix,
                                           idList):
        """Select all objects in the hierarchy of ancestor_id"""
        idStr = ','.join(str(i) for i in idList)
        cmd = self.selectCmd("name LIKE '%s.%%' AND id NOT IN (%s) "
                             % (ancestor_namePrefix, idStr))
        self.executeCommand(cmd)
        return self._results(iterate=False)

    def deleteMissingObjectsByAncestor(self, ancestor_namePrefix, idList):
        """Select all objects in the hierarchy of ancestor_id"""
        idStr = ','.join(str(i) for i in idList)
        cmd = "%s name LIKE '%s.%%' AND id NOT IN (%s) " % (self.DELETE,
                                                            ancestor_namePrefix,
                                                            idStr)
        self.executeCommand(cmd)

    def deleteAll(self):
        """ Delete all objects from the db. """
        self.executeCommand(self.DELETE + "1")
        self.executeCommand(self.DELETE_SEQUENCE) # restart the count of ids
        
    def selectRelationChilds(self, relName, object_parent_id):
        self.executeCommand(self.SELECT_RELATION % ('child', 'parent'), 
                            (relName, object_parent_id))
        return self._results()
        
    def selectRelationParents(self, relName, object_child_id):
        self.executeCommand(self.SELECT_RELATION % ('parent', 'child'), 
                            (relName, object_child_id))
        return self._results()
    
    def selectRelationsByCreator(self, parent_id):
        self.executeCommand(self.SELECT_RELATIONS + "parent_id=?", (parent_id,))
        return self._results()
     
    def selectRelationsByName(self, relationName):
        self.executeCommand(self.SELECT_RELATIONS + "name=?", (relationName,))
        return self._results()
       
    def deleteRelationsByCreator(self, parent_id):
        self.executeCommand("DELETE FROM Relations where parent_id=?", (parent_id,))


class SqliteFlatMapper(Mapper):
    """Specific Flat Mapper implementation using Sqlite database"""
    def __init__(self, dbName, dictClasses=None, tablePrefix=''):
        Mapper.__init__(self, dictClasses)
        self._objTemplate = None
        try:
            self.db = SqliteFlatDb(dbName, tablePrefix)
            self.doCreateTables = self.db.missingTables()
            
            if not self.doCreateTables:
                self.__loadObjDict()
        except Exception, ex:
            raise Exception('Error creating SqliteFlatMapper, dbName: %s, tablePrefix: %s\n error: %s' % (dbName, tablePrefix, ex))
    
    def commit(self):
        self.db.commit()
        
    def close(self):
        self.db.close()
        
    def insert(self, obj):
        if self.doCreateTables:
            self.db.createTables(obj.getObjDict(includeClass=True))
            self.doCreateTables = False
        """Insert a new object into the system, the id will be set"""
        self.db.insertObject(obj.getObjId(), obj.isEnabled(), obj.getObjLabel(), obj.getObjComment(), 
                             *obj.getObjDict().values())
        
    def enableAppend(self):
        """ This will allow to append items to existing db. 
        This is by default not allow, since most sets are not 
        modified after creation.
        """
        if not self.doCreateTables:
            obj = self.selectFirst()
            if obj is not None:
                self.db.setupCommands(obj.getObjDict(includeClass=True))
        
    def clear(self):
        self.db.clear()
        self.doCreateTables = True
    
    def deleteAll(self):
        """ Delete all objects stored """
        self.db.deleteAll()
                
    def delete(self, obj):
        """Delete an object and all its childs"""
        self.db.deleteObject(obj.getObjId())
    
    def updateTo(self, obj, level=1):
        """ Update database entry with new object values. """ 
        if self.db.INSERT_OBJECT is None:
            self.db.setupCommands(obj.getObjDict(includeClass=True))
        args = list(obj.getObjDict().values())
        args.append(obj.getObjId())
        self.db.updateObject(obj.isEnabled(), obj.getObjLabel(), obj.getObjComment(), *args)
            
    def selectById(self, objId):
        """Build the object which id is objId"""
        objRow = self.db.selectObjectById(objId)
        if objRow is None:
            obj = None
        else:
            obj = self.__objFromRow(objRow)
        return obj

    def __loadObjDict(self):
        """ Load object properties and classes from db. """
        # Create a template object for retrieving stored ones
        columnList = []
        rows = self.db.getClassRows()
        attrClasses = {}
        self._objBuildList = []

        for r in rows:
            label = r['label_property']

            if label == SELF:
                self._objClassName = r['class_name']
                self._objTemplate = self._buildObject(self._objClassName)
            else:
                # Lets update the database column mapping
                self.db._columnsMapping[label] = r['column_name']
                columnList.append(label)
                attrClasses[label] = r['class_name']
                attrParts = label.split('.')
                attrJoin = ''
                o = self._objTemplate
                for a in attrParts:
                    attrJoin += a
                    attr = getattr(o, a, None)
                    if attr is None:
                        className = attrClasses[attrJoin]
                        self._objBuildList.append((className, attrJoin.split('.')))
                        attr = self._buildObject(className)
                        setattr(o, a, attr)
                    o = attr
                    attrJoin += '.'
        basicRows = 5
        n = len(rows) + basicRows - 1
        self._objColumns = zip(range(basicRows, n), columnList)
         
    def __buildAndFillObj(self):
        obj = self._buildObject(self._objClassName)
        
        for className, attrParts in self._objBuildList:
            o = obj
            for a in attrParts:
                attr = getattr(o, a, None)
                if not attr:
                    setattr(o, a, self._buildObject(className))
                    break
                o = attr
        return obj
        
    def __objFromRow(self, objRow):
        if self._objTemplate is None:
            self.__loadObjDict()
            
        obj = self._objTemplate #self.__buildAndFillObj()
        obj.setObjId(objRow['id'])
        obj.setObjLabel(self._getStrValue(objRow['label']))
        obj.setObjComment(self._getStrValue(objRow['comment']))
        
        try:
            obj.setEnabled(objRow['enabled'])
            obj.setObjCreation(self._getStrValue(objRow['creation']))
        except Exception:
            # THIS SHOULD NOT HAPPENS
            print "WARNING: 'creation' column not found in object: %s" % obj.getObjId()
            print "         db: %s" % self.db.getDbName()
            print "         objRow: ", dict(objRow)
        
        for c, attrName in self._objColumns:
            obj.setAttributeValue(attrName, objRow[c])

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
            return [obj.clone() for obj in self.__iterObjectsFromRows(objRows, objectFilter)]
        else:
            return self.__iterObjectsFromRows(objRows, objectFilter)
         
    def selectBy(self, iterate=False, objectFilter=None, **args):
        """Select object meetings some criteria"""
        objRows = self.db.selectObjectsBy(**args)
        return self.__objectsFromRows(objRows, iterate, objectFilter)
    
    def selectAll(self, iterate=True
                      , objectFilter=None
                      , orderBy='id'
                      , direction='ASC'
                      , where='1'):
        # Just a sanity check for emtpy sets, that doesn't contains 'Properties' table
        if not self.db.hasTable('Properties'):
            return iter([]) if iterate else []
            
        if self._objTemplate is None:
            self.__loadObjDict()
        objRows = self.db.selectAll(orderBy=orderBy,
                                    direction=direction,
                                    where=where)
        
        return self.__objectsFromRows(objRows, iterate, objectFilter) 

    def aggregate(self, operations, operationLabel, groupByLabels=None):
        rows = self.db.aggregate(operations, operationLabel, groupByLabels)
        results = []
        for row in rows:
            values = {}
            for label in operations:
                values[label] = row[label]
            if groupByLabels is not None:
                for label in groupByLabels:
                    values[label] = row[label]
            results.append(values)

        return results
        #convert row to dictionary

    def count(self):
        if self.doCreateTables:
            return 0
        return self.db.count()   
    
    def __objectsFromIds(self, objIds):
        """Return a list of objects, given a list of id's
        """
        return [self.selectById(rowId['id']) for rowId in objIds]
    
    def hasProperty(self, key):
        return self.db.hasProperty(key)
        
    def getProperty(self, key, defaultValue=None):
        return self.db.getProperty(key, defaultValue)
        
    def setProperty(self, key, value):
        return self.db.setProperty(key, value)
    
    def deleteProperty(self, key):
        return self.db.deleteProperty(key)
    
    def getPropertyKeys(self):
        return self.db.getPropertyKeys()
        

SELF = 'self'

class SqliteFlatDb(SqliteDb):
    """Class to handle a Sqlite database.
    It will create connection, execute queries and commands"""
    # Maintain the current version of the DB schema
    # useful for future updates and backward compatibility
    # version should be an integer number
    VERSION = 1
    
    CLASS_MAP = {'Integer': 'INTEGER',
                 'Float': 'REAL',
                 'Boolean': 'INTEGER'
                 }

    def __init__(self, dbName, tablePrefix='', timeout=1000):
        SqliteDb.__init__(self)
        tablePrefix = tablePrefix.strip()
        if tablePrefix and not tablePrefix.endswith('_'): # Avoid having _ for empty prefix
            tablePrefix += '_'
        #NOTE (Jose Miguel, 2014/01/02
        # Reusing connections is a bit dangerous, since it have lead to
        # unexpected and hard to trace errors due to using an out-of-date
        # reused connection. That's why we are changing now the default to False
        # and only setting to True when the tablePrefix is non-empty, which is the
        # case for classes that are different tables in the same db and it logical
        # to reuse the connection.
        if tablePrefix:
            self._reuseConnections = True
        else:
            self._reuseConnections = False#True
        
        self.CHECK_TABLES = "SELECT name FROM sqlite_master WHERE type='table' AND name='%sObjects';" % tablePrefix
        self.SELECT = "SELECT * FROM %sObjects WHERE " % tablePrefix
        self.FROM   = "FROM %sObjects" % tablePrefix
        self.DELETE = "DELETE FROM %sObjects WHERE " % tablePrefix
        self.INSERT_CLASS = "INSERT INTO %sClasses (label_property, column_name, class_name) VALUES (?, ?, ?)" % tablePrefix
        self.SELECT_CLASS = "SELECT * FROM %sClasses;" % tablePrefix
        self.tablePrefix = tablePrefix
        self._createConnection(dbName, timeout)
        self.INSERT_OBJECT = None
        self.UPDATE_OBJECT = None
        self._columnsMapping = {}

        self.INSERT_PROPERTY = "INSERT INTO Properties (key, value) VALUES (?, ?)"
        self.DELETE_PROPERTY = "DELETE FROM Properties WHERE key=?"
        self.UPDATE_PROPERTY = "UPDATE Properties SET value=? WHERE key=?"
        self.SELECT_PROPERTY = "SELECT value FROM Properties WHERE key=?"
        self.SELECT_PROPERTY_KEYS = "SELECT key FROM Properties"

    def hasProperty(self, key):
        """ Return true if a property with this value is registered. """
        # The database not will not have the 'Properties' table when
        # there is not item inserted (ie an empty set)
        if not self.hasTable('Properties'):
            return False
        self.executeCommand(self.SELECT_PROPERTY, (key,))
        result = self.cursor.fetchone()
        return (result is not None)

    def getProperty(self, key, defaultValue=None):
        """ Return the value of a given property with this key.
        If not found, the defaultValue will be returned.
        """
        # The database not will not have the 'Properties' table when
        # there is not item inserted (ie an empty set)
        if not self.hasTable('Properties'):
            return defaultValue
        
        self.executeCommand(self.SELECT_PROPERTY, (key,))
        result = self.cursor.fetchone()

        if result:
            return result['value']
        else:
            return defaultValue

    def setProperty(self, key, value):
        """ Insert or update the property with a value. """
        # Just ignore the set property for empty sets
        if not self.hasTable('Properties'):
            return

        # All properties are stored as string, except for None type
        value = str(value) if value is not None else None

        if self.hasProperty(key):
            self.executeCommand(self.UPDATE_PROPERTY, (value, key))
        else:
            self.executeCommand(self.INSERT_PROPERTY, (key, value))
            
    def getPropertyKeys(self):
        """ Return all properties stored of this object. """
        self.executeCommand(self.SELECT_PROPERTY_KEYS)
        keys = [r[0] for r in self.cursor.fetchall()]
        return keys        

    def deleteProperty(self, key):
        self.executeCommand(self.DELETE_PROPERTY, (key,))

    def selectCmd(self, whereStr, orderByStr=' ORDER BY id'):
        return self.SELECT + whereStr + orderByStr

    def missingTables(self):
        """ Return True is the needed Objects and Classes table are not created yet. """
        self.executeCommand(self.CHECK_TABLES)
        result = self.cursor.fetchone()

        return result is None

    def clear(self):
        self.executeCommand("DROP TABLE IF EXISTS Properties;")
        self.executeCommand("DROP TABLE IF EXISTS %sClasses;" % self.tablePrefix)
        self.executeCommand("DROP TABLE IF EXISTS %sObjects;" % self.tablePrefix)

    def createTables(self, objDict):
        """Create the Classes and Object table to store items of a Set.
        Each object will be stored in a single row.
        Each nested property of the object will be stored as a column value.
        """
        self.setVersion(self.VERSION)
        # Create a general Properties table to store some needed values
        self.executeCommand("""CREATE TABLE IF NOT EXISTS Properties
                     (key       TEXT UNIQUE, -- property key                 
                      value     TEXT  DEFAULT NULL -- property value
                      )""")
        # Create the Classes table to store each column name and type
        self.executeCommand("""CREATE TABLE IF NOT EXISTS %sClasses
                     (id        INTEGER PRIMARY KEY AUTOINCREMENT,
                      label_property      TEXT UNIQUE, --object label                 
                      column_name TEXT UNIQUE,
                      class_name TEXT DEFAULT NULL  -- relation's class name
                      )""" % self.tablePrefix)
        CREATE_OBJECT_TABLE = """CREATE TABLE IF NOT EXISTS %sObjects
                     (id        INTEGER PRIMARY KEY,
                      enabled   INTEGER DEFAULT 1,   -- used to selected/deselect items from a set
                      label     TEXT DEFAULT NULL,   -- object label, text used for display
                      comment   TEXT DEFAULT NULL,   -- object comment, text used for annotations
                      creation  DATE                 -- creation date and time of the object
                      """ % self.tablePrefix

        c = 0
        for k, v in objDict.iteritems():
            colName = 'c%02d' % c
            className = v[0]
            c += 1
            self.executeCommand(self.INSERT_CLASS, (k, colName, className))
            if k != SELF:
                CREATE_OBJECT_TABLE += ',%s  %s DEFAULT NULL' % (colName, self.CLASS_MAP.get(className, 'TEXT'))

        CREATE_OBJECT_TABLE += ')'
        # Create the Objects table
        self.executeCommand(CREATE_OBJECT_TABLE)
        self.commit()
        # Prepare the INSERT and UPDATE commands
        self.setupCommands(objDict)

    def setupCommands(self, objDict):
        """ Setup the INSERT and UPDATE commands base on the object dictionary. """
        self.INSERT_OBJECT = "INSERT INTO %sObjects (id, enabled, label, comment, creation" % self.tablePrefix
        self.UPDATE_OBJECT = "UPDATE %sObjects SET enabled=?, label=?, comment=?" % self.tablePrefix
        c = 0
        for k in objDict:
            colName = 'c%02d' % c
            self._columnsMapping[k] = colName
            c += 1
            if k != SELF:
                self.INSERT_OBJECT += ',%s' % colName
                self.UPDATE_OBJECT += ', %s=?' % colName

        self.INSERT_OBJECT += ") VALUES (?,?,?,?, datetime('now')" + ',?' * (c-1) + ')'
        self.UPDATE_OBJECT += ' WHERE id=?'

    def getClassRows(self):
        """ Create a dictionary with names of the attributes
        of the colums. """
        self.executeCommand(self.SELECT_CLASS)
        return self._results(iterate=False)

    def getSelfClassName(self):
        """ Return the class name of the attribute named 'self'.
        This is the class of the items stored in a Set.
        """
        self.executeCommand(self.SELECT_CLASS)

        for classRow in self._iterResults():
            if classRow['label_property'] == SELF:
                return classRow['class_name']
        raise Exception("Row '%s' was not found in Classes table. " % SELF)

    def insertObject(self, *args):
        """Insert a new object as a row.
        *args: id, label, comment, ...
        where ... is the values of the objDict from which the tables where created.        
        """
        self.executeCommand(self.INSERT_OBJECT, args)

    def updateObject(self, *args):
        """Update object data """
        self.executeCommand(self.UPDATE_OBJECT, args)

    def selectObjectById(self, objId):
        """Select an object give its id"""
        self.executeCommand(self.selectCmd("id=?"), (objId,))
        return self.cursor.fetchone()

    def selectAll(self, iterate=True, orderBy='id', direction='ASC', where='1'):
        # Handle the specials orderBy values of 'id' and 'RANDOM()'
        # other columns names should be mapped to table column
        # such as: _micId -> c04

        def _getRealCol(colName):
            """ Transform the column name taking into account
             special columns such as: id or RANDOM(), and
             getting the mapping translation otherwise.
            """
            if colName in ['id', 'RANDOM()']:
                return colName
            else:
                return self._columnsMapping[colName]

        if isinstance(orderBy, basestring):
            orderByCol = _getRealCol(orderBy)
        elif isinstance(orderBy, list):
            orderByCol = ','.join([_getRealCol(c) for c in orderBy])
        else:
            raise Exception('Invalid type for orderBy: %s' % type(orderBy))

        # Parse the where string to replace the colunm name with
        # the real table column name ( for example: _micId -> c01 )
        # Right now we are asuming a simple where string in the form
        # colName=VALUE
        if '=' in where:
            whereCol = where.split('=')[0]
            whereRealCol = _getRealCol(whereCol)
            whereStr = where.replace(whereCol, whereRealCol)
        else:
            whereStr = where

        cmd = self.selectCmd(whereStr,
                             orderByStr=' ORDER BY %s %s' % (orderByCol, direction))
        self.executeCommand(cmd)
        return self._results(iterate)

    def aggregate(self, operations, operationLabel, groupByLabels=None):
        #let us count for testing
        selectStr = 'SELECT '
        separator = ' '
        #This cannot be like the following line should be expresed in terms
        #of C1, C2 etc....
        for operation in operations:
            selectStr += "%s %s(%s) AS %s"%(separator
                                            , operation
                                            , self._columnsMapping[operationLabel]
                                            , operation)
            separator = ', '
        if groupByLabels is not None:
            groupByStr = 'GROUP BY '
            separator = ' '
            for groupByLabel in groupByLabels:
                groupByCol = self._columnsMapping[groupByLabel]
                selectStr += ", %(groupByCol)s as %(groupByLabel)s" % locals()
                groupByStr += "%s %s" % (separator, groupByCol)
                separator = ', '
        else:
            groupByStr = ' '
        sqlCommand = selectStr +"\n" + self.FROM + "\n" + groupByStr
        self.executeCommand(sqlCommand)
        return self._results(iterate=False)

    def count(self):
        """ Return the number of element in the table. """
        self.executeCommand(self.selectCmd('1').replace('*', 'COUNT(*)'))
        return self.cursor.fetchone()[0]

    # FIXME: Seems to be duplicated and a subset of selectAll
    def selectObjectsBy(self, iterate=False, **args):
        """More flexible select where the constrains can be passed
        as a dictionary, the concatenation is done by an AND"""
        whereList = ['%s=?' % k for k in args.keys()]
        whereStr = ' AND '.join(whereList)
        whereTuple = tuple(args.values())
        self.executeCommand(self.selectCmd(whereStr), whereTuple)
        return self._results(iterate)

    # FIXME: Seems to be duplicated and a subset of selectAll
    # Moreover, it does not translate between "user colums" and
    # "internal" Objects table columns
    def selectObjectsWhere(self, whereStr, iterate=False):
        self.executeCommand(self.selectCmd(whereStr))
        return self._results(iterate)

    def deleteObject(self, objId):
        """Delete an existing object"""
        self.executeCommand(self.DELETE + "id=?", (objId,))

    def deleteAll(self):
        """ Delete all objects from the db. """
        if not self.missingTables():
            self.executeCommand(self.DELETE + "1")
        
