# **************************************************************************
# *
# * Authors:     Jesus Cuenca-Alba (jcuenca@cnb.csic.es)
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

"""
Map objects to PostgreSql backend
"""

import mapper
from pyworkflow.utils.path import replaceExt, joinExt


FIELD_CLASSNAME="classname"

class PostgresqlMapper(mapper.Mapper):
    """Specific Mapper implementation using postgresql database"""
    def __init__(self ,dbSettingsFile, dictClasses=None,database=None):
        if dictClasses == None:
            dictClasses=self.__getClassesDictionary()
        mapper.Mapper.__init__(self,dictClasses)
        self.__initObjDict()
        if database== None:
            try:
                self.db=PostgresqlDb(dbSettingsFile)
            except Exception, ex:
                raise Exception('Error creating PostgresqlMapper, settings file: %s\n error: %s' % (dbSettingsFile, ex))
        else:
            self.db=database
            
    def __getClassesDictionary(self):
        """ Return a dictionary with all the relevant classes inheriting from Object"""
        # It might produce a warning, but it's not a problem, @see http://stackoverflow.com/questions/3571514/python-why-should-from-module-import-be-prohibited
        from pyworkflow.object import *
        from pyworkflow.em.data import *
        return locals()


    def __getNamePrefix(self, obj):
        if len(obj._objName) > 0 and '.' in obj._objName:
            return replaceExt(obj._objName, obj.strId())
        return obj.strId()


    def __getObjectValue(self, obj):
        if obj.isPointer() and obj.hasValue():
            return obj.get().strId() # For pointers store the id of referenced object
        return obj.getObjValue()

    # insert methods

    def insert(self,obj):
        return(self.__insert(obj))

    def __insert(self, obj, namePrefix=None):
        obj._objId = self.db.insertObject(obj._objName, obj.getClassName(),
                                          self.__getObjectValue(obj), obj._objParentId)
        sid = obj.strId()
        if namePrefix is None:
            namePrefix = sid
        else:
            namePrefix = joinExt(namePrefix, sid)
        self.insertChilds(obj, namePrefix)
        return(obj._objId)


    def insertChild(self, obj, key, attr, namePrefix=None):
        if namePrefix is None:
            namePrefix = self.__getNamePrefix(obj)
        attr._objName = joinExt(namePrefix, key)
        attr._objParentId = obj._objId
        self.__insert(attr, namePrefix)


    # !!!! rename/refactor as insertChildren
    def insertChilds(self,obj,namePrefix):
        """ Inserts children of an object. """

        # This is also done in insertChild, but getting the prefix here avoids 
        # doing it for every child element
        if namePrefix is None:
            namePrefix = self.__getNamePrefix(obj)

        for key, attr in obj.getAttributesToStore():
            self.insertChild(obj, key, attr, namePrefix)


    def commit(self):
        self.db.commit()


    def updateFrom(self, obj):
        objectData = self.db.selectObjectById(obj._objId)
        self.unpack(objectData, obj)


    # select methods

    def selectAll(self, asIterator=False, objectFilter=None):
        """ Return all the objects matching the filter, as a list or as an iterator"""
        self.__initObjDict()
        objectsData = self.db.selectObjectsByParent(parent_id=None)
        return self.__unpackObjects(objectsData, asIterator, objectFilter)


    def selectById(self, objId):
        """Return the object which id in the database is objId"""
        if objId in self.objDict:
            obj = self.objDict[objId]
        else:
            objectData = self.db.selectObjectById(objId)
            if objectData is None:
                obj = None
            else:
                obj= self.unpack(objectData)
        return obj

    
    def getParent(self, obj):
        """ Retrieve obj's  parent object """
        return self.selectById(obj._objParentId)


    def __initObjDict(self):
        """ Clear the object cache dictionary"""        
        self.objDict = {}
        

    # def __objectsFromRows(self, objRows, asIterator=False, objectFilter=None):
    def __unpackObjects(self, objectsData, asIterator=False, objectFilter=None):
        """Returns a group of objects from a set of rows, either as an iterator or as a list.
           You can apply an optional filter to the objects group"""
        if asIterator:
            return self.__iterUnpackObjects(objectsData, objectFilter)
        else:
            return [obj for obj in self.__iterUnpackObjects(objectsData, objectFilter)]

    # def __iterObjectsFromRows(self, objRows, objectFilter=None):
    def __iterUnpackObjects(self, objectsData, objectFilter=None):
        for objectData in objectsData:
            obj = self.unpack(objectData)

            if objectFilter is None or objectFilter(obj):
                yield obj


    # def fillObject(self, obj, objRow):
    # other name would be restoreObjectAndAttributes
    def unpack(self, objectData, obj=None):
        """Restore the object from objectData, and restore all of its attributes (children objects)
           from the database."""
        if obj == None:
            obj = self._buildObject(objectData[FIELD_CLASSNAME])

        self.restoreObject(obj, objectData)

        namePrefix = self.__getNamePrefix(obj)
        childrenData = self.db.selectObjectsByAncestor(namePrefix)
        
        for childData in childrenData:
            childHierarchy = childData['name'].split('.')
            childName = childHierarchy[-1]
            parentId = int(childHierarchy[-2])

            # Here we are assuming that the parent will always  have
            # been processed first, so it will be in the dictiorary
            parentObj = self.objDict[parentId]
            
            childObj = getattr(parentObj, childName, None)
            # setattr(parentObj, childName, childObj) - it is set in restoreObject...
            self.restoreObject(childObj, childData)  
        return obj



    # def fillObjectWithRow(self, obj, objectData):              
    def restoreObject(self, obj, objectData):
        """Restore only the object itself from objectData (no attributes restoration)"""
        obj._objId = objectData['id']

        self.objDict[obj._objId] = obj

        name = objectData['name']
        if name is None or len(name) <= 0:
            obj._objName = ''#obj.strId()
        else:
            obj._objName = name

        objValue = objectData['value']

        obj._objParentId = objectData['parent_id']

        if obj.isPointer():
            if objValue is not None:
                objValue = self.selectById(int(objValue))

        obj.set(objValue)


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

    # !!!! update methods


# There are some python libraries for interfacing Postgres, see
# http://wiki.postgresql.org/wiki/Python
# Initial election: psycopg (http://initd.org/psycopg/)
# http://initd.org/psycopg/docs/usage.html

import psycopg2
import psycopg2.extras
import xml.etree.ElementTree as ET


# getting connection parameters:
# 1) explicit parameters
# 2) get from config file. It's simple and reasonable for the API 

class PostgresqlDb():
    """PostgreSql internals handling"""

    OBJECT_FIELDS="id, parent_id, name, " + FIELD_CLASSNAME + ", value"
    SELECT = "SELECT " + OBJECT_FIELDS + " FROM Objects WHERE "
    DELETE = "DELETE FROM Objects WHERE "


    def __init__(self,configFile=None):
        if configFile:
            self.connectUsing(configFile)
            self.createTables()

    def connect(self, database,user,password,host,port=5432):
        self.connection = psycopg2.connect(database=database, user=user, password=password, host=host, port=port)
        self.cursor = self.connection.cursor(cursor_factory=psycopg2.extras.DictCursor)


    def close(self):
        # auto-closing could be implemented with atexit, @see
        # http://stackoverflow.com/questions/974813/cleaning-up-an-internal-pysqlite-connection-on-object-destruction
        self.connection.close()


    # !!!! refactor - parse the config file using ConfigMapper or XMLmapper
    def connectUsing(self, configFile):
        """configFile is XML. @see settings/postresql.xml"""

        # print "Using %s" % configFile
        tree = ET.parse(configFile)
        root = tree.getroot()
        user=root.find("PostgresqlConfig/user").text
        host=root.find("PostgresqlConfig/host").text
        database=root.find("PostgresqlConfig/database").text
        password=root.find("PostgresqlConfig/password").text
        port = None
        if root.find("PostgresqlConfig/port"):
            port=root.find("PostgresqlConfig/port").text

        self.connect(database,user,password,host,port)


    def insertObject(self, name, classname, value, parent_id):
        """Execute command to insert a new object. Return the inserted object id"""
        self.cursor.execute("INSERT INTO Objects (id, parent_id, name," + FIELD_CLASSNAME + """, value)
                               VALUES (DEFAULT,%s, %s, %s, %s) RETURNING id""", 
                           (parent_id, name, classname, value))
        insertedObjectId=self.cursor.fetchone()[0]
        self.commit()
        return insertedObjectId


    def createTables(self):
        """Create required tables if they don't exist"""
        self.cursor.execute("""CREATE TABLE IF NOT EXISTS Objects
                     (id        SERIAL PRIMARY KEY,
                      parent_id INTEGER REFERENCES Objects(id),
                      name      TEXT,               -- object name 
                      """ + FIELD_CLASSNAME + """ TEXT,               -- object's class name
                      value     TEXT DEFAULT NULL   -- object value, used for Scalars
                      )""")
        self.commit()


    def commit(self):
        self.connection.commit()


    def lastId(self):
        """ Useful for testing"""
        self.cursor.execute("select max(id) from objects")
        row=self.cursor.fetchone()
        print "last id" + str(row)
        return(row[0])

    
    def selectCmd(self, whereStr, orderByStr=' ORDER BY id'):
        return self.SELECT + whereStr + orderByStr


    def executeSelect(self, whereStr, valueTuple=None, orderByStr=' ORDER BY id'):
        # !!!! handle ProgrammingError, OperationalError...
        self.cursor.execute(self.selectCmd(whereStr,orderByStr),valueTuple)


    def selectObjectsByParent(self, parent_id=None, iterate=False):
        """Select object with a parent matching parent_id, or no parent (if parent_id is None)"""
        if parent_id is None:
            self.executeSelect("parent_id is NULL")
        else:
            self.executeSelect("parent_id=%s", (parent_id,))
        return self._results(iterate)  

    # !!!! remove %, use psycopg filter (@see  deleteChildObjects)
    def selectObjectsByAncestor(self, ancestor_namePrefix, asIterator=False):
        """Select all objects in the hierachy of ancestor"""
        self.executeSelect("name LIKE '%s.%%'" % ancestor_namePrefix)
        return self._results(asIterator)          


    def selectObjectById(self, objId):
        """Select an object give its id"""
        self.executeSelect("id=%s", valueTuple=(objId,))  
        return self.cursor.fetchone()


    def selectObjectsBy(self, asIterator=False, **args):     
        """Select based on a constrain dictionary, where all of them must be fulfilled"""
        whereList = ['%s=%%s' % k for k in args.keys()]
        whereStr = ' AND '.join(whereList)
        whereTuple = tuple(args.values())
        self.executeSelect(whereStr, valueTuple=whereTuple)
        return self._results(asIterator)

    # Risk of SQL injection??
    def selectObjectsWhere(self, whereStr, asIterator=False):
        self.executeSelect(whereStr)
        return self._results(asIterator)   


    def _results(self, asIterator=False):
        """ Return the results of a previous query, either as a list or as an iterator"""
        if asIterator:
            return self._iterResults()
        else:
            return self.cursor.fetchall()


    def _iterResults(self):
        row = self.cursor.fetchone()
        while row is not None:
            yield row
            row = self.cursor.fetchone()

    def executeDelete(self, whereStr, valueTuple=None):
        self.cursor.execute(self.DELETE + whereStr,valueTuple)
        self.commit()

    def deleteObject(self, objId):
        """Delete object matching objId"""
        self.executeDelete("id=%s", (objId,))

        
    def deleteChildObjects(self, ancestor_namePrefix):
        """ Delete all children of ancestor, which share the same starting prefix"""
        self.executeDelete( "name LIKE %s", ("%s.%%" % ancestor_namePrefix,) )

        
    def deleteAll(self):
        """ Delete all Objects from the db. """
        self.executeDelete("true")
