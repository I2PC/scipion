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

from mapper import Mapper

class PostgresqlMapper(Mapper):
    """Specific Mapper implementation using postgresql database"""
    def __init__(self ,dbSettingsFile, dictClasses=None):
        Mapper.__init__(self,dictClasses)
        self.db=PostgresqlDb(dbSettingsFile)



    def __getNamePrefix(self, obj):
        if len(obj._objName) > 0 and '.' in obj._objName:
            return replaceExt(obj._objName, obj.strId())
        return obj.strId()


    def __getObjectValue(self, obj):
        if obj.isPointer() and obj.hasValue():
            return obj.get().strId() # For pointers store the id of referenced object
        return obj.getObjValue()


    # !!!! insert
    def insert(self,obj):
        self.__insert(obj)

    def __insert(self, obj, namePrefix=None):
        obj._objId = self.db.insertObject(obj._objName, obj.getClassName(),
                                          self.__getObjectValue(obj), obj._objParentId)
        sid = obj.strId()
        if namePrefix is None:
            namePrefix = sid
        else:
            namePrefix = joinExt(namePrefix, sid)
        self.insertChilds(obj, namePrefix)


    # !!!! refactor namePrefix?


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

    # !!!! commit
    def commit(self):
        pass

    # !!!! @current selectAll
    def selectAll(self):
        


# There are some python libraries for interfacing Postgres, see
# http://wiki.postgresql.org/wiki/Python
# Initial election: psycopg (http://initd.org/psycopg/)
# http://initd.org/psycopg/docs/usage.html

import psycopg2
import xml.etree.ElementTree as ET


# getting connection parameters:
# 1) explicit parameters
# 2) get from config file. It's simple and reasonable for the API 

class PostgresqlDb():
    """PostgreSql internals handling"""

    def __init__(self,configFile=None):
        if configFile:
            self.connectUsing(configFile)
            self.createTables()

    def connect(self, database,user,password,host,port=5432):
        self.connection = psycopg2.connect(database=database, user=user, password=password, host=host, port=port)
        self.cursor = self.connection.cursor()


    def close(self):
        # auto-closing could be implemented with atexit, @see
        # http://stackoverflow.com/questions/974813/cleaning-up-an-internal-pysqlite-connection-on-object-destruction
        self.connection.close()


    # !!!! refactor - parse the config file using ConfigMapper or XMLmapper
    def connectUsing(self, configFile):
        """configFile is XML. @see settings/postresql.xml"""

        print "Using %s" % configFile
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

    # !!!! current insertObject
    def insertObject(self, name, classname, value, parent_id):
        """Execute command to insert a new object. Return the inserted object id"""
        self.cursor.execute("""INSERT INTO Objects (id, parent_id, name, classname, value)
                               VALUES (DEFAULT,%s, %s, %s, %s) RETURNING id""", 
                           (parent_id, name, classname, value))
        return self.cursor.fetchone()[0]
        # !!!! commit?


    def createTables(self):
        """Create required tables if they don't exist"""
        self.cursor.execute("""CREATE TABLE IF NOT EXISTS Objects
                     (id        SERIAL PRIMARY KEY,
                      parent_id INTEGER REFERENCES Objects(id),
                      name      TEXT,               -- object name 
                      classname TEXT,               -- object's class name
                      value     TEXT DEFAULT NULL   -- object value, used for Scalars
                      )""")
        self.commit()


    def commit(self):
        print "commit"
        self.connection.commit()
