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
# !!!! insert
    def insert(self,obj):
        pass

# !!!! commit
    def commit(self):
        pass

# !!!! selectAll
    def selectAll(self):
        pass

# There are some python libraries for interfacing Postgres, see
# http://wiki.postgresql.org/wiki/Python
# Initial election: psycopg (http://initd.org/psycopg/)

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

    def connect(self, database,user,password,host,port=5432):
        self.connection = psycopg2.connect(database=database, user=user, password=password, host=host, port=port)

    # auto-closing could be implemented with atexit, @see
    # http://stackoverflow.com/questions/974813/cleaning-up-an-internal-pysqlite-connection-on-object-destruction
    def close(self):
        self.connection.close()

    # !!!! parse the config file using ConfigMapper or XMLmapper
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
        
