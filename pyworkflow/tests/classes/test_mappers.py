#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_mappers -v
# To run a single test,
# python -m unittest -v test_mappers.TestMappers.test_connectUsing

import os
import os.path
import unittest
import pyworkflow.mapper.postgresql
from pyworkflow.object import *
from pyworkflow.em.data import Microscope

# @see test_object.TestPyworkflow.test_SqliteMapper
class TestMappers(unittest.TestCase):

    def setUp(self):
        self.db=None

    # !!!! add some asserts to the tests

    def getScipionHome(self):
        if "SCIPION_HOME" not in os.environ:
            raise Exception("SCIPION_HOME is not defined as environment variable")

        return os.environ["SCIPION_HOME"]


    def  getConnection(self):
        if self.db == None:
            self.db= pyworkflow.mapper.postgresql.PostgresqlDb()
            dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
            self.db.connectUsing(dbconfig)
        return self.db

    def getLastId(self):
        return self.getConnection().lastId()

    def test_PostgresqlMapper(self):
        # Note: general-purpose exception-handling is handled by Pyunit
        dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
        mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)

        i = Integer(4)

        mapper.insert(i)

        mapper.commit()

        objects = mapper.selectAll()
        self.assertEqual(len(objects),1)

    def test_connectUsing(self):
        db=self.getConnection()
        return db

    def test_createTables(self):
        db=self.getConnection()
        db.createTables()

    def test_insert(self):
       dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
       mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)
       i = Integer(4)
       mapper.insert(i)


    def test_insertChildren(self):
        """ Test mapper insertion of an object with attributes (children)"""
        dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
        mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)
        micro=Microscope()
        micro.voltage=Float(200.0)
        objectId=mapper.insert(micro)
        mapper.commit()
        object = mapper.selectById(objectId)
        object.printAll()
        self.assertEqual(object.voltage.get(),200.0)


    # !!!! actually select some object by its parent id
    def test_selectObjectsByParent(self):
        db=self.getConnection()
        objects=db.selectObjectsByParent()
        print objects

    # !!!! not tested yet
    def test_selectObjectsByAncestor(self):
        db=self.getConnection()
        objects=db.selectObjectsByAncestor("aa")
        print objects


    def test_selectById(self):
       dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
       mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)
       object = mapper.selectById(self.getLastId())
       object.printAll()
        

    def test_selectBy(self):
        db=self.getConnection()
        objects=db.selectObjectsBy(id= 3, value="4")
        print objects

    def test_selectWhere(self):
        db=self.getConnection()
        objects=db.selectObjectsWhere("id= 3 AND value='4'")
        print objects


    def test_selectAll(self):
       dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
       mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)
       for object in mapper.selectAll():
           object.printAll()

    # !!!! assert that the object was deleted indeed
    def test_DeleteObject(self):
        db=self.getConnection()
        db.deleteObject(self.getLastId())


    # !!!! assert that the object was deleted indeed
    def test_DeleteChildObjects(self):
        db=self.getConnection()
        db.deleteChildObjects("qq")

    # !!!! assert that all was deleted
    # This test is dangerous if run against a production DB ;-)
    # Hence, if you really want to run it, add "test_" to the function name
    def DeleteAll(self):
        print "delete All"
        db=self.getConnection()
        db.deleteAll()



