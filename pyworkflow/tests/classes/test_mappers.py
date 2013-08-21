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
        self.mapper=None

    # !!!! add some asserts to the tests

    def getScipionHome(self):
        if "SCIPION_HOME" not in os.environ:
            raise Exception("SCIPION_HOME is not defined as environment variable")
        return os.environ["SCIPION_HOME"]


    def  getConnection(self):
        if self.db == None:
            try:
                if self.mapper != None and self.mapper.db != None:
                    self.db=self.mapper.db
                else:
                    dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
                    if os.path.isfile(dbconfig):
                        self.db= pyworkflow.mapper.postgresql.PostgresqlDb()
                        self.db.connectUsing(dbconfig)
                    else:
                        print "Config file %s not found" % dbconfig
                        return None
            except Exception as e:
                print str(e)
        return self.db

    def getMapper(self):
        if self.mapper == None:
            try:
                dbconfig= os.path.join(self.getScipionHome() , "postgresql.xml")
                if os.path.isfile(dbconfig):
                    self.mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)
                else:
                    print "Config file %s not found" % dbconfig
                    return None
            except Exception as e:
                print str(e)
        return self.mapper

    def getLastId(self):
        return self.getConnection().lastId()

    def test_PostgresqlMapper(self):
        # Note: general-purpose exception-handling is handled by Pyunit
        mapper = self.getMapper()
        if mapper != None:
            i = Integer(4)
            mapper.insert(i)
            mapper.commit()
            object = mapper.selectById(self.getLastId())
            self.assertEqual(object.get(),4)

    def test_connectUsing(self):
        db=self.getConnection()
        return db

    def test_createTables(self):
        db=self.getConnection()
        if db != None:
            db.createTables()

    def test_insert(self,intValue=22):
       mapper=self.getMapper()
       if mapper != None:
           i = Integer(intValue)
           objectId=mapper.insert(i)
           object = mapper.selectById(objectId)
           self.assertEqual(object.get(),intValue)
           return objectId


    def test_insertChildren(self):
        """ Test mapper insertion of an object with attributes (children)"""
        mapper=self.getMapper()
        if mapper != None:
            micro=Microscope()
            micro.voltage=Float(200.0)
            parentId=mapper.insert(micro)
            mapper.commit()
            object = mapper.selectById(parentId)
            object.printAll()
            self.assertEqual(object.voltage.get(),200.0)
            return parentId

    def allChildrenBelongToParent(self,childrenList,parentId):
            return reduce(lambda x,y:  x and y, map(lambda rowDict: parentId in rowDict.values(), childrenList))


    def test_selectObjectsByParent(self):
        db=self.getConnection()
        if db != None:
            parentId=self.test_insertChildren()
            childrenList=db.selectObjectsByParent(parentId)
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))


    def test_selectObjectsByAncestor(self):
        db=self.getConnection()
        if db != None:        
            parentId=self.test_insertChildren()
            childrenList=db.selectObjectsByAncestor(str(parentId))
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))


    def test_selectById(self):
       mapper=self.getMapper()
       if mapper != None:
           object = mapper.selectById(self.getLastId())
           object.printAll()
        

    def test_selectBy(self):
        db=self.getConnection()
        if db != None:
            objects=db.selectObjectsBy(id= 3, value="4")
            print objects

    def test_selectWhere(self):
        db=self.getConnection()
        if db != None:
            objects=db.selectObjectsWhere("id= 3 AND value='4'")
            print objects


    def test_selectAll(self):
       mapper=self.getMapper()
       if mapper != None:
           for object in mapper.selectAll():
               object.printAll()


    def test_DeleteObject(self):
        db=self.getConnection()
        if db != None:
            id=self.test_insert(24)
            # No need to check if the object was inserted, it's already checked inside test_insert
            print "Deleting %s" % str(id)
            db.deleteObject(id)
            row = db.selectObjectById(id)
            self.assertIsNone(row)

    # !!!! assert that the object was deleted indeed
    def test_DeleteChildObjects(self):
        db=self.getConnection()
        if db != None:
            mapper = self.getMapper()
            # test_insertChildren checks that things are really inserted
            id=self.test_insertChildren()
            print str(id)
            db.deleteChildObjects(str(id))
            obj_list=self.selectObjectsByParent(str(id))
            print len(obj_list)


    # !!!! assert that all was deleted
    # This test is dangerous if run against a production DB ;-)
    # Hence, if you really want to run it, add "test_" to the function name
    def DeleteAll(self):
        print "delete All"
        db=self.getConnection()
        if db != None:
            db.deleteAll()



