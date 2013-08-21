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
class TestPostgreSqlMapper(unittest.TestCase):
    def setUp(self):
        self.db=None
        self.mapper=None


    def getScipionHome(self):
        if "SCIPION_HOME" not in os.environ:
            raise Exception("SCIPION_HOME is not defined as environment variable")
        return os.environ["SCIPION_HOME"]


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


    def test_insert(self,intValue=22):
       """Test mapper insertion and selection by Id"""
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
            self.assertEqual(object.voltage.get(),200.0)
            return parentId


    def test_selectAll(self):
       mapper=self.getMapper()
       if mapper != None:
           allObjects= mapper.selectAll()
           self.assertNotEqual(len(allObjects),0)



    # !!!! delete test
    def test_delete(self):
        pass


    # !!!! deleteChilds test
    def test_deleteChilds(self):
        pass

    # !!!! delete test
    def test_delete(self):
        pass



    # !!!! updateFrom test
    def test_updateFrom(self):
        pass


class TestPostgreSqlDb(unittest.TestCase):
    def setUp(self):
        self.db=None
        self.mapper=None
  
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


    def getLastId(self):
        return self.getConnection().lastId()


    def test_connectUsing(self):
        db=self.getConnection()
        return db

    def test_createTables(self):
        db=self.getConnection()
        if db != None:
            db.createTables()


    def allChildrenBelongToParent(self,childrenList,parentId):
            return reduce(lambda x,y:  x and y, map(lambda rowDict: parentId in rowDict.values(), childrenList))


    def test_selectObjectsByParent(self):
        db=self.getConnection()
        if db != None:
            parentId=TestPostgreSqlMapper.test_insertChildren()
            childrenList=db.selectObjectsByParent(parentId)
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))


    def test_selectObjectsByAncestor(self):
        db=self.getConnection()
        if db != None:        
            parentId=self.test_insertChildren()
            childrenList=db.selectObjectsByAncestor(str(parentId))
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))

        

    def test_selectBy(self):
        db=self.getConnection()
        if db != None:
            id=self.test_insert(33)
            objects=db.selectObjectsBy(id= id, value="33")
            self.assertEqual(len(objects),1)


    def test_selectWhere(self):
        db=self.getConnection()
        if db != None:
            id=self.test_insert(44)
            objects=db.selectObjectsWhere("id= %s AND value='%d'" %(id,44))
            self.assertEqual(len(objects),1)




    # This test is dangerous if run against a production DB ;-)
    # Hence, if you really want to run it, add "test_" to the function name
    # !!!! check empty table without using mapper
    def DeleteAll(self):
        print "DELETING ALL..."
        mapper=self.getMapper()
        if mapper != None:
            db=self.getConnection()
            db.deleteAll()
            allObjects= mapper.selectAll()
            self.assertEqual(len(allObjects),0)

    def test_dbDeleteObject(self):
        db=self.getConnection()
        if db != None:
            id=self.test_insert(24)
            # No need to check if the object was inserted, it's already checked inside test_insert
            print "Deleting %s" % str(id)
            db.deleteObject(id)
            row = db.selectObjectById(id)
            self.assertIsNone(row)



    def test_dbDeleteChildObjects(self):
        db=self.getConnection()
        if db != None:
            # test_insertChildren checks that things are really inserted
            id=self.test_insertChildren()
            print str(id)
            db.deleteChildObjects(str(id))
            obj_list=db.selectObjectsByParent(id)
            self.assertEqual(len(obj_list),0)
