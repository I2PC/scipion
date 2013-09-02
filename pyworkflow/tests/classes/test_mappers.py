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
    mapper=None

    @classmethod
    def setUpClass(cls):
        cls.setupMapper()

    @classmethod
    def getScipionHome(self):
        if "SCIPION_HOME" not in os.environ:
            raise Exception("SCIPION_HOME is not defined as environment variable")
        return os.environ["SCIPION_HOME"]

    @classmethod
    def setupMapper(cls):
        try:
            dbconfig= os.path.join(cls.getScipionHome() , "postgresql.xml")
            if os.path.isfile(dbconfig):
                cls.mapper = pyworkflow.mapper.postgresql.PostgresqlMapper(dbconfig)
            else:
                print "Config file %s not found" % (dbconfig,)
                return None
        except Exception as e:
                print str(e)


    def test_insert(self,intValue=22):
       """Test mapper insertion and selection by Id"""
       mapper=TestPostgreSqlMapper.mapper
       if mapper != None:
           i = Integer(intValue)
           objectId=mapper.insert(i)
           object = mapper.selectById(objectId)
           self.assertEqual(object.get(),intValue)
           return objectId


    def test_insertChildren(self):
        """ Test mapper insertion of an object with attributes (children)"""
        mapper=TestPostgreSqlMapper.mapper
        if mapper != None:
            micro=Microscope()
            micro.voltage=Float(200.0)
            parentId=mapper.insert(micro)
            mapper.commit()
            object = mapper.selectById(parentId)
            self.assertEqual(object.voltage.get(),200.0)
            return parentId

    def test_selectById(self):
        # the method is tested in test_insert
        pass

    def test_selectAll(self):
       mapper=TestPostgreSqlMapper.mapper
       if mapper != None:
           allObjects= mapper.selectAll()
           self.assertNotEqual(len(allObjects),0)



    # !!!! delete test
    def test_delete(self):
        pass


    # !!!! deleteChilds test
    def test_deleteChilds(self):
        pass


    # !!!! deleteAll test
    def test_deleteAll(self):
        pass


    def test_updates(self):
        mapper=TestPostgreSqlMapper.mapper
        if mapper != None:
            micro=Microscope()
            micro.voltage=Float(180.0)
            parentId=mapper.insert(micro)
            mapper.commit()
            object = mapper.selectById(parentId)
            self.assertEqual(object.voltage.get(),180.0)

            micro.voltage=Float(160.0)
            # !!!! To test the for in updateTo, we would need an object that allows to add attributes dynamically
            mapper.updateTo(micro)
            object = mapper.selectById(parentId)
            self.assertEqual(object.voltage.get(),160.0)
            micro.voltage=Float(150.0)
            mapper.updateFrom(micro)
            self.assertEqual(micro.voltage.get(),160.0)


class TestPostgreSqlDb(unittest.TestCase):
    database=None
    mapper=None

    @classmethod
    def setUpClass(cls):
        cls.setupDatabase()


    @classmethod
    def getScipionHome(self):
        if "SCIPION_HOME" not in os.environ:
            raise Exception("SCIPION_HOME is not defined as environment variable")
        return os.environ["SCIPION_HOME"]

    @classmethod
    def  setupDatabase(cls):
        try:
            dbconfig= os.path.join(cls.getScipionHome() , "postgresql.xml")
            if os.path.isfile(dbconfig):
                cls.database= pyworkflow.mapper.postgresql.PostgresqlDb()
                cls.database.connectUsing(dbconfig)
            else:
                print "Config file %s not found" % dbconfig
                return None
        except Exception as e:
            print str(e)



    def getMapper(self):
        if self.mapper == None and TestPostgreSqlDb.database != None :
            try:
                self.mapper = pyworkflow.mapper.postgresql.PostgresqlMapper("",database=TestPostgreSqlDb.database)
            except Exception as e:
                print str(e)
        return self.mapper



    def getLastId(self):
        return TestPostgreSqlDb.database.lastId()


    def test_createTables(self):
        db=TestPostgreSqlDb.database
        if db != None:
            db.createTables()


    def insertMicroscope(self):
        mapper=self.getMapper()
        if mapper != None:
            micro=Microscope()
            micro.voltage=Float(200.0)
            parentId=mapper.insert(micro)
            mapper.commit()
            object = mapper.selectById(parentId)
            return parentId



    def insertInteger(self,intValue=22):
       mapper=self.getMapper()
       if mapper != None:
           i = Integer(intValue)
           objectId=mapper.insert(i)
           object = mapper.selectById(objectId)
           return objectId

    def test_insertObject(self):
        # the method is tested in mapper.insert test
        pass



    def allChildrenBelongToParent(self,childrenList,parentId):
            return reduce(lambda x,y:  x and y, map(lambda rowDict: parentId in rowDict.values(), childrenList))


    def test_selectObjectById(self):
        # the method is tested in insertInteger
        pass



    def test_selectObjectsByParent(self):
        if TestPostgreSqlDb.database != None:
            parentId=self.insertMicroscope()
            childrenList=TestPostgreSqlDb.database.selectObjectsByParent(parentId)
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))


    def test_selectObjectsByAncestor(self):
        if TestPostgreSqlDb.database != None:        
            parentId=self.insertMicroscope()
            childrenList=TestPostgreSqlDb.database.selectObjectsByAncestor(str(parentId))
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))

        

    def test_selectBy(self):
        if TestPostgreSqlDb.database != None:
            id=self.insertInteger(33)
            objects=TestPostgreSqlDb.database.selectObjectsBy(id= id, value="33")
            self.assertEqual(len(objects),1)


    def test_selectWhere(self):
        if TestPostgreSqlDb.database != None:
            id=self.insertInteger(44)
            objects=TestPostgreSqlDb.database.selectObjectsWhere("id= %s AND value='%d'" %(id,44))
            self.assertEqual(len(objects),1)



    def test_deleteObject(self):
        if TestPostgreSqlDb.database != None:
            id=self.insertInteger(24)
            print "Deleting %s" % str(id)
            TestPostgreSqlDb.database.deleteObject(id)
            row = TestPostgreSqlDb.database.selectObjectById(id)
            self.assertIsNone(row)



    def test_deleteChildObjects(self):
          if TestPostgreSqlDb.database != None:
            id=self.insertMicroscope()
            print str(id)
            TestPostgreSqlDb.database.deleteChildObjects(str(id))
            obj_list=TestPostgreSqlDb.database.selectObjectsByParent(id)
            self.assertEqual(len(obj_list),0)


    # This test is dangerous if run against a production DB ;-)
    # Hence, if you really want to run it, add "test_" to the function name
    # !!!! check empty table without using mapper
    def DeleteAll(self):
        print "DELETING ALL..."
        mapper=self.getMapper()
        if mapper != None:
            TestPostgreSqlDb.database.deleteAll()
            allObjects= mapper.selectAll()
            self.assertEqual(len(allObjects),0)


    def test_updateObject(self):
        if TestPostgreSqlDb.database != None:
            mapper=self.getMapper()
            if mapper != None:
                i = Integer(66)
                objectId=mapper.insert(i)
                TestPostgreSqlDb.database.updateObject(objectId, i.getName(), i.getClassName(), 67, i.getAttributeValue("parent_id"))
                object = mapper.selectById(objectId)
                self.assertTrue(object.get() == 67)
