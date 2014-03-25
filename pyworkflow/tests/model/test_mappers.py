#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_mappers -v
# To run a single test,
# python -m unittest -v test_mappers.TestMappers.test_connectUsing

import os
import os.path
import unittest
from pyworkflow.mapper import *
from pyworkflow.object import *
from pyworkflow.em.data import Acquisition
from pyworkflow.tests import *
import pyworkflow.dataset as ds



# @see test_object.TestPyworkflow.test_SqliteMapper
class TestPostgreSqlMapper(unittest.TestCase):
    mapper=None

    @classmethod
    def setUpClass(cls):
        cls.setupMapper()

    # this class methods allow for the sharing of the connection by all tests

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
                cls.mapper = postgresql.PostgresqlMapper(dbconfig)
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
            micro=Acquisition()
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



    def test_delete(self):
       # deleteChilds test is implicit in delete, for complex objects
       mapper=TestPostgreSqlMapper.mapper
       if mapper != None:
           micro=Acquisition()
           micro.voltage=Float(200.0)
           parentId=mapper.insert(micro)
           voltageId=micro.voltage.getObjId()
           print parentId,voltageId
           mapper.delete(micro)
           object = mapper.selectById(parentId)
           self.assertIsNone(object)
           object = mapper.selectById(voltageId)
           self.assertIsNone(object)



    # This test is dangerous if run against a production DB ;-)
    # Hence, if you really want to run it, add "test_" to the function name
    def deleteAll(self):
        mapper=TestPostgreSqlMapper.mapper
        if mapper != None:
            print "DELETING ALL..."
            mapper.deleteAll()
            allObjects= mapper.selectAll()
            self.assertEqual(len(allObjects),0)


    def test_updates(self):
        mapper=TestPostgreSqlMapper.mapper
        if mapper != None:
            micro=Acquisition()
            micro.voltage=Float(180.0)
            parentId=mapper.insert(micro)
            mapper.commit()
            object = mapper.selectById(parentId)
            self.assertEqual(object.voltage.get(),180.0)

            micro.voltage=Float(160.0)
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
                cls.database= postgresql.PostgresqlDb()
                cls.database.connectUsing(dbconfig)
            else:
                print "Config file %s not found" % dbconfig
                return None
        except Exception as e:
            print str(e)



    def getMapper(self):
        if self.mapper == None and TestPostgreSqlDb.database != None :
            try:
                self.mapper = postgresql.PostgresqlMapper("",database=TestPostgreSqlDb.database)
            except Exception as e:
                print str(e)
        return self.mapper



    def getLastId(self):
        return TestPostgreSqlDb.database.lastId()


    def test_createTables(self):
        db=TestPostgreSqlDb.database
        if db != None:
            db.createTables()


    def insertAcquisition(self):
        mapper=self.getMapper()
        if mapper != None:
            micro=Acquisition()
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
            parentId=self.insertAcquisition()
            childrenList=TestPostgreSqlDb.database.selectObjectsByParent(parentId)
            self.assertTrue(self.allChildrenBelongToParent(childrenList,parentId))


    def test_selectObjectsByAncestor(self):
        if TestPostgreSqlDb.database != None:        
            parentId=self.insertAcquisition()
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
            id=self.insertAcquisition()
            print str(id)
            TestPostgreSqlDb.database.deleteChildObjects(str(id))
            obj_list=TestPostgreSqlDb.database.selectObjectsByParent(id)
            self.assertEqual(len(obj_list),0)


    # This test is dangerous if run against a production DB ;-)
    # Hence, if you really want to run it, add "test_" to the function name
    def DeleteAll(self):
        if TestPostgreSqlDb.database != None:
            print "DELETING ALL..."
            allObjects= TestPostgreSqlDb.database.selectObjectsWhere(None)
            print "Before: %d" % len(allObjects)
            TestPostgreSqlDb.database.deleteAll()
            allObjects= TestPostgreSqlDb.database.selectObjectsWhere(None)
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





class TestSqliteMapper(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')  
        cls.modelGoldSqlite = cls.dataset.getFile( 'modelGoldSqlite')


    def test_SqliteMapper(self):
        fn = self.getOutputPath("basic.sqlite")
        mapper = SqliteMapper(fn)
        # Insert a Complex
        c = Complex.createComplex()
        mapper.insert(c)
        # Insert an Integer
        i = Integer(1)
        mapper.insert(i)
        # Insert two Boolean
        b = Boolean(False)
        b2 = Boolean(True)
        mapper.insert(b)
        mapper.insert(b2)
        #Test storing pointers
        p = Pointer()
        p.set(c)
        mapper.insert(p)
        
        
        # Store list
        strList = ['1', '2', '3']
        csv = CsvList()
        csv += strList
        mapper.insert(csv)

        # Test to add relations
        relName = 'testRelation'
        creator = c
        mapper.insertRelation(relName, creator, i, b)
        mapper.insertRelation(relName, creator, i, b2)
        
        mapper.insertRelation(relName, creator, b, p)
        mapper.insertRelation(relName, creator, b2, p)        
        
        # Save changes to file
        mapper.commit()

        # Reading test
        fnGold = self.modelGoldSqlite
        mapper2 = SqliteMapper(fnGold, globals())
        
        l = mapper2.selectByClass('Integer')[0]
        self.assertEqual(l.get(), 1)
        
        c2 = mapper2.selectByClass('Complex')[0]
        self.assertTrue(c.equalAttributes(c2))
        
        b = mapper2.selectByClass('Boolean')[0]
        self.assertTrue(not b.get())
        
        p = mapper2.selectByClass('Pointer')[0]
        self.assertEqual(c, p.get())
        
        csv2 = mapper2.selectByClass('CsvList')[0]
        self.assertTrue(list.__eq__(csv2, strList))
        
        # Update a CsvList
#        lc = ListContainer()
#        mapper.store(lc)
#        mapper.commit()
#        
#        lc.csv.append('4')
#        lc.csv.append('3')
#        mapper.store(lc)
#        mapper.commit()
#        
#        mapper3 = SqliteMapper(fn, globals())
#        lc3 = mapper3.selectByClass('ListContainer')[0]
#        print 'csv3: ', lc3.csv
        
        # Iterate over all objects
        allObj = mapper2.selectAll()
        iterAllObj = mapper2.selectAll(iterate=True)
        
        for a1, a2 in zip(allObj, iterAllObj):
            self.assertEqual(a1, a2)
            
        # Test relations
        childs = mapper2.getRelationChilds(relName, i)
        parents = mapper2.getRelationParents(relName, p)
        # In this case both childs and parent should be the same
        for c, p in zip(childs, parents):
            self.assertEqual(c, p, "Childs of object i, should be the parents of object p")

        relations = mapper2.getRelationsByCreator(creator)
        for row in relations:
            print row
        
        
    def test_Protocol(self):
        """Test the list with several Complex"""
        fn = self.getOutputPath("protocol.sqlite")   
        mapper = SqliteMapper(fn, globals())
        prot = MyProtocol(mapper=mapper, n=2, workingDir=self.getOutputPath(''))
        prot._stepsExecutor = StepExecutor(hostConfig=None)
        prot.run()
        
        self.assertEqual(prot._steps[0].status, STATUS_FINISHED)
        
        mapper2 = SqliteMapper(fn, globals())
        prot2 = mapper2.selectById(prot.getObjId())
        
        self.assertEqual(prot.endTime, prot2.endTime)
        self.assertEqual(prot._steps[1].status, prot2._steps[1].status)
        
        
        
class TestXmlMapper(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')  
        cls.modelGoldXml = cls.dataset.getFile( 'modelGoldXml')
        
    def test_XMLMapper(self):
        fn = self.getOutputPath("model.xml")
        c = Complex.createComplex()
        mapper = XmlMapper(fn)
        mapper.insert(c)
        #write file
        mapper.commit()

        fnGold = self.modelGoldXml
        #self.assertTrue(filecmp.cmp(fnGold, fn))
        #read file
        mapper2 = XmlMapper(fnGold, globals())
        c2 = mapper2.selectFirst()
        self.assertEquals(c.imag.get(), c2.imag.get())
        
#    def test_zStep(self):
#        fn = self.getTmpPath(self.sqliteFile)
#        s = MyStep()
#        s.x.set(7)
#        s.y.set(3.0)
#        s.status = "KKK"
#        mapper = SqliteMapper(fn, globals())
#        mapper.insert(s)
#        #write file
#        mapper.commit()
#        
#        s2 = mapper.selectByClass('MyStep')[0]
#        self.assertTrue(s.equalAttributes(s2))
        
#    def test_List(self):
#        """Test the list with several Complex"""
#        n = 10
#        l1 = List()
#        for i in range(n):
#            c = Complex(3., 3.)
#            l1.append(c)
#        fn = self.getTmpPath(self.sqliteFile)        
#        mapper = SqliteMapper(fn, globals())
#        mapper.store(l1)
#        mapper.commit()
#        
#        mapper2 = XmlMapper('kk.xml', globals())
#        mapper2.setClassTag('Complex.Float', 'attribute')
#        mapper2.setClassTag('List.ALL', 'class_name')
#        mapper2.setClassTag('MyStep.ALL', 'attribute')
#        mapper2.setClassTag('MyStep.Boolean', 'name_only')
#        step = MyStep()
#        step.b.set('false')
#        step.status = "running"
#        step.inittime = "now"
#        l1.append(step)
#        mapper2.insert(l1)
#        mapper2.commit()
#        
#        mapper3 = SqliteMapper('kk.sqlite', globals())
#        mapper3.insert(l1)
#        mapper3.commit()

class TestDataSet(BaseTest):
    """ Some tests for DataSet implementation. """

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        
    def test_Table(self):
        table = ds.Table(ds.Column('x', int, 5),
                         ds.Column('y', float, 0.0),
                         ds.Column('name', str))
        
        # Add a row to the table
        table.addRow(1, x=12, y=11.0, name='jose')
        table.addRow(2, x=22, y=21.0, name='juan')
        table.addRow(3, x=32, y=31.0, name='pedro')
        # Expect an exception, since name is not provided and have not default
        self.assertRaises(Exception, table.addRow, 100, y=3.0)
        
        row = table.getRow(1)
        print row
        
        self.assertEqual(table.getSize(), 3, "Bad table size")
        
        # Update a value of a row
        table.updateRow(1, name='pepe')        
        row = table.getRow(1)
        print row
        self.assertEqual(row.name, 'pepe', "Error updating name in row")

        print "Table:"
        print table

#    def test_XmippDataSet(self):
#        """ Create a table from a metadata. """
#        from pyworkflow.em.packages.xmipp3 import XmippDataSet
#        import xmipp
#        mdPath = getInputPath('showj', 'tux_vol.xmd')
#
#        xds = XmippDataSet(mdPath)
#        
#        tableNames = xds.listTables()
#        print '\ntables: ', tableNames
#        
#        tn = tableNames[0]
#        
#        table = xds.getTable(tn)
#        
#        print "\nTable '%s':" % tn
#        print table
#        
#        md = xds._convertTableToMd(table)
#        print md