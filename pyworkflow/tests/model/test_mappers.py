#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os
import os.path
import unittest
from pyworkflow.mapper import *
from pyworkflow.object import *
from pyworkflow.config import *
from pyworkflow.em.data import Acquisition, SetOfImages, Image
from pyworkflow.tests import *
import pyworkflow.dataset as ds
import pyworkflow.utils as pwutils
from pyworkflow.mapper.sqlite import SqliteFlatMapper
from pyworkflow.mapper.sqlite_db import SqliteDb




class TestSqliteMapper(BaseTest):
    _labels = [SMALL]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')  
        cls.modelGoldSqlite = cls.dataset.getFile('modelGoldSqlite')

    def test_SqliteDb(self):
        """ Test the SqliteDb class that is used by the sqlite mappers. """
        db = SqliteDb()
        db._createConnection(self.modelGoldSqlite, timeout=1000)
        
        tables = ['Objects', 'Relations']
        self.assertEqual(tables, db.getTables())
        
        # Test getting the version, for the gold file it should be 0
        self.assertEqual(0, db.getVersion())
        
        db.close()
        
    def test_SqliteMapper(self):
        fn = self.getOutputPath("basic.sqlite")
        fnGoldCopy = self.getOutputPath('gold.sqlite')
        fnGold = self.modelGoldSqlite
        
        print ">>> Using db: ", fn
        print "        gold: ", fnGold
        print "   gold copy: ", fnGoldCopy
        
        pwutils.copyFile(fnGold, fnGoldCopy)

        mapper = SqliteMapper(fn)
        # Insert a Complex
        c = Complex.createComplex() # real = 1, imag = 1
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
        
        # Store csv list
        strList = ['1', '2', '3']
        csv = CsvList()
        csv += strList
        mapper.insert(csv)
        
        # Test normal List
        iList = List()
        mapper.insert(iList) # Insert the list when empty        
        
        i1 = Integer(4)
        i2 = Integer(3)
        iList.append(i1)
        iList.append(i2)
        
        mapper.update(iList) # now update with some items inside
        
        pList = PointerList()
        p1 = Pointer()
        p1.set(b)
        p2 = Pointer()
        p2.set(b2)
        pList.append(p1)
        pList.append(p2)
        
        mapper.store(pList)
        

        # Test to add relations
        relName = 'testRelation'
        creator = c
        mapper.insertRelation(relName, creator, i, b)
        mapper.insertRelation(relName, creator, i, b2)
        
        mapper.insertRelation(relName, creator, b, p)
        mapper.insertRelation(relName, creator, b2, p)        
        
        # Save changes to file
        mapper.commit()
        self.assertEqual(1, mapper.db.getVersion())

        # Intentionally keep gold.sqlite as version 0 to check
        # backward compatibility
        db = SqliteDb()
        db._createConnection(fnGold, timeout=1000)
        print "Checking old version is properly read"
        self.assertEqual(0, db.getVersion())
        colNamesGold = [u'id', u'parent_id', u'name', u'classname', 
                        u'value', u'label', u'comment', u'object_parent_id', 
                        u'object_child_id', u'creation']
        colNames = [col[1] for col in db.getTableColumns('Relations')]
        self.assertEqual(colNamesGold, colNames)
        
        # Reading test
        mapper2 = SqliteMapper(fnGoldCopy, globals())
        print "Checking that Relations table is updated and version to 1"
        self.assertEqual(1, mapper2.db.getVersion())
        # Check that the new column is properly added after updated to version 1
        colNamesGold += [u'object_parent_extended', u'object_child_extended']
        colNames = [col[1] for col in mapper2.db.getTableColumns('Relations')]
        self.assertEqual(colNamesGold, colNames)
        
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
        
        # Iterate over all objects
        allObj = mapper2.selectAll()
        iterAllObj = mapper2.selectAll(iterate=True)

        for a1, a2 in zip(allObj, iterAllObj):
            # Note compare the scalar objects, which have a well-defined comparison
            if isinstance(a1, Scalar):
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
            
    def test_StorePointers(self):
        """ Check that pointers are correctly stored. """
        fn = self.getOutputPath("pointers.sqlite")
        
        print ">>> Using db: ", fn

        mapper = SqliteMapper(fn)
        # Insert a Complex
        c = Complex.createComplex() # real = 1, imag = 1
        mapper.insert(c)
        # Insert an Integer
        p1 = Pointer()
        p1.set(c)
        p1.setExtended('real')
        
        mapper.store(c)
        mapper.store(p1)
        
        self.assertAlmostEqual(c.real.get(), p1.get().get())
        
        p1.set(None) # Reset value and check that is stored properly
        
        self.assertIsNone(p1._extended.get())
        mapper.store(p1)
        mapper.commit()
        
        mapper2 = SqliteMapper(fn, globals())
        p2 = mapper2.selectByClass('Pointer')[0]
        
        #Check the mapper was properly stored when
        # set to None and the _extended property cleanned
        self.assertIsNone(p2.get())

    def test_removeFromLists(self):
        """ Check that lists are properly stored after removing some elements. """
        fn = self.getOutputPath("lists.sqlite")

        print ">>> Using db: ", fn

        # Let's create a Mapper to store a simple List containing two integers
        mapper = SqliteMapper(fn, globals())
        iList = List()
        i1 = Integer(4)
        i2 = Integer(3)
        iList.append(i1)
        iList.append(i2)
        # Store the list and commit changes to db, then close db.
        mapper.store(iList)
        mapper.commit()
        mapper.close()

        # Now let's open again the db with a different connection
        # and load the previously stored list
        mapper2 = SqliteMapper(fn, globals())
        iList2 = mapper2.selectByClass('List')[0]
        # Let's do some basic checks
        self.assertEqual(iList2.getSize(), 2)
        self.assertTrue(Integer(4) in iList2)
        self.assertTrue(Integer(3) in iList2)

        # Now remove one of the integers in the list
        # check consistency in the list elements
        iList2.remove(Integer(4))
        self.assertEqual(iList2.getSize(), 1)
        self.assertTrue(Integer(4) not in iList2)
        self.assertTrue(Integer(3) in iList2)
        # Store once again the new list with one element
        mapper2.store(iList2)
        mapper2.commit()
        mapper2.close()

        # Open the db and load the list once again
        mapper3 = SqliteMapper(fn, globals())
        iList3 = mapper3.selectByClass('List')[0]
        # Check the same consistency before it was stored
        self.assertEqual(iList3.getSize(), 1)
        self.assertTrue(Integer(4) not in iList3)
        self.assertTrue(Integer(3) in iList3)

        
class TestSqliteFlatMapper(BaseTest):
    """ Some tests for DataSet implementation. """
    _labels = [SMALL]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.modelGoldSqlite = cls.dataset.getFile( 'micsGoldSqlite')
        
    def test_SqliteFlatDb(self):
        """ Create a SqliteDataset """
        from pyworkflow.mapper.sqlite import SqliteFlatDb
        print ">>> test_SqliteFlatDb: dbName = '%s'" % self.modelGoldSqlite
        db = SqliteFlatDb(self.modelGoldSqlite)
        # Old db version 0
        self.assertEqual(0, db.getVersion())
        # Test the 'self' class name is correctly retrieved
        self.assertEqual('Micrograph', db.getSelfClassName())  
        # Check the count is equal to 3
        self.assertEqual(3, db.count()) 
        db.close()
        
    def test_insertObjects(self):
        dbName = self.getOutputPath('images.sqlite')
        print ">>> test_insertObjects: dbName = '%s'" % dbName
        mapper = SqliteFlatMapper(dbName, globals())
        n = 10
        
        for i in range(n):
            img = Image()
            img.setLocation(i+1, 'images.stk')
            mapper.store(img)
            
        mapper.setProperty('samplingRate', '3.0')
        mapper.setProperty('defocusU', 1000)
        mapper.setProperty('defocusV', 1000)
        mapper.setProperty('defocusU', 2000) # Test update a property value
        mapper.deleteProperty('defocusV') # Test delete a property
        mapper.commit()
        self.assertEqual(1, mapper.db.getVersion())
        mapper.close()
        
        # Test that values where stored properly
        mapper2 = SqliteFlatMapper(dbName, globals())
        
        self.assertTrue(mapper2.hasProperty('samplingRate'))
        self.assertTrue(mapper2.hasProperty('defocusU'))
        self.assertFalse(mapper2.hasProperty('defocusV'))
        
        self.assertEqual(mapper2.getProperty('samplingRate'), '3.0')
        self.assertEqual(mapper2.getProperty('defocusU'), '2000')
        
    def test_downloads(self):
        dbName = self.getOutputPath('downloads.sqlite')
        #dbName = '/tmp/downloads.sqlite'
        
        print ">>> test_downloads: dbName = '%s'" % dbName
        mapper = SqliteFlatMapper(dbName, globals())
        mapper.enableAppend()
        
        n = 10
        
        for i in range(n):
            download = DownloadRecord(fullName='Paco Perez', 
                                      organization='kkkk')
            mapper.store(download)
            
        mapper.commit()
        mapper.close()
         
    def test_emtpySet(self):
        dbName = self.getOutputPath('empty.sqlite')
        #dbName = '/tmp/downloads.sqlite'
        
        print ">>> test empty set: dbName = '%s'" % dbName
        # Check that writing an emtpy set do not fail
        objSet = Set(filename=dbName)      
        objSet.write()
        objSet.close()       
        
        # Now let's try to open an empty set
        objSet = Set(filename=dbName)
        self.assertEqual(objSet.getSize(), 0)
        items = [obj.clone() for obj in objSet]
        self.assertEqual(len(items), 0)
        
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
        

class TestDataSet(BaseTest):
    """ Some tests for DataSet implementation. """
    _labels = [SMALL]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.modelGoldSqlite = cls.dataset.getFile( 'micsGoldSqlite')
        
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
        
