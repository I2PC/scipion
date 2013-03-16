#!/usr/bin/env python

import os
from os.path import join, dirname, exists
import unittest
import filecmp

from pyworkflow.mapper.sqlite import SqliteMapper
from pyworkflow.object import Integer, Float, String
from pyworkflow.object import Object
from pyworkflow.mapper.xmlmapper import XmlMapper

class Complex(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.imag = Float()
        self.real = Float()
        
    def __str__(self):
        return '(%s, %s)' % (self.imag, self.real)
    
    def __eq__(self, other):
        return self.imag == other.imag and \
            self.real == other.real

class TestPyworkflow(unittest.TestCase):

    def setUp(self):
        #Get the tester.py path
        self.path = dirname(__file__)
        
        self.seq = range(10)
        # Create reference complex values
        self.cGold = complex(1.0, 1.0)
        # Some filenames:
        self.sqliteFile = 'SQLMapper.sqlite'
        self.xmlFile = 'XMLMapper.xml'
        
    def getTestPath(self, filename):
        """Return the path to the pyworkflow/tests dir
        joined with filename"""
        return join(self.path, filename)
    
    def getTmpPath(self, filename):
        """Return the filename in /tmp/ folder.
        If the file exists, it will be deleted"""
        path = join('/tmp', filename)
        if exists(path):
            os.remove(path)        
        return path

    def createComplex(self):
        """Create a Complex object and set
        values with self.cGold standard"""
        c = Complex() # Create Complex object and set values
        c.imag.set(self.cGold.imag)
        c.real.set(self.cGold.real)
        
        return c
            
    def test_Object(self):
        value = 2
        i = Integer(value)
        # make sure the shuffled sequence does not lose any elements
        self.assertEqual(value, i.get())
        value = 2.
        f = Float(value)
        # make sure the shuffled sequence does not lose any elements
        self.assertEqual(value, f.get())
        value = 'thisisanstring'
        s = String(value)
        self.assertEqual(value, s.get())
        self.assertEqual(s.hasValue(), True)
        a = Integer()
        self.assertEqual(a.hasValue(), False)
        c = self.createComplex()
        # Check values are correct
        self.assertTrue(c.imag.get(), self.cGold.imag)
        self.assertTrue(c.real.get(), self.cGold.real)
        
    def test_SqliteMapper(self):
        fn = self.getTmpPath(self.sqliteFile)
        c = self.createComplex()
        mapper = SqliteMapper(fn)
        mapper.insert(c)
        cid = c.id
        i = Integer(1)
        mapper.insert(i)
        #write file
        mapper.commit()

        fnGold = self.getTestPath(self.sqliteFile)

        #self.assertTrue(filecmp.cmp(fnGold, fn))
        #read file
        #TODO
        mapper2 = SqliteMapper(fnGold, globals())
        l = mapper2.select(classname='Integer')[0]
        self.assertTrue(l.get(), 1)
        
        c2 = mapper2.get(cid)
        self.assertTrue(c, c2)

        
    def test_XML(self):
        fn = self.getTmpPath(self.xmlFile)
        c = self.createComplex()
        mapper = XmlMapper(None)
        mapper.insert(c)
        #write file
        mapper.write(fn)

        fnGold = self.getTestPath(self.xmlFile)
        #print goldStandard, fileName
        #self.assertTrue(filecmp.cmp(fnGold, fn))
        #read file
        mapper2 = XmlMapper(globals())
        mapper2.read(fnGold)
        c2 = mapper2.getAll()[0]
        self.assertTrue(c.imag.get(), c2.imag.get())
        #TODO
#        l = mapper2.select(classname='Integer')[0]
#        self.assertEqual(l.get(),1)
#        if os.path.exists(fileName):
#            pass#os.remove(fileName)
        
if __name__ == '__main__':
    unittest.main()