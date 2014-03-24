#!/usr/bin/env python

import os
from os.path import join, dirname, exists
import unittest
import filecmp

from pyworkflow.object import *
from pyworkflow.tests import *

    
class ListContainer(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.csv = CsvList() 


    
    
class TestObject(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')  
        cls.modelGoldSqlite = cls.dataset.getFile( 'modelGoldSqlite')
        cls.modelGoldXml = cls.dataset.getFile( 'modelGoldXml')

            
    def test_Object(self):
        value = 2
        i = Integer(value)
        # make sure the shuffled sequence does not lose any elements
        self.assertEqual(value, i.get())
        # compare objects
        i2 = Integer(value)
        self.assertEqual(i, i2)
        
        value = 2.
        f = Float(value)
        # make sure the shuffled sequence does not lose any elements
        self.assertEqual(value, f.get())
        value = 'thisisanstring'
        s = String(value)
        self.assertEqual(value, s.get())
        self.assertEqual(s.hasValue(), True)
        
        s2 = String()
        # None value is considered empty
        self.assertTrue(s2.empty(), "s2 string should be empty if None")
        s2.set(' ')
        # Only spaces is also empty
        self.assertTrue(s2.empty(), "s2 string should be empty if only spaces")
        s2.set('something')
        # No empty after some value
        self.assertFalse(s2.empty(), "s2 string should not be empty after value")
        
        a = Integer()
        self.assertEqual(a.hasValue(), False)
        c = Complex.createComplex()
        # Check values are correct
        self.assertEqual(c.imag.get(), Complex.cGold.imag)
        self.assertEqual(c.real.get(), Complex.cGold.real)
        
        # Test Boolean logic
        b = Boolean(False)
        self.assertTrue(not b.get())

        b.set('True')
        self.assertTrue(b.get())
        
        b = Boolean()
        b.set(False)
        self.assertTrue(not b.get())
        
        # CsvList should be empty if set to ''
        l = CsvList()
        l.set('')
        self.assertEqual(len(l), 0)
        


        
        
if __name__ == '__main__':
    unittest.main()
    
