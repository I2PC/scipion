#!/usr/bin/env python

import os
from os.path import join, dirname, exists
import unittest
import filecmp

from pyworkflow.object import *
from pyworkflow.em.data import *
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
        

    def test_String(self):
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
        
    def test_Pointer(self):
        c = Complex.createComplex()
        p = Pointer()
        p.set(c)
        p.setExtendedAttribute('Name')
        c.Name = String('Paquito')
        
        self.assertEqual(p.get(), 'Paquito')
        stackFn = "images.stk"
        mrcsFn = "images.mrcs"
        fn = self.getOutputPath('test_images.sqlite')
        imgSet = SetOfImages(filename=fn)
        imgSet.setSamplingRate(1.0)
        for i in range(10):
            img = Image()
            img.setLocation(i+1, stackFn)
            imgSet.append(img)
         
        imgSet.write()
        
        # Test that image number 7 is correctly retrieved
        # from the set
        img7 = imgSet[7]
        self.assertEqual(img7.getFileName(), stackFn)
        
        # Modify some properties of image 7 to test update
        img7.setFileName(mrcsFn)
        img7.setSamplingRate(2.0)
        imgSet.update(img7)
        # Write changes after the image 7 update
        imgSet.write()
        
        # Read again the set to be able to retrieve elements
        imgSet = SetOfImages(filename=fn)

        # Validate that image7 was properly updated        
        img7 = imgSet[7]
        self.assertEqual(img7.getFileName(), mrcsFn)
            
        o = OrderedObject()
        
        o.pointer = Pointer()
        o.pointer.set(imgSet)
        
        o.refC = o.pointer.get()
        attrNames = [k for k, a in o.getAttributes()]
        # Check that 'refC' should not appear in attributes
        # since it is only an "alias" to an existing pointed value
        self.assertNotIn('refC', attrNames)
        
        o.pointer.setExtendedItemId(7)
        # Check that the Item 7 of the set is properly
        # retrieved by the pointer after setting the extendedItemId to 7
        self.assertEqual(imgSet[7], o.pointer.get())
        
    
class TestUtils(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
            
    def test_ListsFunctions(self):
        """ Test of some methods that retrieve lists from string. """
        from pyworkflow.utils import getListFromValues, getFloatListFromValues, getBoolListFromValues
        
        results = [('2x1 2x2 4 5', None, getListFromValues, ['1', '1', '2', '2', '4', '5']),
                   ('2x1 2x2 4 5', None, getFloatListFromValues, [1., 1., 2., 2., 4., 5.]),
                   ('1 2 3x3 0.5', 8, getFloatListFromValues, [1., 2., 3., 3., 3., 0.5, 0.5, 0.5]),
                   ('3x1 3x0 1', 8, getBoolListFromValues, [True, True, True, False, False, False, True, True]),
                   ]
        for s, n, func, goldList in results:
            l = func(s, length=n)#)
            self.assertAlmostEqual(l, goldList)
            if n:
                self.assertEqual(n, len(l))
        
        
            
if __name__ == '__main__':
    unittest.main()
    
