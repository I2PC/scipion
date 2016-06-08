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
        self.assertEqual(value, i.get())
        # compare objects
        i2 = Integer(value)
        self.assertEqual(i, i2)
        
        value = 2.
        f = Float(value)
        self.assertAlmostEqual(value, f.get())
        
        f.multiply(5)
        self.assertAlmostEqual(value*5, f.get())
        
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
        p.setExtended('Name')
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
        
        self.assertFalse(o.pointer.hasExtended(), 
                         'o.pointer should not have extended at this point')
        
        o.pointer.setExtended(7)
        
        self.assertTrue(o.pointer.hasExtended())
        self.assertTrue(o.pointer.hasExtended())
        self.assertEqual(o.pointer.getExtended(), "7")
        
        # Check that the Item 7 of the set is properly
        # retrieved by the pointer after setting the extended to 7
        self.assertEqual(imgSet[7], o.pointer.get())
        
        # Test the keyword arguments of Pointer contructor
        # repeat above tests with new pointer
        ptr = Pointer(value=imgSet, extended=7)
        self.assertTrue(ptr.hasExtended())
        self.assertTrue(ptr.hasExtended())
        self.assertEqual(ptr.getExtended(), "7")
        
        # Check that the Item 7 of the set is properly
        # retrieved by the pointer after setting the extended to 7
        self.assertEqual(imgSet[7], ptr.get())
        
        o2 = OrderedObject()
        o2.outputImages = imgSet
        ptr2 = Pointer()
        ptr2.set(o2)
        # Test nested extended attributes
        ptr2.setExtended('outputImages.7')
        self.assertEqual(imgSet[7], ptr2.get())
        
        # Same as ptr2, but setting extended in constructor
        ptr3 = Pointer(value=o2, extended='outputImages.7')
        self.assertEqual(imgSet[7], ptr3.get())
    
    
    def test_copyAttributes(self):
        """ Check that after copyAttributes, the values
        were properly copied.
        """
        c1 = Complex(imag=10., real=11.)
        c2 = Complex(imag=0., real=1.0001)
        
        # Float values are different, should not be equal
        self.assertFalse(c1.equalAttributes(c2))
        c2.copyAttributes(c1, 'imag', 'real')
        
        self.assertTrue(c1.equalAttributes(c2), 
                        'Complex c1 and c2 have not equal attributes\nc1: %s\nc2: %s\n' % (c1, c2))
        
    def test_equalAttributes(self):
        """ Check that equal attributes function behaves well
        to compare floats with a given precision.
        """
        c1 = Complex(imag=0., real=1.)
        c2 = Complex(imag=0., real=1.0001)
        
        # Since Float precision is 0.001, now c1 and c2
        # should have equal attributes
        self.assertTrue(c1.equalAttributes(c2))
        # Now if we set a more restrictive precision
        # c1 and c2 are not longer equals
        Float.setPrecision(0.0000001)
        self.assertFalse(c1.equalAttributes(c2))

    def test_aggregate(self):
        """test aggregate sql-like functions
        """
        fnDefocusGroups=self.getOutputPath("SetOfDefocusGroups.sqlite")
        setOfDefocus = SetOfDefocusGroup(filename=fnDefocusGroups)
        df1 = DefocusGroup()
        df1.addCTF(CTFModel(defocusU=2000, defocusV=2000))
        df1.addCTF(CTFModel(defocusU=2400, defocusV=2400))
        # At this point defocus Avg should be 2200
        self.assertAlmostEqual(df1.getDefocusAvg(), 2200.)

        df2 = DefocusGroup()
        df2.addCTF(CTFModel(defocusU=3000, defocusV=3000))
        df2.addCTF(CTFModel(defocusU=5000, defocusV=5000))
        self.assertAlmostEqual(df2.getDefocusAvg(), 4000.)

        #TODO: create another test for aggregate
        #operations      = ['min', 'max'] #This should be an enum -> aggregation function
        #operationLabel  = '_defocusMin' #argument of aggregation function
        #groupByLabels   = ['_defocusMin'] # absolute minimum
        
        #print setOfDefocus.aggregate(operations, operationLabel, groupByLabels)
    
    def test_formatString(self):
        """ Test that Scalar objects behave well
        when using string formating such as: %f or %d
        """
        i = Integer(10)
        f = Float(3.345)
        
        s1 = "i = %d, f = %0.3f" % (i, f)
        
        self.assertEqual(s1, "i = 10, f = 3.345")

    def test_Dict(self):
        d = Dict(default='missing')
        d.update({1: 'one', 2: 'two'})

        # Return default value for any non-present key
        self.assertEqual('missing', d[10])

        # Return true for any 'contains' query
        self.assertTrue(100 in d)
        

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
                
    def test_Environ(self):
        """ Test the Environ class with its utilities. """
        from pyworkflow.utils import Environ
        env = Environ({'PATH': '/usr/bin:/usr/local/bin',
                       'LD_LIBRARY_PATH': '/usr/lib:/usr/lib64'
                       })
        env1 = Environ(env)
        env1.set('PATH', '/usr/local/xmipp')
        self.assertEqual(env1['PATH'],'/usr/local/xmipp')
        self.assertEqual(env1['LD_LIBRARY_PATH'], env['LD_LIBRARY_PATH'])
        
        env2 = Environ(env)
        env2.set('PATH', '/usr/local/xmipp', position=Environ.BEGIN)
        self.assertEqual(env2['PATH'],'/usr/local/xmipp' + os.pathsep + env['PATH'])
        self.assertEqual(env2['LD_LIBRARY_PATH'], env['LD_LIBRARY_PATH'])
        
        env3 = Environ(env)
        env3.update({'PATH': '/usr/local/xmipp', 
                     'LD_LIBRARY_PATH': '/usr/local/xmipp/lib'},
                    position=Environ.END)
        self.assertEqual(env3['PATH'],env['PATH'] + os.pathsep +  '/usr/local/xmipp')
        self.assertEqual(env3['LD_LIBRARY_PATH'],env['LD_LIBRARY_PATH'] + os.pathsep +  '/usr/local/xmipp/lib')      
        
        
        
        
            
if __name__ == '__main__':
    unittest.main()
    
