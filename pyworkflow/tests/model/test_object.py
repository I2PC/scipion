#!/usr/bin/env python

import os
from os.path import join, dirname, exists
import unittest
import filecmp

from pyworkflow.object import *
from pyworkflow.protocol import *
from pyworkflow.mapper import *
from pyworkflow.utils.log import *
from pyworkflow.utils.utils import getLineInFile, isInFile
from pyworkflow.tests import *

    
class ListContainer(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.csv = CsvList() 

#Protocol for tests, runs in resume mode, and sleeps for??
class MyProtocol(Protocol):
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        self.name = String(args.get('name', None))
        self.numberOfSleeps = Integer(args.get('n', 1))
        self.runMode = Integer(MODE_RESUME)
        
    def sleep(self, t, s):
        log = self._getPath("step_%02d.txt" % t)
        import time 
        time.sleep(t)
        f = open(log, 'w+')
        f.write("Slept: %d seconds\n" % t)
        f.close()
        return [log]
        
    def _insertAllSteps(self):
        for i in range(self.numberOfSleeps.get()):
            self._insertFunctionStep('sleep', i+1, 'sleeping %d'%i)
    
    
class TestPyworkflow(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')  
        cls.modelGoldSqlite = cls.dataset.getFile( 'modelGoldSqlite')
        cls.modelGoldXml = cls.dataset.getFile( 'modelGoldXml')

    def setUp(self):
        #Get the tester.py path
        self.path = dirname(__file__)
        
        self.seq = range(10)
        
        
   
   

            
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
        
        

        
        
if __name__ == '__main__':
    unittest.main()
    
