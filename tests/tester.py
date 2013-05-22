#!/usr/bin/env python

import os
from os.path import join, dirname, exists
import unittest
import filecmp

from pyworkflow.object import *
from pyworkflow.protocol import *
from pyworkflow.mapper import *

class Complex(Object):
    def __init__(self, imag=0., real=0., **args):
        Object.__init__(self, **args)
        self.imag = Float(imag)
        self.real = Float(real)
        
    def __str__(self):
        return '(%s, %s)' % (self.imag, self.real)
    
    def __eq__(self, other):
        return (self.imag == other.imag and 
                self.real == other.real)
            
    def hasValue(self):
        return True

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
        
    def _defineSteps(self):
        for i in range(self.numberOfSleeps.get()):
            self._insertFunctionStep('sleep', i, 'sleeping %d'%i)
    
    
class TestPyworkflow(unittest.TestCase):

    def setUp(self):
        #Get the tester.py path
        self.path = dirname(__file__)
        
        self.seq = range(10)
        # Create reference complex values
        self.cGold = complex(1.0, 1.0)
        
    def getTestPath(self, *filenames):
        """Return the path to the pyworkflow/tests dir
        joined with filename"""
        return join(self.path, *filenames)
    
    def getTmpPath(self, *filenames):
        """Return the filename in /tmp/ folder.
        If the file exists, it will be deleted"""
        path = self.getTestPath('tmp', *filenames)
        if os.path.exists(path) and not os.path.isdir(path):
            os.remove(path)
        return path
    
    def getGoldPath(self, *filenames):
        return self.getTestPath('gold', *filenames)

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
        a = Integer()
        self.assertEqual(a.hasValue(), False)
        c = self.createComplex()
        # Check values are correct
        self.assertEqual(c.imag.get(), self.cGold.imag)
        self.assertEqual(c.real.get(), self.cGold.real)
        
        # Test Boolean logic
        b = Boolean(False)
        self.assertTrue(not b.get())

        b.set('True')
        self.assertTrue(b.get())
        
        b = Boolean()
        b.set(False)
        self.assertTrue(not b.get())
        
    def test_SqliteMapper(self):
        fn = self.getTmpPath("basic.sqlite")
        c = self.createComplex()
        mapper = SqliteMapper(fn)
        mapper.insert(c)
        cid = c.getObjId()
        i = Integer(1)
        mapper.insert(i)
        b = Boolean(False)
        b2 = Boolean(True)
        mapper.insert(b)
        mapper.insert(b2)
        #Test storing pointers
        p = Pointer()
        p.set(c)
        mapper.insert(p)
        #write file
        mapper.commit()
        
        strList = ['1', '2', '3']
        csv = CsvList()
        csv += strList
        mapper.insert(csv)
        mapper.commit()

        # Reading test
        fnGold = self.getGoldPath("basic.sqlite")
        mapper2 = SqliteMapper(fnGold, globals())
        
        l = mapper2.selectByClass('Integer')[0]
        self.assertEqual(l.get(), 1)
        
        c2 = mapper2.selectByClass('Complex')[0]
        self.assertTrue(c.equalAttributes(c2))
        
        #mapper3 = SqliteMapper(fn, globals())
        b = mapper2.selectByClass('Boolean')[0]
        self.assertTrue(not b.get())
        
        p = mapper2.selectByClass('Pointer')[0]
        self.assertEqual(c, p.get())
        
        csv = mapper2.selectByClass('CsvList')[0]
        self.assertTrue(list.__eq__(csv, strList))
        
        # Iterate over all objects
        allObj = mapper2.selectAll()
        iterAllObj = mapper2.selectAll(iterate=True)
        
        for a1, a2 in zip(allObj, iterAllObj):
            self.assertEqual(a1, a2)

        
    def test_XML(self):
        fn = self.getTmpPath("basic.xml")
        c = self.createComplex()
        mapper = XmlMapper(fn)
        mapper.insert(c)
        #write file
        mapper.commit()

        fnGold = self.getGoldPath("basic.xml")
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

    def test_Protocol(self):
        """Test the list with several Complex"""
        fn = self.getTmpPath("protocol.sqlite")   
        mapper = SqliteMapper(fn, globals())
        prot = MyProtocol(mapper=mapper, n=2, workingDir=self.getTmpPath(''))
        prot._stepsExecutor = StepExecutor()
        prot.run()
        
        self.assertEqual(prot._steps[0].status, STATUS_FINISHED)
        
        mapper2 = SqliteMapper(fn, globals())
        prot2 = mapper2.selectById(prot.getObjId())
        
        self.assertEqual(prot.endTime, prot2.endTime)
        self.assertEqual(prot._steps[1].status, prot2._steps[1].status)
        
if __name__ == '__main__':
    unittest.main()
    