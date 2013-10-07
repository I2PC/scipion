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
    
class ListContainer(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.csv = CsvList() 


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
            self._insertFunctionStep('sleep', i+1, 'sleeping %d'%i)
    
    
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
        
        # CsvList should be empty if set to ''
        l = CsvList()
        l.set('')
        self.assertEqual(len(l), 0)
        
    def test_SqliteMapper(self):
        fn = self.getTmpPath("basic.sqlite")
        mapper = SqliteMapper(fn)
        # Insert a Complex
        c = self.createComplex()
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
        fnGold = self.getGoldPath("basic.sqlite")
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

        relations = mapper2.getRelations(creator)
        for row in relations:
            print row
            row['object_child_id'] = 1
        
    def test_XMLMapper(self):
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
        prot._stepsExecutor = StepExecutor(hostConfig=None)
        prot.run()
        
        self.assertEqual(prot._steps[0].status, STATUS_FINISHED)
        
        mapper2 = SqliteMapper(fn, globals())
        prot2 = mapper2.selectById(prot.getObjId())
        
        self.assertEqual(prot.endTime, prot2.endTime)
        self.assertEqual(prot._steps[1].status, prot2._steps[1].status)
        
        
    def testSimpleFileLog(self):
        import random
        logTestCode = random.randint(1, 100000)
        
        genLogFn = logPath
        log = getGeneralLogger('pyworkflow.test.log.test_scipon_log')        
        genInfoTest = 'General info [' + str(logTestCode) + ']'
        genDebugTest = 'General debug [' + str(logTestCode) + ']'
        genWarningTest = 'General warning [' + str(logTestCode) + ']'
        genErrorTest = 'General error [' + str(logTestCode) + ']'        
        log.info(genInfoTest)
        log.debug(genDebugTest)
        log.warning(genWarningTest)
        
        logFn = self.getTmpPath('fileLog.log')
        log = getFileLogger(logFn)
        fileInfoTest = 'File info [' + str(logTestCode) + ']'
        fileDebugTest = 'File debug [' + str(logTestCode) + ']'
        fileWarningTest = 'File warning [' + str(logTestCode) + ']'
        fileErrorTest = 'File error [' + str(logTestCode) + ']'        
        log.info(fileInfoTest)
        log.debug(fileDebugTest)
        log.warning(fileWarningTest)
        
        log = getGeneralLogger('pyworkflow.test.log.test_scipon_log')
        log.error(genErrorTest)
        
        log = getFileLogger(logFn)
        log.error(fileErrorTest)
        
        # Check general logs
        lineGenInfoTest = getLineInFile(genInfoTest, genLogFn)
        lineGenWarningTest = getLineInFile(genWarningTest, genLogFn)
        lineGenErrorTest = getLineInFile(genErrorTest, genLogFn)
        
        isFileInfoTest = isInFile(fileInfoTest, genLogFn)
        isFileWarningTest = isInFile(fileWarningTest, genLogFn)
        isFileErrorTest = isInFile(fileErrorTest, genLogFn)     
        
        genLoggerChecked = True
        if lineGenInfoTest is None:
            print ('General info log failed!!!')
            genLoggerChecked = False
        if lineGenWarningTest is None:
            print ('General warning log failed!!!')
            genLoggerChecked = False
        if lineGenErrorTest is None:
            print ('General error log failed!!!')
            genLoggerChecked = False
        
        if not((lineGenInfoTest<lineGenWarningTest) & (lineGenWarningTest<lineGenErrorTest)):
            print ('General logs have an incorrect order!!!')
            genLoggerChecked = False
        
        if (isFileInfoTest | isFileWarningTest | isFileErrorTest):
            print ('File logs in general log!!!')
            genLoggerChecked = False
        
        # Check file logs
        lineFileInfoTest = getLineInFile(fileInfoTest, logFn)
        lineFileWarningTest = getLineInFile(fileWarningTest, logFn)
        lineFileErrorTest = getLineInFile(fileErrorTest, logFn)
        
        isGenInfoTest = isInFile(genInfoTest, logFn)
        isGenWarningTest = isInFile(genWarningTest, logFn)
        isGenErrorTest = isInFile(genErrorTest, logFn)    
        
        fileLoggerChecked = True
        if lineFileInfoTest is None:
            print ('File info log failed!!!')
            fileLoggerChecked = False
        if lineFileWarningTest is None:
            print ('File warning log failed!!!')
            fileLoggerChecked = False
        if lineFileErrorTest is None:
            print ('File error log failed!!!')
            fileLoggerChecked = False
        
        if not((lineFileInfoTest<lineFileWarningTest) & (lineFileWarningTest<lineFileErrorTest)):
            print ('File logs have an incorrect order!!!')
            fileLoggerChecked = False
        
        if (isGenInfoTest | isGenWarningTest | isGenErrorTest):
            print ('General logs in file log!!!')
            fileLoggerChecked = False 
        
        self.assertTrue(genLoggerChecked & fileLoggerChecked)  
        
if __name__ == '__main__':
    unittest.main()
    
