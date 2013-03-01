#!/usr/bin/env python

import random
import unittest
from pyworkflow.mapper.sqlite import SqliteMapper
from pyworkflow.object import Integer, Float, String
from random import randint
import os
from pyworkflow.object import Object
from pyworkflow.mapper.xmlmapper import XmlMapper
import filecmp

class ComplexObject(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.a1 = Float()
        self.a2 = Integer()
    def __str__(self):
        className                = self.getClassName()
        partStr = "%s id = %s" % (className,str(self.id))
        for key, attr in self.getAttributesToStore():
            if attr.hasValue():
                partStr += '\n %s: %s' % (key, attr)
            if attr.pointer:
#                partStr += '\n %s %s' % (key, str(attr.id))
                partStr += '\n %s %s' % (key, str(attr.get().id))
        return partStr

class TestPyworkflow(unittest.TestCase):

    def setUp(self):
        self.seq = range(10)

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
        complexObject = ComplexObject()
        goldStandard = {'a1':1.0,'a2':1}
        complexObject.a1.set(goldStandard['a1'])
        complexObject.a2.set(goldStandard['a2'])
        for key, attr in complexObject.getAttributesToStore():
            self.assertEqual(attr.get(),goldStandard[key])
            #self.assertTrue(element in self.seq)
    def test_SqliteMapper(self):
        import sys, os
        baseFilename = 'SQLMapper.sqlite'
        fileName = os.path.join('/tmp/',baseFilename)
        if os.path.exists(fileName):
            os.remove(fileName)
        
        complexObject = ComplexObject()
        goldStandard = {'a1':1.0,'a2':1}
        complexObject.a1.set(goldStandard['a1'])
        complexObject.a2.set(goldStandard['a2'])
        mapper = SqliteMapper(fileName)
        mapper.insert(complexObject)
        i = Integer(1)
        mapper.insert(i)
        #write file
        mapper.commit()

        scriptDir = os.path.dirname(__file__)
        goldStandard = os.path.join(scriptDir,baseFilename)
        #print goldStandard, fileName
        self.assertTrue(filecmp.cmp(goldStandard, fileName) )
        #read file
        #TODO
#        mapper2 = SqliteMapper('goldStandard', globals())
#        l = mapper2.select(classname='Integer')[0]
#        print l
#        self.assertEqual(l.get(),1)
#        if os.path.exists(fileName):
#            os.remove(fileName)
        
    def test_XML(self):
        import sys, os
        baseFilename = 'XMLMapper.xml'
        fileName = os.path.join('/tmp/',baseFilename)
        if os.path.exists(fileName):
            os.remove(fileName)
        
        complexObject = ComplexObject()
        goldStandard = {'a1':1.0,'a2':1}
        complexObject.a1.set(goldStandard['a1'])
        complexObject.a2.set(goldStandard['a2'])
        mapper = XmlMapper(None)
        mapper.insert(complexObject)
        #write file
        mapper.write(fileName)

        scriptDir = os.path.dirname(__file__)
        goldStandard = os.path.join(scriptDir,baseFilename)
        #print goldStandard, fileName
        self.assertTrue(filecmp.cmp(goldStandard, fileName) )
        #read file
        mapper2 = XmlMapper('goldStandard', globals())
        #TODO
#        l = mapper2.select(classname='Integer')[0]
#        self.assertEqual(l.get(),1)
#        if os.path.exists(fileName):
#            pass#os.remove(fileName)
        
if __name__ == '__main__':
    unittest.main()