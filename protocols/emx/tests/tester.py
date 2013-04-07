#!/usr/bin/env python

import random
import unittest
from emx import *
from emx.emxmapper import *
from random import randint
import os, sys
import filecmp
try:
    import collections
except ImportError:
    sys.stderr.write('Could not import OrderedDict, test not available')

class TestEMX(unittest.TestCase):


    def test_00emxMicrograph(self):
        """ Create micrograph and play a little with it 
        """
        m1 = Emxmicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)
        #reference
        dictPrimaryKeys = collections.OrderedDict([('fileName', 'mic'), 
                                                   ('index', 1)])
        dictAttributes  = collections.OrderedDict([('acceleratingVoltage', 100), 
                        ('defocusU', 1000.0)
                       ,('pixelSpacing__X', 5.6) 
                       ,('pixelSpacing__Y', 5.7)
                       ])
        self.assertEqual(m1.dictPrimaryKeys, dictPrimaryKeys)
        self.assertEqual(m1.dictAttributes, dictAttributes)

    def test_10emxParticle(self):
        """ Create a particle and play a little with it 
        """
        m1 = Emxmicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)

        p1 = Emxparticle('part',1)
        p1.set('defocusU',1000.10)
        p1.set('pixelSpacing__X',5.66)
        p1.set('pixelSpacing__Y',5.77)
        
        p1.setForeignKey(m1)
        #reference
 
        dictPrimaryKeys = collections.OrderedDict([('fileName', 'part'), 
                                                   ('index', 1)]) 
        dictAttributes  = collections.OrderedDict([
          ('defocusU', 1000.1), 
          ('pixelSpacing__X', 5.66), 
          ('pixelSpacing__Y', 5.77)
          ]) 
        
        dictForeignKeys = collections.OrderedDict([('fileName', 'mic'), ('index', 1)])

        #dictForeignKeys
        self.assertEqual(p1.dictPrimaryKeys, dictPrimaryKeys)
        self.assertEqual(p1.dictAttributes, dictAttributes)
        self.assertEqual(p1.dictForeignKeys[MICROGRAPH].dictPrimaryKeys, dictForeignKeys)

    def test_20clear(self):
        """Test clear function
        """
        m1 = Emxmicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.27)
        m1.clear()
        
        m2 = Emxmicrograph('mic',1)
        m2.set('activeFlag',None)
        
        self.assertEqual(m1, m2)
        self.assertFalse(m1.strongEq(m2))

    def test_30read(self):
        #fileName = "massive_million.xml"
        #fileName = "massive_100000.xml"
        fileName = 'EMXread.emx'
        emxData=EmxData()
        xmlMapper = XmlMapper(emxData)
        dataPath = os.path.dirname(os.path.realpath(__file__))
        dataPath = os.path.join(dataPath,fileName)
        xmlMapper.readEMXFile(dataPath)

        m1 = Emxmicrograph('mic',1)
        m1.set('acceleratingVoltage',100.)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)

        m2 = Emxmicrograph('mic',2)
        m2.set('acceleratingVoltage',200.)
        m2.set('defocusUAngle',135.)

        p1 = Emxparticle('parti',1)
        p1.set('boxSize__X',11)
        p1.set('boxSize__Y',33)
        p1.set('defocusU',1000.)
        p1.set('pixelSpacing__X',55.6)
        p1.set('pixelSpacing__Y',55.7)

        p1.set('transformationMatrix__t11',11.1)
        p1.set('transformationMatrix__t12',12.1)
        p1.set('transformationMatrix__t13',13.1)
        p1.set('transformationMatrix__t14',14.1)
                                           
        p1.set('transformationMatrix__t21',21.1)
        p1.set('transformationMatrix__t24',24.1)

        p1.set('transformationMatrix__t31',31.1)
        p1.set('transformationMatrix__t32',32.1)
        p1.set('transformationMatrix__t33',33.1)
        p1.set('transformationMatrix__t34',34.1)

        p1.setForeignKey(m1)

        self.assertTrue(emxData.objLists[MICROGRAPH][0].strongEq(m1))
        self.assertTrue(emxData.objLists[MICROGRAPH][1].strongEq(m2))
        self.assertTrue(emxData.objLists[PARTICLE][0].strongEq(p1))


    def test_35size(self):
        emxData=EmxData()

        p1 = Emxparticle('part',1)
        p1.set('defocusU',1000.10)
        p1.set('pixelSpacing__X',5.66)
        p1.set('pixelSpacing__Y',5.77)
        emxData.addObject(p1)


        m1 = Emxmicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)
        

        for i in range (1,10):
            m1.set('index',i)
            emxData.addObject(m1)
        self.assertEqual(emxData.size(),10)

    def test_40write(self):
        fileName = 'EMXwrite.emx'
        emxData=EmxData()
        xmlMapper = XmlMapper(emxData)
        dataPath = os.path.dirname(os.path.realpath(__file__))
        dataPath = os.path.join(dataPath,fileName)
        for i in range (1,3):
            m1 = Emxmicrograph('mic',i)
            m1.set('acceleratingVoltage',100)
            m1.set('defocusU',1000.)
            m1.set('pixelSpacing__X',5.6)
            m1.set('pixelSpacing__Y',5.7)
            emxData.addObject(m1)

        p1 = Emxparticle('part',1)
        p1.set('defocusU',1000.10)
        p1.set('pixelSpacing__X',5.66)
        p1.set('pixelSpacing__Y',5.77)

        p1.set('transformationMatrix__t11',11)
        p1.set('transformationMatrix__t12',12)
        p1.set('transformationMatrix__t13',13)
        p1.set('transformationMatrix__t14',14)

        p1.set('transformationMatrix__t21',21)
        p1.set('transformationMatrix__t22',22)
        p1.set('transformationMatrix__t23',23)
        p1.set('transformationMatrix__t24',24)

        p1.set('transformationMatrix__t31',31)
        p1.set('transformationMatrix__t32',32)
        p1.set('transformationMatrix__t33',33)
        p1.set('transformationMatrix__t34',34)
        
        p1.addForeignKey(MICROGRAPH,m1)
        emxData.addObject(p1)
        fileName= os.path.join('/tmp','EMXwrite.emx')
        xmlMapper.writeEMXFile(fileName)

        dataPath = os.path.dirname(os.path.realpath(__file__))
        dataPath = os.path.join(dataPath,'EMXwrite.emx')
        #print "kdiff3", fileName,dataPath
        self.assertTrue(filecmp.cmp(fileName,dataPath))
        if os.path.exists(fileName):
            os.remove(fileName)

    def test_50MasiveReadWrite(self):
        """This is not really a test
        but an auxiliary function
        """
        
        fileName = 'EMXMasiveWrite.emx'
        dataPath = os.path.dirname(os.path.realpath(__file__))
        dataPath = os.path.join(dataPath,fileName)
        emxDataW=EmxData()
        emxDataR=EmxData()
        xmlMapperW = XmlMapper(emxDataW)
        xmlMapperR = XmlMapper(emxDataR)

        numberMic=3
        numberPartPerMic = 2
#        numberMic=100
#        numberPartPerMic = 1000
        for i in range (1,numberMic):
            m1 = Emxmicrograph('mic',i)
            m1.set('acceleratingVoltage',100)
            m1.set('defocusU',1000.)
            m1.set('pixelSpacing__X',5.6)
            m1.set('pixelSpacing__Y',5.7)
            emxDataW.addObject(m1)
            for p in range (1,numberPartPerMic):
                p1 = Emxparticle('part',p)
                p1.set('defocusU',1000.10)
                p1.set('pixelSpacing__X',5.66)
                p1.set('pixelSpacing__Y',5.77)
                p1.addForeignKey(MICROGRAPH,m1)
                emxDataW.addObject(p1)
        fileName= os.path.join('/tmp',fileName)
        xmlMapperW.writeEMXFile(fileName)
        xmlMapperR.readEMXFile(fileName)
        xmlMapperR.writeEMXFile(fileName+"2")
        self.assertTrue(filecmp.cmp(fileName,fileName+'2'))
        #print "kdiff3", fileName,fileName+"2"
        if os.path.exists(fileName):
            os.remove(fileName)
        if os.path.exists(ileName+'2'):
            os.remove(ileName+'2')
#        import time
#        start = time.time()
#        for i in range (1,1000):
#            if i%100 == 0:
#                end = time.time()
#                print i, "entries", end - start, "s"


if __name__ == '__main__':
    unittest.main()
