#!/usr/bin/env python
'''
/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
MODIFICATION ADVICE:

Please, do not generate or distribute 
a modified version of this file under its original name. 
 '''
import unittest, os, sys

from unittest import TestResult, _TextTestResult
from protlib_filesystem import getXmippPath
#from test.test_support import unlink
try:
   from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
   from unittest import _WritelnDecorator # Python <2.6
try:
 from xmipp import *
except ImportError:
 print "Running outside xmipp, redefining getXmippPath"
 def getXmippPath(directory, fileName):
     return fileName
 
import random
import unittest
from emx import *
from os.path import join
import filecmp
try:
    import collections
except ImportError:
    sys.stderr.write('Could not import OrderedDict, test not available')

class TestEMX(unittest.TestCase):
    
    _labels = [WEEKLY]
    
    testsPath = getXmippPath("resources", "test")
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        os.chdir(self.testsPath)
        
    def test_00emxMicrograph(self):
        """ Create micrograph and play a little with it 
        """
        m1 = EmxMicrograph('mic',1)
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
        m1 = EmxMicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)

        p1 = EmxParticle('part',1)
        p1.set('defocusU',1000.10)
        p1.set('pixelSpacing__X',5.66)
        p1.set('pixelSpacing__Y',5.77)
        
        p1.setMicrograph(m1)
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
        m1 = EmxMicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.27)
        m1.clear()
        
        m2 = EmxMicrograph('mic',1)
        m2.set('activeFlag',None)
        
        self.assertEqual(m1, m2)
        self.assertFalse(m1.strongEq(m2))

    def test_30read(self):
        #fileName = "massive_million.xml"
        #fileName = "massive_100000.xml"
        fileName = join(self.testsPath,'EMX/EMXread.emx')
        emxData = EmxData()
        emxData.read(fileName)

        m1 = EmxMicrograph('mic',1)
        m1.set('acceleratingVoltage',100.)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)

        m2 = EmxMicrograph('mic',2)
        m2.set('acceleratingVoltage',200.)
        m2.set('defocusUAngle',135.)

        p1 = EmxParticle('parti',1)
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

        p1.setMicrograph(m1)

        self.assertTrue(emxData.objLists[MICROGRAPH][0].strongEq(m1))
        self.assertTrue(emxData.objLists[MICROGRAPH][1].strongEq(m2))
        self.assertTrue(emxData.objLists[PARTICLE][0].strongEq(p1))
        
        
    def test_31_firstObject(self):
        #fileName = "massive_million.xml"
        #fileName = "massive_100000.xml"
        fileName = join(self.testsPath,'EMX/EMXread.emx')
        emxData = EmxData()
        mic = emxData.readFirstObject(MICROGRAPH, fileName)

        m1 = EmxMicrograph('mic', 1)
        self.assertEqual(mic, m1, "first micrograph differ from expected")    
 
        particle = emxData.readFirstObject(PARTICLE, fileName)  
        p1 = EmxParticle('parti', 1)
        self.assertEqual(particle, p1, "first particle differ from expected")  
            
        
    def test_35size(self):
        emxData = EmxData()

        p1 = EmxParticle('part',1)
        p1.set('defocusU',1000.10)
        p1.set('pixelSpacing__X',5.66)
        p1.set('pixelSpacing__Y',5.77)
        emxData.addObject(p1)

        m1 = EmxMicrograph('mic',1)
        m1.set('acceleratingVoltage',100)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)

        for i in range (1,10):
            m1.set('index', i)
            emxData.addObject(m1)
        self.assertEqual(emxData.size(), 10)

    def test_40write(self):
        fileName = join(self.testsPath,'EMX/EMXwrite.emx')
        emxData = EmxData()
        
        for i in range (1,3):
            m1 = EmxMicrograph('mic',i)
            m1.set('acceleratingVoltage',100)
            m1.set('defocusU',1000.)
            m1.set('pixelSpacing__X',5.6)
            m1.set('pixelSpacing__Y',5.7)
            emxData.addObject(m1)

        p1 = EmxParticle('part',1)
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
        
        p1.setMicrograph(m1)
        #TODO: check the following line
        #p1.addForeignKey(MICROGRAPH,m1)
        emxData.addObject(p1)
        fileName2= os.path.join('/tmp','EMXwrite.emx')
        emxData.write(fileName2)
        #print "kdiff3", fileName,fileName2
        self.assertTrue(filecmp.cmp(fileName,fileName2))
        if os.path.exists(fileName2):
            os.remove(fileName2)

    def test_50MasiveReadWrite(self):
        """This is not really a test
        but an auxiliary function
        """
        fileName2 = join(self.testsPath,'EMX/EMXMasiveWrite.emx')
        emxDataW = EmxData()
        emxDataR = EmxData()

        numberMic=3
        numberPartPerMic = 2
#        numberMic=100
#        numberPartPerMic = 1000
        for i in range (1,numberMic):
            m1 = EmxMicrograph('mic',i)
            m1.set('acceleratingVoltage',100)
            m1.set('defocusU',1000.)
            m1.set('pixelSpacing__X',5.6)
            m1.set('pixelSpacing__Y',5.7)
            emxDataW.addObject(m1)
            for p in range (1,numberPartPerMic):
                p1 = EmxParticle('part',p)
                p1.set('defocusU',1000.10)
                p1.set('pixelSpacing__X',5.66)
                p1.set('pixelSpacing__Y',5.77)
                p1.setMicrograph(m1)
                emxDataW.addObject(p1)
        fn = os.path.join('/tmp','EMXMasiveWrite.emx')
        fn2 = fn.replace('.emx', '2.emx')
        emxDataW.write(fn)
        emxDataR.read(fn)
        emxDataR.write(fn2)
        self.assertTrue(filecmp.cmp(fn, fn2))
        #print "kdiff3", fn,fn+"2"
        if os.path.exists(fn):
            os.remove(fn)
        if os.path.exists(fn2):
            os.remove(fn2)

    def test_60iterate(self):
        #fileName = "massive_million.xml"
        #fileName = "massive_100000.xml"
        fileName = join(self.testsPath,'EMX/EMXread.emx')
        emxData = EmxData()
        emxData.read(fileName)
        _list = []
        
        m1 = EmxMicrograph('mic',1)
        m1.set('acceleratingVoltage',100.)
        m1.set('defocusU',1000.)
        m1.set('pixelSpacing__X',5.6)
        m1.set('pixelSpacing__Y',5.7)
        _list.append(m1)

        m2 = EmxMicrograph('mic',2)
        m2.set('acceleratingVoltage',200.)
        m2.set('defocusUAngle',135.)
        _list.append(m2)

        p1 = EmxParticle('parti',1)
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

    def test_70schema(self):
        xmlFile    = join(self.testsPath,'EMX/EMXwrite.emx')
        (code,_out,_err) = validateSchema(xmlFile)
        self.assertEqual(code,0)
        
        xmlFile    = join(self.testsPath,'EMX/EMXwrite_badly_formed.emx')
        try:
            validateSchema(xmlFile)
        except ValidateError, v:
            print "EXCEPTION TESTING STARTS HERE: an Error message should appear. It is OK disregard it."
            print "Validate Error"
            print v.getCode(), v.getMessage()
            print "EXCEPTION TESTING ENDS HERE."
        except Exception, e:
            print "ERROR; we should have never arrive here:", e
            self.assertEqual(1,0)

from  XmippPythonTestResult import XmippPythonTestResult

                                        
if __name__ == '__main__':
    #unittest.main()   
    argc = len(sys.argv)      
    if  argc > 1:  
        xmlFile = sys.argv[1]
    else: 
        xmlFile = '/dev/null'

    suite = unittest.TestLoader().loadTestsFromTestCase(TestEMX)
    result = XmippPythonTestResult()
    result.openXmlReport("TestXmippPythonInterface", xmlFile)    
    suite(result)
    result.closeXmlReport()
    
    if result.testFailed != 0:
       result = unittest.TextTestRunner(verbosity=2).run(suite)    
