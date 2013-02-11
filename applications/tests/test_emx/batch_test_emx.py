#!/usr/bin/env xmipp_python
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
from test.test_support import unlink
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
 
import sys
import os
from os.path import join
from emx_data_model import *
from emx_writer import EmxXmlWriter
from emx_reader import EmxXmlReader

class TestEMX(unittest.TestCase):
    testsPath = getXmippPath("resources", "test")
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        os.chdir(self.testsPath)
        
    def test_findObjectType(self):
        inFileName = join(self.testsPath,'EMX/emxTrivialFileRead.xml')
        io = EmxXmlReader()
        emxData = EmxData()
        io.read(inFileName,emxData)
        self.assertTrue(emxData.findObjectType('mic0017.mrc'),MICROGRAPH)
        self.assertTrue(emxData.findObjectType('mic0018.mrc'),MICROGRAPH)
        self.assertTrue(emxData.findObjectType('par0017.mrc'),PARTICLE)
        self.assertFalse(emxData.findObjectType('par0017w.mrc'),PARTICLE)
        

    def test_findObject(self):
        emxData     = EmxData()

        emxMicrograph1 = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph1.set('acceleratingVoltage',100)
        emxMicrograph1.set('pixelSpacing__X',5.)
        emxMicrograph1.set('pixelSpacing__Y',6.)
        emxData.appendMicrograph(emxMicrograph1)
        
        emxMicrograph2 = EmxMicrograph('mic0018.mrc',100)
        emxMicrograph2.set('acceleratingVoltage',101)
        emxMicrograph2.set('pixelSpacing__X',5.5)
        emxMicrograph2.set('pixelSpacing__Y',6.6)
        emxData.appendMicrograph(emxMicrograph2)
        
        valuesCmp={FILENAME:'mic0017.mrc',INDEX:1}
        self.assertTrue(emxData.findObject(emxData.listMicrographs,**valuesCmp) == emxMicrograph1)
        valuesCmp={FILENAME:'mic0017.mrc',INDEX:1}
        self.assertFalse(emxData.findObject(emxData.listMicrographs,**valuesCmp) == emxMicrograph2)
        valuesCmp={FILENAME:'mic0018.mrc',INDEX:100}
        self.assertTrue(emxData.findObject(emxData.listMicrographs,**valuesCmp) == emxMicrograph2)

    def test_EMX_Eq(self):
        #first mic
        emxMicrograph = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph.set('acceleratingVoltage',100.)
        emxMicrograph.set('pixelSpacing__X',5.)
        emxMicrograph.set('pixelSpacing__Y',6.)
        #second mic
        emxMicrograph2 = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph2.set('acceleratingVoltage',100.)
        emxMicrograph2.set('pixelSpacing__X',5.)
        emxMicrograph2.set('pixelSpacing__Y',6.)
        self.assertTrue(emxMicrograph == emxMicrograph2)
        #third mic
        emxMicrograph3 = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph3.set('acceleratingVoltage',100.)
        emxMicrograph3.set('pixelSpacing__X',5.)
        emxMicrograph3.set('pixelSpacing__Y',66.)
        #compare
        self.assertFalse(emxMicrograph3 == emxMicrograph2)

        emxParticle = EmxParticle('par0017.mrc',1)
        emxParticle.set('boxSize__X',11)
        emxParticle.set('boxSize__Y',22)
        emxParticle.set('pixelSpacing__X',6.6)
        emxParticle.set('pixelSpacing__Y',66.6)
        #
        emxParticle2 = EmxParticle('par0017.mrc',1)
        emxParticle2.set('boxSize__X',11)
        emxParticle2.set('boxSize__Y',22)
        emxParticle2.set('pixelSpacing__X',6.6)
        emxParticle2.set('pixelSpacing__Y',66.6)
        self.assertTrue(emxParticle == emxParticle2)
        #
        emxParticle3 = EmxParticle('par0017.mrc',1)
        emxParticle3.set('boxSize__X',11)
        emxParticle3.set('boxSize__Y',22)
        emxParticle3.set('pixelSpacing__X',6.6)
        emxParticle3.set('pixelSpacing__Y',666)
        self.assertFalse(emxParticle3 == emxParticle2)

        emxData  = EmxData()
        emxData2 = EmxData()
        emxData.appendMicrograph(emxMicrograph)
        emxData.appendMicrograph(emxMicrograph2)
        emxData.appendParticle(emxParticle)
        emxData.appendParticle(emxParticle2)
        emxData2.appendMicrograph(emxMicrograph)
        emxData2.appendMicrograph(emxMicrograph2)
        emxData2.appendParticle(emxParticle)
        emxData2.appendParticle(emxParticle2)
        self.assertTrue(emxData == emxData2)
        

    def test_EMX_write(self):
        emxData     = EmxData()
        outFileName = "/tmp/emxTrivialFileWrite.xml"

        emxMicrograph = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph.set('acceleratingVoltage',100)
        emxMicrograph.set('pixelSpacing__X',5.)
        emxMicrograph.set('pixelSpacing__Y',6.)
        emxData.appendMicrograph(emxMicrograph)
        
        emxMicrograph = EmxMicrograph('mic0018.mrc',100)
        emxMicrograph.set('acceleratingVoltage',101)
        emxMicrograph.set('pixelSpacing__X',5.5)
        emxMicrograph.set('pixelSpacing__Y',6.6)
        emxData.appendMicrograph(emxMicrograph)
        
        emxParticle = EmxParticle('par0017.mrc',1,emxMicrograph)
        emxParticle.set('boxSize__X',11)
        emxParticle.set('boxSize__Y',22)
        emxParticle.set('pixelSpacing__X',6.6)
        emxParticle.set('pixelSpacing__Y',66.6)
        emxParticle.set('transformationMatrix__t11',11.)
        emxParticle.set('transformationMatrix__t12',12.)
        emxParticle.set('transformationMatrix__t13',13.)
        emxParticle.set('transformationMatrix__t14',14.)
        emxParticle.set('transformationMatrix__t21',21.)
        emxData.appendParticle(emxParticle)

        emxParticle2 = EmxParticle('par0018.mrc',1)
        emxData.appendParticle(emxParticle2)

        io = EmxXmlWriter()
        io.write(outFileName,emxData)

        from filecmp import cmp
        self.assertTrue(cmp(outFileName,join(self.testsPath,'EMX/emxTrivialFile.xml'),False))
        unlink(outFileName)

    def test_EMX_read(self):
        emxData = EmxData()
        inFileName = join(self.testsPath,'EMX/emxTrivialFileRead.xml')
                          
        emxMicrograph = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph.set('acceleratingVoltage',100)
        emxMicrograph.set('cs',2.2)
        emxMicrograph.set('pixelSpacing__X',5.)
        emxMicrograph.set('pixelSpacing__Y',6.)
        emxData.appendMicrograph(emxMicrograph)
        
        emxMicrograph = EmxMicrograph('mic0018.mrc',100)
        emxMicrograph.set('acceleratingVoltage',101)
        emxMicrograph.set('cs',2.2)
        emxMicrograph.set('pixelSpacing__X',5.5)
        emxMicrograph.set('pixelSpacing__Y',6.6)
        emxData.appendMicrograph(emxMicrograph)
        
        emxParticle = EmxParticle('par0017.mrc',1,emxMicrograph)
        emxParticle.set('boxSize__X',11)
        emxParticle.set('boxSize__Y',22)
        emxParticle.set('pixelSpacing__X',6.6)
        emxParticle.set('pixelSpacing__Y',66.6)
        emxParticle.set('transformationMatrix__t11',11.)
        emxParticle.set('transformationMatrix__t12',12.)
        emxParticle.set('transformationMatrix__t13',13.)
        emxParticle.set('transformationMatrix__t14',14.)
        emxParticle.set('transformationMatrix__t21',21.)
        emxData.appendParticle(emxParticle)

        emxParticle = EmxParticle('par0018.mrc',1)
        emxData.appendParticle(emxParticle)

        io = EmxXmlReader()
        emxData2 = EmxData()
        io.read(inFileName,emxData2)
#        print emxData2,emxData
        self.assertTrue(emxData2 == emxData)
           
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
