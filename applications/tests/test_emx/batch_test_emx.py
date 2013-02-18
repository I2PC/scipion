#!/usr/bin/env xmipp_python
import unittest, os, sys
"""
@summary: This pyUnit test module defines the unit tests for the Xmipp Python Interface
"""
###
### remove xmipp dependency
###

from unittest import TestResult, _TextTestResult
from protlib_filesystem import getXmippPath
from test.test_support import unlink
try:
   from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
   from unittest import _WritelnDecorator # Python <2.6
#scriptdir = getXmippPath('lib')
#sys.path.append(scriptdir) # add default search path
#scriptdir = getXmippPath()))[0] + '/protocols'
#sys.path.append(scriptdir)
from xmipp import *
#from test.test_array import NumberTest
#from json.tests.test_fail import TestFail

import sys
import os
from os.path import join
from protlib_emx import *


class TestEMX(unittest.TestCase):
    testsPath = getXmippPath("resources", "test")
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        os.chdir(self.testsPath)
        
    def test_EMX_write(self):
        emxData = EmxData()
        outFileName = "/tmp/emxTrivialFile.xml"

        emxMicrograph = EmxMicrograph('mic0017.mrc',1)
        emxMicrograph.dictMicrograph['acceleratingVoltage']=100
        emxMicrograph.dictMicrograph['pixelSpacingX']=5.
        emxMicrograph.dictMicrograph['pixelSpacingY']=6.
        emxData.appendObject(emxMicrograph)
        
        emxMicrograph = EmxMicrograph('mic0018.mrc',100)
        emxMicrograph.dictMicrograph['acceleratingVoltage']=101
        emxMicrograph.dictMicrograph['pixelSpacingX']=5.5
        emxMicrograph.dictMicrograph['pixelSpacingY']=6.6
        emxData.appendObject(emxMicrograph)
        
        emxParticle = EmxParticle('par0017.mrc',1)
        emxParticle.dictParticle['boxSizeX']=11
        emxParticle.dictParticle['boxSizeY']=22
        emxParticle.dictParticle['pixelSpacingX']=6.6
        emxParticle.dictParticle['pixelSpacingY']=66.6
        emxParticle.dictParticle['micrographFileName']='mic0001.mrc'
        emxParticle.dictParticle['micrographIndex']=3
        emxParticle.dictParticle['transformationMatrix11']=11.
        emxParticle.dictParticle['transformationMatrix12']=12.
        emxParticle.dictParticle['transformationMatrix13']=13.
        emxParticle.dictParticle['transformationMatrix14']=14.
        emxParticle.dictParticle['transformationMatrix21']=21.
        emxData.appendObject(emxParticle)
        
        emxParticle = EmxParticle('par0018.mrc',1)
        emxData.appendObject(emxParticle)
        
        #_emxData.setEmxClass('emxMicrograph')
        #_emxData.setEmxClass('emxParticle')
        
        emxData.write(outFileName)
        from filecmp import cmp
        print 'ssss',join(self.testsPath,'EMX/emxTrivialFile.xml')
        self.assertTrue(cmp(outFileName,join(self.testsPath,'EMX/emxTrivialFile.xml'),False))
#        self.assertAlmostEqual(a.all(), b.all(), 2)
        #unlink(outFileName)
           
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
