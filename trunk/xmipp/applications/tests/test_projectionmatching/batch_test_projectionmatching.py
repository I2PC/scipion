#!/usr/bin/env python
import unittest, os, sys
"""
@summary: This pyUnit test module defines the unit tests for the Xmipp Python Interface
"""
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
sys.path.append(scriptdir)
from xmipp import *
#from test.test_array import NumberTest
#from json.tests.test_fail import TestFail

class TestProjMatching(unittest.TestCase):
    testsPath = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/applications/tests'
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        pass
        
    def test_Image_compare(self):
        imgPath = os.path.join(self.testsPath, "test_image", "singleImage.spi")
        img1 = Image()
        img1.read(imgPath)
        img2 = Image(imgPath)
        # Test that image is equal to itself
        self.assertEqual(img1, img2)
        # Test different images
        imgPath = "1@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img2.read(imgPath)
        self.assertNotEqual(img1, img2)
           
from  XmippPythonTestResult import XmippPythonTestResult

                                        
if __name__ == '__main__':
    #unittest.main()   
    argc = len(sys.argv)      
    if  argc > 1:  
        xmlFile = sys.argv[1]
    else: 
        xmlFile = '/dev/null'

    suite = unittest.TestLoader().loadTestsFromTestCase(TestProjMatching)
    result = XmippPythonTestResult()
    result.openXmlReport("TestProjMatching", xmlFile)    
    suite(result)
    result.closeXmlReport()
    
    if result.testFailed != 0:
       result = unittest.TextTestRunner(verbosity=2).run(suite)    

