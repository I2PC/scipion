#!/usr/bin/env python
import unittest, os, sys
"""
@summary: This pyUnit test module defines the unit tests for the Xmipp Python Interface
"""
from xmipp import *
from unittest import TestResult, _TextTestResult, _WritelnDecorator
#from test.test_array import NumberTest
#from json.tests.test_fail import TestFail

class TestXmippPythonInterface(unittest.TestCase):
    testsPath = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/applications/tests'
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        pass
    
    def test_compareImage(self):
        imgPath = os.path.join(self.testsPath, "test_image", "singleImage.spi")
        img1 = Image()
        img1.read(imgPath)
        img2 = Image()
        img2.read(imgPath)
        # Test that image is equal to itself
        self.assertEqual(img1, img2)
        # Test different images
        imgPath = "1@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img2.read(imgPath)
        self.assertNotEqual(img1, img2)
        
    def test_compose(self):
         fn1 = FileName("kk000001.xmp")
         fn2 = FileName("")
         fn2.compose("kk",1,"xmp")
         self.assertEqual(str(fn1),str(fn2))
         self.assertNotEqual(str(fn1)+'kk',str(fn2))
    def test_isInStack(self):
         fn1 = FileName("1@.xmp")
         fn2 = FileName("1.xmp")
         self.assertTrue (fn1.isInStack())
         self.assertFalse(fn2.isInStack())
         
         
         
         
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
           
class XmippPythonTestResult(TestResult):
    xml = None
    testFailed  = 0
    numberTests = 0
    def openXmlReport(self, classname, filename):
        self.xml = open(filename, 'w')
        self.xml.write('<testsuite name="%s">\n' % classname)
        
    def closeXmlReport(self):
         self.xml.write('</testsuite>\n')
         self.xml.close()
         sys.stderr.write("%s[==========]%s run %d tests\n" % (bcolors.OKGREEN, bcolors.ENDC, self.numberTests))
         if (self.testFailed):
             sys.stderr.write("%s[  FAILED  ] %s %d tests\n" % (bcolors.FAIL, bcolors.ENDC, self.testFailed))
         else:
             sys.stderr.write("%s[  PASSED  ] %s %d tests\n" % (bcolors.OKGREEN, bcolors.ENDC, self.numberTests-self.testFailed))
             
    
    def getTestNames(self, test):
        parts = str(test).split()
        name = parts[0];
        parts = parts[1].split('.')
        classname = parts[1].replace(")", "")
        self.numberTests += 1 
        return (name, classname)
    
    def addSuccess(self, test):
        name, classname = self.getTestNames(test)
        self.xml.write('   <testcase name="%s" classname="%s"/>\n' % (name, classname))
        sys.stderr.write("%s[ RUN      ]%s %s.%s\n" % (bcolors.OKGREEN, bcolors.ENDC, classname, name))
        sys.stderr.write("%s[      OK  ]%s %s.%s\n" % (bcolors.OKGREEN, bcolors.ENDC, classname, name))
    
    def reportError(self, test, err):
        name, classname = self.getTestNames(test)
        self.xml.write('   <testcase name="%s" classname="%s">\n' % (name, classname))
        self.xml.write('      <failure message=" "/>\n')
        self.xml.write('   </testcase>\n')
        sys.stderr.write("%s[   FAILED ]%s %s.%s\n" % (bcolors.FAIL, bcolors.ENDC, classname, name))
        self.testFailed += 1
                
    def addError(self, test, err):
        self.reportError(test, err)
        
    def addFailure(self, test, err):
        self.reportError(test, err)
                                        
if __name__ == '__main__':
    #unittest.main()      
    xmlFile = sys.argv[1]
    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippPythonInterface)
    result = XmippPythonTestResult()
    result.openXmlReport("TestXmippPythonInterface", xmlFile)    
    suite(result)
    result.closeXmlReport()
    
    if result.testFailed != 0:
       result = unittest.TextTestRunner(verbosity=2).run(suite)    

