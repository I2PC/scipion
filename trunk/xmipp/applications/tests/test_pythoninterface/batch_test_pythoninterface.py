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
        
    def test_Image_add(self):
        path1 = "1@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img1 = Image(path1)
        path2 = "2@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img2 = Image(path2)
        img1.add(img2)
        path2 = os.path.join(self.testsPath, "test_image_generic", "add.spi")
        img2.read(path2)
        self.assertEqual(img1, img2)
        img1.add(img2)        
        self.assertNotEqual(img1, img2)   
             
    def test_Image_minus(self):
        pathSum = os.path.join(self.testsPath, "test_image_generic", "add.spi")
        imgAdd =Image(pathSum)
        path1 = "1@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img1 = Image(path1)
        path2 = "2@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img2 = Image(path2)
        
        imgAdd.minus(img2)
        self.assertEqual(img1, imgAdd)   
             
    def test_FileName_compose(self):
         fn1 = FileName("kk000001.xmp")
         fn2 = FileName("")
         fn2.compose("kk", 1, "xmp")
         self.assertEqual(str(fn1), str(fn2))
         self.assertNotEqual(str(fn1) + 'kk', str(fn2))

    def test_FileName_isInStack(self):
         fn1 = FileName("1@.xmp")
         fn2 = FileName("1.xmp")
         self.assertTrue (fn1.isInStack())
         self.assertFalse(fn2.isInStack())
         
    def test_FileName_isMetaData(self):
         imgPath = os.path.join(self.testsPath, "test_image", "smallStack.stk")
         fn1 = FileName(imgPath)
         self.assertFalse(fn1.isMetaData())
         imgPath = os.path.join(self.testsPath, "test_pythoninterface", "test.xmd")
         fn2 = FileName(imgPath)
         self.assertTrue (fn2.isMetaData())

    def test_Metadata_iter(self):
         mdPath = os.path.join(self.testsPath, "test_pythoninterface", "test.xmd")
         mD = MetaData(mdPath)
         i = 1
         for id in mD:
             img = mD.getValue(MDL_IMAGE, id)
             expImg = '00000%d@Images/proj_ctf_1.stk' % i
             self.assertEqual(img, expImg)
             i += 1
                 
    def test_Metadata_getValue(self):
        '''MetaData_GetValues'''
        mdPath = os.path.join(self.testsPath, "test_pythoninterface", "test.xmd")
        mD = MetaData(mdPath)
        
        img = mD.getValue(MDL_IMAGE, 1L)
        self.assertEqual(img, '000001@Images/proj_ctf_1.stk')
        defocus = mD.getValue(MDL_CTF_DEFOCUSU, 2L)
        self.assertAlmostEqual(defocus, 200.0)
        count = mD.getValue(MDL_COUNT, 3L)
        self.assertEqual(count, 30)
        list = mD.getValue(MDL_ANGLE_COMPARISON, 1L)
        self.assertEqual(list, [1.0, 2.0, 3.0])
        ref = mD.getValue(MDL_REF3D, 2L)
        self.assertEqual(ref, 2)
        
    def test_Metadata_setValue(self):
        '''MetaData_setValues'''
        '''This test should produce the following metadata, which is the same of 'test.xmd'
        data_
        loop_
         _image
         _CTFModel
         _CTF_Defocus_U
         _count
         _angleComparison
         _ref3d
         000001@Images/proj_ctf_1.stk  CTFs/10.ctfparam  -100.000000         10 [     1.000000     2.000000     3.000000 ]         -1
         000002@Images/proj_ctf_1.stk  CTFs/10.ctfparam   200.000000         20 [     1.000000     4.000000     9.000000 ]          2
         000003@Images/proj_ctf_1.stk  CTFs/10.ctfparam  -300.000000         30 [     1.000000     8.000000    27.000000 ]         -3
        ''' 
        mdPath = os.path.join(self.testsPath, "test_pythoninterface", "test.xmd")
        mdRef = MetaData(mdPath)
        md = MetaData()
        ii = -1
        listOrig = [1.0, 2.0, 3.0]
        for i in range(1, 4):
            id = md.addObject() 
            img = '00000%i@Images/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTFMODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_CTF_DEFOCUSU, (i * ii * 100.0), id)
            md.setValue(MDL_COUNT, (i * 10L), id)
            list = [x**i for x in listOrig]
            md.setValue(MDL_ANGLE_COMPARISON, list, id)
            md.setValue(MDL_REF3D, (i * ii), id)
            ii *= -1
            
            
        #print mdRef, md
        self.assertEqual(mdRef, md)
        self.assertRaises(XmippError, md.setValue, MDL_COUNT, 5.5, 1L)
       
   
         
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
           
class XmippPythonTestResult(TestResult):
    xml = None
    testFailed = 0
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
             sys.stderr.write("%s[  PASSED  ] %s %d tests\n" % (bcolors.OKGREEN, bcolors.ENDC, self.numberTests - self.testFailed))
             
    
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

