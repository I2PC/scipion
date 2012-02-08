#!/usr/bin/env xmipp_python
import unittest, os, sys
"""
@summary: This pyUnit test module defines the unit tests for the Xmipp Python Interface
"""
from unittest import TestResult, _TextTestResult
try:
   from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
   from unittest import _WritelnDecorator # Python <2.6
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
sys.path.append(scriptdir)
from xmipp import *
#from test.test_array import NumberTest
#from json.tests.test_fail import TestFail

import sys

def binaryFileComparison(nameo, namet):
    ## open files
    try: file1 = open(nameo, "rb")
    except: print "cannot open file:", nameo ; exit()
    try: file2 = open(namet, "rb")
    except: print "cannot open file:", namet ; file1.close() ; exit()
    ## read 1b from each file until one file reaches eof or bytes don't match
    x = 1; y = 1
    while (x == 1) and (y == 1):
        a = file1.read(1)
        b = file2.read(1)
        x = len(a)
        y = len(b)
        if (a != b) or (x != y):
            return False
    file2.close()
    file1.close()
    return True


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
 
    def test_Image_initConstant(self):
        imgPath = os.path.join(self.testsPath, "test_pythoninterface", "tinyImage.spi")
        img = Image(imgPath)
        img.initConstant(1.) 
        img.write("/tmp/kk.spi") 
        for i in range(0, 3):
            for j in range (0, 3):
                p = img.getPixel(i, j)
                self.assertAlmostEquals(p, 1.0)
    
    def test_Image_initRandom(self):
        imgPath = os.path.join(self.testsPath, "test_pythoninterface", "tinyImage.spi")
        img = Image(imgPath)
        img.resize(1024, 1024) 
        img.initRandom(0., 1., XMIPP_RND_GAUSSIAN)
        mean, dev, min, max = img.computeStats()
        self.assertAlmostEqual(mean, 0., 2)
        self.assertAlmostEqual(dev, 1., 2)
                      
    def test_Image_add(self):
        stackPath = os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img1 = Image("1@" + stackPath)
        img2 = Image("2@" + stackPath)
        sum = img1 + img2
        sumRef = Image(os.path.join(self.testsPath, "test_image_generic", "sum.spi"))
        self.assertEqual(sum, sumRef)
        img1 += img2        
        self.assertEqual(sum, img1)
        img1 += img2
        self.assertNotEqual(sum, img1)   
             
    def test_Image_computeStatistics(self):
        stackPath = os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img1 = Image("1@" + stackPath)
        mean, dev, min, max = img1.computeStats()
        self.assertAlmostEqual(mean, -0.000360, 5)
        self.assertAlmostEqual(dev, 0.105687, 5)
        self.assertAlmostEqual(min, -0.415921, 5)
        self.assertAlmostEqual(max, 0.637052, 5)
             
    def test_Image_minus(self):
        pathSum = os.path.join(self.testsPath, "test_image_generic", "sum.spi")
        imgAdd = Image(pathSum)
        path1 = "1@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img1 = Image(path1)
        path2 = "2@" + os.path.join(self.testsPath, "test_image", "smallStack.stk")
        img2 = Image(path2)
        
        imgAdd -= img2
        self.assertEqual(img1, imgAdd)   
        
    def test_Image_read(self):
        imgPath = os.path.join(self.testsPath, "test_pythoninterface", "tinyImage.spi")
        img = Image(imgPath)        
        count = 0.
        for i in range(0, 3):
            for j in range (0, 3):
                p = img.getPixel(i, j)
                self.assertAlmostEquals(p, count)
                count += 1.
                
    def test_Image_read_header(self):
        
        imgPath = os.path.join(self.testsPath, "test_pythoninterface", "tinyImage.spi")
        img = Image()
        img.read(imgPath, HEADER)        
        
        (x, y, z, n) = img.getDimensions()
        
        self.assertEqual(x, 3)   
        self.assertEqual(y, 3)   
        self.assertEqual(z, 1)   
        self.assertEqual(n, 1)   

        
                
    def test_Image_readApplyGeo(self):
        imgPath = os.path.join(self.testsPath, "test_pythoninterface", "tinyRotated.spi")
        img = Image(imgPath)  
        imgPath = os.path.join(self.testsPath, "test_pythoninterface", "tinyImage.spi")
        md = MetaData()
        id = md.addObject()
        md.setValue(MDL_IMAGE, imgPath, id)
        md.setValue(MDL_ANGLEPSI, 90., id)
        img2 = Image()
        img2.readApplyGeo(md, id)
        self.assertEquals(img, img2)
        
             
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
        
    def test_Metadata_importObjects(self):
        '''import metadata subset'''
        mdPath = os.path.join(self.testsPath, "test_pythoninterface", "test.xmd")
        mD = MetaData(mdPath)
        mDout = MetaData()
        mDout.importObjects(mD, MDValueEQ(MDL_REF3D, -1))
        mdPath = os.path.join(self.testsPath, "test_pythoninterface", "importObject.xmd")
        mD = MetaData(mdPath)
        self.assertEqual(mD, mDout)
        
        
    def test_Metadata_read(self):
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
            list = [x ** i for x in listOrig]
            md.setValue(MDL_ANGLE_COMPARISON, list, id)
            md.setValue(MDL_REF3D, (i * ii), id)
            ii *= -1
            
        tmpFileName = '/tmp/test_pythoninterface_read_tmp.xmd'
        md.write(tmpFileName)
        
        md2 = MetaData()
        md2.read(tmpFileName)
        os.remove(tmpFileName)
        
        self.assertEqual(md, md2)

    def test_Metadata_read_with_labels(self):
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
            list = [x ** i for x in listOrig]
            md.setValue(MDL_ANGLE_COMPARISON, list, id)
            md.setValue(MDL_REF3D, (i * ii), id)
            ii *= -1
            
        tmpFileName = '/tmp/test_pythoninterface_read_tmp.xmd'
        md.write(tmpFileName)
        
        mdList = [MDL_IMAGE, MDL_COUNT]
        md.read(tmpFileName, mdList)
        os.remove(tmpFileName)
        
        md2 = MetaData()
        for i in range(1, 4):
            id = md2.addObject() 
            img = '00000%i@Images/proj_ctf_1.stk' % i
            md2.setValue(MDL_IMAGE, img, id)
            md2.setValue(MDL_COUNT, (i * 10L), id)
        self.assertEqual(md, md2)

        
        
    def test_Metadata_setColumnFormat(self):
        '''MetaData_setValues'''
        '''This test should produce the following metadata, which is the same of 'test_row.xmd'
        data_
         _image 000001@Images/proj_ctf_1.stk
         _CTFModel CTFs/10.ctfparam
        ''' 
        mdPath = os.path.join(self.testsPath, "test_pythoninterface", "test_row.xmd")
        mdRef = MetaData(mdPath)
        md = MetaData()
        ii = -1
        listOrig = [1.0, 2.0, 3.0]
        id = md.addObject() 
        img = '000001@Images/proj_ctf_1.stk'
        md.setValue(MDL_IMAGE, img, id)
        md.setValue(MDL_CTFMODEL, 'CTFs/10.ctfparam', id)
        
        md.setColumnFormat(False)      
        rowFileName = '/tmp/test_row_tmp.xmd'
        md.write(rowFileName)
             
        equalBool = binaryFileComparison(mdPath, rowFileName)
        self.assertEqual(equalBool, True)
        os.remove(rowFileName)

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
            list = [x ** i for x in listOrig]
            md.setValue(MDL_ANGLE_COMPARISON, list, id)
            md.setValue(MDL_REF3D, (i * ii), id)
            ii *= -1
            
            
        #print mdRef, md
        self.assertEqual(mdRef, md)
        self.assertRaises(XmippError, md.setValue, MDL_COUNT, 5.5, 1L)
       
from  XmippPythonTestResult import XmippPythonTestResult

                                        
if __name__ == '__main__':
    #unittest.main()   
    argc = len(sys.argv)      
    if  argc > 1:  
        xmlFile = sys.argv[1]
    else: 
        xmlFile = '/dev/null'

    suite = unittest.TestLoader().loadTestsFromTestCase(TestXmippPythonInterface)
    result = XmippPythonTestResult()
    result.openXmlReport("TestXmippPythonInterface", xmlFile)    
    suite(result)
    result.closeXmlReport()
    
    if result.testFailed != 0:
       result = unittest.TextTestRunner(verbosity=2).run(suite)    

