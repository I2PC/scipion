#!/usr/bin/env xmipp_python
import unittest, os, sys
"""
@summary: This pyUnit test module defines the unit tests for the Xmipp Python Interface
"""
from unittest import TestResult, _TextTestResult
from protlib_filesystem import getXmippPath
from tempfile import NamedTemporaryFile
from time import time
from os.path import  exists
from random import randint
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

def testFile(filename):
    return join("pythoninterface", filename)

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
    testsPath = getXmippPath("resources", "test")
    # Temporary filenames used
    
    def getTmpName(self, suffix='.xmd'):
        """ Get temporary filenames in /tmp' that will be removed 
        after tests execution. 
        """
        ntf = NamedTemporaryFile(dir="/tmp/", suffix=suffix)
        self.tempNames.append(ntf)        
        return ntf.name
    
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        os.chdir(self.testsPath)
        self.tempNames = []
        
    def tearDown(self):
        for ntf in self.tempNames:
            ntf.close()
        
    def test_xmipp_Euler_angles2matrix(self):
        from numpy  import array

        a = array([[ 0.70710678, 0.70710678, -0.        ],
           [-0.70710678, 0.70710678, 0.        ],
           [ 0., 0., 1.        ]])
        rot = 45.
        tilt = 0.
        psi = 0.
        b = Euler_angles2matrix(rot, tilt, psi)
        self.assertAlmostEqual(a.all(), b.all(), 2)

    def test_xmipp_Euler_matrix2angles(self):
        from numpy  import array

        a = array([[ 0.70710678, 0.70710678, -0.        ],
           [-0.70710678, 0.70710678, 0.        ],
           [ 0., 0., 1.        ]])
        rot = 45.
        tilt = 0.
        psi = 0.
        rot1,tilt1,psi1 = Euler_matrix2angles(a)
        if psi1 > 1:
            rot1 = psi1
            psi1=0.
        self.assertAlmostEqual(rot, rot1, 4)
        self.assertAlmostEqual(tilt, tilt1, 4)
        self.assertAlmostEqual(psi, psi1, 4)

    def test_xmipp_activateMathExtensions(self):
        '''activateMathExtensions'''
        md1 = MetaData()
        activateMathExtensions()
        id = md1.addObject()
        md1.setValue(MDL_ANGLE_ROT, 4., id)
        md1.setValue(MDL_ANGLE_TILT, 3., id)
        id = md1.addObject()
        md1.setValue(MDL_ANGLE_ROT, 9., id)
        md1.setValue(MDL_ANGLE_TILT, 5., id)
        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 2., id)
        md2.setValue(MDL_ANGLE_TILT, 3., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 3., id)
        md2.setValue(MDL_ANGLE_TILT, 5., id)
        angleRotLabel = label2Str(MDL_ANGLE_ROT)
        operateString = angleRotLabel + "= sqrt(" + angleRotLabel+")"
        md1.operate(operateString)
        self.assertEqual(md1, md2)

    def test_xmipp_errorBetween2CTFs(self):
        '''activateMathExtensions'''
        md1 = MetaData()
        id = md1.addObject()
#         md1.setValue(MDL_CTF_SAMPLING_RATE, 1., id)
#         md1.setValue(MDL_CTF_VOLTAGE, 300., id);
#         md1.setValue(MDL_CTF_DEFOCUSU, 5000., id);
#         md1.setValue(MDL_CTF_DEFOCUSV, 7500., id);
#         md1.setValue(MDL_CTF_DEFOCUS_ANGLE, -45., id);
#         md1.setValue(MDL_CTF_CS, 2., id);
#         md1.setValue(MDL_CTF_Q0, 0.1, id);
# 
#         md2 = MetaData()
#         id = md2.addObject()
#         md2.setValue(MDL_CTF_SAMPLING_RATE, 1., id)
#         md2.setValue(MDL_CTF_VOLTAGE, 300., id);
#         md2.setValue(MDL_CTF_DEFOCUSU, 10000., id);
#         md2.setValue(MDL_CTF_DEFOCUSV, 10000., id);
#         md2.setValue(MDL_CTF_DEFOCUS_ANGLE, 45., id);
#         md2.setValue(MDL_CTF_CS, 2., id);
#         md2.setValue(MDL_CTF_Q0, 0.1, id);

        md1.setValue(MDL_CTF_SAMPLING_RATE, 1., id)
        md1.setValue(MDL_CTF_VOLTAGE, 200., id);
        md1.setValue(MDL_CTF_DEFOCUSU, 18306.250000, id);
        md1.setValue(MDL_CTF_DEFOCUSV, 16786.470000, id);
        md1.setValue(MDL_CTF_DEFOCUS_ANGLE, 30.100000, id);
        md1.setValue(MDL_CTF_CS, 2., id);
        md1.setValue(MDL_CTF_Q0, 0.07, id);
 
        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_CTF_SAMPLING_RATE, 1., id)
        md2.setValue(MDL_CTF_VOLTAGE, 200., id);
        md2.setValue(MDL_CTF_DEFOCUSU, 17932.700000, id);
        md2.setValue(MDL_CTF_DEFOCUSV, 16930.300000, id);
        md2.setValue(MDL_CTF_DEFOCUS_ANGLE, 45., id);
        md2.setValue(MDL_CTF_CS, 2., id);
        md2.setValue(MDL_CTF_Q0, 0.07, id);
        
        
        error = errorBetween2CTFs(md1,md2, 256, 0.05,0.25)

        self.assertAlmostEqual(error, 10441.1,0)

    def test_xmipp_errorMaxFreqCTFs(self):
        '''activateMathExtensions'''
        md1 = MetaData()
        id = md1.addObject()
        md1.setValue(MDL_CTF_SAMPLING_RATE, 2., id)
        md1.setValue(MDL_CTF_VOLTAGE, 300., id);
        md1.setValue(MDL_CTF_DEFOCUSU, 7500., id);
        md1.setValue(MDL_CTF_DEFOCUSV, 7500., id);
        md1.setValue(MDL_CTF_DEFOCUS_ANGLE, 45., id);
        md1.setValue(MDL_CTF_CS, 2., id);
        md1.setValue(MDL_CTF_Q0, 0.1, id);

        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_CTF_SAMPLING_RATE, 2., id)
        md2.setValue(MDL_CTF_VOLTAGE, 300., id);
        md2.setValue(MDL_CTF_DEFOCUSU, 5000., id);
        md2.setValue(MDL_CTF_DEFOCUSV, 5000., id);
        md2.setValue(MDL_CTF_DEFOCUS_ANGLE, 45., id);
        md2.setValue(MDL_CTF_CS, 2., id);
        md2.setValue(MDL_CTF_Q0, 0.1, id);
        resolution = errorMaxFreqCTFs(md1,md2)

        self.assertAlmostEqual(resolution, 7.0156283,2)

    def test_FileName_compose(self):
         fn1 = FileName("kk000001.xmp")
         fn2 = FileName("")
         fn2.compose("kk", 1, "xmp")
         self.assertEqual(str(fn1), str(fn2))
         self.assertNotEqual(str(fn1) + 'kk', str(fn2))
         
         fn1 = FileName("000001@kk.xmp")
         fn2 = FileName("")
         fn2.compose(1, "kk.xmp")
         self.assertEqual(str(fn1), str(fn2))
         
         fn1 = FileName("jj@kk.xmp")
         fn2 = FileName("")
         fn2.compose("jj", "kk.xmp")
         self.assertEqual(str(fn1), str(fn2))

    def test_FileName_isInStack(self):
         fn1 = FileName("1@.xmp")
         fn2 = FileName("1.xmp")
         self.assertTrue (fn1.isInStack())
         self.assertFalse(fn2.isInStack())
         
    def test_FileName_isMetaData(self):
         imgPath = testFile("smallStack.stk")
         fn1 = FileName(imgPath)
         try: 
             result = fn1.isMetaData()
         except: 
             print "cannot open file:", fn1 ; exit()
         self.assertFalse(result)
         
         imgPath = testFile("test.xmd")
         fn2 = FileName(imgPath)
         self.assertTrue (fn2.isMetaData())

    def test_Image_add(self):
        stackPath = testFile("smallStack.stk")
        img1 = Image("1@" + stackPath)
        img2 = Image("2@" + stackPath)
        sum = img1 + img2
        sumRef = Image(testFile("sum.spi"))
        self.assertEqual(sum, sumRef)
        img1 += img2     
        self.assertEqual(sum, img1)
        img1 += img2
        self.assertNotEqual(sum, img1)   
    
    def test_Image_compare(self):
        imgPath = testFile("singleImage.spi")
        img1 = Image()
        img1.read(imgPath)
        img2 = Image(imgPath)
        # Test that image is equal to itself
        self.assertEqual(img1, img2)
        # Test different images
        imgPath = "1@" + testFile("smallStack.stk")
        img2.read(imgPath)
        self.assertNotEqual(img1, img2)
             
    def test_Image_computeStatistics(self):
        stackPath = testFile("smallStack.stk")
        img1 = Image("1@" + stackPath)
        mean, dev, min, max = img1.computeStats()
        self.assertAlmostEqual(mean, -0.000360, 5)
        self.assertAlmostEqual(dev, 0.105687, 5)
        self.assertAlmostEqual(min, -0.415921, 5)
        self.assertAlmostEqual(max, 0.637052, 5)

    def test_Image_equal(self):
        img = Image()
        img2 = Image()
        img.setDataType(DT_FLOAT)
        img2.setDataType(DT_FLOAT)
        img.resize(3, 3)
        img2.resize(3, 3)
        img.setPixel(1, 1, 1.)
        img2.setPixel(1, 1, 1.01)

        self.assertFalse(img.equal(img2, 0.005))
        self.assertTrue(img.equal(img2, 0.1))
        
    def test_Image_getData(self):
        from numpy import array
        img1 = Image(testFile("singleImage.spi"))
        Z = img1.getData()
        Zref = array([[-0.27623099,-0.13268562, 0.32305956],
                      [ 1.07257104, 0.73135602, 0.49813408],
                      [ 0.90717429, 0.6812411, -0.09380955]])
        self.assertEqual(Z.all(), Zref.all())
 
    def test_Image_initConstant(self):
        imgPath = testFile("tinyImage.spi")
        img = Image(imgPath)
        img.initConstant(1.) 
        for i in range(0, 3):
            for j in range (0, 3):
                p = img.getPixel(i, j)
                self.assertAlmostEquals(p, 1.0)
    
    def test_Image_initRandom(self):
        imgPath = testFile("tinyImage.spi")
        img = Image(imgPath)
        img.resize(1024, 1024) 
        img.initRandom(0., 1., XMIPP_RND_GAUSSIAN)
        mean, dev, min, max = img.computeStats()
        self.assertAlmostEqual(mean, 0., 2)
        self.assertAlmostEqual(dev, 1., 2)
              
    def test_Image_minus(self):
        pathSum = testFile("sum.spi")
        imgAdd = Image(pathSum)
        path1 = "1@" + testFile("smallStack.stk")
        img1 = Image(path1)
        path2 = "2@" + testFile("smallStack.stk")
        img2 = Image(path2)
        
        imgAdd -= img2
        self.assertEqual(img1, imgAdd)   

    def test_Image_project(self):
        vol=Image(testFile('progVol.vol'))
        vol.convert2DataType(DT_DOUBLE)
        proj=vol.projectVolumeDouble(0.,0.,0.)
        self.assertEqual(1,1)
                
    def test_Image_read(self):
        imgPath = testFile("tinyImage.spi")
        img = Image(imgPath)        
        count = 0.
        for i in range(0, 3):
            for j in range (0, 3):
                p = img.getPixel(i, j)
                self.assertAlmostEquals(p, count)
                count += 1.
                
    def test_Image_read_header(self):
        imgPath = testFile("tinyImage.spi")
        img = Image()
        img.read(imgPath, HEADER)        
        
        (x, y, z, n) = img.getDimensions()
        
        self.assertEqual(x, 3)   
        self.assertEqual(y, 3)   
        self.assertEqual(z, 1)   
        self.assertEqual(n, 1)   
                
    def test_Image_readApplyGeo(self):
        imgPath = testFile("tinyRotated.spi")
        img = Image(imgPath)  
        imgPath = testFile("tinyImage.spi")
        md = MetaData()
        id = md.addObject()
        md.setValue(MDL_IMAGE, imgPath, id)
        md.setValue(MDL_ANGLE_PSI, 90., id)
        img2 = Image()
        img2.readApplyGeo(md, id)
        self.assertEquals(img, img2)
        
    def test_Image_setDataType(self):
        img = Image()
        img.setDataType(DT_FLOAT)
        img.resize(3, 3)
        img.setPixel(1, 1, 1.)
        img.computeStats()
        mean, dev, min, max = img.computeStats()
        self.assertAlmostEqual(mean, 0.111111, 5)
        self.assertAlmostEqual(dev, 0.314270, 5)
        self.assertAlmostEqual(min, 0., 5)
        self.assertAlmostEqual(max, 1., 5)   
        
    def test_Image_mirrorY(self):
        img = Image()
        img.setDataType(DT_FLOAT)
        imgY = Image()
        imgY.setDataType(DT_FLOAT)
        dim=3
        img.resize(dim, dim)
        imgY.resize(dim, dim)
        for i in range(0,dim):
            for j in range(0,dim):
                img.setPixel(i, j, 1.*dim*i+j)
        img.mirrorY();
        for i in range(0,dim):
            for j in range(0,dim):
                imgY.setPixel(dim -1 -i, j, 1.*dim*i+j)
        self.assertEquals(img, imgY)
        
    def test_Image_applyTransforMatScipion(self):
            img = Image()
            img2 = Image()
            img.setDataType(DT_FLOAT)
            img2.setDataType(DT_FLOAT)
            dim=3
            img.resize(dim, dim)
            img2.resize(dim, dim)
            for i in range(0,dim):
                for j in range(0,dim):
                    img.setPixel(dim -1 -i, j, 1.*dim*i+j)
            A=[0.,-1.,0.,0.,
               1.,0.,0.,0.,
               0.,0.,1.,1.]
            img2.setPixel(0,0,0)
            img2.setPixel(0,1,3.)
            img2.setPixel(0,2,6.)

            img2.setPixel(1,0,1.)
            img2.setPixel(1,1,4.)
            img2.setPixel(1,2,7.)

            img2.setPixel(2,0,2.)
            img2.setPixel(2,1,5.)
            img2.setPixel(2,2,8.)

            img.applyTransforMatScipion(A)
            #print img.getData(),"\n", img2.getData()
            for i in range(0, dim):
                for j in range (0, dim):
                    p1 = img.getPixel(i, j)
                    p2 = img2.getPixel(i, j)
                    self.assertAlmostEquals(p1, p2)
    def test_Metadata_getValue(self):
        '''MetaData_GetValues'''
        mdPath = testFile("test.xmd")
        mdPath2 = testFile("proj_ctf_1.stk")
        mD = MetaData(mdPath)
        img = mD.getValue(MDL_IMAGE, 1L)
        self.assertEqual(img, '000001@' + mdPath2)
        defocus = mD.getValue(MDL_CTF_DEFOCUSU, 2L)
        self.assertAlmostEqual(defocus, 200)
        count = mD.getValue(MDL_COUNT, 3L)
        self.assertEqual(count, 30)
        list = mD.getValue(MDL_CLASSIFICATION_DATA, 1L)
        self.assertEqual(list, [1.0, 2.0, 3.0])
        ref = mD.getValue(MDL_REF3D, 2L)
        self.assertEqual(ref, 2)
        
    def test_Metadata_importObjects(self):
        '''import metadata subset'''
        mdPath = testFile("test.xmd")
        mD = MetaData(mdPath)
        mDout = MetaData()
        mDout.importObjects(mD, MDValueEQ(MDL_REF3D, -1))
        mdPath = testFile("importObject.xmd")
        mD = MetaData(mdPath)
        self.assertEqual(mD, mDout)
        
    def test_Metadata_iter(self):
         mdPath = testFile("test.xmd")
         mD = MetaData(mdPath)
         i = 1
         for id in mD:
             img = mD.getValue(MDL_IMAGE, id)
             expImg = testFile("test.xmd")
             expImg = ("00000%d@" % i) + testFile("proj_ctf_1.stk")
             #expImg = '00000%d@Images/proj_ctf_1.stk' % i
             self.assertEqual(img, expImg)
             i += 1
        
    def test_Metadata_existsBlockInMetaDataFile(self):
         mdPath = "b2@"+testFile("testBlock.xmd")
         self.assertTrue(existsBlockInMetaDataFile(mdPath))
         mdPath = "nonexisting@"+testFile("testBlock.xmd")
         self.assertFalse(existsBlockInMetaDataFile(mdPath))
         mdPath = "gggg@"+testFile("testBlock.xmd")
         self.assertRaises(Exception,existsBlockInMetaDataFile(mdPath))
         
    def test_Metadata_operate(self):
        md = MetaData()
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 1., id)
        md.setValue(MDL_ANGLE_TILT, 2., id)
        md.setValue(MDL_ANGLE_PSI, 3., id)
        md.setValue(MDL_ANGLE_ROT2, 6., id)
        md.setValue(MDL_ANGLE_TILT2, 5., id)
        md.setValue(MDL_ANGLE_PSI2, 4., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 11., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        md.setValue(MDL_ANGLE_PSI, 13., id)
        md.setValue(MDL_ANGLE_ROT2, 16., id)
        md.setValue(MDL_ANGLE_TILT2, 15., id)
        md.setValue(MDL_ANGLE_PSI2, 14., id)
        
        md2 = MetaData(md)
        anglePsiLabel = label2Str(MDL_ANGLE_PSI)
        angleRotLabel = label2Str(MDL_ANGLE_ROT)
        operateString = angleRotLabel + "=3*" + angleRotLabel
        operateString += "," + anglePsiLabel + "=2*" + anglePsiLabel
        md.operate(operateString)
        for id in md2:
            md2.setValue(MDL_ANGLE_ROT, md2.getValue(MDL_ANGLE_ROT, id) * 3., id);
            md2.setValue(MDL_ANGLE_PSI, md2.getValue(MDL_ANGLE_PSI, id) * 2., id);
        self.assertEqual(md, md2)

    def test_xmipp_addLabelAlias(self):
        tmpFileName = '/tmp/test_pythoninterface_readPlain1.xmd'
        line1="""# XMIPP_STAR_1 * 
# 
data_images

loop_
_rlnVoltage #1
_rlnDefocusU #2
  200.0 10000.0
  201.0 10000.0"""
        #createTestImage
        f = open(tmpFileName,'w')
        f.write(line1)
        f.flush()
        f.close()
        addLabelAlias(MDL_CTF_VOLTAGE,"rlnVoltage")
        addLabelAlias(MDL_CTF_DEFOCUSU,"rlnDefocusU")
        md=MetaData(tmpFileName)
        md2=MetaData()
        id = md2.addObject()
        md2.setValue(MDL_CTF_VOLTAGE, 200.0, id)
        md2.setValue(MDL_CTF_DEFOCUSU, 10000.0, id)
        id = md2.addObject()
        md2.setValue(MDL_CTF_VOLTAGE, 201.0, id)
        md2.setValue(MDL_CTF_DEFOCUSU, 10000.0, id)
        self.assertEqual(md, md2)
        os.remove(tmpFileName)

    def test_xmipp_exportReliontoMetadataFile(self):
        from protlib_import import exportReliontoMetadataFile

        tmpFileNameRe = '/tmp/exportReliontoMetadataFile.star'
        tmpFileNameXm = '/tmp/exportReliontoMetadataFile.xmd'
        line1="""data_images

loop_
_rlnVoltage #1
_rlnDefocusU #2
  200.0 10000.0
  201.0 10000.0"""
        #createR
        f = open(tmpFileNameRe,'w')
        f.write(line1)
        f.flush()
        f.close()
        exportReliontoMetadataFile(tmpFileNameRe,tmpFileNameXm)
        md = MetaData(tmpFileNameXm)
        md2=MetaData()
        id = md2.addObject()
        md2.setValue(MDL_CTF_VOLTAGE, 200.0, id)
        md2.setValue(MDL_CTF_DEFOCUSU, 10000.0, id)
        id = md2.addObject()
        md2.setValue(MDL_CTF_VOLTAGE, 201.0, id)
        md2.setValue(MDL_CTF_DEFOCUSU, 10000.0, id)

        self.assertEqual(md, md2)
        os.remove(tmpFileNameRe)
        os.remove(tmpFileNameXm)
        #print tmpFileNameRe,tmpFileNameXm
        
    def test_Metadata_renameColumn(self):
        '''test_Metadata_renameColumn'''
        md = MetaData()
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 1., id)
        md.setValue(MDL_ANGLE_TILT, 2., id)
        md.setValue(MDL_ANGLE_PSI, 3., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 11., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        md.setValue(MDL_ANGLE_PSI, 13., id)
        
        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT2, 1., id)
        md2.setValue(MDL_ANGLE_TILT2, 2., id)
        md2.setValue(MDL_ANGLE_PSI2, 3., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT2, 11., id)
        md2.setValue(MDL_ANGLE_TILT2, 12., id)
        md2.setValue(MDL_ANGLE_PSI2, 13., id)

        md.renameColumn(MDL_ANGLE_ROT,MDL_ANGLE_ROT2)        
        md.renameColumn(MDL_ANGLE_TILT,MDL_ANGLE_TILT2)        
        md.renameColumn(MDL_ANGLE_PSI,MDL_ANGLE_PSI2)        
        self.assertEqual(md, md2)

        
        md.clear()
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 1., id)
        md.setValue(MDL_ANGLE_TILT, 2., id)
        md.setValue(MDL_ANGLE_PSI, 3., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 11., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        md.setValue(MDL_ANGLE_PSI, 13., id)
        oldLabel=[MDL_ANGLE_ROT,MDL_ANGLE_TILT,MDL_ANGLE_PSI]
        newLabel=[MDL_ANGLE_ROT2,MDL_ANGLE_TILT2,MDL_ANGLE_PSI2]
        md.renameColumn(oldLabel,newLabel)        
        self.assertEqual(md, md2)

    def test_Metadata_sort(self):
        '''test_Metadata_sort'''
        md = MetaData()
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 1., id)
        md.setValue(MDL_ANGLE_TILT, 2., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 11., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 5., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 6., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        
        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 11., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 6., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 5., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 1., id)
        md2.setValue(MDL_ANGLE_TILT, 2., id)

        mdAux=MetaData(md)
        mdAux.sort(MDL_ANGLE_ROT,False)
        self.assertEqual(mdAux, md2)

        mdAux=MetaData(md)
        mdAux.sort(MDL_ANGLE_ROT,False,2)#limit
        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 11., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 6., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        self.assertEqual(mdAux, md2)

        mdAux=MetaData(md)
        mdAux.sort(MDL_ANGLE_ROT,False,2,1)#limit && offset
        md2 = MetaData()
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 6., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        id = md2.addObject()
        md2.setValue(MDL_ANGLE_ROT, 5., id)
        md2.setValue(MDL_ANGLE_TILT, 12., id)
        self.assertEqual(mdAux, md2)


    def test_Metadata_join(self):
         #create metadta
        md = MetaData() 
        md2 = MetaData()
        mdout = MetaData()
        listOrig = [1.0, 2.0, 3.0]
        for i in range(1, 4):
            id = md.addObject() 
            img = '%06d@pythoninterface/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_COUNT, (i * 10L), id)
        for i in range(1, 3):
            id = md2.addObject() 
            img = '%06d@pythoninterface/proj_ctf_1.stk' % i
            md2.setValue(MDL_IMAGE, img, id)
            md2.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md2.setValue(MDL_ANGLE_PSI, 1., id)
        mdout.join (md, md2, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)

        md.clear()
        for i in range(1, 3):
            id = md.addObject() 
            img = '%06d@pythoninterface/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_COUNT, (i * 10L), id)
            md.setValue(MDL_ANGLE_PSI, 1., id)

        self.assertEqual(mdout, md)
                 
    def test_Metadata_intersect(self):
         #create metadta
        md = MetaData()
        md2 = MetaData()
        mdout = MetaData()
        listOrig = [1.0, 2.0, 3.0]
        for i in range(1, 4):
            id = md.addObject() 
            img = '%06d@pythoninterface/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_COUNT, (i * 10L), id)
            md.setValue(MDL_ANGLE_PSI, 1., id)
        for i in range(1, 3):
            id = md2.addObject() 
            img = '%06d@pythoninterface/proj_ctf_1.stk' % i
            md2.setValue(MDL_IMAGE, img, id)
        md.intersection (md2, MDL_IMAGE)

        for i in range(1, 3):
            id = mdout.addObject() 
            img = '%06d@pythoninterface/proj_ctf_1.stk' % i
            mdout.setValue(MDL_IMAGE, img, id)
            mdout.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            mdout.setValue(MDL_COUNT, (i * 10L), id)
            mdout.setValue(MDL_ANGLE_PSI, 1., id)


        self.assertEqual(mdout, md)
                 
    def test_Metadata_read(self):
        '''MetaData_setValues'''
        '''This test should produce the following metadata, which is the same of 'test.xmd'
        data_
        loop_
         _image
         _CTFModel
         _CTF_Defocus_U
         _count
         _classificationData
         _ref3d
         000001@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam  -100.000000         10 ' 1.000000     2.000000     3.000000 '         -1
         000002@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam   200.000000         20 ' 1.000000     4.000000     9.000000 '          2
         000003@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam  -300.000000         30 ' 1.000000     8.000000    27.000000 '         -3
        ''' 
        mdPath = testFile("test.xmd")
        mdRef = MetaData(mdPath)
        md = MetaData()
        ii = -1
        listOrig = [1.0, 2.0, 3.0]
        for i in range(1, 4):
            id = md.addObject() 
            img = '00000%i@pythoninterface/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_CTF_DEFOCUSU, (i * ii * 100.0), id)
            md.setValue(MDL_COUNT, (i * 10L), id)
            list = [x ** i for x in listOrig]
            md.setValue(MDL_CLASSIFICATION_DATA, list, id)
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
         000001@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam  -100.000000         10 ' 1.000000     2.000000     3.000000 '         -1
         000002@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam   200.000000         20 ' 1.000000     4.000000     9.000000 '          2
         000003@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam  -300.000000         30 ' 1.000000     8.000000    27.000000 '         -3
        ''' 
        mdPath = testFile("test.xmd")
        mdRef = MetaData(mdPath)
        md = MetaData()
        ii = -1
        listOrig = [1.0, 2.0, 3.0]
        for i in range(1, 4):
            id = md.addObject() 
            img = '00000%i@pythoninterface/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_CTF_DEFOCUSU, (i * ii * 100.0), id)
            md.setValue(MDL_COUNT, (i * 10L), id)
            list = [x ** i for x in listOrig]
            md.setValue(MDL_CLASSIFICATION_DATA, list, id)
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
            img = '00000%i@pythoninterface/proj_ctf_1.stk' % i
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
        mdPath = testFile("test_row.xmd")
        mdRef = MetaData(mdPath)
        md = MetaData()
        ii = -1
        listOrig = [1.0, 2.0, 3.0]
        id = md.addObject() 
        img = '000001@Images/proj_ctf_1.stk'
        md.setValue(MDL_IMAGE, img, id)
        md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
        
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
         000001@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam  -100.000000         10 [     1.000000     2.000000     3.000000 ]         -1
         000002@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam   200.000000         20 [     1.000000     4.000000     9.000000 ]          2
         000003@pythoninterface/proj_ctf_1.stk  CTFs/10.ctfparam  -300.000000         30 [     1.000000     8.000000    27.000000 ]         -3
        ''' 
        mdPath = testFile("test.xmd")
        mdRef = MetaData(mdPath)
        md = MetaData()
        ii = -1
        listOrig = [1.0, 2.0, 3.0]
        for i in range(1, 4):
            id = md.addObject() 
            img = '00000%i@pythoninterface/proj_ctf_1.stk' % i
            md.setValue(MDL_IMAGE, img, id)
            md.setValue(MDL_CTF_MODEL, 'CTFs/10.ctfparam', id)
            md.setValue(MDL_CTF_DEFOCUSU, (i * ii * 100.0), id)
            md.setValue(MDL_COUNT, (i * 10L), id)
            list = [x ** i for x in listOrig]
            md.setValue(MDL_CLASSIFICATION_DATA, list, id)
            md.setValue(MDL_REF3D, (i * ii), id)
            ii *= -1
            
            
        self.assertEqual(mdRef, md)
        print "THIS IS NOT AN ERROR: we are testing exceptions: an error message should appear regarding count and DOUBLE"
        self.assertRaises(XmippError, md.setValue, MDL_COUNT, 5.5, 1L)
   
        
    def test_Metadata_compareTwoMetadataFiles(self):
        
        try:
            from tempfile import NamedTemporaryFile
            
            fn = self.getTmpName()
            fn2 = self.getTmpName()
            fn3 = self.getTmpName()
        
            auxMd = MetaData()
            auxMd2 = MetaData()
            
            auxMd.setValue(MDL_IMAGE, "image_1.xmp", auxMd.addObject())
            auxMd.setValue(MDL_IMAGE, "image_2.xmp", auxMd.addObject())
        
            auxMd.write("block_000000@" + fn, MD_OVERWRITE)
            auxMd.clear()
            auxMd.setValue(MDL_IMAGE, "image_data_1_1.xmp", auxMd.addObject())
            auxMd.setValue(MDL_IMAGE, "image_data_1_2.xmp", auxMd.addObject())
            auxMd.write("block_000001@" + fn, MD_APPEND)
            auxMd.clear()
            auxMd.setValue(MDL_IMAGE, "image_data_2_1.xmp", auxMd.addObject())
            auxMd.setValue(MDL_IMAGE, "image_data_2_2.xmp", auxMd.addObject())
            auxMd.write("block_000000@" + fn2, MD_OVERWRITE)
            auxMd.clear()
            auxMd.setValue(MDL_IMAGE, "image_data_A_1.xmp", auxMd.addObject())
            auxMd.setValue(MDL_IMAGE, "image_data_A_2.xmp", auxMd.addObject())
            auxMd.write("block_000001@" + fn2, MD_APPEND)
            auxMd.clear()
        
            self.assertFalse(compareTwoMetadataFiles(fn, fn2))
            self.assertTrue(compareTwoMetadataFiles(fn, fn))

            sfnFN = FileName(fn)
            sfn2FN = FileName(fn2)
            self.assertFalse(compareTwoMetadataFiles(sfnFN, sfn2FN))
            self.assertTrue(compareTwoMetadataFiles(sfnFN, sfnFN))
                    
            auxMd.setValue(MDL_IMAGE, "image_1.xmpSPACE", auxMd.addObject())
            auxMd.setValue(MDL_IMAGE, "image_2.xmp", auxMd.addObject())
            auxMd.write("block_000000@" + fn2, MD_OVERWRITE)
            auxMd.clear()
            auxMd.setValue(MDL_IMAGE, "image_data_1_1.xmp", auxMd.addObject())
            auxMd.setValue(MDL_IMAGE, "image_data_1_2.xmp", auxMd.addObject())
            auxMd.write("block_000001@" + fn2, MD_APPEND)
        
            command = "sed 's/SPACE/ /g' " + fn2 + ">" + fn3
            os.system(command)
        
            self.assertTrue(compareTwoMetadataFiles(fn, fn3))
            
        except Exception, e:
            print str(e)
            
    def test_Metadata_agregate(self):
        mdOut = MetaData()
        md    = MetaData()
        mdResult = MetaData()
        for linea in range(1, 101):
            id = md.addObject()
            md.setValue(MDL_IMAGE, '%06d@proj.stk' % linea, id)
            md.setValue(MDL_XCOOR, int(linea%3), id)
        mdOut.aggregate(md, AGGR_COUNT, MDL_XCOOR, MDL_XCOOR, MDL_COUNT)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 0, id)
        mdResult.setValue(MDL_COUNT, 33L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 1, id)
        mdResult.setValue(MDL_COUNT, 34L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 2, id)
        mdResult.setValue(MDL_COUNT, 33L, id)
        
        self.assertEqual(mdOut, mdResult)
            
    def test_Metadata_aggregateMdGroupBy(self):
        mdOut = MetaData()
        md    = MetaData()
        mdResult = MetaData()
        for linea in range(1, 101):
            id = md.addObject()
            md.setValue(MDL_IMAGE, '%06d@proj.stk' % linea, id)
            md.setValue(MDL_XCOOR, int(linea%3), id)
            md.setValue(MDL_YCOOR, int(linea%2), id)
        mdOut.aggregateMdGroupBy(md, AGGR_COUNT, [MDL_XCOOR,MDL_YCOOR], MDL_XCOOR, MDL_COUNT)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 0, id)
        mdResult.setValue(MDL_YCOOR, 0, id)
        mdResult.setValue(MDL_COUNT, 16L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 0, id)
        mdResult.setValue(MDL_YCOOR, 1, id)
        mdResult.setValue(MDL_COUNT, 17L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 1, id)
        mdResult.setValue(MDL_YCOOR, 0, id)
        mdResult.setValue(MDL_COUNT, 17L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 1, id)
        mdResult.setValue(MDL_YCOOR, 1, id)
        mdResult.setValue(MDL_COUNT, 17L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 2, id)
        mdResult.setValue(MDL_YCOOR, 0, id)
        mdResult.setValue(MDL_COUNT, 17L, id)
        id = mdResult.addObject()
        mdResult.setValue(MDL_XCOOR, 2, id)
        mdResult.setValue(MDL_YCOOR, 1, id)
        mdResult.setValue(MDL_COUNT, 16L, id)
        
        self.assertEqual(mdOut, mdResult)
            
    def test_Metadata_fillExpand(self):
        '''MetaData_fillExpand'''
        ctf = MetaData()
        objId = ctf.addObject()
        ctf.setColumnFormat(False)
        # Create ctf param metadata 1
        u1, v1 = 10000.0, 10001.0
        u2, v2 = 20000.0, 20001.0
        
        ctf.setValue(MDL_CTF_DEFOCUSU, u1, objId)
        ctf.setValue(MDL_CTF_DEFOCUSV, v1, objId)
        ctfFn1 = self.getTmpName(suffix='_ctf1.param')
        ctf.write(ctfFn1)
        # Create ctf param metadata 2
        ctf.setValue(MDL_CTF_DEFOCUSU, u2, objId)
        ctf.setValue(MDL_CTF_DEFOCUSV, v2, objId)
        ctfFn2 = self.getTmpName(suffix='_ctf2.param')
        ctf.write(ctfFn2)
        
        # True values
        ctfs = [ctfFn1, ctfFn1, ctfFn2]
        defocusU = [u1, u1, u2]
        defocusV = [v1, v1, v2]
        # Create image metadata
        md1 = MetaData()
        for i in range(3):
            objId = md1.addObject() 
            img = '00000%i@pythoninterface/proj_ctf_1.stk' % (i + 1)
            md1.setValue(MDL_IMAGE, img, objId)
            md1.setValue(MDL_CTF_MODEL, ctfs[i], objId)
            
        # Create the expected metadata after call fillExpand
        md2 = MetaData()
        for i in range(3):
            objId = md2.addObject() 
            img = '00000%i@pythoninterface/proj_ctf_1.stk' % (i + 1)
            md2.setValue(MDL_IMAGE, img, objId)
            md2.setValue(MDL_CTF_MODEL, ctfs[i], objId)
            md2.setValue(MDL_CTF_DEFOCUSU, defocusU[i], objId)
            md2.setValue(MDL_CTF_DEFOCUSV, defocusV[i], objId)            
        
        md1.fillExpand(MDL_CTF_MODEL)
        self.assertEqual(md1, md2)
        
#RESULTS test_metadata_stress
#no blocks:1000, no lines=100
#BUFFER FILLING took 11.34 seconds
#write(append) sqlite took 113.89 seconds
#write(append) star took 11.15 seconds
#write(random) sqlite took 215.04 seconds
#write(random) star took 13.18 seconds
#read(random) sqlite took 22.55 seconds
#read(random) star took 9.16 seconds
#
#no blocks:1, no lines=1000000
#BUFFER FILLING took 53.60 seconds
#write(append) sqlite took 35.81 seconds
#write(append) star took 53.84 seconds
#write(random) sqlite took 35.70 seconds
#write(random) star took 53.57 seconds
#read(random) sqlite took 2.34 seconds
#read(random) star took 28.91 seconds
#
#
#    def test_metadata_stress(self):
#        #create set and save 1000 block with 100 lines each
#        numberMetadatas=1000+1
#        numberLines=100+1
#        md = MetaData()
#        fnStar = self.getTmpName('.xmd')
#        fnSqlite = self.getTmpName('.sqlite')
#        outFileName=""
#        #fill buffers so comparison is fair
#        #if you do not fill buffers first meassurements will be muchfaster
#        print "fill buffers so comparisons are fair"
#        tmpFn="/tmp/deleteme.xmd"
#        timestamp1 = time()
#        if exists(tmpFn):
#            os.remove(tmpFn)
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (met,tmpFn)
#            for lin in range(1, numberLines):
#                id = md.addObject() 
#                
#                md.setValue(MDL_IMAGE, '%06d@proj.stk' % met, id)
#                md.setValue(MDL_XCOOR,   lin+met*numberLines, id)
#                md.setValue(MDL_YCOOR, 1+lin+met*numberLines, id)
#            md.write(outFileName,MD_APPEND)
#            md.clear()
#        timestamp2 = time()
#        print "\nBUFFER FILLING took %.2f seconds" % (timestamp2 - timestamp1)
#        
#
#        #write sqlite lineal
#        timestamp1 = time()
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (met,fnSqlite)
#            for lin in range(1, numberLines):
#                id = md.addObject() 
#                
#                md.setValue(MDL_IMAGE, '%06d@proj.stk' % met, id)
#                md.setValue(MDL_XCOOR,   lin+met*numberLines, id)
#                md.setValue(MDL_YCOOR, 1+lin+met*numberLines, id)
#            md.write(outFileName,MD_APPEND) 
#            md.clear()
#        timestamp2 = time()
#        print "\nwrite(append) sqlite took %.2f seconds" % (timestamp2 - timestamp1)
#        #write star lineal
#        timestamp1 = time()
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (met,fnStar)
#            for lin in range(1, numberLines):
#                id = md.addObject() 
#                
#                md.setValue(MDL_IMAGE, '%06d@proj.stk' % met, id)
#                md.setValue(MDL_XCOOR,   lin+met*numberLines, id)
#                md.setValue(MDL_YCOOR, 1+lin+met*numberLines, id)
#            md.write(outFileName,MD_APPEND) 
#            md.clear()
#        timestamp2 = time()
#        print "\nwrite(append) star took %.2f seconds" % (timestamp2 - timestamp1)
#        self.assertEqual(1, 1)
#        #sqlite random print
#        timestamp1 = time()
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (randint(1,numberMetadatas),fnSqlite)
#            for lin in range(1, numberLines):
#                id = md.addObject() 
#                
#                md.setValue(MDL_IMAGE, '%06d@proj.stk' % met, id)
#                md.setValue(MDL_XCOOR,   lin+met*numberLines, id)
#                md.setValue(MDL_YCOOR, 1+lin+met*numberLines, id)
#            md.write(outFileName,MD_APPEND) 
#            md.clear()
#        timestamp2 = time()
#        print "\nwrite(random) sqlite took %.2f seconds" % (timestamp2 - timestamp1)
#        #star random
#        timestamp1 = time()
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (randint(1,numberMetadatas),fnStar)
#            for lin in range(1, numberLines):
#                id = md.addObject() 
#                
#                md.setValue(MDL_IMAGE, '%06d@proj.stk' % met, id)
#                md.setValue(MDL_XCOOR,   lin+met*numberLines, id)
#                md.setValue(MDL_YCOOR, 1+lin+met*numberLines, id)
#            md.write(outFileName,MD_APPEND) 
#            md.clear()
#        timestamp2 = time()
#        print "\nwrite(random) star took %.2f seconds" % (timestamp2 - timestamp1)
#
#        #read random 
#        timestamp1 = time()
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (randint(1,numberMetadatas),fnSqlite)
#            md.read(outFileName) 
#            md.clear()
#        timestamp2 = time()
#        print "\nread(random) sqlite took %.2f seconds" % (timestamp2 - timestamp1)
#
#        timestamp1 = time()
#        for met in range(1, numberMetadatas):
#            outFileName = "b%06d@%s" % (randint(1,numberMetadatas),fnStar)
#            md.read(outFileName) 
#            md.clear()
#        timestamp2 = time()
#        print "\nread(random) star took %.2f seconds" % (timestamp2 - timestamp1)
#        
    def test_SymList_readSymmetryFile(self):
        '''readSymmetryFile'''
        a = SymList()
        a.readSymmetryFile("i3")
        self.assertEqual(True, True)
        
    def test_SymList_computeDistance(self):
        '''computeDistance'''
        SL = SymList()
        SL.readSymmetryFile("i3h")
        md = MetaData()
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 1., id)
        md.setValue(MDL_ANGLE_TILT, 2., id)
        md.setValue(MDL_ANGLE_PSI, 3., id)
        md.setValue(MDL_ANGLE_ROT2, 6., id)
        md.setValue(MDL_ANGLE_TILT2, 5., id)
        md.setValue(MDL_ANGLE_PSI2, 4., id)
        id = md.addObject()
        md.setValue(MDL_ANGLE_ROT, 11., id)
        md.setValue(MDL_ANGLE_TILT, 12., id)
        md.setValue(MDL_ANGLE_PSI, 13., id)
        md.setValue(MDL_ANGLE_ROT2, 16., id)
        md.setValue(MDL_ANGLE_TILT2, 15., id)
        md.setValue(MDL_ANGLE_PSI2, 14., id)
        mdOut = MetaData(md)
        SL.computeDistance(md, False, False, False)
        id = md.firstObject()
        total = md.getValue(MDL_ANGLE_DIFF, id)
        self.assertAlmostEqual(total, 5.23652, 4)
#    //EXPECT_NEAR (total, 5.23652,0.00001);
#    XMIPP_CATCH
#}

       
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
