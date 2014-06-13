'''
Created on May 20, 2013

@author: laura
'''

from glob import glob, iglob
import unittest
import os, traceback
from itertools import izip
from pyworkflow.tests import *
from pyworkflow.em.data import *
from pyworkflow.utils.path import makePath
from pyworkflow.em.packages.xmipp3.convert import *
import pyworkflow.em.packages.eman2.convert as e2convert
from pyworkflow.em.protocol import EMProtocol



class TestImage(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.mic1 = cls.dataset.getFile( 'mic1')
    
    def testLocation(self):
        fn = self.mic1
        mic = Micrograph()
        mic.setFileName(fn)

        # Check that setFileName-getFileName is working properly        
        self.assertEqual(fn, mic.getFileName())
        
        
class TestSetOfMicrographs(BaseTest):
    
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.dbGold = cls.dataset.getFile( 'micsGoldSqlite')
        cls.micsPattern = cls.dataset.getFile('allMics')
        cls.dbFn = cls.getOutputPath('micrographs.sqlite')
        cls.acquisition = Acquisition(magnification=60000, voltage=300, sphericalAberration=2, amplitudeContrast=0.1)
        
        #cls.mics = glob(cls.micsPattern)
        cls.mics = []
        for mic in iglob(cls.micsPattern):
            cls.mics.append(cls.getRelPath(cls.dataset.getPath(), mic))
        
        if len(cls.mics) == 0:
            raise Exception('There are not micrographs matching pattern')
        cls.mics.sort()                 
    
        
    def checkSet(self, micSet):
        idCount = 1
    
        for mic1, fn in izip(micSet, self.mics):  
            #traceback.print_stack(file=sys.stdout)
            micFn = mic1.getFileName()
            self.assertEqual(fn, micFn, 
                             "micrograph NAME in the set is wrong, \n   expected: '%s'\n        got: '%s'" 
                             % (fn, micFn))
            self.assertEqual(idCount, mic1.getObjId(), 
                             "micrograph ID in the set is wrong, \n   expected: '%s'\n        got: '%s'" 
                             % (idCount, mic1.getObjId()))
            mic2 = micSet[idCount] # Test getitem
            self.assertEqual(mic1.getObjId(), mic2.getObjId(), "micrograph got from ID is wrong")
            idCount += 1           
             
    def testCreate(self):
        cwd = os.getcwd()
        # Change to test path
        os.chdir(self.dataset.getPath())
        """ Create a SetOfMicrographs from a list of micrographs """
        micSet = SetOfMicrographs(filename=self.dbFn)
        
        micSet.setAcquisition(self.acquisition)
        micSet.setSamplingRate(1.2)
        for fn in self.mics:
            mic = Micrograph()
            mic.setFileName(fn)
            micSet.append(mic)
        micSet.write()    
        
        self.checkSet(micSet)
        os.chdir(cwd)
        
    def testRead(self):
        """ Read micrographs from a SetOfMicrographs """
        micSet = SetOfMicrographs(filename=self.dbGold)
        self.checkSet(micSet)
        
#    def testXmippConvert(self):
#        """ Test the convertion of a SetOfMicrographs to Xmipp"""
#        micSet = SetOfMicrographs()
#        micSet.setFileName(self.dbGold)
#        mdFn = getOutputPath('test_data', 'micrographs.xmd')
#        
#        writeSetOfMicrographs(micSet, mdFn)
#        
#        # Test reading a set of coordinates
#        posDir = getInputPath('Picking_XmippBPV3_Down3')
#        print "reading pos from :", posDir
#        coordSet = SetOfCoordinates()
#        fn = getOutputPath('test_data', 'coordinates.sqlite')
#        coordSet.setFileName(fn)
#        
#        readSetOfCoordinates(posDir, micSet, coordSet)
#        coordSet.write()
#        
#        
#        cwd = os.getcwd()
#        # Change to test path
#        os.chdir(getPath())
#        
#        # Test writing micrgraphs to an hdf        
#        filename = getOutputPath('test_data', 'micrographs.hdf')
#        e2convert.writeSetOfParticles(micSet, filename)
#        
#        # Test writing a set of particles
#        #partSet = SetOfParticles()
#        #readSetOfParticles(fnImages, imgSet)
#        
#        os.chdir(cwd)


class TestSetOfParticles(BaseTest):
    """ Check if the information of the images is copied to another image when a new SetOfParticles is created"""
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
#        
#    def testCreateFromOther(self):
#        inImgSet = SetOfParticles(filename=self.dbGold)
#        inImgSet.setHasCTF(True)
#        outImgFn = self.getOutputPath("particles.sqlite")
#        outImgSet = SetOfParticles(filename=outImgFn)
#        outImgSet.copyInfo(inImgSet)
#        
#        print "inputs particles has CTF?", inImgSet.hasCTF()
#        teststk = self.getOutputPath('test.stk')
#        for i, img in enumerate(inImgSet):
#            j = i + 1
#            img.setLocation(j, teststk)
#            outImgSet.append(img)
#        
#        outImgSet.write()
#        
#        if outImgSet.hasCTF():
#            print "everything OK!"
#        else:
#            print "The info of the particles was not copied"
#        
#        cleanPath(outImgFn)
        
    def test_readStack(self):
        """ Read an stack of 29 particles from .hdf file.
        Particles should be of 500x500 pixels.
        """
        size = 29
        xdim = 500        
        inStack = self.dataset.getFile( 'particles1')
        outFn = self.getOutputPath('particles.sqlite')
        
        imgSet = SetOfParticles(filename=outFn)
        imgSet.setSamplingRate(1.0)
        imgSet.readStack(inStack) # This should add 29 new items to the set
        
        self.assertEquals(size, imgSet.getSize()) # Check same size
        self.assertEquals(xdim, imgSet.getDim()[0]) # Check same dimensions
        
        print "writing particles to: ", outFn
        imgSet.write()
        
    def test_hugeSet(self):
        """ Create a set of a big number of particles to measure
        creation time with sqlite operations. 
        """
        # Allow what huge means to be defined with environment var
        n = int(os.environ.get('SCIPION_TEST_HUGE', 10000))
        print ">>>> Creating a set of %d particles." % n
        
        dbFn = self.getOutputPath('huge_set.sqlite')
        #dbFn = ':memory:'
        
        img = Particle()
        imgSet = SetOfParticles(filename=dbFn)
        imgSet.setSamplingRate(1.0)
        
        for i in range(1, n+1):
            img.setLocation(i, "images.stk")
            
            imgSet.append(img)
            img.cleanObjId()
            
        imgSet.write()
        
            
        

class TestSetOfClasses2D(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  
        cls.outputParticles = cls.getOutputPath('output_particles.sqlite')
        #cls.dbGold = cls.dataset.getFile('particlesGold.sqlite')
        
#    def testCreateFromOther(self):
#        inImgSet = SetOfParticles(filename=self.dbGold)
#        inImgSet.setHasCTF(True)
#        outImgFn = self.outputPath + "_particles.sqlite"
#        outImgSet = SetOfParticles(filename=outImgFn)
#        outImgSet.copyInfo(inImgSet)
#        teststk = self.getOutputPath('test.stk')
#        print "inputs particles has CTF?", inImgSet.hasCTF()
#        for i, img in enumerate(inImgSet):
#            j = i + 1
#            img.setLocation(j, teststk)
#            outImgSet.append(img)
#        
#        outImgSet.write()
#        
#        if outImgSet.hasCTF():
#            print "everything OK!"
#        else:
#            print "The info of the particles was not copied"
#        
#        cleanPath(outImgFn)

if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippCTFModel.testConvertXmippCtf')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()
