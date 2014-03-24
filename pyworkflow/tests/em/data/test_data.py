'''
Created on May 20, 2013

@author: laura
'''

from glob import glob, iglob
import unittest
from pyworkflow.tests import *
from pyworkflow.em.data import *
from pyworkflow.utils.path import makePath
from pyworkflow.em.packages.xmipp3.convert import *
import pyworkflow.em.packages.eman2.convert as e2convert
from pyworkflow.em.protocol import EMProtocol



class TestImage(unittest.TestCase):
        
    def setUp(self):
        pass
    
    def testLocation(self):
        fn = 'mic0001.mrc'
        mic = Micrograph()
        mic.setFileName(fn)

        # Check that setFileName-getFileName is working properly        
        self.assertEqual(fn, mic.getFileName())
        
        
class TestSetOfMicrographs(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data')
        
        cls.dbGold = getGoldPath('Micrographs_BPV3', 'micrographs_gold.sqlite')
        
        cls.micsPattern = getInputPath('Micrographs_BPV3', '*.mrc')
        
        cls.dbFn = getOutputPath(cls.outputPath, 'micrographs.sqlite')
        
        #cls.mics = glob(cls.micsPattern)
        cls.mics = []
        for mic in iglob(cls.micsPattern):
            cls.mics.append(getRelPath(mic))
        
        if len(cls.mics) == 0:
            raise Exception('There are not micrographs matching pattern')
        cls.mics.sort()
                  
        cleanPath(cls.outputPath)
        makePath(cls.outputPath)
    
        
    def checkSet(self, micSet):
        idCount = 1
        micSet.loadIfEmpty()
        
        for fn, mic in zip(self.mics, micSet):            
            micFn = mic.getFileName()
            self.assertEqual(fn, micFn, 
                             "micrograph NAME in the set is wrong, \n   expected: '%s'\n        got: '%s'" 
                             % (fn, micFn))
            self.assertEqual(idCount, mic.getObjId(), 
                             "micrograph ID in the set is wrong, \n   expected: '%s'\n        got: '%s'" 
                             % (idCount, mic.getObjId()))
            mic2 = micSet[idCount] # Test getitem
            self.assertEqual(mic.getObjId(), mic2.getObjId(), "micrograph got from ID is wrong")
            idCount += 1            
        
    def testCreate(self):
        """ Create a SetOfMicrographs from a list of micrographs """
        micSet = SetOfMicrographs()
        micSet.setFileName(self.dbFn)
        micSet.setSamplingRate(1.2)
        for fn in self.mics:
            mic = Micrograph()
            mic.setFileName(fn)
            micSet.append(mic)
            
        micSet.write()    
        self.checkSet(micSet)
        
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


class TestSetOfParticles(unittest.TestCase):
    """ Check if the information of the images is copied to another image when a new SetOfParticles is created"""
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data')
        cleanPath(cls.outputPath)
        makePath(cls.outputPath)
        
        cls.outputParticles = getOutputPath('test_data', 'output_particles.sqlite')
        
        cls.dbGold = getGoldPath('SetOfParticles', 'input_particles.sqlite')
        
    def aaatestCreateFromOther(self):
        inImgSet = SetOfParticles(filename=self.dbGold)
        inImgSet.setHasCTF(True)
        outImgFn = self.outputPath + "_particles.sqlite"
        outImgSet = SetOfParticles(filename=outImgFn)
        outImgSet.copyInfo(inImgSet)
        
        print "inputs particles has CTF?", inImgSet.hasCTF()
        for i, img in enumerate(inImgSet):
            j = i + 1
            img.setLocation(j, "test.stk")
            outImgSet.append(img)
        
        outImgSet.write()
        
        if outImgSet.hasCTF():
            print "everything OK!"
        else:
            print "The info of the particles was not copied"
        
        cleanPath(outImgFn)
        
    def test_str(self):
        """ Test the string representation of a SetOfParticles. """
        stackFn = getInputPath('images_LTA.stk')
        imgSet = SetOfParticles(filename=self.outputParticles)
        imgSet.setSamplingRate(1.0)
        imgSet.readStack(stackFn)
        imgSet.printAll()
        
        print imgSet
        
#        for img in imgSet:
#            print img.getLocation()
        
        

class TestSetOfClasses2D(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data')
        
        cls.dbGold = getGoldPath('SetOfParticles', 'input_particles.sqlite')
        
    def testCreateFromOther(self):
        inImgSet = SetOfParticles(filename=self.dbGold)
        inImgSet.setHasCTF(True)
        outImgFn = self.outputPath + "_particles.sqlite"
        outImgSet = SetOfParticles(filename=outImgFn)
        outImgSet.copyInfo(inImgSet)
        
        print "inputs particles has CTF?", inImgSet.hasCTF()
        for i, img in enumerate(inImgSet):
            j = i + 1
            img.setLocation(j, "test.stk")
            outImgSet.append(img)
        
        outImgSet.write()
        
        if outImgSet.hasCTF():
            print "everything OK!"
        else:
            print "The info of the particles was not copied"
        
        cleanPath(outImgFn)

if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippCTFModel.testConvertXmippCtf')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()