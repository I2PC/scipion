'''
Created on May 20, 2013

@author: laura
'''

from glob import glob, iglob
import unittest
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.tests import *


class TestXmippSetOfMicrographs(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        cls.outputPath = getOutputPath('test_data_xmipp')
#        cleanPath(cls.outputPath)
#        makePath(cls.outputPath)
        
        cls.mdGold = getGoldPath('Micrographs_TiltedPhantom', 'micrographs_gold.xmd')
        cls.dbGold = getGoldPath('Micrographs_TiltedPhantom', 'micrographs_gold.sqlite')
        
        cls.micsPattern = getInputPath('Micrographs_TiltedPhantom', '*.mrc')
        
        cls.mdFn = getOutputPath(cls.outputPath, 'micrographs.xmd')
        cls.dbFn = getOutputPath(cls.outputPath, 'micrographs.sqlite')
        
        #cls.mics = glob(cls.micsPattern)
        cls.mics = []
        for mic in iglob(cls.micsPattern):
            cls.mics.append(getRelPath(mic))
        
        if len(cls.mics) == 0:
            raise Exception('There are not micrographs matching pattern')
        cls.mics.sort()
    
        """ Create tilted pairs """
        cls.tiltedDict = {}
        for i, fn in enumerate(cls.mics):
            if i%2==0:
                fn_u = fn
            else:
                cls.tiltedDict[fn] = fn_u  
                
    def setUp(self):
        cleanPath(self.outputPath)
        makePath(self.outputPath)
                    
    def testWriteMd(self):
        """ Test creating and writing a XmippSetOfMicrographs from a list of micrographs """
        
        xmippSet = XmippSetOfMicrographs(self.mdFn, tiltPairs=True)
        xmippSet.setSamplingRate(1.2)
        for fn in self.mics:
            mic = XmippMicrograph(fn)
            xmippSet.append(mic)
            
        for u, t in self.tiltedDict.iteritems():
            xmippSet.appendPair(u, t)  
        
        # Write the metadata
        xmippSet.write()

        self.assertTrue(self.checkMicrographsMetaData(xmippSet), "micrographs metadata does not exist")
        
    def testReadMd(self):
        """ Test reading an XmippSetOfMicrographs from an existing  metadata """
        xmippSet = XmippSetOfMicrographs(self.mdGold, tiltPairs=True)
        
        #Check that micrographs on metadata corresponds to the input ones (same order too)
        for i, mic in enumerate(xmippSet):
            self.assertEqual(self.mics[i], mic.getFileName(), "Micrograph %d in set is different from expected" % i)

        for (iU, iT) in xmippSet.iterTiltPairs():
            self.assertEqual(self.tiltedDict[iU], iT, "")
            
    def testConvert(self):
        """ Test converting a SetOfMicrographs to a XmippSetOfMicrographs """
        setMics = self.createSetOfMicrographs()
                
        xmippSet = XmippSetOfMicrographs.convert(setMics, self.mdFn)
        
        self.assertTrue(self.checkMicrographsMetaData(xmippSet), "micrographs metadata does not exist")
        
    def testCopy(self):
        """ Test copying from a SetOfMicrographs to a XmippSetOFMicrographs """
        setMics = self.createSetOfMicrographs()
        
        xmippSet = XmippSetOfMicrographs(self.mdFn)
        
        mapsId = {}
        
        for mic in setMics:
            xmic = XmippMicrograph(mic.getFileName())
            xmippSet.append(xmic)
            mapsId[mic.getId()] = xmic.getId()
            
        xmippSet.copyInfo(setMics)
        
        # If input micrographs have tilt pairs copy the relation
        if setMics.hasTiltPairs():
            xmippSet.copyTiltPairs(setMics, mapsId.get)
            
        xmippSet.write()
        
    def testReadBd(self):
        """ Read micrographs from a SetOfMicrographs """
        setMics = SetOfMicrographs(self.dbGold)

        for i, mic in enumerate(setMics):
            self.assertEqual(self.mics[i], mic.getFileName(), "Micrograph %d in set is different from expected" % i)
        
    def createSetOfMicrographs(self):
        """ Create a SetOfMicrographs from a list of micrographs """
        setMics = SetOfMicrographs(self.dbFn, tiltPairs=True)
        setMics.setSamplingRate(1.2)
        for fn in self.mics:
            mic = Micrograph(fn)
            setMics.append(mic)
            if fn in self.tiltedDict.values():
                mic_t = mic
            else:
                setMics.appendPair(mic.getObjId(), mic_t.getObjId()) 
            
        setMics.write()
        
        self.assertTrue(self.checkMicrographsDb(setMics), "micrographs database does not exist")
        
        return setMics
        
            
    def checkMicrographsMetaData(self, xmippSet):
        """ Check that a metadata is equal to the gold one """
        #TODO: Implement how to check that two metadatas are equal
        md_gold = xmipp.FileName(self.mdGold)
        md_fn = xmipp.FileName(xmippSet.getFileName())
        return (sys.getsizeof(md_fn) == sys.getsizeof(md_gold))

    def checkMicrographsDb(self, setMics):
        """ Check that a database is equal to the gold one """
        #TODO: Implement how to check that two databases are equal
        return (os.path.getsize(setMics.getFileName()) == os.path.getsize(self.dbGold))    
            
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_data_xmipp.TestXmippSetOfMicrographs.testReadBd')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()