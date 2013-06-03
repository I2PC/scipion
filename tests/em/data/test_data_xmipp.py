'''
Created on May 20, 2013

@author: laura
'''

import os, sys
from glob import glob
import unittest
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.tests import *


class TestXmippSetOfMicrographs(unittest.TestCase):
    
    def setUp(self):
        self.outputPath = getOutputPath('test_data_xmipp')
        cleanPath(self.outputPath)
        makePath(self.outputPath)
        
        self.mdGold = getGoldPath('micrographs_gold.xmd')
        self.dbGold = getGoldPath('micrographs_gold.sqlite')
        
        self.micsPattern = getInputPath('MicrographsTilted', '*.mrc')
        
        self.mdFn = getOutputPath(self.outputPath, 'micrographs.xmd')
        self.dbFn = getOutputPath(self.outputPath, 'micrographs.sqlite')
        
        self.mics = glob(self.micsPattern)
        
        if len(self.mics) == 0:
            raise Exception('There is not mics matching pattern')
        self.mics.sort()
                    
        self.tiltedDict = self.createTiltedDict()
                
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
        """ Test reading an XmippSetOfMicrographs from a metadata """
        xmippSet = XmippSetOfMicrographs(self.mdFn, tiltPairs=True)
        
        #Check that micrographs on metadata corresponds to the input ones (same order too)
        #TODO: Check against gold metadata????
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
        
    def createSetOfMicrographs(self):
        setMics = SetOfMicrographs(self.dbFn, tiltPairs=True)
        setMics.setSamplingRate(1.2)
        for fn in self.mics:
            mic = Micrograph(fn)
            setMics.append(mic)
            if fn in self.tiltedDict:
                mic_u = mic
            else:
                setMics.appendPair(mic_u.getObjId(), mic.getObjId()) 
            
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
    
    def createTiltedDict(self):
        """ Create tilted pairs """
        tiltedDict = {}
        for i, fn in enumerate(self.mics):
            if i%2==0:
                fn_u = fn
            else:
                tiltedDict[fn_u] = fn  
        return tiltedDict
            
if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromName('test_xmipp_data.TestXmippSetOfMicrographs.testMerge')
#    unittest.TextTestRunner(verbosity=2).run(suite)
    
    unittest.main()