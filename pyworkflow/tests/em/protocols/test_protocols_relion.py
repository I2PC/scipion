# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *


# Some utility functions to import micrographs that are used
# in several tests.
class TestRelionBase(unittest.TestCase):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        protImport = ProtImportParticles(pattern=pattern, samplingRate=samplingRate, checkStack=checkStack)
        cls.proj.launchProtocol(protImport, wait=True)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return protImport
    
    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = ProtImportVolumes(pattern=pattern, samplingRate=samplingRate)
        cls.proj.launchProtocol(protImport, wait=True)
        return protImport


class TestRelionClassify2D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
            
    def testRelion2D(self):
        print "Run relion2D"
        prot2D = ProtRelionClassify2D(doCTF=False, maskRadiusA=170, numberOfMpi=4, numberOfThreads=1)
        prot2D.numberOfClasses.set(4)
        prot2D.numberOfIterations.set(3)
        prot2D.inputParticles.set(self.protImport.outputParticles)
        self.proj.launchProtocol(prot2D, wait=True)        
        
        self.assertIsNotNone(getattr(prot2D, 'outputClasses', None), 
                             "There was a problem with Relion 2D:\n" + prot2D.getErrorMessage()) 


class TestRelionClassify3D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionClassify3D(self):
        print "Run ProtRelionClassify3D"
        relion3DClass = ProtRelionClassify3D(numberOfClasses=3, numberOfIterations=4, doCTF=False, runMode=1, 
                                 numberOfMpi=2, numberOfThreads=2)
        relion3DClass.inputParticles.set(self.protImport.outputParticles)
        relion3DClass.referenceVolume.set(self.protImportVol.outputVolume)
        self.proj.launchProtocol(relion3DClass, wait=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        className = sys.argv[1]
        cls = globals().get(className, None)
        if cls:
            suite = unittest.TestLoader().loadTestsFromTestCase(cls)
            unittest.TextTestRunner(verbosity=2).run(suite)
        else:
            print "Test: '%s' not found." % className
    else:
        unittest.main()
