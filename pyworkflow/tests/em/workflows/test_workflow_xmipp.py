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
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.em.packages.xmipp3.constants import OTHER 
from pyworkflow.em.packages.xmipp3.protocol_projmatch import XmippProtProjMatch
from test_workflow import TestWorkflow
from pyworkflow.em.protocol import ProtImportCoordinates
   
       
class TestXmippWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
    
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.micsFn, samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(), "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        #Import a set of volumes        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         filesPath=self.vol1, samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)        

        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs,
                                            doDownsample=True, downFactor=5, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling)
        
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = self.newProtocol(XmippProtCTFMicrographs,
                                   numberOfThreads=4, minDefocus=2.2, maxDefocus=2.5)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
        # After CTF estimation, the output micrograph should have CTF info
        self.validateFiles('protCTF', protCTF)
        

        protPP = self.newProtocol(ProtImportCoordinates,
                                  importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                  filesPath=self.allCrdsDir,
                                  filesPattern='*.pos', boxSize=110)

        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
        self.protDict['protPicking'] = protPP
        self.assertIsNotNone(protPP.outputCoordinates,
                             "There was a problem with the import of coordinates")

        print "Run extract particles with other downsampling factor"
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=64, downsampleType=OTHER,
                                       doFlip=True,
                                       runMode=1, doInvert=True)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)
        
        print "Run Screen Particles"
        protScreen = self.newProtocol(XmippProtScreenParticles,
                                      autoParRejection=XmippProtScreenParticles.REJ_MAXZSCORE, maxZscore=3.0)
        protScreen.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protScreen)
        self.assertIsNotNone(protScreen.outputParticles, "There was a problem with Screen Particles")
        
        print "Run ML2D"
        protML2D = self.newProtocol(XmippProtML2D, numberOfClasses=1, maxIters=4, doMlf=True,
                                 numberOfMpi=2, numberOfThreads=1)
        protML2D.inputParticles.set(protScreen.outputParticles)
        self.launchProtocol(protML2D)        
        self.assertIsNotNone(protML2D.outputClasses, "There was a problem with ML2D") 
        self.validateFiles('protML2D', protML2D)
        
        print "Run CL2D"
        protCL2D = self.newProtocol(XmippProtCL2D, numberOfClasses=2, numberOfInitialClasses=1,
                                 numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protCL2D)        
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")
        self.validateFiles('protCL2D', protCL2D) 

        print "Run Only Align2d"
        protOnlyAlign = self.newProtocol(XmippProtCL2DAlign,
                                         maximumShift=5, numberOfIterations=2,
                                         numberOfMpi=2, numberOfThreads=1, useReferenceImage=False)
        protOnlyAlign.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protOnlyAlign)        
        self.assertIsNotNone(protOnlyAlign.outputParticles, "There was a problem with Only align2d")  
        self.validateFiles('protOnlyAlign', protOnlyAlign)
        
        print "Run kerdensom"
        ProtKerdensom = self.newProtocol(XmippProtKerdensom, useMask=False, SomXdim=2, SomYdim=2,
                                 SomReg0=800, SomReg1=400, SomSteps=2)
        ProtKerdensom.inputParticles.set(protOnlyAlign.outputParticles)
        self.launchProtocol(ProtKerdensom)        
        self.assertIsNotNone(ProtKerdensom.outputClasses, "There was a problem with kerdensom")  
        #self.validateFiles('ProtKerdensom', ProtKerdensom)
        
        print "Run Rotational Spectra"
        xmippProtRotSpectra = self.newProtocol(XmippProtRotSpectra, 
                                               SomXdim=2, SomYdim=2,
                                               spectraInnerRadius=4,
                                               spectraOuterRadius=24)
        xmippProtRotSpectra.inputParticles.set(protOnlyAlign.outputParticles)
        self.launchProtocol(xmippProtRotSpectra)        
        self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")
        


if __name__ == "__main__":
    unittest.main()
