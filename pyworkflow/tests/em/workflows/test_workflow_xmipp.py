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
   
       
class TestXmippWorkflow(TestWorkflow):
#     GOLD_FILES = {'protImport': ['protImport/BPV_1388.mrc',
#                     'protImport/micrographs.sqlite', 
#                     'protImport/BPV_1387.mrc',
#                     'protImport/BPV_1386.mrc'],
#               'protDownsampling': [
#                     'protImport/BPV_1386.mrc', 
#                     'protImport/BPV_1388.mrc', 
#                     'protImport/BPV_1387.mrc',
#                     'protImport/micrographs.sqlite', 
#                     'protDownsampling/BPV_1386.mrc', 
#                     'protDownsampling/BPV_1387.mrc', 
#                     'protDownsampling/BPV_1388.mrc', 
#                     'protDownsampling/micrographs.sqlite', 
#                     'protDownsampling/logs/run.log',
#                     'protDownsampling/logs/run.db', 
#                     ],
#               'protCTF': ['protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_quadrant.xmp', 
#                     'protCTF/logs/run.log', 
#                     'protCTF/logs/run.db', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.psd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf_enhanced_psd.xmp', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.psd', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_halfplane.xmp', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf_enhanced_psd.xmp', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf_ctfmodel_halfplane.xmp', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf_enhanced_psd.xmp', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_quadrant.xmp', 
#                     'protDownsampling/BPV_1388.mrc', 
#                     'protDownsampling/BPV_1387.mrc', 
#                     'protDownsampling/micrographs.sqlite', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.psd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf_ctfmodel_halfplane.xmp', 
#                     'protCTF/micrographs.sqlite',
#                     'protCTF/micrographs.xmd',  
#                     'protDownsampling/BPV_1386.mrc', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf_ctfmodel_quadrant.xmp'],
#               'protExtract':[
#                     'protImport/BPV_1386.mrc',
#                     'protImport/BPV_1388.mrc',
#                     'protImport/BPV_1387.mrc',
#                     'protImport/micrographs.sqlite',
#                     'protPicking/coordinates.sqlite',
#                     'protExtract/tmp/BPV_1388_flipped.xmp', 
#                     'protExtract/tmp/BPV_1387_flipped.xmp', 
#                     'protExtract/tmp/BPV_1386_noDust.xmp', 
#                     'protExtract/extra/BPV_1386.xmd', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protExtract/images.xmd', 
#                     'protExtract/particles.sqlite',
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protExtract/extra/BPV_1388.xmd', 
#                     'protExtract/extra/BPV_1387.stk', 
#                     'protExtract/tmp/BPV_1387_noDust.xmp', 
#                     'protExtract/tmp/BPV_1388_noDust.xmp', 
#                     'protExtract/extra/BPV_1387.xmd',
#                     'protExtract/extra/BPV_1386.pos',
#                     'protExtract/extra/BPV_1387.pos',
#                     'protExtract/extra/BPV_1388.pos',  
#                     'protExtract/tmp/BPV_1386_flipped.xmp',
#                     'protExtract/tmp/BPV_1387_downsampled.xmp',
#                     'protExtract/tmp/BPV_1386_downsampled.xmp',
#                     'protExtract/tmp/BPV_1388_downsampled.xmp',
#                     'protExtract/tmp/BPV_1386.ctfParam',
#                     'protExtract/tmp/BPV_1387.ctfParam',
#                     'protExtract/tmp/BPV_1388.ctfParam',
#                     'protExtract/logs/run.log',
#                     'protExtract/logs/run.db',
#                     ],
#               'protML2D': [
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protExtract/extra/BPV_1387.stk', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protExtract/images.xmd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam',
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam',
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam',
#                     'protML2D/mlf2d_extra/iter002/iter_classes.xmd', 
#                     'protML2D/mlf2d_extra/iter001/iter_classes.xmd', 
#                     'protML2D/mlf2d_extra/iter002/iter_images.xmd', 
#                     'protML2D/mlf2d_classes.stk', 
#                     'protML2D/mlf2d_extra/iter003/iter_images.xmd', 
#                     'protML2D/mlf2d_images.xmd', 
#                     'protML2D/mlf2d_extra/iter004/iter_classes.stk', 
#                     'protML2D/mlf2d_extra/iter001/iter_images.xmd', 
#                     'protML2D/mlf2d_extra/iter002/iter_classes.stk', 
#                     'protML2D/mlf2d_extra/iter004/iter_classes.xmd', 
#                     'protML2D/mlf2d_extra/iter004/iter_images.xmd', 
#                     'protML2D/mlf2d_extra/iter003/iter_classes.xmd', 
#                     'protML2D/mlf2d_extra/iter001/iter_noise.xmd',
#                     'protML2D/mlf2d_extra/iter002/iter_noise.xmd',
#                     'protML2D/mlf2d_extra/iter002/iter_ssnr.xmd',
#                     'protML2D/mlf2d_extra/iter000/iter_classes.xmd',
#                     'protML2D/mlf2d_extra/iter001/iter_ssnr.xmd',
#                     'protML2D/mlf2d_extra/iter003/iter_noise.xmd',
#                     'protML2D/mlf2d_extra/iter000/iter_noise.xmd',
#                     'protML2D/mlf2d_extra/iter004/iter_ssnr.xmd',
#                     'protML2D/mlf2d_noise.xmd',
#                     'protML2D/mlf2d_extra/iter000/iter_classes.stk',
#                     'protML2D/mlf2d_extra/iter003/iter_ssnr.xmd',
#                     'protML2D/mlf2d_extra/cref_classes.stk',
#                     'protML2D/mlf2d_extra/iter000/iter_ssnr.xmd',
#                     'protML2D/mlf2d_extra/cref_classes.xmd',
#                     'protML2D/mlf2d_extra/iter004/iter_noise.xmd',
#                     'protML2D/logs/run.log', 
#                     'protML2D/logs/run.db',
#                     'protML2D/mlf2d_classes.xmd', 
#                     'protML2D/mlf2d_extra/iter001/iter_classes.stk', 
#                     'protML2D/mlf2d_extra/iter003/iter_classes.stk'],
#               'protCL2D': ['protCL2D/extra/classes_core_hierarchy.txt', 
#                     'protCL2D/extra/level_01/level_classes_core.xmd', 
#                     'protCL2D/extra/level_01/level_classes.stk', 
#                     'protCL2D/extra/level_01/classes_sorted.stk', 
#                     'protCL2D/extra/level_00/level_classes.xmd', 
#                     'protCL2D/extra/level_00/classes_sorted.xmd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
#                     'protCL2D/extra/level_00/level_classes_core.stk', 
#                     'protCL2D/logs/run.log', 
#                     'protCL2D/logs/run.db', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
#                     'protCL2D/extra/level_00/classes_core_sorted.xmd', 
#                     'protCL2D/extra/level_01/classes_core_sorted.stk', 
#                     'protExtract/images.xmd', 
#                     'protCL2D/extra/level_00/classes_sorted.stk', 
#                     'protCL2D/extra/level_01/classes_sorted.xmd', 
#                     'protCL2D/extra/classes_hierarchy.txt', 
#                     'protCL2D/extra/images.xmd', 
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protCL2D/extra/level_01/classes_core_sorted.xmd', 
#                     'protCL2D/extra/level_00/classes_core_sorted.stk', 
#                     'protExtract/extra/BPV_1387.stk', 
#                     'protCL2D/extra/level_00/level_classes_core.xmd', 
#                     'protCL2D/extra/level_01/level_classes.xmd', 
#                     'protCL2D/extra/level_00/level_classes.stk', 
#                     'protCL2D/extra/level_01/level_classes_core.stk'],
#               'protOnlyAlign': ['protOnlyAlign/extra/images.xmd', 
#                     'protOnlyAlign/extra/level_00/class_classes.stk', 
#                     'protOnlyAlign/extra/level_00/class_classes.xmd', 
#                     'protCTF/extra/BPV_1386/xmipp_ctf.ctfparam', 
#                     'protCTF/extra/BPV_1388/xmipp_ctf.ctfparam', 
#                     'protOnlyAlign/logs/run.log', 
#                     'protOnlyAlign/logs/run.db', 
#                     'protExtract/extra/BPV_1388.stk', 
#                     'protCTF/extra/BPV_1387/xmipp_ctf.ctfparam', 
#                     'protExtract/images.xmd', 
#                     'protExtract/extra/BPV_1386.stk', 
#                     'protExtract/extra/BPV_1387.stk']
#               }

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')
    
    def testXmippWorkflow(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs, filesPath=self.micsFn, samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(), "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        #Import a set of volumes        
        print "Import Volume"
        protImportVol = self.newProtocol(ProtImportVolumes, filesPath=self.vol1, samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
#        self.validateFiles('protImportVol', protImportVol)        

        # Perform a downsampling on the micrographs
        print "Downsampling..."
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs, doDownsample=True, downFactor=5, doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs, "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling)
        
        # Now estimate CTF on the downsampled micrographs 
        print "Performing CTF..."   
        protCTF = self.newProtocol(XmippProtCTFMicrographs, numberOfThreads=4, minDefocus=2.2, maxDefocus=2.5)                
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF, "There was a problem with the CTF estimation")
        # After CTF estimation, the output micrograph should have CTF info
        self.validateFiles('protCTF', protCTF)
        
        print "Running fake particle picking..."   
        protPP = self.newProtocol(XmippProtParticlePicking, importFolder=self.allCrdsDir)                
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)        
        self.launchProtocol(protPP)
        self.protDict['protPicking'] = protPP
        self.assertIsNotNone(protPP.outputCoordinates, "There was a problem with the faked picking")
            
        print "Run extract particles with other downsampling factor"
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=64, downsampleType=OTHER, doFlip=True, downFactor=8, runMode=1, doInvert=True)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles, "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)
        
        print "Run Screen Particles"
        protScreen = self.newProtocol(XmippProtScreenParticles, autoParRejection=XmippProtScreenParticles.REJ_MAXZSCORE, maxZscore=3.0)
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
        protOnlyAlign = self.newProtocol(XmippProtCL2DAlign, maximumShift=5, numberOfIterations=2, 
                                 numberOfMpi=2, numberOfThreads=1, useReferenceImage=False)
        protOnlyAlign.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protOnlyAlign)        
        self.assertIsNotNone(protOnlyAlign.outputParticles, "There was a problem with Only align2d")  
        self.validateFiles('protOnlyAlign', protOnlyAlign)
        
        print "Run kerdensom"
        ProtKerdensom = self.newProtocol(XmippProtKerdensom, useMask=False, SomXdim=2, SomYdim=2,
                                 SomReg0=800, SomReg1=400, SomSteps=2)
        ProtKerdensom.inputImages.set(protOnlyAlign.outputParticles)
        self.launchProtocol(ProtKerdensom)        
        self.assertIsNotNone(ProtKerdensom.outputClasses, "There was a problem with kerdensom")  
        #self.validateFiles('ProtKerdensom', ProtKerdensom)
        
        print "Run Rotational Spectra"
        xmippProtRotSpectra = self.newProtocol(XmippProtRotSpectra, SomXdim=2, SomYdim=2)
        xmippProtRotSpectra.inputImages.set(protOnlyAlign.outputParticles)
        self.launchProtocol(xmippProtRotSpectra)        
        self.assertIsNotNone(xmippProtRotSpectra.outputClasses, "There was a problem with Rotational Spectra")
        


if __name__ == "__main__":
    unittest.main()
