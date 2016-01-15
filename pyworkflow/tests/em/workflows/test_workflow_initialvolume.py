# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

"""
Test the execution of the initial volume service with different datasets.
Execute:
1) Import averages
2) Run initial volume algorithms (ransac, eman, significant, prime)
3) Run align volumes

Datasets:
1) Groel
2) Ribosome
3) BPV
"""

import os

import pyworkflow.tests as tests
import pyworkflow.em as em
import pyworkflow.em.packages.xmipp3 as xmipp3
import pyworkflow.em.packages.eman2 as eman2
  
   
class TestGroel(tests.BaseTest):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        tests.setupTestProject(cls)
        cls.ds = tests.DataSet.getDataSet('initial_volume')
    
    def test_groel(self):
        """ Run an Import particles protocol. """
        cpus = os.environ.get('SCIPION_TEST_CPU', 4)
        # 1. Run import of averages
        groelAvg = self.ds.getFile('groel')
        sym = 'd7'
        protImport = self.newProtocol(em.ProtImportAverages, 
                                      objLabel='import averages (groel)',
                                      filesPath=groelAvg, 
                                      samplingRate=1)
        self.launchProtocol(protImport)
        
        # 2. Run initial models
        # 2a. Ransac 
        protRansac = self.newProtocol(xmipp3.XmippProtRansac,
                                      objLabel='xmipp - ransac',
                                      symmetryGroup=sym,
                                      numberOfMpi=1,
                                      numberOfThreads=cpus
                                      )
        protRansac.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protRansac)
        
        # 2b. Eman 
        protEmanInitVol = self.newProtocol(eman2.EmanProtInitModel,
                                           objLabel='eman - initial vol',
                                           symmetry=sym,
                                           numberOfThreads=cpus)
        protEmanInitVol.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protEmanInitVol)
        
        # 3. Significant
        protSignificant = self.newProtocol(xmipp3.XmippProtReconstructSignificant,
                                      objLabel='xmipp - significant',
                                      symmetryGroup=sym,
                                      numberOfMpi=cpus,
                                      numberOfThreads=1,
                                      iter=15,
                                      alpha0=95
                                      )
        protSignificant.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protSignificant)
        
        # 4. Align all volumes
        protAlign = self.newProtocol(xmipp3.XmippProtAlignVolumeForWeb,
                                      objLabel='xmipp - align volumes',
                                      numberOfMpi=1,
                                      numberOfThreads=cpus
                                      )
        protAlign.inputReference.set(protSignificant.outputVolume)
        protAlign.inputVolumes.append(protRansac.outputVolumes)
        protAlign.inputVolumes.append(protEmanInitVol.outputVolumes)
        protAlign.inputVolumes.append(protSignificant.outputVolume)
        self.launchProtocol(protAlign)        
                
      
      
class TestBPV(tests.BaseTest):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        tests.setupTestProject(cls)
        cls.ds = tests.DataSet.getDataSet('initial_volume')
    
    def test_bpv(self):
        """ Run an Import particles protocol. """
        cpus = os.environ.get('SCIPION_TEST_CPU', 4)
        # 1. Run import of averages
        groelAvg = self.ds.getFile('bpv')
        sym = 'i1'
        protImport = self.newProtocol(em.ProtImportAverages, 
                                      objLabel='import averages (bpv)',
                                      filesPath=groelAvg, 
                                      samplingRate=1)
        self.launchProtocol(protImport)
        
        # 2. Run initial models
        # 2a. Ransac 
        protRansac = self.newProtocol(xmipp3.XmippProtRansac,
                                      objLabel='xmipp - ransac',
                                      objComment='Since there are only 8 projections, a dimensionality reduction cannot be safely done. In this case, it is better to take only 3 images in every RANSAC iterations and lower the inlier threshold to 0.65 so that more images can have the chances of being considered during the reconstruction process.',
                                      symmetryGroup=sym,
                                      dimRed=False,
                                      numSamples=3, # less than 8 classes provided
                                      corrThresh=0.65, 
                                      numberOfMpi=1,
                                      numberOfThreads=cpus
                                      )
        protRansac.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protRansac)
        
        # 2b. Eman 
        protEmanInitVol = self.newProtocol(eman2.EmanProtInitModel,
                                           objLabel='eman - initial vol',
                                           symmetry='icos',
                                           numberOfThreads=cpus)
        protEmanInitVol.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protEmanInitVol)  
         
        # 3. Significant
        protSignificant = self.newProtocol(xmipp3.XmippProtReconstructSignificant,
                                      objLabel='xmipp - significant',
                                      symmetryGroup=sym,
                                      numberOfMpi=cpus,
                                      numberOfThreads=1,
                                      iter=15,
                                      alpha0=99.0
                                      )
        protSignificant.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protSignificant)
        
        # 4. Align all volumes
        protAlign = self.newProtocol(xmipp3.XmippProtAlignVolumeForWeb,
                                      objLabel='xmipp - align volumes',
                                      numberOfMpi=1,
                                      numberOfThreads=cpus
                                      )
        protAlign.inputReference.set(protSignificant.outputVolume)
        protAlign.inputVolumes.append(protRansac.outputVolumes)
        protAlign.inputVolumes.append(protEmanInitVol.outputVolumes)
        protAlign.inputVolumes.append(protSignificant.outputVolume)
        self.launchProtocol(protAlign)        
        
        
        
             
class TestRibosome(tests.BaseTest):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        tests.setupTestProject(cls)
        cls.ds = tests.DataSet.getDataSet('initial_volume')
    
    def test_ribosome(self):
        """ Run an Import particles protocol. """
        cpus = os.environ.get('SCIPION_TEST_CPU', 4)
        # 1. Run import of averages
        groelAvg = self.ds.getFile('ribosome')
        sym = 'c1'
        protImport = self.newProtocol(em.ProtImportAverages, 
                                      objLabel='import averages (ribosome)',
                                      filesPath=groelAvg, 
                                      samplingRate=1)
        self.launchProtocol(protImport)
        
        # 2. Run initial models
        # 2a. Ransac 
        protRansac = self.newProtocol(xmipp3.XmippProtRansac,
                                      objLabel='xmipp - ransac',
                                      symmetryGroup=sym,
                                      numberOfMpi=1,
                                      numberOfThreads=cpus
                                      )
        protRansac.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protRansac)
        
        # 2b. Eman 
        protEmanInitVol = self.newProtocol(eman2.EmanProtInitModel,
                                           objLabel='eman - initial vol',
                                           symmetry=sym,
                                           numberOfThreads=cpus)
        protEmanInitVol.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protEmanInitVol)  
         
        # 3. Significant
        protSignificant = self.newProtocol(xmipp3.XmippProtReconstructSignificant,
                                      objLabel='xmipp - significant',
                                      symmetryGroup=sym,
                                      numberOfMpi=cpus,
                                      numberOfThreads=1,
                                      iter=15
                                      )
        protSignificant.inputSet.set(protImport.outputAverages)
        self.launchProtocol(protSignificant)
        
        # 4. Align all volumes
        protAlign = self.newProtocol(xmipp3.XmippProtAlignVolumeForWeb,
                                      objLabel='xmipp - align volumes',
                                      numberOfMpi=1,
                                      numberOfThreads=cpus
                                      )
        protAlign.inputReference.set(protSignificant.outputVolume)
        protAlign.inputVolumes.append(protRansac.outputVolumes)
        protAlign.inputVolumes.append(protEmanInitVol.outputVolumes)
        protAlign.inputVolumes.append(protSignificant.outputVolume)
        self.launchProtocol(protAlign) 
                
             
class TestSignificant(tests.BaseTest):
    """ Test only significant execution with BPV virus. """
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        tests.setupTestProject(cls)
        cls.ds = tests.DataSet.getDataSet('initial_volume')
    
    def _runSignificant(self, inputSet, args):
        myargs = dict(args)
        prot1 = self.newProtocol(xmipp3.XmippProtReconstructSignificant,
            objLabel='significant d7', 
            **myargs
            )
        #prot1.inputClasses.set(inputSet)
        self.launchProtocol(prot1)
        
        output = prot1.outputVolume
        
        myargs['thereisRefVolume'] = True
         
        prot2 = self.newProtocol(xmipp3.XmippProtReconstructSignificant,
            objLabel='significant d7 (with ref)', 
            **myargs
            )
        #prot2.inputClasses.set(inputSet)
        prot2.refVolume.set(output)
        self.launchProtocol(prot2)
        
        
    def test_significant(self):
        """ Run an Import particles protocol. """
        cpus = os.environ.get('SCIPION_TEST_CPU', 4)
        # 1. Run import of averages
        avg = self.ds.getFile('groel')
        
        protImport = self.newProtocol(em.ProtImportAverages, 
                                      filesPath=avg, 
                                      samplingRate=1)
        self.launchProtocol(protImport)
        
        args = {'symmetryGroup': 'd7',
                'iter': 3,
                'numberOfMpi': cpus,
                'inputSet': protImport.outputAverages,
                'alpha0':98.0
                }
        # Run significant with one volume
        self._runSignificant(protImport.outputAverages, args)
        