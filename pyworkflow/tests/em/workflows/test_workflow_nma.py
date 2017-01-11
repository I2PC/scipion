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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.em.protocol import ProtImportPdb, ProtImportVolumes, ProtImportParticles
from pyworkflow.tests import setupTestProject, DataSet
from test_workflow import TestWorkflow  
from pyworkflow.em.packages.xmipp3.nma import (XmippProtNMA, XmippProtAlignmentNMA, 
                                               XmippProtDimredNMA, NMA_CUTOFF_ABS)
from pyworkflow.em.packages.xmipp3 import XmippProtConvertToPseudoAtoms
   
   
   
class TestNMA(TestWorkflow):
    """ Check the images are converted properly to spider format. """
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('nma')
    
    def test_nma1(self):
        """ Run NMA simple workflow for both Atomic and Pseudoatoms. """
        
        #------------------------------------------------
        # Case 1. Import a Pdb -> NMA
        #------------------------------------------------
        
        # Import a PDB
        protImportPdb = self.newProtocol(ProtImportPdb, inputPdbData=1,
                                      pdbFile=self.ds.getFile('pdb'))
        self.launchProtocol(protImportPdb)
        
        # Launch NMA for PDB imported
        protNMA1 = self.newProtocol(XmippProtNMA,
                                    cutoffMode=NMA_CUTOFF_ABS)
        protNMA1.inputStructure.set(protImportPdb.outputPdb)
        self.launchProtocol(protNMA1)
        
        # Import the set of particles 
        # (in this order just to be in the middle in the tree)
        protImportParts = self.newProtocol(ProtImportParticles,
                                           filesPath=self.ds.getFile('particles'),
                                           samplingRate=1.0)
        self.launchProtocol(protImportParts) 

        # Launch NMA alignment, but just reading result from a previous metadata
        protAlignment = self.newProtocol(XmippProtAlignmentNMA,
                                         modeList='7-9',
                                         copyDeformations=self.ds.getFile('gold/pseudo_run1_images.xmd'))
        protAlignment.inputModes.set(protNMA1.outputModes)
        protAlignment.inputParticles.set(protImportParts.outputParticles)
        self.launchProtocol(protAlignment)       
        
        # Launch Dimred after NMA alignment 
        protDimRed = self.newProtocol(XmippProtDimredNMA,
                                      dimredMethod=0, # PCA
                                      reducedDim=2)
        protDimRed.inputNMA.set(protAlignment)
        self.launchProtocol(protDimRed)     
    
    def test_nma2(self):
        #------------------------------------------------        
        # Case 2. Import Vol -> Pdb -> NMA
        #------------------------------------------------
        # Import the set of particles 
        # (in this order just to be in the middle in the tree)
        protImportParts = self.newProtocol(ProtImportParticles,
                                           filesPath=self.ds.getFile('particles'),
                                           samplingRate=1.0)
        self.launchProtocol(protImportParts) 
        
        # Import a Volume
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         filesPath=self.ds.getFile('vol'),
                                         samplingRate=1.0)
        self.launchProtocol(protImportVol)
        
        # Convert the Volume to Pdb
        from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import NMA_MASK_THRE
        protConvertVol = self.newProtocol(XmippProtConvertToPseudoAtoms)
        protConvertVol.inputStructure.set(protImportVol.outputVolume)
        protConvertVol.maskMode.set(NMA_MASK_THRE)
        protConvertVol.maskThreshold.set(0.2)
        protConvertVol.pseudoAtomRadius.set(2.5)
        self.launchProtocol(protConvertVol)
        
        # Launch NMA with Pseudoatoms
        protNMA2 = self.newProtocol(XmippProtNMA,
                                    cutoffMode=NMA_CUTOFF_ABS)
        protNMA2.inputStructure.set(protConvertVol.outputPdb)
        self.launchProtocol(protNMA2)
                                          
        # Launch NMA alignment, but just reading result from a previous metadata
        protAlignment = self.newProtocol(XmippProtAlignmentNMA,
                                         modeList='7-9',
                                         copyDeformations=self.ds.getFile('gold/pseudo_run1_images.xmd'))
        protAlignment.inputModes.set(protNMA2.outputModes)
        protAlignment.inputParticles.set(protImportParts.outputParticles)
        self.launchProtocol(protAlignment)       
        
        # Launch Dimred after NMA alignment 
        protDimRed = self.newProtocol(XmippProtDimredNMA,
                                      dimredMethod=0, # PCA
                                      reducedDim=2)
        protDimRed.inputNMA.set(protAlignment)
        self.launchProtocol(protDimRed)     
                                      
