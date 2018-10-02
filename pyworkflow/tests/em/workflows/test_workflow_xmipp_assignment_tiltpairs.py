# **************************************************************************
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

from pyworkflow.tests import *
from test_workflow import TestWorkflow
from pyworkflow.utils import importFromPlugin
from pyworkflow.em import ProtImportMicrographsTiltPairs, ProtImportMicrographs, \
    ProtImportCoordinates

XmippProtAssignmentTiltPair = importFromPlugin('xmipp3.protocols',
                                               'XmippProtAssignmentTiltPair', doRaise=True)


# update this test when RCT workflow are implemented
class TestXmippAssignmentTiltPairsWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('rct')
        cls.allCrdsDir = cls.dataset.getFile('positions')
        cls.micsUFn = cls.dataset.getFile('untilted')
        cls.micsTFn = cls.dataset.getFile('tilted')
        
    def test1(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographsTiltPairs, 
                                      patternUntilted=self.micsUFn,
                                      patternTilted=self.micsTFn,
                                      samplingRate=2.28, voltage=100,
                                      sphericalAberration=2.9)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographsTiltPair,
                             "There was a problem with the import")
        
        protImportCoorU = self.newProtocol(ProtImportCoordinates,
                         importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                         filesPath=self.allCrdsDir,
                         filesPattern='F_rct_u_*.pos', boxSize=100)
        uMics = protImport.outputMicrographsTiltPair.getUntilted()
        protImportCoorU.inputMicrographs.set(uMics)
        self.launchProtocol(protImportCoorU)

        tMics = protImport.outputMicrographsTiltPair.getTilted()
        protImportCoorT = self.newProtocol(ProtImportCoordinates,
                         importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                         filesPath=self.allCrdsDir,
                         filesPattern='F_rct_t_*.pos', boxSize=100)
        protImportCoorT.inputMicrographs.set(tMics)
        self.launchProtocol(protImportCoorT)
                
        # Then simulate a particle picking
        print "Running tilt pairs assignment..."   
        protAssigning = self.newProtocol(XmippProtAssignmentTiltPair)
        micsTiltPair = protImport.outputMicrographsTiltPair
        protAssigning.inputMicrographsTiltedPair.set(micsTiltPair)
        print self.micsUFn
        print self.micsTFn
        protAssigning.untiltedSet.set(protImportCoorU.outputCoordinates)
        protAssigning.tiltedSet.set(protImportCoorT.outputCoordinates)
        self.launchProtocol(protAssigning)
        self.assertIsNotNone(protAssigning.outputCoordinatesTiltPair,
                             "There was a problem with the protocol assignment tilt pairs")
        print '-----------------------------------------------------------'
        num_particles = protAssigning.outputCoordinatesTiltPair.getUntilted().getSize()
        print num_particles
        if (num_particles>1000):
            out_ = True
        else:
            out_ = None
        self.assertIsNotNone(out_, "There was a problem with the protocol assignment tilt pairs") 
            
        #self.validateFiles('protPicking', protPicking)  
