# ***************************************************************************
# * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************/

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes
from pyworkflow.em.packages.xmipp3.protocol_apply_transformation_matrix import XmippProtApplyTransformationMatrix
from pyworkflow.em.packages.xmipp3.protocol_align_volume import XmippProtAlignVolume, ALIGN_ALGORITHM_LOCAL
    
class TestApplyTransformationMatrix(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.protImport1 = cls.newProtocol(ProtImportVolumes,
                                         filesPath=cls.dsRelion.getFile('volumes/reference_rotated.vol'), 
                                         samplingRate=1.0)
        cls.launchProtocol(cls.protImport1)        
        cls.protImport2 = cls.newProtocol(ProtImportVolumes,
                                         filesPath=cls.dsRelion.getFile('volumes/reference.mrc'), 
                                         samplingRate=1.0)
        cls.launchProtocol(cls.protImport2)
    
    def importFromRelionRefine3D(self):
        """ Import aligned Particles
        """
        prot = self.newProtocol(ProtImportParticles,
                                 objLabel='particles from relion (auto-refine 3d)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=self.dsRelion.getFile('import/classify3d/extra/relion_it015_data.star'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot)        
        self.assertIsNotNone(prot.outputParticles.getFileName(), 
                             "There was a problem with the import")
        return prot
                 
    def testAlignVolumeLocal(self):
        protAlign = self.newProtocol(XmippProtAlignVolume,
                                     inputReference=self.protImport1.outputVolume,
                                     alignmentAlgorithm=ALIGN_ALGORITHM_LOCAL,
                                     initialRotAngle=0,
                                     initialTiltAngle=0,
                                     initialInplaneAngle=0,
                                     initialShiftX=0,
                                     initialShiftY=0,
                                     initialShiftZ=0,
                                     initialScale=1,
                                     optimizeScale=True,
                                     numberOfMpi=1, numberOfThreads=1                               
                                     )
        protAlign.inputVolumes.append(self.protImport2.outputVolume)
        self.launchProtocol(protAlign)    
        return protAlign
    
    def testApplyTransformationMatrix(self):
        """ 
            Test protocol apply transformation matrix
        """
        protImportVols = self.testAlignVolumeLocal()
        protImportParts = self.importFromRelionRefine3D()        
        protApplyTransformationMatrix = self.newProtocol(XmippProtApplyTransformationMatrix,
                                objLabel='apply transformation matrix')
        protApplyTransformationMatrix.inputParticles.set(protImportParts.outputParticles)
        protApplyTransformationMatrix.inputVolume.set(protImportVols.outputVolume)
        self.launchProtocol(protApplyTransformationMatrix)