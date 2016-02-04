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

#import os
#from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes
from pyworkflow.em.packages.xmipp3.protocol_validate_overfitting import XmippProtValidateOverfitting

        
    
class TestValidateOverfitting(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        
    def importVolumes(self):
        """ Import several Volumes
        """
        '''args = {'filesPath': self.dsRelion.getFile('import/case2/relion_volumes.mrc'),
                'filesPattern': '',
                'samplingRate': 7.08
                }
        
        prot = self.newProtocol(ProtImportVolumes, **args)
        prot.setObjLabel('volumes from mrc stack')
        self.launchProtocol(prot)
        # Check the number of output volumes and dimensions
        self.assertEqual(60, prot.outputVolumes.getDim()[0])
        self.assertIsNotNone(prot.outputVolumes.getFileName(), "There was a problem with the import")
        '''
        
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': '*_64.vol',
                'samplingRate': 1
                }
        
        prot = self.newProtocol(ProtImportVolumes, **args)
        prot.setObjLabel('import volume')
        self.launchProtocol(prot)
        
        return prot
    
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
        
        self.assertIsNotNone(prot.outputParticles.getFileName(), "There was a problem with the import")
        return prot
         
    def test_validateOverfitting(self):
        """ Test protocol validate overfitting
        """
        protImportVols = self.importVolumes()
        protImportPars = self.importFromRelionRefine3D()
        
        protValidateOverfitting = self.newProtocol(XmippProtValidateOverfitting,
                                objLabel='validate overfitting',
                                numberOfMpi=1, numberOfThreads=4)
        protValidateOverfitting.inputParticles.set(protImportPars.outputParticles)
        protValidateOverfitting.input3DReference.set(protImportVols.outputVolume)
        protValidateOverfitting.numberOfParticles.set("10 20 50")
        protValidateOverfitting.numberOfIterations.set(3)
        protValidateOverfitting.symmetryGroup.set('c1')
        
        self.launchProtocol(protValidateOverfitting)