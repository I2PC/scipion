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
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes, ProtSubSet
from pyworkflow.em.packages.xmipp3.protocol_volume_homogenizer import XmippProtVolumeHomogenizer
    
class TestVolumeHomogenizer(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        
    def importVolumes(self):
        """ Import several Volumes
        """                
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
        self.assertIsNotNone(prot.outputParticles.getFileName(), 
                             "There was a problem with the import")
        
        protSubset = self.newProtocol(ProtSubSet, 
                                      objLabel='subset of particles',
                                      chooseAtRandom=True,
                                      nElements=200)
        protSubset.inputFullSet.set(prot.outputParticles)        
        self.launchProtocol(protSubset)
        self.assertIsNotNone(protSubset.outputParticles.getFileName(), 
                             "There was a problem with the protSubset")        
                          
        return protSubset
         
    def test_volumeHomogenizer(self):
        """ Test protocol volume homogenizer
        """
        protImportVols = self.importVolumes()
        protImportParts = self.importFromRelionRefine3D()        
        protVolumeHomogenizer = self.newProtocol(XmippProtVolumeHomogenizer,
                                objLabel='volume homogenizer')
        protVolumeHomogenizer.referenceVolume.set(protImportVols.outputVolume)
        protVolumeHomogenizer.referenceParticles.set(protImportParts.outputParticles)
        protVolumeHomogenizer.inputVolume.set(protImportVols.outputVolume)
        protVolumeHomogenizer.inputParticles.set(protImportParts.outputParticles)
        protVolumeHomogenizer.doAlignment.set(True)  
        self.launchProtocol(protVolumeHomogenizer)