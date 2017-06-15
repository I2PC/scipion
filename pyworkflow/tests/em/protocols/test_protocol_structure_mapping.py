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
# *  e-mail address 'scipion@cnb.csic.es'
# ***************************************************************************/

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.em.packages.xmipp3.nma.protocol_structure_mapping import XmippProtStructureMapping

    
class TestStructureMapping(BaseTest):
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
    
    def testStructureMapping(self):
        """ 
            Test protocol structure mapping
        """             
        protStrucMap = self.newProtocol(XmippProtStructureMapping,
                                objLabel='structure mapping',
                                numberOfMpi=1, numberOfThreads=4)
        protStrucMap.inputVolumes.append(self.protImport1.outputVolume)
        protStrucMap.inputVolumes.append(self.protImport2.outputVolume)
        protStrucMap.numberOfModes.set(20)
        protStrucMap.pseudoAtomTarget.set(2.0)
        protStrucMap.rcPercentage.set(97.0)
        protStrucMap.collectivityThreshold.set(0.2)
        self.launchProtocol(protStrucMap)