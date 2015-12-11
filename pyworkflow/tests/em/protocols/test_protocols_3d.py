# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import pyworkflow.em as em
import pyworkflow.tests as pwtests
from pyworkflow.em.protocol.protocol_create_volume_set import ProtCreateVolumeSet



class TestProtCreateVolumeSet(pwtests.BaseTest):
    
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('xmipp_tutorial')

        cls.protImport1 = cls.newProtocol(em.ProtImportVolumes,
                                          objLabel='import set(1)',
                                          filesPath=cls.ds.getFile('volumes/volume_?_iter*mrc'), 
                                          samplingRate=1.0)
        cls.launchProtocol(cls.protImport1)
        
        cls.protImport2 = cls.newProtocol(em.ProtImportVolumes,
                                          objLabel='import vol(2)',
                                          filesPath=cls.ds.getFile('volumes/BPV_*64.vol'), 
                                          samplingRate=1.0)
        cls.launchProtocol(cls.protImport2)
        
        cls.protImport3 = cls.newProtocol(em.ProtImportVolumes,
                                          objLabel='import vol(3)',
                                          filesPath=cls.ds.getFile('volumes/align_test/bpv64.vol'), 
                                          samplingRate=1.0)
        cls.launchProtocol(cls.protImport3)
        
        # Run a copy of import1 just to have another set to test
        cls.protImport4 = cls.proj.copyProtocol(cls.protImport1)
        cls.protImport4.setObjLabel('import set(4)')
        cls.launchProtocol(cls.protImport4)
        
    def runCreateSet(self, *inputs):
        prot = self.newProtocol(ProtCreateVolumeSet)
        for obj in inputs:
            prot.inputVolumes.append(obj)
        self.launchProtocol(prot)
        return prot
        
    def test_individualInputs(self):
        self.runCreateSet(self.protImport2.outputVolume,
                          self.protImport3.outputVolume)
        
        
    def test_mixedInputs(self):
        self.runCreateSet(self.protImport1.outputVolumes,
                          self.protImport2.outputVolume)
        
        
    def test_allInputs(self):
        self.runCreateSet(self.protImport1.outputVolumes,
                          self.protImport2.outputVolume,
                          self.protImport3.outputVolume)

        
