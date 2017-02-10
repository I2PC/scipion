# ***************************************************************************
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os
#from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles
from pyworkflow.em.packages.emxlib.dataimport import ProtEmxExport


class TestExportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsEmx = DataSet.getDataSet('emx')
        cls.dsMda = DataSet.getDataSet('mda')
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        
    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and
        the condition is valid. 
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o 
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
        
    
class TestExportParticlesEMX(TestExportBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        # Import only once for all tests
        cls.protImport = cls.runImportFromScipion()
        
    @classmethod
    def runImportFromScipion(cls):
        """ Import the particles with 3D projection directions for reconstruct
        a volume. Actually this test use a similar .sqlite file than
        the output of
        test_fromRelionRefine3D
        """
        prot1 = cls.newProtocol(ProtImportParticles,
                                 objLabel='from scipion (to-reconstruct)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=cls.dsRelion.getFile('import/case2/particles.sqlite'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        cls.launchProtocol(prot1)
        return prot1
    
    def test_exportOneStack(self):
        protExport = self.newProtocol(ProtEmxExport,
                                      objLabel='to single stack',)
        protExport.inputSet.set(self.protImport.outputParticles)
        
        self.launchProtocol(protExport)

        # Check the files were generated properly        
        self.assertTrue(os.path.exists(protExport._getPath('emxData', 'data.emx')))
        self.assertTrue(os.path.exists(protExport._getPath('emxData', 'data.mrc')))
        
    def test_exportMultipleStacks(self):
        protExport = self.newProtocol(ProtEmxExport, 
                                      objLabel='to multiple stacks',
                                      outputStack=ProtEmxExport.STACK_MICS)
        protExport.inputSet.set(self.protImport.outputParticles)
        
        self.launchProtocol(protExport)

        # Check the files were generated properly
        self.assertTrue(os.path.exists(protExport._getPath('emxData', 'data.emx')))
        for i in range(1, 21):
            self.assertTrue(os.path.exists(protExport._getPath('emxData', 'particles_%05d.mrc' % i)))
