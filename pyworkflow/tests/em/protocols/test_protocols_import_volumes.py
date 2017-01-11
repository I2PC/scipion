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


from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportVolumes


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsEmx = DataSet.getDataSet('emx')
        cls.dsMda = DataSet.getDataSet('mda')
        
    
class TestImportVolumes(TestImportBase):
    
    def test_pattern(self):
        """ Import several Particles from a given pattern.
        """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_*mrc',
                'samplingRate': 2.1
                }
        
        
        # Id's should be set increasing from 1 if ### is not in the 
        # pattern
        prot1 = self.newProtocol(ProtImportVolumes, **args)
        prot1.setObjLabel('import mrc')
        self.launchProtocol(prot1)
        
        # Id's should be taken from filename   
        args['filesPath'] = self.dsRelion.getFile('import/case2/relion_volumes.mrc') 
        args['filesPattern'] = ''
        prot2 = self.newProtocol(ProtImportVolumes, **args)
        prot2.setObjLabel('from mrc stack')
        self.launchProtocol(prot2)
        # Check the number of output volumes and dimensions
        self.assertEqual(3, prot2.outputVolumes.getSize())
        self.assertEqual(60, prot2.outputVolumes.getDim()[0])

        # Id's should be taken from filename   
        args['filesPath'] = self.dsRelion.getFile('import/case2/relion_volumes.stk') 
        args['filesPattern'] = ''
        prot3 = self.newProtocol(ProtImportVolumes, **args)
        prot3.setObjLabel('from spider stack')
        self.launchProtocol(prot3)
        # Check the number of output volumes and dimensions
        self.assertEqual(3, prot3.outputVolumes.getSize())
        self.assertEqual(60, prot3.outputVolumes.getDim()[0])