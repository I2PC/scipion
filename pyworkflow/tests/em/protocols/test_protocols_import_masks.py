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
from pyworkflow.em.protocol import ProtImportMask


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')

    
class TestImportMasks(TestImportBase):
    
    def test_import_mask2d(self):
        """ Import a mask 2d.
        """
        args = {'maskPath': self.dsXmipp.getFile('mask2d'),
                'samplingRate': 2.1
                }

        prot = self.newProtocol(ProtImportMask, **args)
        prot.setObjLabel('import mask 2d')
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputMask, "There was a problem when importing a 2d mask.")

    def test_import_mask3d(self):
        """ Import a mask 3d.
        """
        args = {'maskPath': self.dsXmipp.getFile('mask3d'),
                'samplingRate': 2.1
                }

        prot = self.newProtocol(ProtImportMask, **args)
        prot.setObjLabel('import mask 3d')
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputMask, "There was a problem when importing a 3d mask.")

