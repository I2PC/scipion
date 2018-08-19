# ***************************************************************************
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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
from pyworkflow.em.protocol import ProtImportCoordinates, ProtImportMicrographs



class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        #cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsGroel = DataSet.getDataSet('groel')
        
    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and
        the condition is valid. 
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o 
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
        
    
class TestImportCoordinates(TestImportBase):

    def testImportCoordinates(self):
        #First, import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs, filesPath=self.dsXmipp.getFile('allMics'), samplingRate=1.237, voltage=300)
        protImport.setObjLabel('import micrographs from xmipp tutorial ')
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(), "There was a problem with the import")

        prot1 = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                 filesPath=self.dsXmipp.getFile('pickingXmipp'),
                                 filesPattern='*.pos', boxSize=550,
                                 scale=3.,
                                 invertX=False,
                                 invertY=False
                                 )
        prot1.inputMicrographs.set(protImport.outputMicrographs)
        prot1.setObjLabel('import coords from xmipp ')
        self.launchProtocol(prot1)
        
        # Make sure that all 264 coordinates where correctly imported
        self.assertTrue(prot1.outputCoordinates.getSize()==264)

        # prot2 = self.newProtocol(ProtImportCoordinates,
        #                          importFrom=ProtImportCoordinates.IMPORT_FROM_RELION,
        #                          filesPath=self.dsXmipp.getFile('boxingDir'),#no dataset with picking
        #                          filesPattern='info/*_info.json',
        #                          boxSize=110,
        #                          scale=2,
        #                          invertX=False,
        #                          invertY=False
        #                          )
        # prot2.inputMicrographs.set(protImport.outputMicrographs)
        # prot2.setObjLabel('import coords from relion ')

        prot3 = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_EMAN,
                                 filesPath=self.dsXmipp.getFile('boxingDir'),
                                 filesPattern='*_info.json', boxSize=550,
                                 scale=5.,
                                 invertX=False,
                                 invertY=False)
        prot3.inputMicrographs.set(protImport.outputMicrographs)
        prot3.setObjLabel('import coords from eman ')

        self.launchProtocol(prot3)

        #First, import a set of micrographs
        protImportGroel = self.newProtocol(ProtImportMicrographs, filesPath=self.dsGroel.getFile('mic1'), samplingRate=1)
        protImportGroel.setObjLabel('import micrographs from groel')
        self.launchProtocol(protImportGroel)
        self.assertIsNotNone(protImportGroel.outputMicrographs.getFileName(), "There was a problem with the import")

        protPickGroel = self.newProtocol(ProtImportCoordinates,
                                 importFrom=ProtImportCoordinates.IMPORT_FROM_DOGPICKER,
                                 filesPath=self.dsGroel.getFile('pickingDogpicker'),
                                 filesPattern='*.txt', boxSize=10,
                                 threshold=0.7)
        protPickGroel.inputMicrographs.set(protImportGroel.outputMicrographs)
        protPickGroel.setObjLabel('import coords from dogpicker ')

        self.launchProtocol(protPickGroel)

