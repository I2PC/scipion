# ***************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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
# ***************************************************************************/

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.protocol import ProtImportMicrographs, ProtCreateStreamData
from pyworkflow.utils import importFromPlugin

ProtCTFFind = importFromPlugin('grigoriefflab.protocols', 'ProtCTFFind', doRaise=True)
XmippProtCTFMicrographs = importFromPlugin('xmipp3.protocols',
                                           'XmippProtCTFMicrographs', doRaise=True)
XmippProtCTFConsensus = importFromPlugin('xmipp3.protocols',
                                         'XmippProtCTFConsensus', doRaise=True)


class TestCtfConsensus(BaseTest):
    """ Check if the Xmipp-CTFconsensus rejects CTFs (and the coorrespondig mics)
        when two CTF estimations give different results,
        and accept when the two estimations give similar results.
    """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.micsFn = cls.dataset.getFile('allMics')

    def checkCTFs(self, protConsensus, refCTFs, refMics, label=''):
        outputCTF = getattr(protConsensus, "outputCTF"+label)
        outputMicrographs = getattr(protConsensus, "outputMicrographs"+label)

        self.assertIsNotNone(outputCTF,
                             "There was a problem with the CTF-Consensus. "
                             "No outputCTF is created.")
        self.assertIsNotNone(outputMicrographs,
                             "There was a problem with the CTF-Consensus. "
                             "No outputMicrographs is created.")

        self.assertEqual(outputCTF.getSize(), refCTFs.getSize(),
                         "The outputCTF size is wrong.")
        self.assertEqual(outputMicrographs.getSize(), refCTFs.getSize(),
                         "The outputMicrographs size is wrong.")
        self.assertTupleEqual(outputMicrographs.getDim(), refMics.getDim(),
                              "The outputMicrographs dimension is wrong.")

        firstCTF = outputCTF.getFirstItem()
        refCTF = refCTFs.getFirstItem()
        self.assertTrue(firstCTF.equalAttributes(refCTF),
                        "The outputCTF has different attributes than the input.")

    def test1(self):
        # Import a set of micrographs
        protImport = self.newProtocol(ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")

        # Create pure noise micrographs (to force a Discarded consensus set)
        protStream = self.newProtocol(ProtCreateStreamData,
                                      xDim=9216,
                                      yDim=9441,
                                      nDim=3,
                                      samplingRate=1.237,
                                      setof=2,  # 2 -> SetOfRandomMicrographs
                                      creationInterval=1)
        self.proj.launchProtocol(protStream, wait=False)

        # Computes the CTF with Xmipp
        protCTF1 = self.newProtocol(XmippProtCTFMicrographs)
        protCTF1.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF1, wait=False)

        # Computes the CTF with CTFFind4
        protCTF2 = self.newProtocol(ProtCTFFind)
        protCTF2.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF2, wait=False)

        self._waitOutput(protStream, "outputMicrographs")
        # Computes the CTF with CTFFind4 for the noise mics
        protCTF3 = self.newProtocol(ProtCTFFind)
        protCTF3.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF3, wait=False)


        # Computes the Consensus of GOOD CTFs
        self._waitOutput(protCTF1, "outputCTF")
        self._waitOutput(protCTF2, "outputCTF")
        protCTFcons = self.newProtocol(XmippProtCTFConsensus,
                                       useDefocus=False,
                                       useAstigmatism=False,
                                       useResolution=False,
                                       calculateConsensus=True)
        protCTFcons.inputCTF.set(protCTF1.outputCTF)
        protCTFcons.inputCTF2.set(protCTF2.outputCTF)
        self.launchProtocol(protCTFcons)

        protCTF1.outputCTF.load()  # Needed to update the setSize
        self.checkCTFs(protCTFcons,
                       refMics=protImport.outputMicrographs,
                       refCTFs=protCTF1.outputCTF)


        # Computes the Consensus comparing a good CTF to a RANDOM one
        self._waitOutput(protCTF3, "outputCTF")
        protCTFcons2 = self.newProtocol(XmippProtCTFConsensus,
                                       useDefocus=False,
                                       useAstigmatism=False,
                                       useResolution=False,
                                       calculateConsensus=True)
        protCTFcons2.inputCTF.set(protCTF1.outputCTF)
        protCTFcons2.inputCTF2.set(protCTF3.outputCTF)
        self.launchProtocol(protCTFcons2)
        self.checkCTFs(protCTFcons2,
                       refMics=protImport.outputMicrographs,
                       refCTFs=protCTF1.outputCTF,
                       label="Discarded")
