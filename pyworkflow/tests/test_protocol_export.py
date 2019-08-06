# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
""" Module to test json exporting options"""
from pyworkflow.em import Pointer, ProtMultiPointerTest
from pyworkflow.object import PointerList, String
from tests import *

MY_OUTPUT = "myOutput"


class TestProtocolExport(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_exportMultiPointerToSets(self):
        """Test how multipointers to are exported and imported"""

        protMpOutput = self.newProtocol(ProtMultiPointerTest,
                                objLabel='multipointer for output')

        protMpOutput._defineOutputs(**{MY_OUTPUT: String("hola!")})
        self.saveProtocol(protMpOutput)

        # Add multiPointers with extended
        plWithExtended = PointerList()
        plWithExtended.append(Pointer(protMpOutput, extended=MY_OUTPUT))

        protMp = self.newProtocol(ProtMultiPointerTest,
                                  objLabel='multipointer with extended export to json',
                                  mpToAttr=plWithExtended)

        self.saveProtocol(protMp)
        # Trigger the refresh of the runsGraph!!
        self.proj._runsGraph = None

        protDict = self.proj.getProtocolsDict()

        # Get the multipointer params items for the second prot
        # Get the second prot
        # This has the shape of :  (33, OrderedDict([('object.className', 'ProtMultiPointerTest'), ...) ])
        ndProtAttrs = protDict[protMp.getObjId()]

        # Look for the mpToAttr
        for key, value in ndProtAttrs.items():

            if key == "mpToAttr":
                self.assertEqual(1, len(value),
                                 "multipointer param items exceeds the "
                                 "expected number of items")
                self.assertEqual("%s.%s" % (protMpOutput.getObjId(), MY_OUTPUT),
                                 value[0],
                                 "Multipointer item value %s seems to be wrong."
                                 % value[0])
