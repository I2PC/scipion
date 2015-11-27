# ***************************************************************************
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from itertools import izip
import urllib

import pyworkflow.tests as tests
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtImportMicrographs, ProtImportParticles
from pyworkflow.em.packages.xmipp3 import XmippProtExtractParticles
from pyworkflow.em.packages.emxlib import ProtEmxExport


class TestEmxWeb(tests.BaseTest):
    """test emx web page
    scipion test tests.em.workflows.test_workflow_emx
    """
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.baseUrl = "http://i2pc.cnb.csic.es/emx/resourcesEmx/Tests/"

    def downloadFile(self, filename):
        outFile  = self.proj.getTmpPath(filename)
        url = self.baseUrl + self.url
        urllib.urlretrieve(url + filename, outFile)
        return outFile

    def test_coordinate1(self):
        """
        Download a micrograph and a set of coordinates in the EMX interchange standard.
        Convert both files to your package format.
        Extract the particles from the micrograph with a box size of 128 pixels. (Remember that coordinates refer to particle centers.)
        Upload the extracted images to the web as a 2D CCP4 stack (standard exchange format).
        Three galleries of images will be displayed: the gold standard, the one just uploaded and the differences between them. The test has been sucessful if the gold standard and the images updated are identical.
        As extra check, the Web Site will make a pixel by pixel comparison between images belonging to both galleries. A green tick will appear if both images are identical and a red cross if any pair of pixels differ more than 10**-2.
        """
        #download data
        self.url = "Coordinates/Test1/"
        micFn = self.downloadFile("micrograph.mrc")
        emxFn = self.downloadFile("coordinates.emx")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                           objLabel='from emx (coordinatesTest1)',
                                           importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                           emxFile=emxFn,
                                           alignType=3,
                                           voltage=100,
                                           magnification=10000,
                                           samplingRate=2.46)
        self.launchProtocol(protEmxImport)

        outputMics = getattr(protEmxImport, 'outputMicrographs', None)
        outputCoords = getattr(protEmxImport, 'outputCoordinates', None)

        self.assertIsNotNone(outputMics)
        self.assertIsNotNone(outputCoords)
        coordsGolds = [(539, 414), (509, 192), (711, 158), (634, 349),
                       (403, 157), (347, 437), (728, 389)]
        for coord, coordGold in izip(outputCoords, coordsGolds):
            self.assertEquals(coord.getPosition(), coordGold)

        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=128,
                                       downsampleType=0,
                                       doFlip=False,
                                       downFactor=1)
        protExtract.inputCoordinates.set(protEmxImport.outputCoordinates)
        protExtract.inputMicrographs.set(protEmxImport.outputMicrographs)

        self.launchProtocol(protExtract)
        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protExtract.outputParticles)
        self.launchProtocol(protEmxExport)
        #TODO: upload result to emx web site. Now it is down

    def test_coordinate2(self):
        """
        as test_coordinate1 but with a 129 box
        """
        #download data
        self.url = "Coordinates/Test2/"
        micFn = self.downloadFile("micrograph.mrc")
        emxFn = self.downloadFile("coordinates.emx")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                           objLabel='from emx (coordinatesTest2)',
                                           importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                           emxFile=emxFn,
                                           alignType=3,
                                           voltage=100,
                                           magnification=10000,
                                           samplingRate=2.46)
        self.launchProtocol(protEmxImport)

        outputMics = getattr(protEmxImport, 'outputMicrographs', None)
        outputCoords = getattr(protEmxImport, 'outputCoordinates', None)

        self.assertIsNotNone(outputMics)
        self.assertIsNotNone(outputCoords)
        coordsGolds = [(539, 414), (509, 192), (711, 158), (634, 349),
                       (403, 157), (347, 437), (728, 389)]
        for coord, coordGold in izip(outputCoords, coordsGolds):
            self.assertEquals(coord.getPosition(), coordGold)

        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=129,
                                       downsampleType=0,
                                       doFlip=False,
                                       downFactor=1)
        protExtract.inputCoordinates.set(protEmxImport.outputCoordinates)
        protExtract.inputMicrographs.set(protEmxImport.outputMicrographs)

        self.launchProtocol(protExtract)
        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protExtract.outputParticles)
        self.launchProtocol(protEmxExport)
        #TODO: upload result to emx web site. Now it is down

    def test_coordinate3(self):
        #download data and create set of particles
        self.url = "Coordinates/Test1/"# the 1 is OK
        micFn = self.downloadFile("micrograph.mrc")
        emxFn = self.downloadFile("coordinates.emx")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                           objLabel='from emx (coordinatesTest3)',
                                           importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                           emxFile=emxFn,
                                           alignType=3,
                                           voltage=100,
                                           magnification=10000,
                                           samplingRate=2.46)
        self.launchProtocol(protEmxImport)

        outputMics = getattr(protEmxImport, 'outputMicrographs', None)
        outputCoords = getattr(protEmxImport, 'outputCoordinates', None)

        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=129,
                                       downsampleType=0,
                                       doFlip=False,
                                       downFactor=1)
        protExtract.inputCoordinates.set(protEmxImport.outputCoordinates)
        protExtract.inputMicrographs.set(protEmxImport.outputMicrographs)

        self.launchProtocol(protExtract)
        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protExtract.outputParticles)
        self.launchProtocol(protEmxExport)

        #create output file
        #TODO: upload result to emx web site. Now it is down