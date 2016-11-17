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
import os

import pyworkflow.tests as tests
from pyworkflow.em.protocol import ProtImportParticles
from pyworkflow.em.packages.xmipp3 import XmippProtExtractParticles
from pyworkflow.em.packages.emxlib import ProtEmxExport
from pyworkflow.em.packages.xmipp3 import XmippProtFilterParticles, XmippProtApplyAlignment, XmippProtReconstructFourier
from pyworkflow.em.protocol.protocol_sets import ProtSplitSet,ProtUnionSet
from pyworkflow.em.convert import ImageHandler



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
        particle_even = self.downloadFile("particle_even.mrcs")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (coordinatesTest1)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=3,
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.0)
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
                                       doRemoveDust=False,
                                       doNormalize = False,
                                       doInvert=False,
                                       doFlip=False)
        protExtract.inputCoordinates.set(protEmxImport.outputCoordinates)
        protExtract.inputMicrographs.set(protEmxImport.outputMicrographs)

        self.launchProtocol(protExtract)
        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protExtract.outputParticles)
        self.launchProtocol(protEmxExport)
        stackFn = os.path.join(protEmxExport._getPath('emxData'),"data.mrc")
        self.assertTrue(ImageHandler().compareData(particle_even, stackFn, tolerance=0.01))
        #There is no way to compare emx fiels, they should bw different
        #stackEMX = os.path.join(protEmxExport._getPath('emxData'),"data.emx")
        #self.assertTrue(ImageHandler().compareData(particle_even, stackFn, tolerance=0.01))

    def test_coordinate2(self):
        """
        as test_coordinate1 but with a 129 box
        """
        #download data
        self.url = "Coordinates/Test2/"
        micFn = self.downloadFile("micrograph.mrc")
        emxFn = self.downloadFile("coordinates.emx")
        particle_odd = self.downloadFile("particle_odd.mrcs")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (coordinatesTest2)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=3,
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.0)
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
                                       doRemoveDust=False,
                                       doNormalize = False,
                                       doInvert=False,
                                       doFlip=False,
                                       doSort=False)
        protExtract.inputCoordinates.set(protEmxImport.outputCoordinates)
        protExtract.inputMicrographs.set(protEmxImport.outputMicrographs)

        self.launchProtocol(protExtract)
        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protExtract.outputParticles)
        self.launchProtocol(protEmxExport)
        stackFn = os.path.join(protEmxExport._getPath('emxData'),"data.mrc")
        # this assert does not work is a compare the stack as a block
        for num in range (1,8):
            self.assertTrue(ImageHandler().compareData("%d@"%num+particle_odd, "%d@"%num+stackFn, tolerance=0.01))

    def test_coordinate3(self):
        #download data and create set of particles
        self.url = "Coordinates/Test1/"# the 1 is OK
        micFn = self.downloadFile("micrograph.mrc")
        emxFn = self.downloadFile("coordinates.emx")
        particle_even = self.downloadFile("particle_even.mrcs")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (coordinatesTest3)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=3,
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.0)
        self.launchProtocol(protEmxImport)

        outputMics = getattr(protEmxImport, 'outputMicrographs', None)
        outputCoords = getattr(protEmxImport, 'outputCoordinates', None)

        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=128,
                                       downsampleType=0,
                                       doRemoveDust=False,
                                       doNormalize = False,
                                       doInvert=False,
                                       doFlip=False,
                                       doSort=False)
        protExtract.inputCoordinates.set(protEmxImport.outputCoordinates)
        protExtract.inputMicrographs.set(protEmxImport.outputMicrographs)

        self.launchProtocol(protExtract)
        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protExtract.outputParticles)
        self.launchProtocol(protEmxExport)

        #create output file  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        stackFn = os.path.join(protEmxExport._getPath('emxData'),"data.mrc")
        self.assertTrue(ImageHandler().compareData(particle_even, stackFn, tolerance=0.01))

    def test_ctf1(self):
        """
        This test requires:

        Download a file with two parametric CTFs in EMX exchange format.
        Download a file with two test images in CCP4 format (these images are zero valued except at the central pixel). The sampling rate of these test image is 2 A/px.
        Convert the CTFs to your package format.
        Convert the test images to your package format.
        Apply each CTFs to the corresponding test image (the resulting images are the microscope PSF). Conceptually, the steps required to apply the CTF to each test image are:
        Create a digital version of the parametric CTF
        Fourier transform the test image
        Multiply the CTF by the Fourier transformed test image
        Perform an inverse Fourier transform on the above image
        Create a 2D CCP4 stack file with the two images obtained from applying the CTF to the test images
        The sampling of the resulting images should be 2 Angstrom/pixel
        Upload the stack file.
        Three galleries of images will be displayed: the gold standard, the one just uploaded and the differences between them. The test has been successful if the gold standard and the images updated are identical.
        As extra check, the Web Site will make a pixel by pixel comparison between images belonging to both galleries. A green tick will appear if both images are identical and a red cross if any pair of pixels differ more than 10-2.

        Note: both CTFs as well as the test image have different sampling rate . The sampling rate of the output image should be 2 A/px.
        """
        #download data
        self.url = "CTF/Test1/"
        deltaParticle = self.downloadFile("delta.mrc")
        emxFn = self.downloadFile("ctf.emx")
        deltaParticleCTF = self.downloadFile("ctf.mrcs")
        #import ctf
        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (ctfTest1)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=3,#none?
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.)
        self.launchProtocol(protEmxImport)

        outputMics = getattr(protEmxImport, 'outputMicrographs', None)
        outputCTFs = getattr(protEmxImport, 'outputCTF', None)
        outputParts = getattr(protEmxImport, 'outputParticles', None)
        self.assertIsNotNone(outputMics)
        self.assertIsNotNone(outputCTFs)
        self.assertIsNotNone(outputParts)
        #split input particles
        protSplit = self.newProtocol(ProtSplitSet,
                                     numberOfSets=2,
                                     )
        protSplit.inputSet.set(outputParts)
        protSplit.randomize.set(False)
        self.proj.launchProtocol(protSplit, wait=True)
        outputParts01 = getattr(protSplit, 'outputParticles01', None)
        outputParts02 = getattr(protSplit, 'outputParticles02', None)

        #filter with CTF1
        protFilter1 = self.newProtocol(XmippProtFilterParticles)
        protFilter1.inputParticles.set(outputParts01)
        protFilter1.inputCTF.set(outputCTFs)
        protFilter1.inputCTF.setExtended(1)
        protFilter1.filterSpace.set(0)#fourier space
        protFilter1.filterModeFourier.set(3)#ctf
        self.launchProtocol(protFilter1)
        #filter with CTF2
        protFilter2 = self.newProtocol(XmippProtFilterParticles)
        protFilter2.inputParticles.set(outputParts02)
        protFilter2.inputCTF.set(outputCTFs)
        protFilter2.inputCTF.setExtended(2)
        protFilter2.filterSpace.set(0)#fourier space
        protFilter2.filterModeFourier.set(3)#ctf
        self.launchProtocol(protFilter2)
        #create a single set of particles
        protUnion = self.newProtocol(ProtUnionSet)
        protUnion.inputSets.append(protFilter1.outputParticles)
        protUnion.inputSets.append(protFilter2.outputParticles)
        self.launchProtocol(protUnion)

        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protUnion.outputSet)
        self.launchProtocol(protEmxExport)
        #TODO: upload result to emx web site. Now it is down
        stackFn = os.path.join(protEmxExport._getPath('emxData'),"data.mrc")
        self.assertTrue(ImageHandler().compareData(deltaParticleCTF, stackFn, tolerance=0.01))

    def test_ctf2(self):
        """
        This test requires:

            Download a file with one parametric CTF in EMX exchange format.
            Download a file with two test images in CCP4 format (these images are zero valued except in a single pixel). Their sampling rate is 2A/px
            Convert them to your package format.
            Apply the CTF to the corresponding image.
            Create a 2D CCP4 stack file with the two images resulting from applying the CTF to the test image
            The sampling of the resulting stack should be the same 2 Angstrom/pixel
            Upload the stack file.
            Three galleries of images will be displayed: the gold standard, the one just uploaded and the differences between them. The test has been successful if the gold standard and the images updated are identical.
            As extra check, the Web Site will make a pixel by pixel comparison between images belonging to both galleries. A green tick will appear if both images are identical and a red cross if any pair of pixels differ more than 10-2.

        Note: both CTFs as well as the test image have different sampling rate . The sampling rate of the output image should be 2 A/px.
        """
        #download data
        self.url = "CTF/Test2/"
        deltaPArticle = self.downloadFile("data.mrc")
        emxFn = self.downloadFile("ctf2.emx")
        deltaParticleCTF = self.downloadFile("ctf.mrcs")
        #import ctf
        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (ctfTest2)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=3,#none?
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.)
        self.launchProtocol(protEmxImport)

        outputMics = getattr(protEmxImport, 'outputMicrographs', None)
        outputCTFs = getattr(protEmxImport, 'outputCTF', None)
        outputParts = getattr(protEmxImport, 'outputParticles', None)
        self.assertIsNotNone(outputMics)
        self.assertIsNotNone(outputCTFs)
        self.assertIsNotNone(outputParts)
        #split input particles
        protSplit = self.newProtocol(ProtSplitSet,
                                     numberOfSets=2,
                                     )
        protSplit.inputSet.set(outputParts)
        protSplit.randomize.set(False)
        self.proj.launchProtocol(protSplit, wait=True)
        outputParts01 = getattr(protSplit, 'outputParticles01', None)
        outputParts02 = getattr(protSplit, 'outputParticles02', None)

        #filter with CTF1
        protFilter1 = self.newProtocol(XmippProtFilterParticles)
        protFilter1.inputParticles.set(outputParts01)
        #protFilter1.inputCTF.set(outputCTFs)
        #protFilter1.inputCTF.setExtended(1)
        protFilter1.filterSpace.set(0)#fourier space
        protFilter1.filterModeFourier.set(3)#ctf
        self.launchProtocol(protFilter1)
        #filter with CTF2
        protFilter2 = self.newProtocol(XmippProtFilterParticles)
        protFilter2.inputParticles.set(outputParts02)
        #protFilter2.inputCTF.set(outputCTFs)
        #ssprotFilter2.inputCTF.setExtended(2)
        protFilter2.filterSpace.set(0)#fourier space
        protFilter2.filterModeFourier.set(3)#ctf
        self.launchProtocol(protFilter2)
        #create a single set of particles
        protUnion = self.newProtocol(ProtUnionSet)
        protUnion.inputSets.append(protFilter1.outputParticles)
        protUnion.inputSets.append(protFilter2.outputParticles)
        self.launchProtocol(protUnion)

        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(protUnion.outputSet)
        self.launchProtocol(protEmxExport)
        #TODO: upload result to emx web site. Now it is down
        stackFn = os.path.join(protEmxExport._getPath('emxData'),"data.mrc")
        self.assertTrue(ImageHandler().compareData(deltaParticleCTF, stackFn, tolerance=0.01))

    def test_orientation1(self):
        """
        This test requires:
            Download a stack of 2D images (particles) and a set of transformation matrices in the EMX interchange standard.
            Convert both files to your package format.
            Apply the transformation matrix to the stack of images.
            Compute an average image of the whole stack
            Upload the average image as a 2D CCP4 image (standard exchange format).
            Two images will be displayed: the gold standard and the one just uploaded. The test has been successful if the gold standard and the image updated are identical.
            Note: test data created using a subunit of GroEL (PDB id 1SS8)
        """
        #download data
        self.url = "Orientation/Test1/"
        imgFn    = self.downloadFile("images.mrc")
        emxFn    = self.downloadFile("images.emx")
        average  = self.downloadFile("average.mrc")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (orientation1)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=0,#2D align
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.46)
        self.launchProtocol(protEmxImport)
        outputParticles = getattr(protEmxImport, 'outputParticles', None)
        #apply alignment
        protApply = self.newProtocol(XmippProtApplyAlignment)
        protApply.inputParticles.set(protEmxImport.outputParticles)
        self.launchProtocol(protApply)
        # We check that protocol generates output
        self.assertIsNotNone(protApply.outputParticles, 
                             "There was a problem generating output particles")
        # Check that output particles do not have alignment information
        self.assertFalse(protApply.outputParticles.hasAlignment(), 
                         "Output particles should not have alignment information")

        average = getattr(protApply, 'outputAverage', None)
        outputParticles = getattr(protApply, 'outputParticles', None)

        #export as emx
        protEmxExport = self.newProtocol(ProtEmxExport)
        protEmxExport.inputSet.set(outputParticles)
        self.launchProtocol(protEmxExport)
        #TODO: upload result to emx web site. Now it is down
        stackFn = os.path.join(protEmxExport._getPath('emxData'),"data.mrc")
        firstImg = outputParticles.getFirstItem()
        self.assertTrue(ImageHandler().compareData(firstImg.getFileName(), 
                                                   stackFn, tolerance=0.01))

    def test_orientation2(self):#really this is test 3
        """
        This test requires:
            Download a stack of 2D images (particles) and a set of transformation matrices in the EMX interchange standard.
            Convert both files to your package format.
            Perform a 3D reconstruction
            Upload the resulting reconstruction as a 3D CCP4 image (standard exchange format).
            Three galleries of images will be displayed: the gold standard, the one just uploaded and the differences between them. The test has been sucessful if the gold standard and the image updated are identical.
            Note: test data created using the large ribosomal subunit (PDB entry 1FFK) (origin is at volume center as described here)
        """
        #download data
        self.url = "Orientation/Test3/"
        imgFn = self.downloadFile("stack2D.mrc")
        emxFn = self.downloadFile("stack2D.emx")

        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='from emx (orientation3)',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=2,#3D proj
                                         voltage=100,
                                         magnification=10000,
                                         samplingRate=2.46)
        self.launchProtocol(protEmxImport)
        outputParticles = getattr(protEmxImport, 'outputParticles', None)
        #reconstruct outputVolume
        protReconstruct = self.newProtocol(XmippProtReconstructFourier)
        protReconstruct.inputParticles.set(protEmxImport.outputParticles)
        self.launchProtocol(protReconstruct)
        # We check that protocol generates output
        self.assertIsNotNone(protReconstruct.outputVolume, "There was a problem generating the volume")

        outputVolume = getattr(protReconstruct, 'outputVolume', None)
