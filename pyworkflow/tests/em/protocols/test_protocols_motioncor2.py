# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *             Grigory Sharov (sharov@igbmc.fr)
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

from pyworkflow.em import ProtImportMovies, ImageHandler
from pyworkflow.tests import *
from pyworkflow.em.packages.motioncorr import ProtMotionCorr


# Some utility functions to import movies that are used
# in several tests.
class TestMotioncor2AlignMovies(BaseTest):

    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize
        #  + microscope magnification
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern,
                  'dosePerFrame': 1.3
                  }
        if 'samplingRate' not in kwargs:
            del params['samplingRate']
            params['samplingRateMode'] = 0
        else:
            params['samplingRateMode'] = 1

        params.update(kwargs)

        protImport = cls.newProtocol(ProtImportMovies, **params)
        cls.launchProtocol(protImport)
        return protImport
    
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        cls.protImport1 = cls.runImportMovies(cls.ds.getFile('qbeta/qbeta.mrc'),
                                              magnification=50000)
        cls.protImport2 = cls.runImportMovies(cls.ds.getFile('cct/cct_1.em'),
                                              magnification=61000)

    def _checkMicrographs(self, protocol):
        self.assertIsNotNone(getattr(protocol, 'outputMicrographs', None),
                             "Output SetOfMicrographs were not created.")

    def _checkAlignment(self, movie, goldRange, goldRoi):
        alignment = movie.getAlignment()
        range = alignment.getRange()
        aliFrames = range[1] - range[0] + 1
        msgRange = "Alignment range must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRange, range, msgRange % (goldRange,
                                                       type(goldRange),
                                                       range,
                                                       type(range)))
        roi = alignment.getRoi()
        shifts = alignment.getShifts()
        zeroShifts = (aliFrames * [0], aliFrames * [0])
        nrShifts = len(shifts[0])
        msgRoi = "Alignment ROI must be %s (%s) and it is %s (%s)"
        msgShifts = "Alignment SHIFTS must be non-zero!"
        self.assertEqual(goldRoi, roi, msgRoi % (goldRoi, type(goldRoi),
                                                 roi, type(roi)))
        self.assertNotEqual(zeroShifts, shifts, msgShifts)
        self.assertEqual(nrShifts, aliFrames, "Number of shifts is not equal"
                                              " number of aligned frames.")

    def test_cct_motioncor2_patch(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='cct - motioncor2 test1 (patch-based)',
                                useMotioncor2=True,
                                patchX=2, patchY=2)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0])

    def test_qbeta_motioncor2_patch(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncor2 test2 (patch-based)',
                                useMotioncor2=True,
                                patchX=2, patchY=2,
                                group=2)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0])

    def test_qbeta_motioncor2_sel(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncor2 test3 (frame range)',
                                useMotioncor2=True,
                                patchX=2, patchY=2,
                                alignFrame0=2,
                                alignFrameN=6,
                                sumFrame0=2,
                                sumFrameN=6)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (2, 6), [0, 0, 0, 0])


# Some utility functions to import movies that are used
# in several tests.
class TestMotioncor2GainFile(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the ScannedPixelSize
        #  + microscope magnification
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern,
                  'dosePerFrame': 1.3,
                  'gainFile': None,
                  'darkFile': None,
                  }
        if 'samplingRate' not in kwargs:
            del params['samplingRate']
            params['samplingRateMode'] = 0
        else:
            params['samplingRateMode'] = 1

        params.update(kwargs)

        protImport = cls.newProtocol(ProtImportMovies, **params)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        mrcGainFn = cls.ds.getFile("mrcGain")
        spiGainFn = cls.ds.getFile("spiGain")
        mrcDarkFn = cls.ds.getFile("mrcDark")
        spiDarkFn = cls.ds.getFile("spiDark")

        cls.protImportNoGainNoDark = cls.runImportMovies(
                                            cls.ds.getFile('falcon2012_1'),
                                            magnification=45000,
                                            objLabel="Import No gain or dark")
        cls.protImportMrcGain = cls.runImportMovies(
                                            cls.ds.getFile('falcon2012_1'),
                                            magnification=45000,
                                            gainFile=mrcGainFn,
                                            darkFile=mrcDarkFn,
                                            objLabel="Import MRC gain and dark")

        cls.protImportSpiGain = cls.runImportMovies(
                                            cls.ds.getFile('falcon2012_1'),
                                            magnification=45000,
                                            gainFile=spiGainFn,
                                            darkFile=spiDarkFn,
                                            objLabel='Import SPI gain and dark')

    def _checkMicrographs(self, protocol):
        self.assertSetSize(protocol.outputMicrographs, size=1)

    def testGainAndDarkCombinations(self):
        protNoGain = self.newProtocol(ProtMotionCorr,
                                      objLabel='NoGainFile - motioncor2',
                                      useMotioncor2=True,
                                      patchX=2, patchY=2,
                                      group=2)

        protNoGain.inputMovies.set(self.protImportNoGainNoDark.outputMovies)
        self.launchProtocol(protNoGain)

        # Check the is a step for converting the gain and dark
        steps = protNoGain.loadSteps()

        # Check new step is there
        gainStep = steps[0]
        self.assertEqual(gainStep.funcName, '_convertInputStep')

        # Check movie steps depends on it
        movieStep = steps[1]
        self.assertIn(str(gainStep.getObjId()), movieStep._prerequisites)

        self._checkMicrographs(protNoGain)

        # Gain and dark in MRC
        protMRC = self.newProtocol(ProtMotionCorr,
                                   objLabel='MrcGainFile - motioncor2',
                                   useMotioncor2=True,
                                   patchX=2, patchY=2)

        protMRC.inputMovies.set(self.protImportMrcGain.outputMovies)
        self.launchProtocol(protMRC)

        gain = self.protImportMrcGain.outputMovies.getGain()
        finalGain = protMRC.getFinalCorrectionImagePath(gain)
        self.assertTrue(os.path.exists(finalGain))
        self.assertEquals(gain, finalGain)

        self._checkMicrographs(protMRC)

        # SPI gain and dark format
        protSPI = self.newProtocol(ProtMotionCorr,
                                   objLabel='SpiGainFile - motioncor2',
                                   useMotioncor2=True,
                                   patchX=2, patchY=2,
                                   group=2)

        protSPI.inputMovies.set(self.protImportSpiGain.outputMovies)

        self.launchProtocol(protSPI)

        self._checkMicrographs(protSPI)

        gain = self.protImportSpiGain.outputMovies.getGain()
        finalGain = protSPI.getFinalCorrectionImagePath(gain)
        self.assertTrue(os.path.exists(finalGain))
        self.assertTrue('.spi' not in finalGain)

        # Get first micrograph from MC
        mic = protSPI.outputMicrographs.getFirstItem()
        micFn = mic.getFileName()

        # Get the mic without gain
        micNoGain = protNoGain.outputMicrographs.getFirstItem()
        micNoGainFn = micNoGain.getFileName()

        # Check images are different
        self.assertFalse(
            ImageHandler().compareData(micFn, micNoGainFn, tolerance=0.00),
            msg="mic with SPI gain and without it are not different.")
