# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
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

from pyworkflow.em import ProtImportMovies
from pyworkflow.tests import *
from pyworkflow.em.packages.motioncorr import ProtMotionCorr
from pyworkflow.em.packages.motioncorr.convert import MOTIONCORR


# Some utility functions to import movies that are used
# in several tests.
class TestMotioncorrAlignMovies(BaseTest):
    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('movies')

    @classmethod
    def runImportMovies(cls, pattern, **kwargs):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or the
        # ScannedPixelSize + microscope magnification
        params = {'samplingRate': 1.14,
                  'voltage': 300,
                  'sphericalAberration': 2.7,
                  'magnification': 50000,
                  'scannedPixelSize': None,
                  'filesPattern': pattern
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

    def _checkAlignment(self, movie, goldRange, goldRoi, goldShifts):
        alignment = movie.getAlignment()
        range = alignment.getRange()
        msgRange = "Alignment range must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRange, range, msgRange % (goldRange,
                                                       type(goldRange),
                                                       range,
                                                       type(range)))
        roi = alignment.getRoi()
        shifts = alignment.getShifts()
        msgRoi = "Alignment ROI must be %s (%s) and it is %s (%s)"
        msgShifts = "Alignment SHIFTS must be %s (%s) and it is %s (%s)"
        self.assertEqual(goldRoi, roi, msgRoi % (goldRoi, type(goldRoi),
                                                 roi, type(roi)))
        self.assertEqual(goldShifts, shifts, msgShifts % (goldShifts, 
                                                          type(goldShifts),
                                                          shifts,
                                                          type(shifts)))

    def test_qbeta(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncorr test1',
                                useMotioncor2=False,)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0],
                             ([0.0, 1.7946, 2.3304, 2.5357, 2.5893, 2.125, 1.7188],
                              [0.0, 0.5848, 0.692, 0.4821, 0.2545, 0.3348, 0.1205]))

    def test_cct(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='cct - motioncorr test',
                                useMotioncor2=False,
                                doSaveMovie=True)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0], ([0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    
    def test_qbeta_SkipCrop(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncorr test2',
                                useMotioncor2=False,
                                alignFrame0=3, alignFrameN=5,
                                cropOffsetX=10, cropOffsetY=10)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)

        expected = ([0.0, 0.3646, 0.2604], [0.0, -0.1354, -0.3958])

        if "7.5" in MOTIONCORR:
            expected = ([0.0, 0.3333, 0.2292], [0.0, -0.2187, -0.4688])
        self._checkAlignment(prot.outputMovies[1],
                             (3, 5), [10, 10, 0, 0], expected)
