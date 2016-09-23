# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *             Josue Gomez Blanco (jgomez@cnb.csic.es)
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

import unittest, sys
from pyworkflow.em import ProtImportMovies
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
        # We have two options: passe the SamplingRate or the ScannedPixelSize + microscope magnification
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

    def test_cct_motioncor2_patch(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='cct - motioncor2 test1 (patch-based)',
                                useMotioncor2=True,
                                patch='2 2',
                                frameDose=1.3)
        prot.inputMovies.set(self.protImport2.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0], ([0.0, -1.31, -2.14, -2.13, -1.87, -1.53, -1.07],
                                                    [0.0, -1.32, -2.22, -2.47, -2.61, -2.6, -2.52]))

    def test_qbeta_motioncor2_patch(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncor2 test2 (patch-based)',
                                useMotioncor2=True,
                                patch='2 2',
                                group=2)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (1, 7), [0, 0, 0, 0], ([0.0, 0.85, 1.7, 2.07, 1.95, 1.84, 1.72],
                                                    [0.0, 0.22, 0.43, 0.46, 0.3, 0.15, -0.01]))

    def test_qbeta_motioncor2_sel(self):
        prot = self.newProtocol(ProtMotionCorr,
                                objLabel='qbeta - motioncor2 test3 (frame range)',
                                useMotioncor2=True,
                                patch='2 2',
                                alignFrame0=2,
                                alignFrameN=6,
                                sumFrame0=2,
                                sumFrameN=6)
        prot.inputMovies.set(self.protImport1.outputMovies)
        self.launchProtocol(prot)

        self._checkMicrographs(prot)
        self._checkAlignment(prot.outputMovies[1],
                             (2, 6), [0, 0, 0, 0], ([0.0, 0.58, 0.73, 0.74, 0.3],
                                                    [0.0, 0.14, -0.1, -0.35, -0.15]))
