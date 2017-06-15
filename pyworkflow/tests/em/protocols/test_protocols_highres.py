# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
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

from pyworkflow.em import exists
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.xmipp3 import XmippProtReconstructHighRes
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportParticles
from pyworkflow.em.protocol.protocol_sets import ProtSubSet
import xmipp


class TestHighres(BaseTest):
    @classmethod
    def runImportVolume(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportParticles(cls):
        """ Import Particles.
        """
        args = {'importFrom': ProtImportParticles.IMPORT_FROM_SCIPION,
                'sqliteFile': cls.particles,
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.0,
                'voltage': 200,
                'samplingRate': 0.99,
                'haveDataBeenPhaseFlipped': True
                }

        # Id's should be set increasing from 1 if ### is not in the 
        # pattern
        protImport = cls.newProtocol(ProtImportParticles, **args)
        protImport.setObjLabel('import particles')
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        # Data
        cls.dataset = DataSet.getDataSet('10010')
        cls.initialVolume = cls.dataset.getFile('initialVolume')
        cls.particles = cls.dataset.getFile('particles')

        cls.protImportVol = cls.runImportVolume(cls.initialVolume, 4.95)
        cls.protImportParts = cls.runImportParticles()

    def test(self):
        subset = self.newProtocol(ProtSubSet,
                                  inputFullSet=self.protImportParts.outputParticles,
                                  chooseAtRandom=True,
                                  nElements=400)
        self.launchProtocol(subset)

        highres = self.newProtocol(XmippProtReconstructHighRes,
                                   inputParticles=subset.outputParticles,
                                   inputVolumes=self.protImportVol.outputVolume,
                                   particleRadius=180,
                                   symmetryGroup="i1",
                                   alignmentMethod=XmippProtReconstructHighRes.AUTOMATIC_ALIGNMENT,
                                   maximumTargetResolution="15 10 7",
                                   numberOfMpi=8)
        self.launchProtocol(highres)
        self.assertIsNotNone(highres.outputParticles,
                             "There was a problem with Highres")

        fnResolution = highres._getExtraPath("Iter005/iterInfo.xmd")
        if not exists(fnResolution):
            self.assertTrue(False, fnResolution + " does not exist")
        else:
            md = xmipp.MetaData("resolution@" + fnResolution)
            R = md.getValue(xmipp.MDL_RESOLUTION_FREQREAL, md.firstObject())
            self.assertTrue(R < 8, "Resolution is not below 8A")
