# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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


from pyworkflow.utils import exists
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.packages.cryoEF import ProtCryoEF
from pyworkflow.em.protocol import ProtImportParticles


class TestCryoEFBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='relion_tutorial'):
        cls.ds = DataSet.getDataSet(dataProject)
        cls.partFn = cls.ds.getFile('import/refine3d/extra/relion_it025_data.star')

    @classmethod
    def runImportParticlesStar(cls, partStar, mag, samplingRate):
        """ Import particles from Relion star file. """
        protImport = cls.newProtocol(ProtImportParticles,
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=partStar,
                                     magnification=mag,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=True)
        cls.launchProtocol(protImport)
        return protImport


class TestCryoEF(TestCryoEFBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestCryoEFBase.setData()
        cls.protImportParts = cls.runImportParticlesStar(cls.partFn, 50000, 7.08)

    def test_cryoEF(self):
        protFsc = self.newProtocol(ProtCryoEF,
                                   inputParticles=self.protImportParts.outputParticles,
                                   diam=300)
        self.launchProtocol(protFsc)
        protFsc._initialize()
        self.assertTrue(exists(protFsc._getFileName('real space PSF')),
                        "cryoEF has failed")
