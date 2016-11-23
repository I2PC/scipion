# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *


class TestMLTomoBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_programs'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.volsFn = cls.dataset.getFile('input/Ml_tomo/Subnadnorm/*.norm')
        cls.refsFn = cls.dataset.getFile('input/Ml_tomo/vir_norm.spi')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        protImport = cls.newProtocol(ProtImportVolumes, 
                                     filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport


class TestXmippMLtomo(TestMLTomoBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestMLTomoBase.setData()
        cls.protImportVol = cls.runImportVolumes(cls.volsFn, 3.5)
        cls.protImportRefVol = cls.runImportVolumes(cls.refsFn, 3.5)
    
    def testProtXmippMLTomo(self):
        prot = self.newProtocol(XmippProtMLTomo,
                                generateRefs=False,
                                numberOfIterations=3,
                                symmetry='i3',
                                missingDataType=0,
                                missingAng='-63.20 64.47',
                                maxCC=True,
                                angSampling=5,
                                dim=16,
                                maxRes=0.45,
                                numberOfMpi=2, numberOfThreads=1)
        prot.inputVols.set(self.protImportVol.outputVolumes)
        prot.inputRefVols.set(self.protImportRefVol.outputVolume)
        self.launchProtocol(prot)
        
        self.assertIsNotNone(getattr(prot, 'outputClasses', None),
                             "There was a problem with ML_tomo:\n" + prot.getErrorMessage())
