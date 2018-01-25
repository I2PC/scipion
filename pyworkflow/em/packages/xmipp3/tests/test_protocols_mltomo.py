# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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
from pyworkflow.em.packages.xmipp3 import ProtImportVolumes, XmippProtMLTomo


class TestMLTomoBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_programs'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.volsFn = cls.dataset.getFile('input/ml_tomo/*em')

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
        cls.protImportVol = cls.runImportVolumes(cls.volsFn, 15.0)

    def testProtXmippMLTomo(self):
        prot = self.newProtocol(XmippProtMLTomo,
                                numberOfReferences=2,
                                numberOfIterations=5,
                                symmetry='c1',
                                missingAng='-60 60',
                                angSampling=30,
                                doPerturb=True,
                                regIni=5.0,
                                numberOfMpi=3, numberOfThreads=1)
        prot.inputVols.set(self.protImportVol.outputVolumes)
        self.launchProtocol(prot)
        
        self.assertIsNotNone(getattr(prot, 'outputClasses', None),
                             "There was a problem when running ml_tomo:\n" + prot.getErrorMessage())
