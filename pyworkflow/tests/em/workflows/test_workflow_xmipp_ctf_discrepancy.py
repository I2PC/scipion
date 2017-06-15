"""
/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'scipion@cnb.csic.es'
 ***************************************************************************/
MODIFICATION ADVICE:

Please,  do not  generate or  distribute 
a modified version of this file under its original name. 
"""

import os
from itertools import izip

from pyworkflow.em.data import SetOfCTF, CTFModel, Micrograph, SetOfMicrographs
from pyworkflow.em.packages.xmipp3 import XmippProtCTFDiscrepancy
from pyworkflow.em.protocol import ProtImportMicrographs, ProtImportCTF
from pyworkflow.object import PointerList, Pointer
from test_workflow import TestWorkflow
import pyworkflow.tests as tests
from pyworkflow.em import ImageHandler
import xmipp

class TestXmippCTFDiscrepancyBase(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        #cls.dataset = tests.DataSet.getDataSet('CTFDiscrepancy')
        
    
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def testCtfdiscrepancyWorkflow(self):
        #create one micrograph set
        fnMicSet = self.proj.getTmpPath("mics.sqlite")
        fnMic    = self.proj.getTmpPath("mic.mrc")
        mic = Micrograph()
        mic.setFileName(fnMic)
        micSet = SetOfMicrographs(filename=fnMicSet)

        # create two CTFsets
        fnCTF1 = self.proj.getTmpPath("ctf1.sqlite")
        fnCTF2 = self.proj.getTmpPath("ctf2.sqlite")
        ctfSet1 = SetOfCTF(filename=fnCTF1)
        ctfSet2 = SetOfCTF(filename=fnCTF2)

        ###create one fake micrographs image
        projSize = 32
        img = xmipp.Image()
        img.setDataType(xmipp.DT_FLOAT)
        img.resize(projSize, projSize)
        img.write(fnMic)

        #fill the sets
        for i in range(1,4):
            mic = Micrograph()
            mic.setFileName(fnMic)
            micSet.append(mic)

            defocusU = 1000+10*i
            defocusV = 1000+i
            defocusAngle = i*10
            psdFile="psd_1%04d"%i
            ctf = self._getCTFModel(defocusU,
                                    defocusV,
                                    defocusAngle,
                                    psdFile)
            ctf.setMicrograph(mic)
            ctfSet1.append(ctf)

            defocusU = 1000+20*i
            defocusV = 1000+i
            defocusAngle = i*20
            psdFile="psd_2%04d"%i
            ctf = self._getCTFModel(defocusU,
                                    defocusV,
                                    defocusAngle,
                                    psdFile)
            ctf.setMicrograph(mic)
            ctfSet2.append(ctf)
        ctfSet1.write()
        ctfSet2.write()
        micSet.write()

        
        #import micrograph set
        args = {'importFrom': ProtImportMicrographs.IMPORT_FROM_SCIPION,
                'sqliteFile': fnMicSet,
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 100,
                'samplingRate': 2.1
                }

        protMicImport = self.newProtocol(ProtImportMicrographs, **args)
        protMicImport.setObjLabel('import micrographs from sqlite ')
        self.launchProtocol(protMicImport)
        
        #import ctfsets
        protCTF1 = self.newProtocol(ProtImportCTF,
                                 importFrom=ProtImportCTF.IMPORT_FROM_SCIPION,
                                 filesPath=fnCTF1)
        protCTF2 = self.newProtocol(ProtImportCTF,
                                 importFrom=ProtImportCTF.IMPORT_FROM_SCIPION,
                                 filesPath=fnCTF2)
        protCTF1.inputMicrographs.set(protMicImport.outputMicrographs)
        protCTF2.inputMicrographs.set(protMicImport.outputMicrographs)
        protCTF1.setObjLabel('import ctfs from scipion_1 ')
        protCTF2.setObjLabel('import ctfs from scipion_2 ')
        self.launchProtocol(protCTF1)
        self.launchProtocol(protCTF2)

        # launch CTF discrepancy protocol
        protCtfDiscrepancy = self.newProtocol(XmippProtCTFDiscrepancy)
        protCtfDiscrepancy.inputCTF1.set(protCTF1.outputCTF)
        protCtfDiscrepancy.inputCTF2.set(protCTF2.outputCTF)
        protCtfDiscrepancy.setObjLabel('ctf discrepancy')
        self.launchProtocol(protCtfDiscrepancy)
        ctf0 = protCtfDiscrepancy.outputCTF.getFirstItem()
        resolution = int(ctf0._discrepancy_resolution.get())
        defocusU = int(ctf0.getDefocusU())
        self.assertEqual(resolution, 2)
        self.assertEqual(defocusU, 1010)

