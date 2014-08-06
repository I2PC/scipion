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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
MODIFICATION ADVICE:

Please,  do not  generate or  distribute 
a modified version of this file under its original name. 
"""

import unittest
from pyworkflow.em.packages.emxlib import ProtEmxImport
from pyworkflow.em.packages.xmipp3 import XmippProtCTFDiscrepancy
import pyworkflow.em as em
import pyworkflow.tests as tests
import os
import itertools
from test_workflow import TestWorkflow

class TestXmippCTFDiscrepancyBase(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = tests.DataSet.getDataSet('CTFDiscrepancy')

    def testCtfdiscrepancyWorkflow(self):
        """ Import  3 EMX files with micrographs and defocus and compare them
        """
        emxFn1 = self.dataset.getFile('emxMicrographCtf1')
        emxFn2 = self.dataset.getFile('emxMicrographCtf2')
        emxFn3 = self.dataset.getFile('emxMicrographCtf3')
        protEmxImport1 = self.newProtocol(ProtEmxImport,
                                          inputEMX=emxFn1
        )
        protEmxImport2 = self.newProtocol(ProtEmxImport,
                                          inputEMX=emxFn2
        )
        protEmxImport3 = self.newProtocol(ProtEmxImport,
                                          inputEMX=emxFn3
        )
        self.launchProtocol(protEmxImport1)
        self.launchProtocol(protEmxImport2)
        self.launchProtocol(protEmxImport3)

        import pdb
        pdb.set_trace()
        _list=[protEmxImport1.outputCTF,protEmxImport2.outputCTF,protEmxImport3.outputCTF]
        protCtfDiscrepancy = self.newProtocol(XmippProtCTFDiscrepancy, inputCTFs=_list)
        self.launchProtocol(protCtfDiscrepancy)


        f = open('/tmp/kk.txt', 'w')
        f.write("protCtfDiscrepancy.inputCTFs=%s\n"% type(protCtfDiscrepancy.inputCTFs))
        f.write("protEmxImport1.outputCTF=%s\n"% type(protEmxImport1.outputCTF))
        protCtfDiscrepancy.inputCTFs.append(protEmxImport1.outputCTF)
        protCtfDiscrepancy.inputCTFs.append(protEmxImport2.outputCTF)
        protCtfDiscrepancy.inputCTFs.append(protEmxImport3.outputCTF)
        f.write("protCtfDiscrepancy.inputCTFs=%s\n"% protCtfDiscrepancy.inputCTFs)
        f.write("type protCtfDiscrepancy.inputCTFs=%s\n"% type(protCtfDiscrepancy.inputCTFs))
        self.launchProtocol(protCtfDiscrepancy)
        f.close()


    # for mic1, mic2 in itertools.izip(mics, em.protEmxImport.outputMicrographs):
    #     # Remove the absolute path in the micrographs to
    #     # really check that the attributes should be equal
    #     mic1.setFileName(os.path.basename(mic1.getFileName()))
    #     mic2.setFileName(os.path.basename(mic2.getFileName()))
    #     self.assertTrue(mic1.equalAttributes(mic2))



if __name__ == "__main__":
    unittest.main()