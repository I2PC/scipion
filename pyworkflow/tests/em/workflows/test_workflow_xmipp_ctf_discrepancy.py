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

from itertools import izip
from pyworkflow.em.data import SetOfCTF
from pyworkflow.em.packages.emxlib import ProtEmxImport
from pyworkflow.em.packages.xmipp3 import XmippProtCTFDiscrepancy
from pyworkflow.object import PointerList
from test_workflow import TestWorkflow
import pyworkflow.tests as tests
import unittest
import os

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
        pl = PointerList()
        pl.append(protEmxImport1.outputCTF)
        pl.append(protEmxImport2.outputCTF)
        pl.append(protEmxImport3.outputCTF)
        protCtfDiscrepancy = self.newProtocol(XmippProtCTFDiscrepancy)
        protCtfDiscrepancy.inputCTFs.set(pl)
        self.launchProtocol(protCtfDiscrepancy)

        ctfsGold = SetOfCTF(filename = self.dataset.getFile('ctfsGold'))
        ctfComputed = protCtfDiscrepancy.outputCTFPair
        for ctf1, ctf2 in izip(ctfComputed, ctfsGold):
            ctf1.getMicrograph().setFileName(os.path.basename(ctf1.getMicrograph().getFileName()))
            ctf2.getMicrograph().setFileName(os.path.basename(ctf2.getMicrograph().getFileName()))
            self.assertTrue(ctf1.equalAttributes(ctf2))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        className = sys.argv[1]
        cls = globals().get(className, None)
        if cls:
            suite = unittest.TestLoader().loadTestsFromTestCase(cls)
            unittest.TextTestRunner(verbosity=2).run(suite)
        else:
            print "Test: '%s' not found." % className
    else:
        unittest.main()
