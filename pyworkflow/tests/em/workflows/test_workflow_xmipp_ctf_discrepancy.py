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

import os
from itertools import izip

from pyworkflow.em.data import SetOfCTF
from pyworkflow.em.packages.xmipp3 import XmippProtCTFDiscrepancy
from pyworkflow.em.protocol import ProtImportMicrographs
from pyworkflow.object import PointerList, Pointer
from test_workflow import TestWorkflow
import pyworkflow.tests as tests



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
        protEmxImport1 = self.newProtocol(ProtImportMicrographs,
                                          importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                          samplingRate=1,
                                          emxFile=emxFn1
        )
        protEmxImport2 = self.newProtocol(ProtImportMicrographs,
                                          importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                          samplingRate=1,
                                          emxFile=emxFn2
        )
        protEmxImport3 = self.newProtocol(ProtImportMicrographs,
                                          importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                          samplingRate=1,
                                          emxFile=emxFn3
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
        ctfSetFn, ctfSetPairFn = protCtfDiscrepancy._getAnalyzeFiles()
        
        ctfComputed = SetOfCTF(filename=ctfSetPairFn)
        for ctf1, ctf2 in izip(ctfComputed, ctfsGold):
            ctf1.getMicrograph().setFileName(os.path.basename(ctf1.getMicrograph().getFileName()))
            ctf2.getMicrograph().setFileName(os.path.basename(ctf2.getMicrograph().getFileName()))
            self.assertTrue(ctf1.equalAttributes(ctf2))
            
            
class TestXmippCTFDiscrepancyBase2(TestWorkflow):
    """ 
    Same test as previous one, but using a different way 
    of setting the inputCTF pointers using the extended
    property of pointers.
    """
    
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
        protEmxImport1 = self.newProtocol(ProtImportMicrographs,
                                          importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                          samplingRate=1,
                                          emxFile=emxFn1
        )
        protEmxImport2 = self.newProtocol(ProtImportMicrographs,
                                          importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                          samplingRate=1,
                                          emxFile=emxFn2
        )
        protEmxImport3 = self.newProtocol(ProtImportMicrographs,
                                          importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                          samplingRate=1,
                                          emxFile=emxFn3
        )
        pl = PointerList([Pointer(value=protEmxImport1, extended='outputCTF'),
                          Pointer(value=protEmxImport2, extended='outputCTF'),
                          Pointer(value=protEmxImport3, extended='outputCTF')
                          ])
        protCtfDiscrepancy = self.newProtocol(XmippProtCTFDiscrepancy)
        protCtfDiscrepancy.inputCTFs.set(pl)

        self.proj.saveProtocol(protEmxImport1)
        self.proj.saveProtocol(protEmxImport2)
        self.proj.saveProtocol(protEmxImport3)

        self.proj.saveProtocol(protCtfDiscrepancy)
    
        self.launchProtocol(protEmxImport1)
        self.launchProtocol(protEmxImport2)
        self.launchProtocol(protEmxImport3)

        self.launchProtocol(protCtfDiscrepancy)
        
        ctfsGold = SetOfCTF(filename = self.dataset.getFile('ctfsGold'))
        ctfSetFn, ctfSetPairFn = protCtfDiscrepancy._getAnalyzeFiles()
        
        ctfComputed = SetOfCTF(filename=ctfSetPairFn)
        for ctf1, ctf2 in izip(ctfComputed, ctfsGold):
            ctf1.getMicrograph().setFileName(os.path.basename(ctf1.getMicrograph().getFileName()))
            ctf2.getMicrograph().setFileName(os.path.basename(ctf2.getMicrograph().getFileName()))
            self.assertTrue(ctf1.equalAttributes(ctf2))
            
                       