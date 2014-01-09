# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.brandeis import *


class TestBrandeisCtffind(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):    
        # Create a new project for the tests
        setupProject(cls)
        
    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage):
        """ Run an Import micrograph protocol. """
        cls.protImport = ProtImportMicrographs(pattern=pattern, samplingRate=samplingRate, voltage=voltage)
        cls.proj.launchProtocol(cls.protImport, wait=True)
        # check that input micrographs have been imported (a better way to do this?)
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. outputMicrographs is None.' % pattern)
        return cls.protImport
    
    def testCtffind(self):
        #First, import a set of micrographs
        pattern = getInputPath('Micrographs_BPV1', '*.mrc')
        protImport = self.runImportMicrograph(pattern, samplingRate=1.273, voltage=300)
        
        protCTF = ProtCTFFind()
        protCTF.inputMicrographs.set(protImport.outputMicrographs)
        self.proj.launchProtocol(protCTF, wait=True)
        
        self.assertIsNotNone(protCTF.outputCTF, "SetOfCTF has not been produced.") 
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBrandeisCtffind)
    #suite = unittest.TestLoader().loadTestsFromName('test_protocols_brandeis.TestBrandeisCtffind.testCtffind')
    unittest.TextTestRunner(verbosity=2).run(suite)