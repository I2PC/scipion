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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from test_protocols_xmipp import TestXmippBase


    
class TestXmippCLTomo(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testCLTomo(self):
        print "Import volumes"
        protImportVol = ProtImportVolumes(pattern=getInputPath('CLTomo', 'subvols*.spi'), samplingRate=1)
        self.proj.launchProtocol(protImportVol, wait=True)
    
        print "Run CLTomo"
        protCLTomo = XmippProtCLTomo(numberOfReferences=1,numberOfIterations=1)
        protCLTomo.volumelist.set(protImportVol.outputVolumes)
        self.proj.launchProtocol(protCLTomo, wait=True)        
        
        self.assertIsNotNone(protCLTomo.outputClasses, "There was a problem with CLTomo output classes")          
        self.assertIsNotNone(protCLTomo.alignedVolumes, "There was a problem with CLTomo output aligned volumes")          

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
