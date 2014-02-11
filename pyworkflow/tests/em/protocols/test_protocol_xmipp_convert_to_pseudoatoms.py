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
    
class TestXmippConvertToPseudotaoms(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testCLTomo(self):
        print "Import volumes"
        protImportVol = ProtImportVolumes(pattern=getInputPath('Volumes_BPV', 'BPV_scale_filtered_windowed_110.vol'), samplingRate=6.5)
        self.proj.launchProtocol(protImportVol, wait=True)
    
        print "Run convert to pseudoatoms"
        prot = XmippProtConvertToPseudoAtoms(pseudoAtomTarget=15)
        prot.inputStructure.set(protImportVol.outputVolume)
        self.proj.launchProtocol(prot, wait=True)        
        
        self.assertIsNotNone(prot.outputPdb, "There was a problem with Convert to pseudoatoms output Pdb")          

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
