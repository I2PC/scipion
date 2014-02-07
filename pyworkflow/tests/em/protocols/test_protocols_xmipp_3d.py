# **************************************************************************
# *
# * Authors:    Josue Gomez Blanco (jgomez@cnb.csic.es)
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


    
class TestXmippCeateMask3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testCreateMask1(self):
        print "Import volume"
        protImportVol = ProtImportVolumes(pattern=getInputPath('Volumes_BPV', 'BPV_scale_filtered_windowed_64.vol'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol, wait=True)
    
        print "Run create mask from volume"
        protMask1 = XmippProtCreateMask3D(source=0, threshold=0.4)
        protMask1.volume.set(protImportVol.outputVolumes)
        self.proj.launchProtocol(protMask1, wait=True)        
        
        self.assertIsNotNone(protMask1.outputMask, "There was a problem with create mask from volume")          

        print "Run create mask from geometry"
        protMask2 = XmippProtCreateMask3D(source=1, size=64, samplingRate=9.89, geo=6, innerRadius=10, outerRadius=25, borderDecay=2)
        self.proj.launchProtocol(protMask2, wait=True)        
        
        self.assertIsNotNone(protMask2.outputMask, "There was a problem with create mask from geometry")          

        print "Run create mask from another mask"
        protMask3 = XmippProtCreateMask3D(source=2, doMorphological=True, elementSize=3)
        protMask3.inputMask.set(protMask1.outputMask)
        self.proj.launchProtocol(protMask3, wait=True)        
         
        self.assertIsNotNone(protMask3.outputMask, "There was a problem with mask from another mask")          

class TestXmippResolution3D(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        setupProject(cls)
        
#         cls.protImport = cls.runImportParticles(pattern=images, samplingRate=1, checkStack=False)
#         cls.iniVol = getInputPath('ml3dData', 'icoFiltered.vol')
    
    def testCreateMask1(self):
        print "Import volume 1"
        protImportVol1 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_1_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol1, wait=True)

        print "Import volume 2"
        protImportVol2 = ProtImportVolumes(pattern=getInputPath('Test_Resolution_Xmipp', 'volume_2_iter_002.mrc'), samplingRate=9.896)
        self.proj.launchProtocol(protImportVol2, wait=True)
    
        print "Run resolution 3D"
        protResol3D = XmippProtResolution3D(doSSNR=False)
        protResol3D.inputVolume.set(protImportVol1.outputVolume)
        protResol3D.referenceVolume.set(protImportVol2.outputVolume)
        self.proj.launchProtocol(protResol3D, wait=True)        
         
        self.assertIsNotNone(protResol3D._defineFscName(), "There was a problem with fsc")
        self.assertIsNotNone(protResol3D._defineStructFactorName(), "There was a problem with structure factor")

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
