# ***************************************************************************
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es) (2017)
#                Carlos Oscar S. Sorzano   (coss@cnb.csic.es)
# *
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
# ***************************************************************************/

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles, ProtSplitSet
from pyworkflow.em.packages.xmipp3 import XmippProtCompareAngles


class TestXmippProtCompareAngles(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        
    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and the condition is valid.
            In addition, returns the size of the set.
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o 
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
        return o.getSize()

    def importParticles(self, numberOfParticles, path):
        """ Import an EMX file with Particles and defocus.
            Returns a Subset of the first numberOfParticles.
        """
        protImport = self.newProtocol(ProtImportParticles,
                                 objLabel='from relion (auto-refine 3d)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=self.dsRelion.getFile(path),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(protImport)
        nParticles = self.checkOutput(protImport, 'outputParticles', 
                         ['outputParticles.hasAlignmentProj()'])

        nOfSets=int(nParticles/numberOfParticles)
        # We are going to make a subset to speed-up following processes
        protSubset = self.newProtocol(ProtSplitSet, 
                                      randomize=False,
                                      numberOfSets=nOfSets)
        protSubset.inputSet.set(protImport.outputParticles)
        self.launchProtocol(protSubset)

        return getattr(protSubset, 'outputParticles01', None)
    
    def test_validate(self):
        pathNoCTF = 'import/refine3d/extra/relion_it001_data.star'
        partNoCTF = self.importParticles(20,pathNoCTF)

        pathCTF = 'import/case2/relion_it015_data.star'
        partCTF = self.importParticles(20,pathCTF)

        
        '''Same particles'''
        protValidateSameParts = self.newProtocol(XmippProtCompareAngles)
        protValidateSameParts.inputParticles1.set(partNoCTF)
        protValidateSameParts.inputParticles2.set(partNoCTF)

        '''Slighly different particles'''
        protValidateDiffParts = self.newProtocol(XmippProtCompareAngles)
        protValidateDiffParts.inputParticles1.set(partCTF)
        protValidateDiffParts.inputParticles2.set(partNoCTF)

        '''With Symmetry'''
        protValidateSymm = self.newProtocol(XmippProtCompareAngles,
                                            symmetryGroup='d2')
        protValidateSymm.inputParticles1.set(partNoCTF)
        protValidateSymm.inputParticles2.set(partCTF)
        

        '''Assertions'''
        self.launchProtocol(protValidateSameParts)      
        self.checkOutput(protValidateSameParts, 'outputParticles') 
        firstPart = protValidateSameParts.outputParticles.getFirstItem()
        self.assertAlmostEqual(float(firstPart._xmipp_shiftDiff), 0, 0.01) 
        self.assertAlmostEqual(float(firstPart._xmipp_angleDiff), 0, 0.01) 

        self.launchProtocol(protValidateDiffParts)      
        self.checkOutput(protValidateDiffParts, 'outputParticles') 
        firstPart = protValidateDiffParts.outputParticles.getFirstItem()
        self.assertAlmostEqual(float(firstPart._xmipp_shiftDiff), 1.96, 1) 
        self.assertAlmostEqual(float(firstPart._xmipp_angleDiff), 58.9, 1) 

        self.launchProtocol(protValidateSymm)      
        self.checkOutput(protValidateSymm, 'outputParticles') 
        firstPart = protValidateSymm.outputParticles.getFirstItem()
        self.assertAlmostEqual(float(firstPart._xmipp_shiftDiff), 1.96, 1) 
        self.assertAlmostEqual(float(firstPart._xmipp_angleDiff), 58.9, 1)
