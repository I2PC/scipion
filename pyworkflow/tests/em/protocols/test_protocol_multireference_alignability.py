# ***************************************************************************
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
#                J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************/

#import os
#from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes, ProtSubSet
from pyworkflow.em.packages.xmipp3 import XmippProtMultiRefAlignability


class TestMultireferenceAlignability(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        
    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and
        the condition is valid. 
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o 
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
            
    def importVolumes(self):
        prot = self.newProtocol(ProtImportVolumes, 
                                objLabel='import volume', 
                                filesPath=self.dsRelion.getFile('volumes/reference.mrc'),
                                samplingRate=7.08)
        self.launchProtocol(prot)
        
        return prot

    def importParticles(self, numberOfParticles):
        """ Import an EMX file with Particles and defocus
        """
        prot = self.newProtocol(ProtImportParticles,
                                 objLabel='from relion (auto-refine 3d)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=self.dsRelion.getFile('import/refine3d/extra/relion_it001_data.star'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot)
        self.checkOutput(prot, 'outputParticles', 
                         ['outputParticles.hasAlignmentProj()',
                          'outputParticles.isPhaseFlipped()'])
        # We are going to make a subset to speed-up following processes
        protSubset = self.newProtocol(ProtSubSet, 
                                      objLabel='subset of particles',
                                      chooseAtRandom=True,
                                      nElements=numberOfParticles)
        protSubset.inputFullSet.set(prot.outputParticles)
        
        self.launchProtocol(protSubset)
        self.checkOutput(protSubset, 'outputParticles', 
                         ['outputParticles.getSize() == %d' 
                          % numberOfParticles])
        
        return protSubset
    
    def test_validate(self):
        protImportVols = self.importVolumes()
        protImportPars = self.importParticles(20)        
        protValidate = self.newProtocol(XmippProtMultiRefAlignability,
                                objLabel='validate alignability',
                                angularSampling=10,
                                numberOfMpi=4, numberOfThreads=1)
        protValidate.inputParticles.set(protImportPars.outputParticles)
        protValidate.inputVolumes.set(protImportVols.outputVolume)
        
        self.launchProtocol(protValidate)
        
        self.checkOutput(protValidate, 'outputParticles')
        self.checkOutput(protValidate, 'outputVolumes')
        
        firstVol = protValidate.outputVolumes.getFirstItem()
        self.assertTrue(firstVol.weightAlignabilityPrecision > 0.7)
        self.assertTrue(firstVol.weightAlignabilityAccuracy > 0.1)

