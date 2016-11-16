# ***************************************************************************
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
from pyworkflow.em.protocol import (ProtImportParticles, ProtImportVolumes,
                                    ProtSubSet)


from pyworkflow.em.packages.xmipp3 import XmippProtSolidAngles


class TestDirectionalClasses(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and
        the condition is valid.
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
        
    def runImportParticles(self):
        """ Import an EMX file with Particles and defocus
        """
        partStar = self.ds.getFile('import/case2/relion_it015_data.star')
        prot = self.newProtocol(ProtImportParticles,
                                 objLabel='from relion (data.star)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=partStar,
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot)
        self.checkOutput(prot, 'outputParticles',
                         ['outputParticles.hasAlignmentProj()',
                          'outputParticles.isPhaseFlipped()'])

        return prot

    def runImportVolume(self):
        volFn = self.ds.getFile('import/case2/volume.mrc')
        prot = self.newProtocol(ProtImportVolumes,
                                filesPath=volFn, filesPattern='',
                                samplingRate=7.08)
        self.launchProtocol(prot)

        return prot

    def test_oneClass(self):
        protImportParts = self.runImportParticles()
        protImportVol = self.runImportVolume()

        # Let's keep a smaller subset of particles to speed-up computations
        protSubset = self.newProtocol(ProtSubSet,
                                      objLabel='subset 1K',
                                      chooseAtRandom=True,
                                      nElements=1000)

        protSubset.inputFullSet.set(protImportParts.outputParticles)
        self.launchProtocol(protSubset)

        # We use a coarse angular sampling of 20 to speed-up test
        prot = self.newProtocol(XmippProtSolidAngles,
                                objLabel='directional classes',
                                angularSampling=20,
                                angularDistance=25,
                                numberOfMpi=4
                                )

        prot.inputVolume.set(protImportVol.outputVolume)
        prot.inputParticles.set(protSubset.outputParticles)

        self.launchProtocol(prot)




        
