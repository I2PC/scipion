# ***************************************************************************
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.tests as tests
from pyworkflow.em.data import Particle, SetOfParticles
from pyworkflow.em.protocol import ProtSplitSet
from pyworkflow.tests import envVarOn
from pyworkflow.protocol import MODE_RESTART
from pyworkflow.em.protocol import ProtImportParticles
from os import remove
from os.path import exists
import xmipp


class TestSetsBase(tests.BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        #cls.dataset = tests.DataSet.getDataSet('emx')

    def test_orderBy(self):
        """ create set of particles and orderby a given attribute
        """
        #create set of particles

        inFileNameMetadata = self.proj.getTmpPath('particlesOrderBy.sqlite')
        inFileNameData = '/tmp/images.stk'

        imgSet = SetOfParticles(filename=inFileNameMetadata)
        imgSet.setSamplingRate(1.5)
        img = Particle()

        for i in range(1, 10):
            img.setLocation(i, inFileNameData)
            img.setMicId(i%3)
            img.setClassId(i%5)
            imgSet.append(img)
            img.cleanObjId()
        imgSet.write()
        imgSet.close()
        #now import the dataset
        #this newProtocols should be a class method?
        prot1 = self.newProtocol(ProtImportParticles,
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=inFileNameMetadata,
                                 magnification=10000,
                                 samplingRate=1.5
                                 )
        prot1.setObjLabel('from sqlite (test-sets)')
        self.launchProtocol(prot1)

        if prot1.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % inFileNameMetadata)
        
        #save it
        protSplitSet   = self.newProtocol(ProtSplitSet,
                                          inputSet=prot1.outputParticles,
                                          numberOfSets=2,
                                          randomize=True)
        self.launchProtocol(protSplitSet)
        protSplitSet.outputParticles01.printAll()
        #assert
        #delete file
        ##remove(inFileNameMetadata)
        ##remove(inFileNameData)
