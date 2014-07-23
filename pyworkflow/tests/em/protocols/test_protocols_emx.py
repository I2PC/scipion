"""
/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *


class TestEmxBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('emx')
    
    def test_coodinatesTest1(self):
        """ Import an EMX file with just one micrograph 
        and a few coordinates.
        """
        emxFn = self.dataset.getFile('coordinatesT1')
        protEmxImport = self.newProtocol(ProtEmxImport, 
                                         inputEMX=emxFn,
                                         amplitudeContrast=0.1,
                                         sphericalAberration=2.,
                                         voltage=100)
        self.launchProtocol(protEmxImport)

        # Reference coordinates
        coords = SetOfCoordinates(filename=self.dataset.getFile('coordinatesGoldT1'))
        BaseTest.compareSets(self, protEmxImport.outputCoordinates, coords)
        
    def test_particleDefocus(self):
        """ Import an EMX file with a stack of particles
        that has defocus
        """
        #SCIPION_TESTS=/home/roberto/Scipion/pyworkflow-code/data/tests
        #SCIPION_USER_DATA=/home/roberto/ScipionUserData

        print "self.folder", self.dataset.folder
        #os.chdir(os.environ['SCIPION_TESTS'])
        print "os.environ['SCIPION_HOME']", os.environ['SCIPION_HOME']
        emxFn = self.dataset.getFile('defocusParticleT2')
        print "emxFn1",emxFn
        protEmxImport = self.newProtocol(ProtEmxImport, 
                                         inputEMX=emxFn
                                         )
        self.launchProtocol(protEmxImport)
        micFn =self.dataset.getFile('micrographsGoldT2')
        mics  = SetOfMicrographs(filename = micFn)
        
        micsTmp1  = protEmxImport._createSetOfMicrographs('kk1')
        micsTmp2  = protEmxImport._createSetOfMicrographs('kk2')
        for mic in mics:
            fn = os.path.basename(mic.getFileName())
            mic.setFileName(String(fn))
            micsTmp1.append(mic)
        for mic in protEmxImport.outputMicrographs:
            fn = os.path.basename(mic.getFileName())
            mic.setFileName(fn)
            micsTmp2.append(mic)

        BaseTest.compareSets(self, micsTmp2, micsTmp1)




