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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************/

import os
from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles
from pyworkflow.em.data import SetOfParticles


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsEmx = DataSet.getDataSet('emx')
        cls.dsMda = DataSet.getDataSet('mda')
        
    
class TestImportParticles(TestImportBase):
    
    def test_pattern(self):
        """ Import several Particles from a given pattern.
        """
        args = {'importFrom': ProtImportParticles.IMPORT_FROM_FILES,
                'filesPath': self.dsXmipp.getFile('particles/'),
                'filesPattern': 'BPV_????_ptcls.hdf',
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 100,
                'samplingRate': 2.1
                }
        
        
        # Id's should be set increasing from 1 if ### is not in the 
        # pattern
        protMicImport = self.newProtocol(ProtImportParticles, **args)
        protMicImport.setObjLabel('from files')
        self.launchProtocol(protMicImport)
        
        # Id's should be taken from filename    
        args['filesPattern'] = 'BPV_####_ptcls.hdf'
        protMicImport = self.newProtocol(ProtImportParticles, **args)
        protMicImport.setObjLabel('from files (with mic id)')
        self.launchProtocol(protMicImport)

        # Import some Particles from EMX
        emxFn = self.dsEmx.getFile('coordinatesT1')
        args['importFrom'] = ProtImportParticles.IMPORT_FROM_EMX
        args['micrographsEMX'] = emxFn
        protEmxImport = self.newProtocol(ProtImportParticles, **args)
        protEmxImport.setObjLabel('from emx (with coords)')
        self.launchProtocol(protEmxImport)
    
    def test_fromEmx(self):
        """ Import an EMX file with Particles and defocus
        """
        protEmxImport = self.newProtocol(ProtImportParticles,
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=self.dsEmx.getFile('particles/particles.emx'),
                                         magnification=10000,
                                         samplingRate=2.46
                                         )
        protEmxImport.setObjLabel('from emx ')
        self.launchProtocol(protEmxImport)
        
    def test_fromXmipp(self):
        """ Import an EMX file with Particles and defocus
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                         importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                         mdFile=self.dsXmipp.getFile('gold/xmipp_ml2d_images.xmd'),
                                         magnification=10000,
                                         samplingRate=1.
                                         )
        prot1.setObjLabel('from xmipp (ml2d)')
        self.launchProtocol(prot1)
        
    def test_fromXmippWithMic(self):
        """ Import an EMX file with Particles and defocus
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                         importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                         mdFile=self.dsXmipp.getFile('gold/images10.xmd'),
                                         magnification=10000,
                                         samplingRate=1.
                                         )
        prot1.setObjLabel('from xmipp (with mic id)')
        self.launchProtocol(prot1)

