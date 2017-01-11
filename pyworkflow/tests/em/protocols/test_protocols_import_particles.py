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

#import os
#from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportParticles


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsEmx = DataSet.getDataSet('emx')
        cls.dsMda = DataSet.getDataSet('mda')
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
                'samplingRate': 2.1,
                'haveDataBeenPhaseFlipped': True
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

    def test_fromEmx(self):
        """ Import an EMX file with Particles and defocus
        """
        args = {'importFrom': ProtImportParticles.IMPORT_FROM_EMX,
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 100,
                'alignType' : 3,
                'samplingRate': 2.46,
                'magnification': 10000,
                'haveDataBeenPhaseFlipped': True
                }
        args['emxFile'] = self.dsEmx.getFile('particles/particles.emx')
        protEmxImport = self.newProtocol(ProtImportParticles, **args)
        protEmxImport.setObjLabel('from emx (particles)')
        self.launchProtocol(protEmxImport)

        # Import some Particles from EMX
        args['emxFile'] = self.dsEmx.getFile('coordinatesT1')
        protEmxImport = self.newProtocol(ProtImportParticles, **args)
        protEmxImport.setObjLabel('from emx (with coords)')
        self.launchProtocol(protEmxImport)

    def test_fromXmipp(self):
        """ Import an Xmipp file with Particles and defocus
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                 importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                 mdFile=self.dsXmipp.getFile('gold/xmipp_ml2d_images.xmd'),
                                 magnification=10000,
                                 samplingRate=1,
                                 haveDataBeenPhaseFlipped=True
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
                                 samplingRate=1,
                                 haveDataBeenPhaseFlipped=True
                                 )
        prot1.setObjLabel('from xmipp (with mic id)')
        self.launchProtocol(prot1)

    def test_fromRelionRefine3D(self):
        """ Import an EMX file with Particles and defocus
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from relion (auto-refine 3d)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=self.dsRelion.getFile('import/refine3d/extra/relion_it001_data.star'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot1)
        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignmentProj()',
                                                    'outputParticles.isPhaseFlipped()'])
        
    def test_fromRelionClassify2D(self):
        """ Import an EMX file with Particles and defocus
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from relion (classify 2d)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=self.dsRelion.getFile('import/classify2d/extra/relion_it015_data.star'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot1)
        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignment2D()',
                                                    'outputParticles.isPhaseFlipped()'])
        self.checkOutput(prot1, 'outputClasses')    
 
    def test_fromRelionClassify3D(self):
        """ Import an EMX file with Particles and defocus
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from relion (classify 3d)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=self.dsRelion.getFile('import/classify3d/extra/relion_it015_data.star'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot1)         
        
    def test_fromScipionReconstruct(self):
        """ Import the particles with 3D projection directions for reconstruct
        a volume. Actually this test use a similar .sqlite file than
        the output of
        test_fromRelionRefine3D
        """
        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from scipion (to-reconstruct)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dsRelion.getFile('import/case2/particles.sqlite'),
                                 magnification=10000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(prot1)
        
        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasCTF()', 
                                                    'outputParticles.hasAlignmentProj()',
                                                    'outputParticles.isPhaseFlipped()'])
        
    def test_fromImagic(self):
        """ Import particles from imagic stack.
        """
        args = {'importFrom': ProtImportParticles.IMPORT_FROM_FILES,
                'filesPath': self.dsMda.getFile('particles/'),
                'filesPattern': 'particles.hed',
                'amplitudConstrast': 0.1,
                'sphericalAberration': 2.,
                'voltage': 100,
                'samplingRate': 2.1,
                }


        # Id's should be set increasing from 1 if ### is not in the 
        # pattern
        protImport = self.newProtocol(ProtImportParticles, **args)
        protImport.setObjLabel('from imagic')
        self.launchProtocol(protImport)
        
        self.checkOutput(protImport, 'outputParticles', [])
