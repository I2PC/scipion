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
from pyworkflow.em.protocol import ProtImportVolumes


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsEmx = DataSet.getDataSet('emx')
        cls.dsMda = DataSet.getDataSet('mda')
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportVolumes(TestImportBase):

    def test_import_volume(self):
        """ 1) Import single volume and set origin in default position (the
            volume will be centered in the center of coordinates)
        """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_1_iter_002.mrc',
                'setOrigCoord': False,
                'samplingRate': 2.1,
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot1 = self.newProtocol(ProtImportVolumes, **args)
        prot1.setObjLabel('import vol,\n default origin,\n chimera show '
                          'axis,\n vol origin in coord origin')
        self.launchProtocol(prot1)
        volume = prot1.outputVolume
        t = volume.getOrigin()
        x, y, z = t.getShifts()
        # x, y, z in Angstroms
        # Chimera will show (x, y, z) divided by the samplingRate
        # in pixels = (32, 32, 32)
        self.assertEqual(-67.2, x)
        self.assertEqual(-67.2, y)
        self.assertEqual(-67.2, z)

        # """ 2) Import single volume and set origin in userDefined position
        # """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_1_iter_002.mrc',
                'setOrigCoord': True,
                'samplingRate': 2.1,
                'x': 16.8,
                'y': 33.6,
                'z': 50.4
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot2 = self.newProtocol(ProtImportVolumes, **args)
        prot2.setObjLabel('import vol,\n user origin,\n chimera show axis')
        self.launchProtocol(prot2)
        volume = prot2.outputVolume
        t = volume.getOrigin()
        x, y, z = t.getShifts()

        # x, y, z in Angstroms
        # in Chimera we will see (x, y, z) divided by the samplingRate
        # in pixels = (8, 16, 24)
        self.assertEqual(-16.8, x)
        self.assertEqual(-33.6, y)
        self.assertEqual(-50.4, z)

        # """ 3) Import single volume and set origin in userDefined position (the
        #     volume will be centered in the center of coordinates because
        #     coordinates are the default ones)
        # """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_1_iter_002.mrc',
                'setOrigCoord': True,
                'samplingRate': 2.1,
                'x': 67.2,
                'y': 67.2,
                'z': 67.2
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot2 = self.newProtocol(ProtImportVolumes, **args)
        prot2.setObjLabel('import vol,\n user origin,\n chimera show axis')
        self.launchProtocol(prot2)
        volume = prot2.outputVolume
        t = volume.getOrigin()
        x, y, z = t.getShifts()

        # x, y, z in Angstroms
        # in Chimera we will see (x, y, z) divided by the samplingRate
        # in pixels = (8, 16, 24)
        self.assertEqual(-67.2, x)
        self.assertEqual(-67.2, y)
        self.assertEqual(-67.2, z)

        # """ 4) Import two volumes and set origin in default position
        # """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_*mrc',
                'setOrigCoord': False,
                'samplingRate': 2.1,
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot3 = self.newProtocol(ProtImportVolumes, **args)
        prot3.setObjLabel('import 2 vols,\n default origin,\n chimera no '
                          'axis')
        self.launchProtocol(prot3)
        for volume in prot3.outputVolumes:
            t = volume.getOrigin()
            x, y, z = t.getShifts()
            self.assertEqual(-67.2, x)
            self.assertEqual(-67.2, y)
            self.assertEqual(-67.2, z)

        # """ 5) Import two volumes and set origin in userDefined position
        # """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_*mrc',
                'setOrigCoord': True,
                'samplingRate': 2.1,
                'x': 16.8,
                'y': 33.6,
                'z': 50.4
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot4 = self.newProtocol(ProtImportVolumes, **args)
        prot4.setObjLabel('import 2 vols,\n user origin,\n chimera no axis')
        self.launchProtocol(prot4)
        for volume in prot4.outputVolumes:
            t = volume.getOrigin()
            x, y, z = t.getShifts()
            self.assertEqual(-16.8, x)
            self.assertEqual(-33.6, y)
            self.assertEqual(-50.4, z)

        # """ 6) Import three volumes (mrc stack) and set origin in userDefined
        # position
        # """
        # Id's should be taken from filename
        args['filesPath'] = self.dsRelion.getFile('import/case2/'
                                                  'relion_volumes.mrc')
        args['filesPattern'] = ''
        prot2 = self.newProtocol(ProtImportVolumes, **args)
        prot2.setObjLabel('import 3 vols from mrc stack,\n user origin,'
                          '\n chimera no axis')
        self.launchProtocol(prot2)
        # Check the number of output volumes and dimensions
        self.assertEqual(3, prot2.outputVolumes.getSize())
        self.assertEqual(60, prot2.outputVolumes.getDim()[0])

        # """ 7) Import three volumes (spider stack) and set origin in
        # userDefined position
        # """
        # Id's should be taken from filename
        args['filesPath'] = self.dsRelion.getFile('import/case2/'
                                                  'relion_volumes.stk')
        args['filesPattern'] = ''
        prot3 = self.newProtocol(ProtImportVolumes, **args)
        prot3.setObjLabel('import 3 vols from spider stack\n user origin,'
                          '\n chimera no axis')
        self.launchProtocol(prot3)
        # Check the number of output volumes and dimensions
        self.assertEqual(3, prot3.outputVolumes.getSize())
        self.assertEqual(60, prot3.outputVolumes.getDim()[0])


        #  """ 8)To test old data where volumes have no origin at all"""
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_1_iter_002.mrc',
                # 'setOrigCoord': False,
                'samplingRate': 2.1,
                }

        prot4 = self.newProtocol(ProtImportVolumes, **args)
        prot4.setObjLabel('import vol,\n no origin at all, legacy data,'
                          '\n chimera show '
                          'axis,\n vol origin in coord origin')
        self.launchProtocol(prot4)
        volume = prot4.outputVolume
        volume.setOrigin(None)
        # The volume has no origin
        t = volume.getOrigin(force=True)
        x, y, z = t.getShifts()
        # x, y, z in Angstroms
        # Chimera will show (x, y, z) divided by the samplingRate
        # in pixels = (32, 32, 32)
        self.assertEqual(-67.2, x)
        self.assertEqual(-67.2, y)
        self.assertEqual(-67.2, z)

        ## TODO: associate origen coordinates from header file
        # args = {'filesPath': self.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
        #         'samplingRate': 1.5,
        #         'setDefaultOrigin': False,
        #         'x': ,
        #         'y': ,
        #         'z':
        #         }
        #
        # protImportVol = self.newProtocol(ProtImportVolumes, **args)
        # protImportVol.setObjLabel('import vol 1ake_4-5A,\n header origin,\n chimera show axis')
        # self.launchProtocol(protImportVol)
        # volume = protImportVol.outputVolume
        #
        # # Id's should be set increasing from 1 if ### is not in the
        # # pattern
        #
        # t = volume.getOrigin()
        # x, y, z = t.getShifts()
        # # x, y, z in Angstroms
        # # in Chimera we will see (x, y, z) divided by the samplingRate
        # self.assertEqual(-11.99, x)
        # self.assertEqual(7.88, y)
        # self.assertEqual(-10.90, z)
