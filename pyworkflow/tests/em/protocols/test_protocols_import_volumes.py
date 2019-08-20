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
import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.em.constants import SYM_I222r, SCIPION_SYM_NAME
from pyworkflow.em.convert.symmetry import Icosahedron
import pyworkflow.utils as pwutils


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
                'setHalfMaps': False,
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
                'setHalfMaps': False,
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

        # Assert path is relative
        self.assertFalse(os.path.isabs(volume.getFileName()),
                        "single volume path is not relative")

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
                'setHalfMaps': False,
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
                'setHalfMaps': False,
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
            # Assert path is relative
            self.assertFalse(os.path.isabs(volume.getFileName()),
                            "volume path is not relative")

        # """ 5) Import two volumes and set origin in userDefined position
        # """
        args = {'filesPath': self.dsXmipp.getFile('volumes/'),
                'filesPattern': 'volume_*mrc',
                'setHalfMaps': False,
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
                'setHalfMaps': False,
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

    def test_import_volume2(self):
        """
        Test to import a full map (Icosahedron) and two maps (half1 and half2)
        to compute the FSC
        """
        volFeatName = '/tmp/Icosahedron_map.txt'
        volMapNamefull = '/tmp/Icosahedron_map_full.mrc'
        volMapNamehalf1 =  '/tmp/Icosahedron_map_half1.mrc'
        volMapNamehalf2 = '/tmp/Icosahedron_map_half2.mrc'
        self.createFeatVolume(volFeatName, volMapNamefull, sym=SYM_I222r)
        self.createFeatVolume(volFeatName, volMapNamehalf1, sym=SYM_I222r)
        self.createFeatVolume(volFeatName, volMapNamehalf2, sym=SYM_I222r)
        _samplingRate = 1.0
        args = {'filesPath': volMapNamefull,
                'filesPattern': '',
                'setHalfMaps': True,
                'half1map': volMapNamehalf1,
                'half2map': volMapNamehalf2,
                'samplingRate': _samplingRate
                }
        prot5 = self.newProtocol(ProtImportVolumes, **args)
        prot5.setObjLabel('import phantom icosahedron,\n half1 and half2')
        self.launchProtocol(prot5)
        volume = prot5.outputVolume
        volume.setOrigin(None)
        # The volume has no origin
        t = volume.getOrigin(force=True)
        x, y, z = t.getShifts()
        # x, y, z in Angstroms
        # Chimera will show (x, y, z) divided by the samplingRate
        # in pixels = ()
        self.assertEqual(-90.0, x)
        self.assertEqual(-90.0, y)
        self.assertEqual(-90.0, z)

    def __runXmippProgram(self, program, args):
        """ Internal shortcut function to launch a Xmipp program.
        If xmipp not available o fails return False, else Tru"""
        try:
            xmipp3 = pwutils.importFromPlugin('xmipp3', doRaise=True)
            xmipp3.Plugin.runXmippProgram(program, args)
        except:
            return False
        return True


    def createFeatVolume(self, volFeatName, volMapName, sym=SYM_I222r):
        f = open(volFeatName, "w")
        f.write("""# Phantom description file, (generated with phantom help)
# General Volume Parameters:
#      Xdim      Ydim      Zdim   Background_Density Scale
       180 180 180 0 1.0
# Feature Parameters:
#Type  +/=  Density X_Center Y_Center Z_Center
""")
        icosahedron = Icosahedron(orientation=SCIPION_SYM_NAME[sym][1:])
        x = 0.;
        y = 0.;
        z = 0.
        f.write("# large sphere at the center\n")
        f.write("sph  + 1. %.3f %.3f %.3f 36.\n" % (x, y, z))
        f.write("# 5-fold\n")

        for i, vertice in enumerate(icosahedron.getVertices()):
            vertice = 55.0 * vertice
            f.write("sph  + 3 %.3f %.3f %.3f 8.25\n" %
                    (vertice[0], vertice[1], vertice[2]))
            if i==0:
                self.pentonDir =  "%.3f, %.3f, %.3f"%(vertice[0], vertice[1], vertice[2])

        # print 3fold points
        f.write("# 3-fold\n")

        for _3fold in icosahedron.get3foldAxis():
            x, y, z = _3fold
            f.write("sph  + 0.8 %.3f %.3f %.3f 6.0\n" % (55.0 * x, 55.0 * y, 55.0 * z))

        # print 2fold points
        f.write("# 2-fold\n")
        for _2fold in icosahedron.get2foldAxis():
            x, y, z = _2fold
            f.write("sph  + 0.7 %.3f %.3f %.3f 3.0\n" %
                    (55.0 * x, 55.0 * y, 55.0 * z))
        f.close()
        #    map
        program = "xmipp_phantom_create"
        args= '-i {featFile} -o {mapFile}'.format(
            featFile=volFeatName, mapFile=volMapName)
        self.__runXmippProgram(program, args)