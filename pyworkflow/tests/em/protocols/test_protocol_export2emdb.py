# **************************************************************************
# *
# * Authors:    Amaya Jimenez (ajimenez@cnb.csic.es)
# *             Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.protocol import ProtImportVolumes, ProtImportPdb
from pyworkflow.em.protocol.protocol_export import ProtExportEMDB
from pyworkflow.utils import importFromPlugin
from pyworkflow.em.constants import SYM_I222r, SCIPION_SYM_NAME
from pyworkflow.em.convert.symmetry import Icosahedron
import pyworkflow.utils as pwutils

from pyworkflow.em import Volume

try:
    XmippProtMultipleFSCs = importFromPlugin('xmipp3.protocols', 'XmippProtMultipleFSCs', doRaise=True)
    XmippProtResolution3D = importFromPlugin('xmipp3.protocols', 'XmippProtResolution3D', doRaise=True)
    XmippProtCreateMask3D = importFromPlugin('xmipp3.protocols', 'XmippProtCreateMask3D', doRaise=True)
    ChimeraProtRigidFit = importFromPlugin('chimera.protocols', 'ChimeraProtRigidFit', doRaise=True)
    ChimeraProtOperate = importFromPlugin('chimera.protocols', 'ChimeraProtOperate', doRaise=True)
except Exception as e:
    print(str(e))
    exit(0)

class TestExport2EMDB(BaseTest):
    VOLFEATNAME = '/tmp/Icosahedron_map.txt'
    VOLMAPNAMEFULL = '/tmp/Icosahedron_map_full.mrc'
    VOLMAPNAMEHALF1 = '/tmp/Icosahedron_map_half1.mrc'
    VOLMAPNAMEHALF2 = '/tmp/Icosahedron_map_half2.mrc'
    SAMPLINGRATE = 1.0

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate, label):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         objLabel=label,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        cls.protImportHalf1  = cls.runImportVolumes(cls.half1, 3.54,
                                                    'import half1')
        cls.protImportHalf2  = cls.runImportVolumes(cls.half2, 3.54,
                                                    'import half2')

    def _importAtomStructCIF(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import atom struct\nmmCIF\n5ni1.cif')
        self.launchProtocol(protImportPDB)
        structure_mmCIF = protImportPDB.outputPdb
        return structure_mmCIF
    
    def _importAtomStructPDB(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import atom struct\nPDB\n5ni1.pdb')
        self.launchProtocol(protImportPDB)
        structure_PDB = protImportPDB.outputPdb
        return structure_PDB


    def _import_volume_icosahedron(self):
        """
        Test to import a full map (Icosahedron) and two maps (half1 and half2)
        to compute the FSC
        """
        self.createFeatVolume(self.VOLFEATNAME, self.VOLMAPNAMEFULL, sym=SYM_I222r)
        self.createFeatVolume(self.VOLFEATNAME, self.VOLMAPNAMEHALF1, sym=SYM_I222r)
        self.createFeatVolume(self.VOLFEATNAME, self.VOLMAPNAMEHALF2, sym=SYM_I222r)

        args = {'filesPath': self.VOLMAPNAMEFULL,
                'filesPattern': '',
                'setHalfMaps': True,
                'half1map': self.VOLMAPNAMEHALF1,
                'half2map': self.VOLMAPNAMEHALF2,
                'samplingRate': self.SAMPLINGRATE
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
        return volume

    def __runXmippProgram(self, program, args):
        """ Internal shortcut function to launch a Xmipp program.
        If xmipp not available o fails return False, else True"""
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
            if i == 0:
                self.pentonDir = "%.3f, %.3f, %.3f" % (vertice[0], vertice[1], vertice[2])

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
        args = '-i {featFile} -o {mapFile}'.format(
            featFile=volFeatName, mapFile=volMapName)
        self.__runXmippProgram(program, args)

    def test_volume1_AtomStructCIF(self):

        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.protImportHalf1.outputVolume)
        prot.referenceVolume.set(self.protImportHalf2.outputVolume)
        self.launchProtocol(prot)

        protExp = self.newProtocol(ProtExportEMDB)
        protExp.exportVolume.set(self.protImportHalf1.outputVolume)
        protExp.exportFSC.set(prot.outputFSC)

        protExp.exportAtomStruct.set(self._importAtomStructCIF())

        protExp.filesPath.set(os.getcwd() + "/dir1")
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        #protExp._createFileNamesTemplates()
        dirName = protExp.filesPath.get()
        nameVolume = os.path.join(dirName, protExp.VOLUMENAME)
        nameFsc = os.path.join(dirName, ProtExportEMDB.FSC % 0)
        nameAtomStruct = os.path.join(dirName, protExp.COORDINATEFILENAME)
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameFsc))
        self.assertTrue(os.path.exists(nameAtomStruct))

        #Chek if the files have the correct data
        orig_list_x, orig_list_y = protExp.exportFSC.get().getData()
        fo = open(nameFsc, "rU")
        saved_x=[]
        orig_x=[]
        count=0
        for line in fo:
            if line[0:3]=='<x>':
                saved_x.append(int(float(line[3:-5])*1000))
                orig_x.append(int(orig_list_x[count]*1000))
                count=count+1

        self.assertListEqual(orig_x, saved_x)

    def test_volume1_AtomStructPDB(self):

        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.protImportHalf1.outputVolume)
        prot.referenceVolume.set(self.protImportHalf2.outputVolume)
        self.launchProtocol(prot)

        protExp = self.newProtocol(ProtExportEMDB)
        protExp.exportVolume.set(self.protImportHalf1.outputVolume)
        protExp.exportFSC.set(prot.outputFSC)

        # run Chimera rigid fit to get a PDB file
        extraCommands = ""
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 0')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': self.protImportHalf1.outputVolume,
                'pdbFileToBeRefined': self._importAtomStructPDB()
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and PDB\n save model')
        self.launchProtocol(protChimera)
        protExp.exportAtomStruct.set(protChimera.outputPdb_01)

        protExp.filesPath.set(os.getcwd() + "/dir2")
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        #protExp._createFileNamesTemplates()
        dirName = protExp.filesPath.get()
        nameVolume = os.path.join(dirName, protExp.VOLUMENAME)
        nameFsc = os.path.join(dirName, ProtExportEMDB.FSC % 0)
        nameAtomStruct = os.path.join(dirName, protExp.COORDINATEFILENAME)
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameFsc))
        self.assertTrue(os.path.exists(nameAtomStruct))

        #Chek if the files have the correct data
        orig_list_x, orig_list_y = protExp.exportFSC.get().getData()
        fo = open(nameFsc, "rU")
        saved_x=[]
        orig_x=[]
        count=0
        for line in fo:
            if line[0:3]=='<x>':
                saved_x.append(int(float(line[3:-5])*1000))
                orig_x.append(int(orig_list_x[count]*1000))
                count=count+1

        self.assertListEqual(orig_x, saved_x)

    def test_volume2_AtomStructPDB(self):
        volume = self._import_volume_icosahedron()
        vol1 = self._import_volume_half1()
        vol2 = self._import_volume_half2()
        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(vol1)
        prot.referenceVolume.set(vol2)
        self.launchProtocol(prot)

        protExp = self.newProtocol(ProtExportEMDB)
        protExp.exportVolume.set(volume)
        protExp.exportFSC.set(prot.outputFSC)

        # run Chimera operate to get a PDB file
        extraCommands = ""
        extraCommands += "runCommand('scipionwrite model #1')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': self._importAtomStructPDB()
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        protChimera.setObjLabel('chimera operate\n volume and PDB\nHUMAN_hemoglobin')
        self.launchProtocol(protChimera)
        protExp.exportAtomStruct.set(protChimera.outputPdb_01)

        args = {
                'inputVolume': volume,
                'threshold': 1,
                'doMorphological': True,
                }
        protMask = self.newProtocol(XmippProtCreateMask3D, **args)
        protMask.setObjLabel('create 3d mask\nicosahedron\n')
        self.launchProtocol(protMask)
        protExp.exportMask.set(protMask.outputMask)

        protExp.filesPath.set(os.getcwd() + "/dir3")
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        #protExp._createFileNamesTemplates()
        dirName = protExp.filesPath.get()
        nameVolume = os.path.join(dirName, protExp.VOLUMENAME)
        nameHalfVolume1 = os.path.join(dirName, protExp.HALF1NAME)
        nameHalfVolume2 = os.path.join(dirName, protExp.HALF2NAME)
        nameFsc         = os.path.join(dirName, ProtExportEMDB.FSC % 0)
        nameAtomStruct  = os.path.join(dirName, protExp.COORDINATEFILENAME)
        nameMask        = os.path.join(dirName, protExp.MASKFILENAME)
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameHalfVolume1))
        self.assertTrue(os.path.exists(nameHalfVolume2))
        self.assertTrue(os.path.exists(nameFsc))
        self.assertTrue(os.path.exists(nameAtomStruct))
        self.assertTrue(os.path.exists(nameMask))

        #Check if the files have the correct data
        orig_list_x, orig_list_y = protExp.exportFSC.get().getData()
        fo = open(nameFsc, "rU")
        saved_x=[]
        orig_x=[]
        count=0
        for line in fo:
            if line[0:3]=='<x>':
                saved_x.append(int(float(line[3:-5])*1000))
                orig_x.append(int(orig_list_x[count]*1000))
                count=count+1
        self.assertListEqual(orig_x, saved_x)

    def _import_volume_half1(self):

        args = {'filesPath': self.VOLMAPNAMEHALF1,
                'setHalfMaps': False,
                'setOrigCoord': False,
                'samplingRate': self.SAMPLINGRATE,
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot1 = self.newProtocol(ProtImportVolumes, **args)
        prot1.setObjLabel('import vol\n icosahedron half1')
        self.launchProtocol(prot1)
        volume = prot1.outputVolume
        return volume

    def _import_volume_half2(self):

        args = {'filesPath': self.VOLMAPNAMEHALF2,
                'setHalfMaps': False,
                'setOrigCoord': False,
                'samplingRate': self.SAMPLINGRATE,
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot2 = self.newProtocol(ProtImportVolumes, **args)
        prot2.setObjLabel('import vol\n icosahedron half2')
        self.launchProtocol(prot2)
        volume = prot2.outputVolume
        return volume

