# **************************************************************************
# *
# * Authors:    Amaya Jimenez (ajimenez@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from chimera.protocols import ChimeraProtRigidFit
from pyworkflow.utils import importFromPlugin

XmippProtMultipleFSCs = importFromPlugin('xmipp3.protocols', 'XmippProtMultipleFSCs', doRaise=True)
XmippProtResolution3D = importFromPlugin('xmipp3.protocols', 'XmippProtResolution3D')

class TestExport2EMDB(BaseTest):
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

    @classmethod
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

    def test_volume1_AtomStructCIF(self):

        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.protImportHalf1.outputVolume)
        prot.referenceVolume.set(self.protImportHalf2.outputVolume)
        self.launchProtocol(prot)

        protExp = self.newProtocol(ProtExportEMDB)
        protExp.exportVolume.set(self.protImportHalf1.outputVolume)
        protExp.exportFSC.set(prot.outputFSC)

        protExp.exportAtomStruct.set(self._importAtomStructCIF())

        protExp.filesPath.set(os.getcwd())
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        #protExp._createFileNamesTemplates()
        nameVolume = protExp.VOLUMENAME
        dirName = protExp.filesPath.get()
        nameFsc = os.path.join(dirName, "fsc_%02d.xml" % 0)
        nameAtomStruct = protExp.COORDINATEFILENAME
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

        # run Chimera rigid fit ti get a PDB file
        extraCommands = ""
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "saverefmodel 0')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': self.protImportHalf1.outputVolume,
                'pdbFileToBeRefined': self._importAtomStructCIF()
                }
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and PDB\n save model')
        self.launchProtocol(protChimera)
        protExp.exportAtomStruct.set(protChimera.outputPdb_01)

        protExp.filesPath.set(os.getcwd())
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        #protExp._createFileNamesTemplates()
        nameVolume = protExp.VOLUMENAME
        dirName = protExp.filesPath.get()
        nameFsc = os.path.join(dirName, "fsc_%02d.xml" % 0)
        nameAtomStruct = protExp.COORDINATEFILENAME
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
