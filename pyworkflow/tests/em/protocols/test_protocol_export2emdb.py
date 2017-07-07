# **************************************************************************
# *
# * Authors:    Amaya Jimenez (ajimenez@cnb.csic.es)
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
from pyworkflow.em.protocol import ProtImportVolumes
from pyworkflow.em.protocol.protocol_export import ProtExportEMDB
from pyworkflow.em.packages.xmipp3 import XmippProtMultipleFSCs, XmippProtResolution3D

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

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        cls.protImportHalf1  = cls.runImportVolumes(cls.half1, 3.54,
                                                    'import half1')
        cls.protImportHalf2  = cls.runImportVolumes(cls.half2, 3.54,
                                                    'import half2')


    def test_volume1(self):

        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.protImportHalf1.outputVolume)
        prot.referenceVolume.set(self.protImportHalf2.outputVolume)
        self.launchProtocol(prot)

        protExp = self.newProtocol(ProtExportEMDB)
        protExp.exportVolume.set(self.protImportHalf1.outputVolume)
        protExp.exportFSC.set(prot.outputFSC)

        protExp.filesPath.set(os.getcwd())
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        protExp._createFileNamesTemplates()
        nameVolume = protExp.getFnPath()
        nameFsc = protExp.getFnPath('fsc')
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameFsc))

        #Chek if the files have the correct data
        orig_list_x, orig_list_y = protExp.exportFSC.get().getData()
        fo = open(protExp.getFnPath('fsc'), "rU")
        saved_x=[]
        orig_x=[]
        count=0
        for line in fo:
            if line[0:3]=='<x>':
                saved_x.append(int(float(line[3:-5])*1000))
                orig_x.append(int(orig_list_x[count]*1000))
                count=count+1

        self.assertListEqual(orig_x, saved_x)
