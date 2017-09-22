# ***************************************************************************
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es)
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
from itertools import izip

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportCTF, ProtImportMicrographs



class TestImportCTFs(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsGrigorieff = DataSet.getDataSet('grigorieff')

        # First, import a set of micrographs that will be used
        # from all ctf test cases
        cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                      filesPath=cls.dsXmipp.getFile('allMics'),
                                      samplingRate=1.237, voltage=300)
        cls.launchProtocol(cls.protImport)

    def testImportCTFFromXmipp(self):
        protCTF = self.newProtocol(ProtImportCTF,
                                 importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                 filesPath=self.dsXmipp.getFile('ctfsDir'),
                                 filesPattern='*.ctfparam')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import ctfs from xmipp ')
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")

    def testImportCtffind3(self):
        protCTF = self.newProtocol(ProtImportCTF,
                                 importFrom=ProtImportCTF.IMPORT_FROM_GRIGORIEFF,
                                 filesPath=self.dsGrigorieff.getFile('ctffind3'),
                                 filesPattern='BPV*/*txt')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import from ctffind3')
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")
        
    def testImportCtffind4(self):
        protCTF = self.newProtocol(ProtImportCTF,
                                 importFrom=ProtImportCTF.IMPORT_FROM_GRIGORIEFF,
                                 filesPath=self.dsGrigorieff.getFile('ctffind4'),
                                 filesPattern='BPV*/*txt')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import from ctffind4')
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")


    def testImportFromScipion(self):
        ctfSqlite =  self.dsGrigorieff.getFile('ctffind3/ctfs.sqlite')

        protCTF = self.newProtocol(ProtImportCTF,
                                   objLabel='import from scipion',
                                   importFrom=ProtImportCTF.IMPORT_FROM_SCIPION,
                                   filesPath=ctfSqlite)

        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")

    def testImportFromCTffind4Conflict(self):
        micsPath = os.path.abspath(self.proj.getTmpPath('micrographs'))
        ctfsPath = os.path.abspath(self.proj.getTmpPath('ctfs'))
        pwutils.makePath(micsPath)
        pwutils.makePath(ctfsPath)

        micDict = {'BPV_1386': 'mic1',
                   'BPV_1387': 'mic10',
                   'BPV_1388': 'mic2'}

        pattern = 'micrographs/%s.mrc'

        for k, v in micDict.iteritems():
            # Create a micrograph link with a given naming convention
            pwutils.createAbsLink(self.dsXmipp.getFile(pattern % k),
                                  self.proj.getTmpPath(pattern % v))

            pwutils.createAbsLink(self.dsGrigorieff.getFile('ctffind4/%s' % k),
                                  self.proj.getTmpPath('ctfs/%s' % v))

        protCTF = self.newProtocol(ProtImportCTF,
                                   importFrom=ProtImportCTF.IMPORT_FROM_GRIGORIEFF,
                                   filesPath=ctfsPath,
                                   filesPattern='mic*/*txt')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import from ctffind3 - empty')
        self.launchProtocol(protCTF)

        # FIXME: After the execution of the above case, some empty sets
        # are created, if there is no matching, there should not be any output

        # Let's import now the correct names, but with the problematic matching
        protImport2 = self.newProtocol(ProtImportMicrographs,
                                       filesPath=micsPath,
                                       filesPattern='mic*mrc',
                                       samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport2)

        protCTF2 = self.newProtocol(ProtImportCTF,
                                   importFrom=ProtImportCTF.IMPORT_FROM_GRIGORIEFF,
                                   filesPath=ctfsPath,
                                   filesPattern='mic*/*txt')
        protCTF2.inputMicrographs.set(protImport2.outputMicrographs)
        protCTF2.setObjLabel('import from ctffind3 - match')
        self.launchProtocol(protCTF2)

        # The whole set (3 items) should find its corresponding CTF
        self.assertEqual(protCTF2.outputCTF.getSize(), 3)
