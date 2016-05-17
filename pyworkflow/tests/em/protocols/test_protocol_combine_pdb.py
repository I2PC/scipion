# **************************************************************************
# *
# * Authors:    Mohsen Kazemi (mkazemi@cnb.csic.es)
# *             
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import *
from pyworkflow.em.protocol.protocol_import import ProtImportPdb
from pyworkflow.em.packages.xmipp3.pdb.protocol_combine_pdb import XmippProtCombinePdb

class TestCombinePdb(BaseTest):     
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('nma')
        cls.pdb = cls.dataset.getFile('pdb')
   
    def testXmippPdbCombineFromDb(self):
        print "Run combine two pdbs from database"
        protCombine = self.newProtocol(XmippProtCombinePdb,                                        
                                       pdbId1="3j3i", pdbId2="2mj4")
        self.launchProtocol(protCombine)
        
    def testXmippPdbCombineFromObj(self):
        print "Run combine two pdbs from import"
        protImport = self.newProtocol(ProtImportPdb, 
                                      inputPdbData=ProtImportPdb.IMPORT_FROM_FILES, 
                                      pdbFile=self.pdb)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputPdb.getFileName(), 
                             "There was a problem with the import")        
        protCombine = self.newProtocol(XmippProtCombinePdb, 
                                       inputPdbData1=XmippProtCombinePdb.IMPORT_OBJ, 
                                       inputPdbData2=XmippProtCombinePdb.IMPORT_OBJ)                                       
        protCombine.pdbObj1.set(protImport.outputPdb)
        protCombine.pdbObj2.set(protImport.outputPdb)
        self.launchProtocol(protCombine)
         
    def testXmippPdbCombineFromFn(self):
        print "Run combine two pdbs from file"
        protCombine = self.newProtocol(XmippProtCombinePdb,
                                       inputPdbData1=XmippProtCombinePdb.IMPORT_FROM_FILES, 
                                       inputPdbData2=XmippProtCombinePdb.IMPORT_FROM_FILES, 
                                       pdbFile1=self.pdb, pdbFile2=self.pdb)
        self.launchProtocol(protCombine)
        
