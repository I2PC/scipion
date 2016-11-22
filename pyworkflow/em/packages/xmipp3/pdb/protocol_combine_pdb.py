# **************************************************************************
# *
# * Authors:  Mohsen Kazemi (mkazemi@cnb.csic.es)
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

import os, ftplib, gzip
import sys

import pyworkflow.protocol.params as params
import pyworkflow.em as em


class XmippProtCombinePdb(em.ProtImportFiles):
    """ Combine two PDB files into a single one.  """
    _label = 'combine PDBs'
    IMPORT_FROM_ID = 0
    IMPORT_OBJ = 1
    IMPORT_FROM_FILES = 2 
       
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        """ Define the parameters that will be input for the Protocol.
        This definition is also used to generate automatically the GUI.
        """
        form.addSection(label='Input')
        form.addParam('inputPdbData1', params.EnumParam,
                      choices=['id', 'object', 'file'],
                      label="Retrieve the 1st PDB from",
                      default=self.IMPORT_FROM_ID,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Retrieve PDB data from server, '
                           'use a pdb Object, or a local file')
        form.addParam('pdbId1', params.StringParam, 
                      condition='inputPdbData1 == IMPORT_FROM_ID',
                      label="The 1st Pdb Id ", allowsNull=True,
                      help='Type a pdb Id (four alphanumeric characters).')
        form.addParam('pdbObj1', params.PointerParam, pointerClass='PdbFile',
                      label="Input the 1st pdb ", 
                      condition='inputPdbData1 == IMPORT_OBJ', allowsNull=True,
                      help='Specify a pdb object.')
        form.addParam('pdbFile1', params.FileParam,
                      label="File path of the 1st pdb", 
                      condition='inputPdbData1 == IMPORT_FROM_FILES', 
                      allowsNull=True,
                      help='Specify a path to desired PDB structure.')        
        form.addParam('inputPdbData2', params.EnumParam, 
                      choices=['id', 'object', 'file'],
                      label="Retrieve the 2nd PDB from", 
                      default=self.IMPORT_FROM_ID,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Retrieve PDB data from server, '
                           'use a pdb Object, or a local file')
        form.addParam('pdbId2', params.StringParam, 
                      condition='inputPdbData2 == IMPORT_FROM_ID',
                      label="The 2nd Pdb Id ", allowsNull=True,
                      help='Type a pdb Id (four alphanumeric characters).')
        form.addParam('pdbObj2', params.PointerParam, pointerClass='PdbFile',
                      label="Input the 2nd pdb ", 
                      condition='inputPdbData2 == IMPORT_OBJ', allowsNull=True,
                      help='Specify a pdb object.')
        form.addParam('pdbFile2', params.FileParam,
                      label="File path of the 2nd pdb", 
                      condition='inputPdbData2 == IMPORT_FROM_FILES', 
                      allowsNull=True,
                      help='Specify a path to desired PDB structure.')
        
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ In this function the steps that are going to be executed should
        be defined. Two of the most used functions are:
        _insertFunctionStep or _insertRunJobStep
        """
        if self.inputPdbData1 == self.IMPORT_FROM_ID:
            self._insertFunctionStep('pdbDownloadStep')
        if self.inputPdbData2 == self.IMPORT_FROM_ID:
            self._insertFunctionStep('pdbDownloadStep')
        self._insertFunctionStep('combinePdbStep')        
    
    #--------------------------- STEPS functions --------------------------------------------
    def pdbDownloadStep(self):
        """Download all pdb files in file_list and unzip them."""
        if self.inputPdbData1 == self.IMPORT_FROM_ID:
            em.downloadPdb(self.pdbId1.get(), self._getPdbFileName1(), 
                           self._log)
        if self.inputPdbData2 == self.IMPORT_FROM_ID:
            em.downloadPdb(self.pdbId2.get(), self._getPdbFileName2(), 
                           self._log)
        
    def combinePdbStep(self):
        """ Although is not mandatory, usually is used by the protocol to
        register the resulting outputs in the database.
        """
        pdbFn1 = self._getPdbFileName1()
        pdbFn2 = self._getPdbFileName2()
        fnOut = self._getPdbOutName()
        
        self.info("The 1st input file: " + pdbFn1)
        self.info("The 2nd input file: " + pdbFn2)
        self.info("Output file: " + fnOut)
        
        fh = open(fnOut, "a")
        for line in open(pdbFn1,"r"): 
            charac = line.split()            
            if charac[0] <> 'END':          
                fh.write(line)
                fh.write("\n")
        for line in open(pdbFn2,"r"):
            charac = line.split()            
            if charac[0] == 'HETATM' or charac[0] == 'ATOM' or charac[0] == 'END':            
                fh.write(line)
                fh.write("\n")          
            
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        """ Even if the full set of parameters is available, 
        this function provides summary information about an specific run.
        """ 
        summary = [ ] 
        # Add some lines of summary information
        if self.inputPdbData1 == self.IMPORT_FROM_ID:
            summary.append("The 1st input PDB ID: %s" % self.pdbId1.get())
        elif self.inputPdbData1 == self.IMPORT_OBJ:
            summary.append("The 1st input PDB File: %s" % self.pdbObj1.get().getFileName())
        else:
            summary.append("The 1st input PDB File: %s" % self.pdbFile1.get())
        if self.inputPdbData2 == self.IMPORT_FROM_ID:
            summary.append("The 2nd input PDB ID: %s" % self.pdbId2.get())
        elif self.inputPdbData2 == self.IMPORT_OBJ:
            summary.append("The 2nd input PDB File: %s" % self.pdbObj2.get().getFileName())
        else:
            summary.append("The 2nd input PDB File: %s" % self.pdbFile2.get())
        
        return summary
      
    def _validate(self):
        """ The function of this hook is to add some validation before 
        the protocol is launched to be executed. It should return a 
        list of errors. If the list is empty the protocol can be executed.
        """
        errors = []
        if self.inputPdbData1 == self.IMPORT_FROM_ID:
            lenStr = len(self.pdbId1.get())
            if lenStr <> 4:
                errors = ["Pdb id is composed only by four alphanumeric characters"]
        if self.inputPdbData2 == self.IMPORT_FROM_ID:
            lenStr = len(self.pdbId2.get())
            if lenStr <> 4:
                errors = ["Pdb id is composed only by four alphanumeric characters"]
        
        return errors
    
    #--------------------------- UTLIS functions --------------------------------------------
    def _getPdbFileName1(self):
        if self.inputPdbData1 == self.IMPORT_FROM_ID:
            return self._getExtraPath('%s.pdb' % self.pdbId1.get())
        elif self.inputPdbData1 == self.IMPORT_OBJ:
            return self.pdbObj1.get().getFileName()
        else:
            return self.pdbFile1.get()
    
    def _getPdbFileName2(self):
        if self.inputPdbData2 == self.IMPORT_FROM_ID:
            return self._getExtraPath('%s.pdb' % self.pdbId2.get())
        elif self.inputPdbData2 == self.IMPORT_OBJ:
            return self.pdbObj2.get().getFileName()
        else:
            return self.pdbFile2.get() 
    
    def _getPdbOutName(self):
        return self._getExtraPath("combined.pdb")