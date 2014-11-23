# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
In this module are protocol base classes related to EM imports of Micrographs, Particles, Volumes...
"""

from os.path import exists, basename
from glob import glob

from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import (PathParam, FloatParam, BooleanParam, FileParam,
                                        EnumParam, IntParam, StringParam, PointerParam,
                                        LabelParam, LEVEL_ADVANCED)
from pyworkflow.utils.path import expandPattern, createLink, copyFile
from pyworkflow.em.data import Volume, PdbFile

from base import ProtImport



class ProtImportVolumes(ProtImport):
    """Protocol to import a set of volumes to the project"""
    _label = 'import volumes'
    
    def __init__(self, **args):
        ProtImport.__init__(self, **args)         
       
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('pattern', PathParam, label=Message.LABEL_PATTERN)
        form.addParam('samplingRate', FloatParam,
                   label=Message.LABEL_SAMP_RATE)
    
    def _insertAllSteps(self):
        self._insertFunctionStep('importVolumes', self.getPattern(), self.samplingRate.get())
        
    def createVolume(self, volumePath):
        """ Copy the volume to WorkingDir and create
        the volumen object.
        """
        dst = self._getPath(basename(volumePath))            
        createLink(volumePath, dst)
        vol = Volume()
        if dst.endswith('.mrc'):
            dst += ':mrc'
        vol.setFileName(dst)
        vol.setSamplingRate(self.samplingRate.get())
        return vol
    
    def importVolumes(self, pattern, samplingRate):
        """ Copy volumes matching the filename pattern
        Register other parameters.
        """
        n = self._getNumberFilePaths(pattern)
        filePaths = self._getFilePaths(pattern)
        
        if n == 0:
            raise Exception(Message.ERROR_IMPORT_VOL)
        elif n == 1:
            volume = self.createVolume(filePaths[0])
            self._defineOutputs(outputVolume=volume)
        else:
            # Create a set of volumes
            volSet = self._createSetOfVolumes()
            volSet.setSamplingRate(self.samplingRate.get())
#             filePaths.sort()
            for f in filePaths:
                volSet.append(self.createVolume(f))
            self._defineOutputs(outputVolumes=volSet)

    def getPattern(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
        """
        pattern = expandPattern(self.pattern.get())    
        return pattern  
        
    def getFiles(self):
        pattern = self.getPattern()
        n = self._getNumberFilePaths(pattern)
        
        if n == 1:
            return self.outputVolume.getFileName()
        else:
            return self.outputVolumes.getFiles()
    
    def _getFilePaths(self, pattern):
        """ Return a sorted list with the paths of files"""
        filePaths = glob(pattern)
        filePaths.sort()
        
        return filePaths
    
    def _getNumberFilePaths(self, pattern):
        """ Return the number of files""" 
        filePaths = self._getFilePaths(pattern)
        n = len(filePaths)
        return n

    def _summary(self):
        summary = []
        pattern = self.getPattern()
        n = self._getNumberFilePaths(pattern)
        
        if n == 1:
            if not hasattr(self, 'outputVolume'):
                summary.append("Output volume not ready yet.") 
            else:
                summary.append("Import of %d volumes from %s" % (1, self.getPattern()))
                summary.append("Sampling rate : %f" % self.samplingRate.get())
        else:
            if not hasattr(self, 'outputVolumes'):
                summary.append("Output volumes not ready yet.") 
            else:
                summary.append("Import of %d volumes from %s" % (n, self.getPattern()))
                summary.append("Sampling rate : %f" % self.samplingRate.get())
        
        return summary
    

class ProtImportPdb(ProtImport):
    """ Protocol to import a set of pdb volumes to the project"""
    _label = 'import pdb volumes'
    
    def __init__(self, **args):
        ProtImport.__init__(self, **args)         
       
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('pdbPath', FileParam, 
                      label="PDB file",
                      help='Specify a path to desired PDB structure.')
         
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep', self.pdbPath.get())
        
    def createOutputStep(self, pdbPath):
        """ Copy the PDB structure and register the output object.
        """
        if not exists(pdbPath):
            raise Exception("PDB not found at *%s*" % pdbPath)
        
        baseName = basename(pdbPath)
        localPath = self._getExtraPath(baseName)
        copyFile(pdbPath, localPath)
        pdb = PdbFile()
        pdb.setFileName(localPath)
        self._defineOutputs(outputPdb=pdb)

    def _summary(self):
        summary = ['PDB file imported from *%s*' % self.pdbPath.get()]

        return summary
    
    def _validate(self):
        errors = []
        if not exists(self.pdbPath.get()):
            errors.append("PDB not found at *%s*" % self.pdbPath.get())
        #TODO: maybe also validate that if exists is a valid PDB file 
        return errors
    
    
