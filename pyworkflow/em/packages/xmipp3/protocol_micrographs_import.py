# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This module contains the protocol to import micrographs from xmipp projects.
"""

import os
import xmipp
from pyworkflow.protocol.params import FileParam, IntParam, FloatParam, String
from pyworkflow.em.protocol import ProtImport
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import copyFile
from convert import xmippToLocation



class ProtXmippMicsImport(ProtImport):
    """    
    Protocol to import micrographs from existing Xmipp runs.
    """
    _label = 'import'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputXmipp', FileParam, 
                      label="Input micrograph metadata file",  
                      help='Select the input metadata XMD file from a Xmipp run.'
                           'Also the *microscope.xmd and *acquisition_info.xmd files '
                           'should be present.')
        form.addParam('magnification', IntParam, default=50000,
           label=Message.LABEL_MAGNI_RATE)
        form.addParam('ampContrast', FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        
        
        #TODO: Add an option to allow the user to decide if copy binary files or not        
            
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        self._insertFunctionStep('createMdTmpFiles', self.inputXmipp.get())
        self._insertFunctionStep('createOutputStep', self.magnification.get(), self.ampContrast.get())
        
    #--------------------------- STEPS functions --------------------------------------------
    def createMdTmpFiles(self, dataFile):
        from convert import locationToXmipp
        path = os.path.dirname(dataFile)
        self.micSetMd = self._getExtraPath("input_micrographs.xmd")
        
        microscopeFn = os.path.join(path, "microscope.xmd")
        self.microscopeMd = self._getExtraPath("microscope.xmd")
        
        acqFn = os.path.join(path, "acquisition_info.xmd")
        newAcqMd = self._getExtraPath("acquisition_info.xmd")
        
        copyFile(dataFile, self.micSetMd)
        copyFile(microscopeFn, self.microscopeMd)
        copyFile(acqFn, newAcqMd)
        
        md = xmipp.MetaData(self.micSetMd)
        self.hasCtf = False
        for objId in md:
            index, oldFn = xmippToLocation(md.getValue(xmipp.MDL_MICROGRAPH, objId))
            newFn = self._getMicPath(path, oldFn)
            micLocation = locationToXmipp(index, newFn)
            md.setValue(xmipp.MDL_MICROGRAPH, micLocation, objId)
            
            if md.getValue(xmipp.MDL_CTF_MODEL, objId):
                self.hasCtf = True
                for attr in [xmipp.MDL_CTF_MODEL, xmipp.MDL_PSD, xmipp.MDL_PSD_ENHANCED, xmipp.MDL_IMAGE1, xmipp.MDL_IMAGE2]:
                    oldAttr = md.getValue(attr, objId)
                    newAttr = self._changePath(path, oldAttr)
                    md.setValue(attr, newAttr, objId)
                
        md.write(self.micSetMd)
    
    
    def createOutputStep(self, magnification, ampContrast):
        from convert import readSetOfMicrographs, readCTFModel
        self.info('Creating the set of micrographs...')
        
        # Create the set of micrographs
        micSet = self._createSetOfMicrographs()
        
        # Setting Acquisition properties
        acquisition = micSet.getAcquisition()
        acMd = xmipp.MetaData(self.microscopeMd)
        
        voltage = acMd.getValue(xmipp.MDL_CTF_VOLTAGE, 1)
        sphericalAberration = acMd.getValue(xmipp.MDL_CTF_CS, 1)
        
        acquisition.setVoltage(voltage)
        acquisition.setSphericalAberration(sphericalAberration)
        acquisition.setAmplitudeContrast(ampContrast)
        acquisition.setMagnification(magnification)
        
        readSetOfMicrographs(self.micSetMd, micSet)
        self._defineOutputs(outputMicrographs=micSet)
        
        md = xmipp.MetaData(self.micSetMd)
        if self.hasCtf:
            ctfSet = self._createSetOfCTF()
            ctfSet.setMicrographs(micSet)

            for mic in micSet:
            
                for objId in md:
                    _, xmippfileName = xmippToLocation(md.getValue(xmipp.MDL_MICROGRAPH, objId))
                    
                    if mic.getFileName() == xmippfileName:
                        ctfparam = md.getValue(xmipp.MDL_CTF_MODEL, objId)
                        ctfModel = readCTFModel(ctfparam, mic)
                        
                        ctfModel._psdFile = String(md.getValue(xmipp.MDL_PSD, objId))
                        ctfModel._xmipp_enhanced_psd = String(md.getValue(xmipp.MDL_PSD_ENHANCED, objId))
                        ctfModel._xmipp_ctfmodel_quadrant = String(md.getValue(xmipp.MDL_IMAGE1, objId))
                        ctfModel._xmipp_ctfmodel_halfplane = String(md.getValue(xmipp.MDL_IMAGE2, objId))
                        
                        break
                ctfSet.append(ctfModel)
            self._defineOutputs(outputCTF=ctfSet)
            self._defineCtfRelation(micSet, ctfSet)

    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
    
    def _getMicPath(self,path, oldFn):
        
        if os.path.isabs(oldFn):
            newFn = oldFn
        else:
            newFn = self._changePath(path, oldFn)
        return newFn
    
    def _changePath(self, path, oldPath):
        
        
        relPathList = path.split(os.sep)[:-3]
        relPath = "/" + os.path.join(*relPathList)
        newPath = os.path.abspath(os.path.join(relPath, oldPath))
        return newPath
#     
#     def _parseCommand(self, optimiserFile):
#         """ Read from the optimiser.star the second line which should 
#         contain the exact command line used to launch relion and 
#         grap the parameters in a dictionary way. """
#         opts = {}
#         self.info('Parsing parameters from optimiser file: %s' % optimiserFile)
#         f = open(optimiserFile)
#         for line in f:
#             if '--angpix' in line:
#                 parts = line.split()
#                 for i, p in enumerate(parts):
#                     if p.startswith('--'): # it an option
#                         opts[p] = parts[i+1] # take what follows the option
#                 break
#         f.close()
#         return opts
#     
#     def _findImagesPath(self, starFile):
#         from convert import findImagesPath
#         self._imagesPath = findImagesPath(starFile)
#         
#         if self._imagesPath:
#             self.info('Images path found in: %s' % self._imagesPath)
#         else:
#             self.warning('Images binaries not found!!!')
#         
#     def _preprocessRow(self, imgRow):
#         from convert import setupCTF, prependToFileName
#         if self._imagesPath is not None:
#             prependToFileName(imgRow, self._imagesPath)
#         
#         setupCTF(imgRow, self.samplingRate.get())