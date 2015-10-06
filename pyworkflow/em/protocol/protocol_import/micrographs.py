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

from os.path import exists

from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import FileParam, FloatParam, EnumParam
from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER

from images import ProtImportImages



class ProtImportMicBase(ProtImportImages):
    """ Just to have a base class to both 
    ProtImportMicrographs and ProtImportMovies
    """
    _checkStacks = False
    
    def _defineParams(self, form):
        ProtImportImages._defineParams(self, form)
        
    def _defineAcquisitionParams(self, form):
        group = ProtImportImages._defineAcquisitionParams(self, form)
        group.addParam('samplingRateMode', EnumParam, 
                       choices=[Message.LABEL_SAMP_MODE_1, Message.LABEL_SAMP_MODE_2],
                       default=SAMPLING_FROM_IMAGE,
                       label=Message.LABEL_SAMP_MODE,
                       help=Message.TEXT_SAMP_MODE)
        group.addParam('samplingRate', FloatParam,  default=1.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE, 
                       label=Message.LABEL_SAMP_RATE,
                       help=Message.TEXT_SAMP_RATE)
        group.addParam('scannedPixelSize', FloatParam, default=7.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER,
                       label=Message.LABEL_SCANNED,
                       help='')
        return group
        
    def setSamplingRate(self, micSet):
        """ Set the sampling rate to the given set. """
        if self.samplingRateMode == SAMPLING_FROM_IMAGE:
            micSet.setSamplingRate(self.samplingRate.get())
        else:
            micSet.setScannedPixelSize(self.scannedPixelSize.get())
        
    
class ProtImportMicrographs(ProtImportMicBase):
    """Protocol to import a set of micrographs to the project"""
    _label = 'import micrographs'
    _outputClassName = 'SetOfMicrographs' 
    
    IMPORT_FROM_EMX = 1
    IMPORT_FROM_XMIPP3 = 2
    IMPORT_FROM_SCIPION = 3

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        choices = ProtImportImages._getImportChoices(self)
        return choices + ['emx', 'xmipp3', 'scipion']
    
    def _defineImportParams(self, form):
        """ Just redefine to put some import parameters
        before the acquisition related parameters.
        """
        form.addParam('emxFile', FileParam,
              condition = '(importFrom == %d)' % self.IMPORT_FROM_EMX,
              label='Input EMX file',
              help="Select the EMX file containing micrographs information.\n"
                   "See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]")
        
        form.addParam('mdFile', FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Micrographs metadata file',
                      help="Select the micrographs Xmipp metadata file.\n"
                           "It is usually a _micrograph.xmd_ file result\n"
                           "from import, preprocess or downsample protocols.")
        
        form.addParam('sqliteFile', FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_SCIPION,
                      label='Micrographs sqlite file',
                      help="Select the micrographs sqlite file.\n")
    
    #--------------------------- INSERT functions ---------------------------------------------------
    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        ci = self.getImportClass()
        
        if ci is None:
            ProtImportMicBase._insertAllSteps(self)
        else:
            self._insertFunctionStep('importMicrographsStep', importFrom,
                                     self.importFilePath)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importMicrographsStep(self, importFrom, *args):
        ci = self.getImportClass()
        ci.importMicrographs()
        
        summary = "Import from *%s* file:\n" % self.getEnumText('importFrom')
        summary += self.importFilePath + '\n'
        
        if self.hasAttribute('outputParticles'):
            particles = self.outputParticles
            summary += '   Particles: *%d* (ctf=%s, alignment=%s)\n' % (particles.getSize(),
                                                                        particles.hasCTF(),
                                                                        particles.getAlignment())
                                                                      
        if self.hasAttribute('outputCoordinates'): # EMX files can contain only Coordinates information
            summary += '   Coordinates: *%d* \n' % (self.outputCoordinates.getSize())
            
        if self.hasAttribute('outputMicrographs'): # EMX files can contain only Coordinates information
            summary += '   Micrographs: *%d* \n' % (self.outputMicrographs.getSize())
        
        if self.copyFiles:
            summary += '\n_WARNING_: Binary files copied into project (extra disk space)'
            
        self.summaryVar.set(summary)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        from pyworkflow.em.convert import ImageHandler
        ci = self.getImportClass()
        if ci is None:
            errors = ProtImportMicBase._validate(self)
            for micFn, _ in self.iterFiles():
                imgh = ImageHandler()
                if imgh.isImageFile(micFn):
                    _, _, z, n = imgh.getDimensions(micFn)
                    if n > 1 or z > 1:
                        errors.append("The protocol not support micrographs stored in stacks. "
                                      "If you want to obtain your micrographs individually, "
                                      "you can run the following command:\n"
                                      "scipion run scipion_directory/scripts/split_stacks.py --files *your files* --ext *extension*")
                # JMRT: only check the first image, for large dataset
                # even reading the header can take a while
                break 
            return errors
            
        else:
            return ci.validateMicrographs()
    
    def _summary(self):
        if self.importFrom == self.IMPORT_FROM_FILES:
            return ProtImportMicBase._summary(self)
        else:
            return [self.summaryVar.get('No summary information.')]
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        if self.importFrom == self.IMPORT_FROM_EMX:
            from pyworkflow.em.packages.emxlib import EmxImport
            self.importFilePath = self.emxFile.get('').strip()
            return EmxImport(self, self.importFilePath)
        elif self.importFrom == self.IMPORT_FROM_XMIPP3:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport
            self.importFilePath = self.mdFile.get('').strip()
            return XmippImport(self, self.mdFile.get())
        elif self.importFrom == self.IMPORT_FROM_SCIPION:
            from dataimport import ScipionImport
            self.importFilePath = self.sqliteFile.get('').strip()
            return ScipionImport(self, self.importFilePath) 
        else:
            self.importFilePath = ''
            return None       
    
    def loadAcquisitionInfo(self):
        ci = self.getImportClass()
        if exists(self.importFilePath):
            return ci.loadAcquisitionInfo()
        else:
            return None


class ProtImportMovies(ProtImportMicBase):
    """Protocol to import a set of movies (from direct detector cameras) to the project"""
    _label = 'import movies'
    _outputClassName = 'SetOfMovies'
        
    def _defineParams(self, form):
        ProtImportMicBase._defineParams(self, form)    
        form.addParam('gainFile', FileParam,  
                      label='Gain image', 
                      help='A gain reference related to a set of movies'
                           ' for gain correction')
        form.addParam('darkFile', FileParam,  
                      label='Dark image', 
                      help='A dark image related to a set of movies')
        
    def setSamplingRate(self, movieSet):
        ProtImportMicBase.setSamplingRate(self, movieSet)
        movieSet.setGain(self.gainFile.get())
        movieSet.setDark(self.darkFile.get())
                    
