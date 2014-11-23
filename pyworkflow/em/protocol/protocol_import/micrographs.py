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
        group = self.acquisitionGroup
        group.addParam('samplingRateMode', EnumParam, 
                       choices=[Message.LABEL_SAMP_MODE_1, Message.LABEL_SAMP_MODE_2],
                       default=SAMPLING_FROM_IMAGE,
                       label=Message.LABEL_SAMP_MODE,
                       help=Message.TEXT_SAMP_MODE)
        group.addParam('samplingRate', FloatParam,  
                       condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE, 
                       label=Message.LABEL_SAMP_RATE,
                       help=Message.TEXT_SAMP_RATE)
        group.addParam('scannedPixelSize', FloatParam, default=7.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER,
                       label=Message.LABEL_SCANNED,
                       help='')
        
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

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        return ['emx', 'xmipp3']
    
    def _defineBasicParams(self, form):
        ProtImportMicBase._defineBasicParams(self, form)
        
        form.addParam('micrographsEMX', FileParam,
              condition = '(importFrom == %d)' % self.IMPORT_FROM_EMX,
              label='Input EMX file',
              help="Select the EMX file containing micrographs information.\n"
                   "See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]")
        
        form.addParam('micrographsMd', FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Micrographs metadata file',
                      help="Select the micrographs Xmipp metadata file.\n"
                           "It is usually a _micrograph.xmd_ file result\n"
                           "from import, preprocess or downsample protocols.")

    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        
        if importFrom == self.IMPORT_FROM_FILES:
            ProtImportMicBase._insertAllSteps(self)
        elif importFrom == self.IMPORT_FROM_EMX:
            self._insertFunctionStep('importMicrographsStep', importFrom,
                                     self.micrographsEMX.get())
        elif importFrom == self.IMPORT_FROM_XMIPP3:
            self._insertFunctionStep('importMicrographsStep', importFrom,
                                     self.micrographsMd.get())
            
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        if self.importFrom == self.IMPORT_FROM_EMX:
            from pyworkflow.em.packages.emxlib import EmxImport
            self.importFilePath = self.micrographsEMX.get('').strip()
            return EmxImport(self, self.importFilePath)
        elif self.importFrom == self.IMPORT_FROM_XMIPP3:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport
            self.importFilePath = self.micrographsMd.get('').strip()
            return XmippImport(self, self.micrographsMd.get())
        else:
            self.importFilePath = ''
            return None       
        
    def importMicrographsStep(self, importFrom, *args):
        ci = self.getImportClass()
        ci.importMicrographs()
        
        size = self.outputMicrographs.getSize()
        sampling = self.outputMicrographs.getSamplingRate()
        summary = "Import of *%d* Micrographs from %s file:\n" % (size, self.getEnumText('importFrom'))
        summary += self.importFilePath + '\n'
        summary += "Sampling rate : *%0.2f* A/px\n\n" % sampling
        if self.copyFiles:
            summary += 'WARNING: Binary files copied into project (extra disk space)'
        self.summaryVar.set(summary)
        
    def loadAcquisitionInfo(self):
        ci = self.getImportClass()
        if exists(self.importFilePath):
            return ci.loadAcquisitionInfo()
        else:
            return None
        
    def _validate(self):
        ci = self.getImportClass()
        if ci is None:
            return ProtImportMicBase._validate(self)
        else:
            return ci.validate()
    
    def _summary(self):
        if self.importFrom == self.IMPORT_FROM_FILES:
            return ProtImportMicBase._summary(self)
        else:
            return [self.summaryVar.get()]
        

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
        
    def setSamplingRate(self, movieSet):
        ProtImportMicBase.setSamplingRate(self, movieSet)
        movieSet.setGain(self.gainFile.get())
                    
