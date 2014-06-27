# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano         (ldelcano@cnb.csic.es)
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
This module contains protocols classes related to Random Conical Tilt.
"""

import sys
from glob import glob
from protocol_import import ProtImport
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import copyFile, basename, createLink, expandPattern
from pyworkflow.protocol.params import PathParam, FloatParam, BooleanParam, EnumParam, IntParam
from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data_tiltpairs import MicrographsTiltPair


class ProtImportMicrographsTiltPairs(ProtImport):
    """Protocol to import a pair of tilted micrographs in the project"""
    
    _className = 'Micrograph'
    _label = 'import tilted micrographs'
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('patternUntilted', PathParam, label=Message.LABEL_PATTERNU,
                      help=Message.TEXT_PATTERN)
        form.addParam('patternTilted', PathParam, label=Message.LABEL_PATTERNT,
                      help=Message.TEXT_PATTERN)
        form.addParam('copyToProj', BooleanParam, label=Message.LABEL_COPYFILES, default=False)
        form.addParam('voltage', FloatParam, default=200,
                   label=Message.LABEL_VOLTAGE)
        form.addParam('sphericalAberration', FloatParam, default=2.26,
                   label=Message.LABEL_SPH_ABERRATION)
        form.addParam('ampContrast', FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        form.addParam('samplingRateMode', EnumParam, default=SAMPLING_FROM_IMAGE,
                   label=Message.LABEL_SAMP_MODE,
                   choices=[Message.LABEL_SAMP_MODE_1, Message.LABEL_SAMP_MODE_2])
        form.addParam('samplingRate', FloatParam, default=1, 
                   label=Message.LABEL_SAMP_RATE,
                   condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE)
        form.addParam('magnification', IntParam, default=50000,
                   label=Message.LABEL_MAGNI_RATE,
                   condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER)
        form.addParam('scannedPixelSize', FloatParam, default=7.0,
                   label=Message.LABEL_SCANNED,
                   condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importMicrographs(self, pattern, suffix, voltage, sphericalAberration, amplitudeContrast):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        from pyworkflow.em import findClass
        filePaths = glob(expandPattern(pattern))
        
        imgSet = self._createSetOfMicrographs(suffix)
        acquisition = imgSet.getAcquisition()
        # Setting Acquisition properties
        acquisition.setVoltage(voltage)
        acquisition.setSphericalAberration(sphericalAberration)
        acquisition.setAmplitudeContrast(amplitudeContrast)
        
        # Call a function that should be implemented by each subclass
        self._setOtherPars(imgSet)
        
        outFiles = [imgSet.getFileName()]
        imgh = ImageHandler()
        img = imgSet.ITEM_TYPE()
        n = 1
        size = len(filePaths)
        
        filePaths.sort()
        
        for i, fn in enumerate(filePaths):
#             ext = os.path.splitext(basename(f))[1]
            dst = self._getExtraPath(basename(fn))
            if self.copyToProj:
                copyFile(fn, dst)
            else:
                createLink(fn, dst)
            
            if n > 1:
                for index in range(1, n+1):
                    img.cleanObjId()
                    img.setFileName(dst)
                    img.setIndex(index)
                    imgSet.append(img)
            else:
                img.cleanObjId()
                img.setFileName(dst)
                imgSet.append(img)
            outFiles.append(dst)
            
            sys.stdout.write("\rImported %d/%d" % (i+1, size))
            sys.stdout.flush()
            
        print "\n"
        
        imgSet.write()
        
        
        return imgSet
    
    def createTiltedPairs(self):
        args = {} 
        micsTiltPair = MicrographsTiltPair()
        micsU = self.importMicrographs(self.patternUntilted.get(), 'untilted',
                                 self.voltage.get(), self.sphericalAberration.get(), self.ampContrast.get())
        micsT = self.importMicrographs(self.patternTilted.get(), 'tilted',
                                 self.voltage.get(), self.sphericalAberration.get(), self.ampContrast.get())
        
        micsTiltPair.setUntilted(micsU)
        micsTiltPair.setTilted(micsT)
        output = self._getOutput()
        args[output] = micsTiltPair
        self._defineOutputs(**args)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        if (not self.patternUntilted.get() or not self.patternTilted.get()):
            errors.append(Message.ERROR_PATTERN_EMPTY)
        else:
            filePathsUntilted = glob(expandPattern(self.patternUntilted.get()))
            filePathsTilted = glob(expandPattern(self.patternTilted.get()))
        
            if (len(filePathsUntilted) == 0 or len(filePathsTilted) == 0):
                errors.append(Message.ERROR_PATTERN_FILES)

        return errors
    
    def _summary(self):
        summary = []
        output = self._getOutput()
        if not hasattr(self, output):
            summary.append("Output " + output + " not ready yet.") 
            if self.copyToProj:
                summary.append("*Warning*: Import step could take a long time due to the images are copying in the project.")
        else:
            summary.append("Import of %d " % getattr(self, output).getTilted().getSize() + self._className + "s tilted from %s" % self.patternTilted.get())
            summary.append("Import of %d " % getattr(self, output).getUntilted().getSize() + self._className + "s untilted from %s" % self.patternUntilted.get())
            summary.append("Sampling rate : %0.2f A/px" % self.samplingRate.get())
        
        return summary
    
    def _methods(self):
        methods = []
        output = self._getOutput()
        if hasattr(self, output):
            methods.append("%d " % getattr(self, output).getTilted().getSize() + self._className + "s tilted has been imported")
            methods.append("%d " % getattr(self, output).getUntilted().getSize() + self._className + "s untilted has been imported")
            methods.append("with a sampling rate : %0.2f A/px" % self.samplingRate.get())
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getFiles(self):
        return getattr(self, self._getOutput().getFiles())
    
    def _getOutput(self):
        return "outputMicrographsTiltPair"

    
    def _insertAllSteps(self):
        self._insertFunctionStep('createTiltedPairs')
        
    def _setOtherPars(self, micSet):
        micSet.getAcquisition().setMagnification(self.magnification.get())
        
        if self.samplingRateMode == SAMPLING_FROM_IMAGE:
            micSet.setSamplingRate(self.samplingRate.get())
        else:
            micSet.setScannedPixelSize(self.scannedPixelSize.get())
            
    