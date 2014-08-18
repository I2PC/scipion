# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
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

import sys
from os.path import exists, basename
from glob import glob

from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import (PathParam, FloatParam, BooleanParam, 
                                        EnumParam, IntParam, StringParam, PointerParam)
from pyworkflow.utils.path import removeBaseExt, expandPattern, createLink, copyFile
from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Micrograph, Volume, Movie, PdbFile, CTFModel

from protocol import EMProtocol
from protocol_micrographs import ProtCTFMicrographs


class ProtImport(EMProtocol):
    """ Base class for other Import protocols. """
    pass


class ProtImportImages(ProtImport):
    """Common protocol to import a set of images in the project"""
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('pattern', PathParam, label=Message.LABEL_PATTERN,
                      help=Message.TEXT_PATTERN)
        form.addParam('checkStack', BooleanParam, label=Message.LABEL_CHECKSTACK, default=False)
        form.addParam('copyToProj', BooleanParam, label=Message.LABEL_COPYFILES, default=False)
        form.addParam('voltage', FloatParam, default=200,
                   label=Message.LABEL_VOLTAGE)
        form.addParam('sphericalAberration', FloatParam, default=2.26,
                   label=Message.LABEL_SPH_ABERRATION)
        form.addParam('ampContrast', FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        form.addParam('magnification', IntParam, default=50000,
                   label=Message.LABEL_MAGNI_RATE)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importImagesStep(self, pattern, checkStack, voltage, sphericalAberration, amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)
        filePaths = glob(pattern)
        
        imgSet = self._createSet()
        acquisition = imgSet.getAcquisition()
        # Setting Acquisition properties
        acquisition.setVoltage(voltage)
        acquisition.setSphericalAberration(sphericalAberration)
        acquisition.setAmplitudeContrast(amplitudeContrast)
        acquisition.setMagnification(magnification)
        
        # Call a function that should be implemented by each subclass
        self._setOtherPars(imgSet)
        
        outFiles = [imgSet.getFileName()]
        imgh = ImageHandler()
        img = imgSet.ITEM_TYPE()
        img.setAcquisition(acquisition)
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
            
            if checkStack:
                _, _, _, n = imgh.getDimensions(dst)
            
            if self._idRegex:
                # Try to match the micrograph id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fn)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
                micId = int(match.group(1))
            else:
                micId = None
                
            if n > 1:
                for index in range(1, n+1):
                    img.cleanObjId()
                    img.setMicId(micId)
                    img.setFileName(dst)
                    img.setIndex(index)
                    imgSet.append(img)
            else:
                img.setObjId(micId)
                img.setFileName(dst)
                imgSet.append(img)
            outFiles.append(dst)
            
            sys.stdout.write("\rImported %d/%d" % (i+1, size))
            sys.stdout.flush()
            
        print "\n"
        
        args = {}
        outputSet = self._getOutputSet(self._className)
        args[outputSet] = imgSet
        self._defineOutputs(**args)
        
        return outFiles
    
    #--------------------------- INFO functions ----------------------------------------------------
    def getPattern(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
        """
        import re
        self._idRegex = None
        pattern = expandPattern(self.pattern.get())
        match = re.match('[^#]*(#+)[^#]*', pattern)
        if match is not None:
            g = match.group(1)
            n = len(g)
            self._idRegex = re.compile(pattern.replace(g, '(%s)' % ('\d'*n)))
            pattern = pattern.replace(g, '[0-9]'*n)
        return pattern   
        
    def _validate(self):
        errors = []
        if not self.getPattern():
            errors.append(Message.ERROR_PATTERN_EMPTY)
        else:
            print "Validating pattern '%s'" % self.getPattern()
            filePaths = glob(self.getPattern())
        
            if len(filePaths) == 0:
                errors.append(Message.ERROR_PATTERN_FILES + " %s" % self.getPattern())

        return errors
    
    def _summary(self):
        summary = []
        outputSet = self._getOutputSet(self._className)
        if not hasattr(self, outputSet):
            summary.append("Output " + self._className + "s not ready yet.") 
            if self.copyToProj:
                summary.append("*Warning*: Import step could be take a long time due to the images are copying in the project.")
        else:
            summary.append("Import of %d " % getattr(self, outputSet).getSize() + self._className + "s from %s" % self.pattern.get())
            summary.append("Sampling rate : %0.2f A/px" % getattr(self, outputSet).getSamplingRate())
        
        return summary
    
    def _methods(self):
        methods = []
        outputSet = self._getOutputSet(self._className)
        if hasattr(self, outputSet):
            methods.append("%d " % getattr(self, outputSet).getSize() + self._className + "s has been imported")
            methods.append("with a sampling rate of %0.2f A/px" % getattr(self, outputSet).getSamplingRate())
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getFiles(self):
        return getattr(self, self._getOutputSet(self._className)).getFiles()
    
    def _getOutputSet(self, setName):
        return "output" + setName + "s"


class ProtImportMicrographs(ProtImportImages):
    """Protocol to import a set of micrographs to the project"""

    _className = 'Micrograph'
    _label = 'import micrographs'
    
    def _defineParams(self, form):
        ProtImportImages._defineParams(self, form)
        form.addParam('samplingRateMode', EnumParam, default=SAMPLING_FROM_IMAGE,
                   label=Message.LABEL_SAMP_MODE,
                   choices=[Message.LABEL_SAMP_MODE_1, Message.LABEL_SAMP_MODE_2])
        form.addParam('samplingRate', FloatParam, default=1, 
                   label=Message.LABEL_SAMP_RATE,
                   condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE)
        form.addParam('scannedPixelSize', FloatParam, default=7.0,
                   label=Message.LABEL_SCANNED,
                   condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER)
    
    def _insertAllSteps(self):
        self._createSet = self._createSetOfMicrographs
        self._insertFunctionStep('importImagesStep', self.getPattern(), self.checkStack.get(), 
                                 self.voltage.get(), self.sphericalAberration.get(), 
                                 self.ampContrast.get(), self.magnification.get()) #, self.samplingRate.get(),
    
    def _validate(self):
        errors = ProtImportImages._validate(self)
#         if self._checkMrcStack():
#             errors.append("The micrographs can't be a mrc stack")
        return errors
    
    def _setOtherPars(self, micSet):
        micSet.getAcquisition().setMagnification(self.magnification.get())
        
        if self.samplingRateMode == SAMPLING_FROM_IMAGE:
            micSet.setSamplingRate(self.samplingRate.get())
        else:
            micSet.setScannedPixelSize(self.scannedPixelSize.get())
    
    def _checkMrcStack(self):
        filePaths = glob(self.getPattern())
        imgh = ImageHandler()
        stack = False
        for f in filePaths:
            ext = removeBaseExt(f)
            _, _, _, n = imgh.getDimensions(f)
            if ext == ".mrc" and n > 1:
                stack = True
                break
        return stack


class ProtImportParticles(ProtImportImages):
    """Protocol to import a set of particles to the project"""
 
    _className = 'Particle'
    _label = 'import particles'
        
    def _defineParams(self, form):
        ProtImportImages._defineParams(self, form)
        form.addParam('samplingRate', FloatParam,
                   label=Message.LABEL_SAMP_RATE)
        
    def _insertAllSteps(self):
        self._createSet = self._createSetOfParticles
        self._insertFunctionStep('importImagesStep', self.getPattern(),
                                self.checkStack.get(), self.voltage.get(), self.sphericalAberration.get(),
                                self.ampContrast.get(), self.magnification.get())
        
    def _setOtherPars(self, imgSet):
        imgSet.setSamplingRate(self.samplingRate.get())
    
    def getFiles(self):
        return self.outputParticles.getFiles()


class ProtImportVolumes(ProtImport):
    """Protocol to import a set of volumes to the project"""
    _label = 'import volumes'
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)         
       
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
    
    def _validate(self):
        errors = []
        if self.getPattern():
            errors.append(Message.ERROR_PATTERN_EMPTY)
        
        if self._getNumberFilePaths(self.pattern.get()) == 0:
                errors.append(Message.ERROR_PATTERN_FILES)

        return errors


class ProtImportPdb(ProtImport):
    """Protocol to import a set of pdb volumes to the project"""
    _label = 'import pdb volumes'
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)         
       
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('path', StringParam, 
                      label="Pattern",
                      help='Specify a path or an url to desired PDB structure.')
         
    def _insertAllSteps(self):
        self._insertFunctionStep('importPdbStep', self.path.get())
        
    def importPdbStep(self, path):
        """ Copy volumes matching the filename pattern
        Register other parameters.
        """
        if not exists(path):
            raise Exception("PDB not found at *%s*" % path)
        pdb = PdbFile()
        pdb.setFileName(path)
        self._defineOutputs(outputPdb=pdb)

    def _summary(self):
        summary = ['PDB file imported from *%s*' % self.path.get()]

        return summary
    
    def _validate(self):
        errors = []
        if not exists(self.path.get()):
            errors.append("PDB not found at *%s*" % self.path.get())
        return errors


class ProtImportMovies(ProtImportImages):
    """Protocol to import a set of movies (from direct detector cameras) to the project"""
    
    _className = 'Movie'
    _label = 'import movies'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImportImages._defineParams(self, form)
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
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importMoviesStep', self.getPattern(), self.voltage.get(), self.magnification.get(),
                                 self.sphericalAberration.get(), self.ampContrast.get())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importMoviesStep(self, pattern, voltage, magnification, sphericalAberration, amplitudeContrast):
        """ Copy movies matching the filename pattern
        Register other parameters.
        """
        filePaths = self._getFilePaths(pattern)
        movSet = self._createSetOfMovies()
        acquisition = movSet.getAcquisition()
        # Setting Acquisition properties
        acquisition.setVoltage(voltage)
        acquisition.setSphericalAberration(sphericalAberration)
        acquisition.setAmplitudeContrast(amplitudeContrast)
        acquisition.setMagnification(magnification)
        
        if self.samplingRateMode == SAMPLING_FROM_IMAGE:
            movSet.setSamplingRate(self.samplingRate.get())
        else:
            movSet.setScannedPixelSize(self.scannedPixelSize.get())
        
        imgh = ImageHandler()
        
        for f in filePaths:
            dst = self._getPath(basename(f))
            if self.copyToProj:
                copyFile(f, dst)
            else:
                createLink(f, dst)
            _, _, _, n = imgh.getDimensions(dst)
            
            mov = Movie()
            mov.setAcquisition(acquisition)
            mov.setSamplingRate(movSet.getSamplingRate())
            movSet.append(mov)
            
            for i in range(1, n + 1):
                mic = Micrograph()
                mic.setAcquisition(acquisition)
                mic.setLocation(i, dst)
                mic.setSamplingRate(mov.getSamplingRate())
                mov.append(mic)
        
        self._defineOutputs(outputMovies=movSet)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getFilePaths(self, pattern):
        """ Return a sorted list with the paths of files"""
        filePaths = glob(expandPattern(pattern))
        filePaths.sort()
        
        return filePaths
    
    
class ProtImportCTF(ProtImport, ProtCTFMicrographs):
    """Protocol to import a set of movies (from direct detector cameras) to the project"""
    _label = 'import ctf'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                          label='Input particles',
                          help='Select the particles that you want to update the CTF parameters.')
        form.addParam('pattern', PathParam, 
                      label='CTF pattern',
                      help='Select files containing the CTF estimation.\n'
                           'You should use #### characters in the pattern\n'
                           'to mark where the micrograph id will be taken from. ')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importCtfStep', 
                                 self.inputMicrographs.get().getObjId(),
                                 self.getPattern())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importCtfStep(self, micsId, pattern):
        """ Copy movies matching the filename pattern
        Register other parameters.
        """
        inputMics = self.inputMicrographs.get()
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(inputMics)
        
        ctfFiles = self._getFilePaths(pattern)
        ctfDict = {}
        for fn in ctfFiles:
            # Try to match the micrograph id from filename
            # this is set by the user by using #### format in the pattern
            match = self._idRegex.match(fn)
            if match is None:
                raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
            ctfId = int(match.group(1))
            ctfDict[ctfId] = fn
            
        from pyworkflow.em.packages.brandeis.convert import parseCtffindOutput

        for mic in inputMics:
            result = parseCtffindOutput(ctfDict[mic.getObjId()])
            defocusU, defocusV, defocusAngle = result
            # save the values of defocus for each micrograph in a list
            ctf = CTFModel()
            ctf.copyObjId(mic)
            ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
            ctf.setMicrograph(mic)
            ctfSet.append(ctf)
        
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(inputMics, ctfSet)            
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getPattern(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
        """
        import re
        self._idRegex = None
        pattern = expandPattern(self.pattern.get())
        match = re.match('[^#]*(#+)[^#]*', pattern)
        if match is not None:
            g = match.group(1)
            n = len(g)
            self._idRegex = re.compile(pattern.replace(g, '(%s)' % ('\d'*n)))
            pattern = pattern.replace(g, '[0-9]'*n)
        else:
            raise Exception("Importing ctf will required match micrographs ids.")
        return pattern   
    
    def _getFilePaths(self, pattern):
        """ Return a sorted list with the paths of files"""
        filePaths = glob(pattern)
        filePaths.sort()
        
        return filePaths
    
    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
