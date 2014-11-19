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
from os.path import exists, basename, join
from glob import glob
import re

from pyworkflow.object import Float, Integer
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import (PathParam, FloatParam, BooleanParam, FileParam,
                                        EnumParam, IntParam, StringParam, PointerParam,
                                        LEVEL_EXPERT)
from pyworkflow.utils.path import expandPattern, createLink, copyFile
from pyworkflow.em.constants import SAMPLING_FROM_IMAGE, SAMPLING_FROM_SCANNER
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Volume, PdbFile, Coordinate, Acquisition

from protocol import EMProtocol
from protocol_particles import ProtParticlePicking





class ProtImport(EMProtocol):
    """ Base class for other Import protocols. 
    All imports protocols will have:
    1) Several options to import from (_getImportOptions function)
    2) First option will always be "from files". (for this option 
      files with a given pattern will be retrieved  and the ### will 
      be used to mark an ID part from the filename.
      - For each file a function to process it will be called (_importFile(fileName, fileId))
    """
    IMPORT_FROM_FILES = 0
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineBasicParams(form)
        
    def _defineBasicParams(self, form):
        importChoices = ['files'] + self._getImportChoices()
        filesCondition = self._getFilesCondition()
        
        form.addSection(label='Import')
        form.addParam('importFrom', EnumParam, 
                      choices=importChoices, default=self.IMPORT_FROM_FILES,
                      label='Import from',
                      help='Select the type of import.')
        form.addParam('filesPath', PathParam, 
                      condition=filesCondition,
                      label="Files path",
                      help="Select the files path from where do you want to import\n"
                           "the files. The path can also contains wildcards to select\n"
                           "files from different folders.\n"
                           "Examples:\n"
                           "data/day??_micrographs/ \n"
                           "~/Particles/")
        form.addParam('filesPattern', StringParam,
                      label='Pattern', 
                      condition=filesCondition,
                      help="Select the pattern of the files to be imported.\n"
                           "The pattern can contains standard wildcards such as:\n"
                           "*, ?, etc ... or special ones as ### to mark some\n"
                           "digits in the filename to be used as ID. ")
        
        form.addParam('copyFiles', BooleanParam, default=False, 
                      expertLevel=LEVEL_EXPERT, 
                      label="Copy files?",
                      help="By default the files are not copied into\n"
                           "the project to avoid data duplication and to save\n"
                           "disk space. Instead of copying, symbolic links are\n"
                           "created pointing to original files. This approach\n"
                           "has the drawback that if the project is moved to\n"
                           "another computer, the links needs to be restored.\n")
        
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        if self.importFrom == self.IMPORT_FROM_FILES:
            if not self.getPattern():
                errors.append("The path and pattern can not be both empty!!!")
            else:
                # Just check the number of files matching the pattern
                self.getMatchFiles()
                if self.numberOfFiles == 0:
                    errors.append("There are no files matching the pattern " + "%s" % self.getPattern())

        return errors
    
    #--------------------------- BASE methods to be overriden ------------------
    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        return []
    
    def _getFilesCondition(self):
        """ Return an string representing the condition
        when to display the files path and pattern to grab
        files.
        """
        return '(importFrom == %d)' % self.IMPORT_FROM_FILES
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getPattern(self):
        """ Expand the pattern using environ vars or username
        and also replacing special character # by digit matching.
        """
        self._idRegex = None
        fullPattern = join(self.filesPath.get(''), self.filesPattern.get(''))
        pattern = expandPattern(fullPattern)
        match = re.match('[^#]*(#+)[^#]*', pattern)
        
        if match is not None:
            g = match.group(1)
            n = len(g)
            self._idRegex = re.compile(pattern.replace(g, '(%s)' % ('\d'*n)))
            pattern = pattern.replace(g, '[0-9]'*n)
        
        return pattern   
    
    def getMatchFiles(self, pattern=None):
        """ Return a sorted list with the paths of files that matched the pattern"""
        if pattern is None:
            pattern = self.getPattern()
        filePaths = glob(pattern)
        filePaths.sort()
        self.numberOfFiles = len(filePaths)
        
        return filePaths
    
    def iterFiles(self):
        """ Iterate throught the files matched with the pattern.
        Provide the fileName and fileId.
        """
        filePaths = self.getMatchFiles()
        # Set a function to copyFile or createLink
        # depending in the user selected option 
        if self.copyFiles:
            self.copyOrLink = copyFile
        else:
            self.copyOrLink = createLink
        
        for fileName in filePaths:
            if self._idRegex:
                # Try to match the file id from filename
                # this is set by the user by using #### format in the pattern
                match = self._idRegex.match(fileName)
                if match is None:
                    raise Exception("File '%s' doesn't match the pattern '%s'" % (fileName, self.getPattern()))
                fileId = int(match.group(1))
            else:
                fileId = None
                
            yield fileName, fileId            


class ProtImportImages(ProtImport):
    """Common protocol to import a set of images into the project"""
    # This label should be set in subclasses
    _label = 'None'
    # The following class property should be set in each import subclass
    # for example, if set to SetOfParticles, this will the output classes
    # It is also assumed that a function with the name _createSetOfParticles
    # exists in the EMProtocol base class
    _outputClassName = 'None'  
    # If set to True, each binary file will be inspected to
    # see if it is a binary stack containing more items
    _checkStacks = True
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImport._defineParams(self, form)
        
        group = form.addGroup('Acquisition info')
        group.addParam('voltage', FloatParam, default=200,
                   label=Message.LABEL_VOLTAGE, 
                   help=Message.TEXT_VOLTAGE)
        group.addParam('sphericalAberration', FloatParam, default=2,
                   label=Message.LABEL_SPH_ABERRATION, 
                   help=Message.TEXT_SPH_ABERRATION)
        group.addParam('amplitudeContrast', FloatParam, default=0.1,
                      label=Message.LABEL_AMPLITUDE,
                      help=Message.TEXT_AMPLITUDE)
        group.addParam('magnification', IntParam, default=50000,
                   label=Message.LABEL_MAGNI_RATE, 
                   help=Message.TEXT_MAGNI_RATE)
        self.acquisitionGroup = group
        
    #--------------------------- INSERT functions ---------------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importImagesStep', self.getPattern(), 
                                 self.voltage.get(), self.sphericalAberration.get(), 
                                 self.amplitudeContrast.get(), self.magnification.get()) #, self.samplingRate.get(),
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def importImagesStep(self, pattern, voltage, sphericalAberration, amplitudeContrast, magnification):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)
        
        createSetFunc = getattr(self, '_create' + self._outputClassName)
        imgSet = createSetFunc()
        acquisition = imgSet.getAcquisition()
        self._fillAcquistion(acquisition)
        # Setting Acquisition properties

        
        # Call a function that should be implemented by each subclass
        self._setSampling(imgSet)
        
        outFiles = [imgSet.getFileName()]
        imgh = ImageHandler()
        img = imgSet.ITEM_TYPE()
        img.setAcquisition(acquisition)
        n = 1
        
        for i, (fileName, fileId) in enumerate(self.iterFiles()):
            dst = self._getExtraPath(basename(fileName))
            self.copyOrLink(fileName, dst)
            
            if self._checkStacks:
                _, _, _, n = imgh.getDimensions(dst)
                
            if n > 1:
                for index in range(1, n+1):
                    img.cleanObjId()
                    img.setMicId(fileId)
                    img.setFileName(dst)
                    img.setIndex(index)
                    imgSet.append(img)
            else:
                img.setObjId(fileId)
                img.setFileName(dst)
                imgSet.append(img)
            outFiles.append(dst)
            
            sys.stdout.write("\rImported %d/%d" % (i+1, self.numberOfFiles))
            sys.stdout.flush()
            
        print "\n"
        
        args = {}
        outputSet = self._getOutputName()
        args[outputSet] = imgSet
        self._defineOutputs(**args)
        
        return outFiles
    
    #--------------------------- INFO functions ----------------------------------------------------
    
    def _summary(self):
        summary = []
        outputSet = self._getOutputSet()
        if outputSet is None:
            summary.append("Output " + self._outputClassName + " not ready yet.") 
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import. ")
        else:
            summary.append("Import of *%d* %s from %s" % (outputSet.getSize(), self._getOutputItemName(), self.getPattern()))
            summary.append("Sampling rate : *%0.2f* A/px" % outputSet.getSamplingRate())
        
        return summary
    
    def _methods(self):
        methods = []
        outputSet = self._getOutputSet()
        if outputSet is not None:
            methods.append("We used *%d* %s" % (outputSet.getSize(), self._getOutputItemName())+\
                           " with a sampling rate of *%0.2f* A/px (microscope voltage %d kV, magnification %dx)" %
                            (outputSet.getSamplingRate(),round(self.voltage.get()),round(self.magnification.get())))
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getFiles(self):
        outputSet = self._getOutputSet()
        if outputSet is not None:
            return self._getOutputSet().getFiles()
        else:
            return []
    
    def _getOutputName(self):
        # We assume that the import output is always a 'SetOfSomething'
        return self._outputClassName.replace('SetOf', 'output')
    
    def _getOutputItemName(self):
        return self._outputClassName.replace('SetOf', '')
    
    def _getOutputSet(self):
        return getattr(self, self._getOutputName(), None)
    
    def _fillAcquistion(self, acquisition):
        """ Fill the acquition object with protocol params. """
        acquisition.setVoltage(self.voltage.get())
        acquisition.setSphericalAberration(self.sphericalAberration.get())
        acquisition.setAmplitudeContrast(self.amplitudeContrast.get())
        acquisition.setMagnification(self.magnification.get())
        
    def getAcquisition(self):
        """ Build and fill an acquisition object. """
        acquisition = Acquisition()
        self._fillAcquistion(acquisition)
        
        return acquisition    


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
        group.addParam('samplingRate', FloatParam,  allowsNull=True,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_IMAGE, 
                       label=Message.LABEL_SAMP_RATE,
                       help=Message.TEXT_SAMP_RATE)
        group.addParam('scannedPixelSize', FloatParam, default=7.0,
                       condition='samplingRateMode==%d' % SAMPLING_FROM_SCANNER,
                       label=Message.LABEL_SCANNED,
                       help='')
        
    def _setSampling(self, micSet):
        if self.samplingRateMode == SAMPLING_FROM_IMAGE:
            micSet.setSamplingRate(self.samplingRate.get())
        else:
            micSet.setScannedPixelSize(self.scannedPixelSize.get())
        
    
class ProtImportMicrographs(ProtImportMicBase):
    """Protocol to import a set of micrographs to the project"""
    _label = 'import micrographs'
    _outputClassName = 'SetOfMicrographs' 
    
    IMPORT_FROM_XMIPP3 = 1
    IMPORT_FROM_EMX = 2

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        return ['xmipp3', 'emx']
    
    def _defineBasicParams(self, form):
        ProtImportMicBase._defineBasicParams(self, form)
        form.addParam('micrographsMd', PathParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Micrographs metadata file',
                      help="Select the micrographs Xmipp metadata file.\n"
                           "It is usually a _micrograph.xmd_ file result\n"
                           "from import, preprocess or downsample protocols.")
        form.addParam('acquisitionFromXmd', BooleanParam, default=True,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Read acquisition from Metadatas?', 
                      help='If you set to _Yes_, some acquistion parameters\n'
                           'will try to be recovered from the _microscope.xmd_\n'
                           'and _acquisition_info.xmd_\n')
        
        form.addParam('micrographsEMX', PathParam,
              condition = '(importFrom == %d)' % self.IMPORT_FROM_EMX,
              label='Input EMX file',
              help="Select the EMX file containing micrographs information.\n"
                   "See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]")
        form.addParam('acquisitionFromEmx', BooleanParam, default=True,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_EMX,
                      label='Read acquisition from EMX?', 
                      help='If you set to _Yes_, the acquistion parameters\n'
                           'will be recovered from the EMX file.')

    def _insertAllSteps(self):
        if self.importFrom == self.IMPORT_FROM_FILES:
            ProtImportMicBase._insertAllSteps(self)
        elif self.importFrom == self.IMPORT_FROM_EMX:
            self._insertFunctionStep('importFromEmxStep', self.micrographsEMX.get())
            
    def importFromEmxStep(self, emxFile):
        from pyworkflow.em.packages.emxlib import EmxImport
        emx = EmxImport()
        emx.importData(self, emxFile)
        
    def importFromXmippStep(self, micrographsMd):
        from pyworkflow.em.packages.xmipp3.convert import XmippImport
        xi = XmippImport(self)
        xi.importMicrographs(micrographsMd)
        
        

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
        
    def _setSampling(self, movieSet):
        ProtImportMicBase._setSampling(self, movieSet)
        movieSet.setGain(self.gainFile.get())
                    

class ProtImportParticles(ProtImportImages):
    """Protocol to import a set of particles to the project"""
    _label = 'import particles'
    _outputClassName = 'SetOfParticles'
        
    def _defineParams(self, form):
        ProtImportImages._defineParams(self, form)
        group = self.acquisitionGroup
        group.addParam('samplingRate', FloatParam,
                   label=Message.LABEL_SAMP_RATE)
        
    def _setSampling(self, imgSet):
        imgSet.setSamplingRate(self.samplingRate.get())
    

class ProtImportAverages(ProtImportParticles):
    """Protocol to import a set of averages to the project"""
    _label = 'import averages'
    _outputClassName = 'SetOfAverages'    
    

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
        EMProtocol.__init__(self, **args)         
       
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
    
    
class ProtImportCoordinates(ProtImport, ProtParticlePicking):
    """ Protocol to import a set of coordinates """
    _label = 'import coordinates'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs', 
                          label='Input micrographs',
                          help='Select the particles that you want to import coordinates.')
        form.addParam('pattern', PathParam, 
                      label='Coordinate files pattern',
                      help='Select files containing the coordinate files in Arachnid format.\n'
                           'You should use #### characters in the pattern\n'
                           'to mark where the micrograph id will be taken from. ')
        form.addParam('boxSize', IntParam, 
                      label='Box size',
                      help='')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importCoordinateStep', 
                                 self.inputMicrographs.get().getObjId(),
                                 self.getPattern())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importCoordinateStep(self, micsId, pattern):
        """ Copy movies matching the filename pattern
        Register other parameters.
        """
        inputMics = self.inputMicrographs.get()
        coordSet = self._createSetOfCoordinates(inputMics)
        coordSet.setBoxSize(self.boxSize.get())
        
        coordFiles = self._getFilePaths(pattern)
        coordDict = {}
        for fn in coordFiles:
            # Try to match the micrograph id from filename
            # this is set by the user by using #### format in the pattern
            match = self._idRegex.match(fn)
            if match is None:
                raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
            ctfId = int(match.group(1))
            coordDict[ctfId] = fn
            
        coordinate = Coordinate()
        coordinate._araId = Integer()
        coordinate._araPeak = Float()
        
        for mic in inputMics:
            micId = mic.getObjId()
            if micId in coordDict:
                araFile = coordDict[micId]
                for araId, araPeak, araX, araY in self._parseCoordinateFile(araFile):
                    coordinate.setX(araX)
                    coordinate.setY(araY)
                    coordinate.setMicId(micId)
                    coordinate._araId.set(araId)
                    coordinate._araPeak.set(araPeak)
                    coordinate.cleanObjId()
                    coordSet.append(coordinate)
            else:
                self.warning("Coordinates for micrograph with id %d were not found." % mic.getObjId())
            
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(inputMics, coordSet)            
    
    def _parseCoordinateFile(self, filename):
        """ Parse filename and return iterator over rows data """
        f = open(filename)
        for line in f:
            l = line.strip()
            if not l.startswith(';'):
                parts = l.split()
                yield parts[-4:]                
        f.close()
        
    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
    
