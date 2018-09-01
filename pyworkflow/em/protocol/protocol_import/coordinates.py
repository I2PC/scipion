# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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
# *
# **************************************************************************

from os.path import join, exists
from glob import glob
from itertools import izip

from pyworkflow.em import metadata as md
from pyworkflow.protocol.params import (IntParam, PointerParam, FloatParam,
                                        BooleanParam, PathParam, EnumParam,
                                        LEVEL_ADVANCED)
from pyworkflow.utils.path import removeBaseExt, expandPattern
from pyworkflow.utils.properties import Message
from pyworkflow.em.data_tiltpairs import CoordinatesTiltPair
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.utils import importFromPlugin
import xmippLib
from base import ProtImportFiles


class ProtImportCoordinates(ProtImportFiles, ProtParticlePicking):
    """ Protocol to import a set of coordinates """
    _label = 'import coordinates'

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_XMIPP = 1
    IMPORT_FROM_RELION = 2
    IMPORT_FROM_EMAN = 3
    IMPORT_FROM_DOGPICKER = 4

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['auto', 'xmipp', 'relion', 'eman', 'dogpicker']

    def _getDefaultChoice(self):
        return self.IMPORT_FROM_AUTO

    def _getFilesCondition(self):
        """ Return an string representing the condition
        when to display the files path and pattern to grab
        files.
        """
        return True

    # ---------------------- DEFINE param functions ---------------------------
    def _defineParams(self, form):


        ProtImportFiles._defineParams(self, form)

    def _defineImportParams(self, form):
        form.addParam('inputMicrographs', PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs',
                      help='Select the micrographs for which you want to '
                           'import coordinates.')

        form.addParam('boxSize', IntParam, label='Box size')
        form.addParam('scale', FloatParam,
                      label='Scale', default=1,
                      help='Factor to scale coordinates')
        form.addParam('invertX', BooleanParam, default=False,
                      label='Invert X')
        form.addParam('invertY', BooleanParam, default=False,
                      label='Invert Y',
                      help='Invert Y for EMAN coordinates taken on dm3 or'
                           ' tif micrographs')

    # ------------------- INSERT steps functions ------------------------------
    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        self._insertFunctionStep('createOutputStep', importFrom,
                                 self.filesPath.get())

    # ------------------ STEPS functions --------------------------------------
    def createOutputStep(self, importFrom, *args):
        inputMics = self.inputMicrographs.get()
        coordsSet = self._createSetOfCoordinates(inputMics)
        coordsSet.setBoxSize(self.boxSize.get())
        ci = self.getImportClass()
        for coordFile, fileId in self.iterFiles():
            mic = self.getMatchingMic(coordFile, fileId)
            if mic is not None:
                def addCoordinate(coord):
                    coord.setMicrograph(mic)
                    self.correctCoordinatePosition(coord)
                    coordsSet.append(coord)
                # Parse the coordinates in the given format for this micrograph
                ci.importCoordinates(coordFile, addCoordinate=addCoordinate)

        self._defineOutputs(outputCoordinates=coordsSet)
        self._defineSourceRelation(self.inputMicrographs, coordsSet)

    # ---------------- INFO functions -----------------------------------------
    def _summary(self):
        summary = []

        if not hasattr(self, 'outputCoordinates'):
            msg = 'Output coordinates not ready yet'
        else:
            msg = "%s  coordinates from " % self.outputCoordinates.getSize()
            msg += "micrographs %s " % self.getObjectTag('inputMicrographs')
            importFrom = self._getImportChoices()[self.getImportFrom()]
            msg += "were imported using %s format." % importFrom
            if self.scale.get() != 1.:
                msg += " Scale factor %0.2f was applied." % self.scale
            if self.invertX.get():
                msg += " X coordinate was inverted."
            if self.invertY.get():
                msg += " Y coordinate was inverted."

            summary.append(msg)
            summary.append("Output coordinates: %s." % self.getObjectTag('outputCoordinates'))
        return summary

    def _methods(self):
        return self._summary()

    # ------------------ UTILS functions --------------------------------------
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        filesPath = self.filesPath.get()
        importFrom = self.getImportFrom()

        if importFrom == self.IMPORT_FROM_XMIPP:
            XmippImport = importFromPlugin('xmipp3.convert', 'XmippImport',
                                           'Xmipp is needed to import .xmd files',
                                           doRaise=True)
            return XmippImport(self, filesPath)

        elif importFrom == self.IMPORT_FROM_RELION:
            RelionImport = importFromPlugin('relion.convert', 'RelionImport',
                                            errorMsg='Relion is needed to import .star files',
                                            doRaise=True)
            return RelionImport(self, filesPath)

        elif importFrom == self.IMPORT_FROM_EMAN:
            EmanImport = importFromPlugin('eman2.convert', 'EmanImport',
                                          errorMsg='Eman is needed to import .json or '
                                                   '.box files',
                                          doRaise=True)
            return EmanImport(self, None)

        elif importFrom == self.IMPORT_FROM_DOGPICKER:
            DogpickerImport = importFromPlugin('appion.convert', 'DogpickerImport',
                                               errorMsg='appion plugin is needed to import '
                                                        'dogpicker files',
                                               doRaise=True)
            return DogpickerImport(self)
        else:
            self.importFilePath = ''
            return None

    def getImportFrom(self):
        importFrom = self.importFrom.get()
        if importFrom == self.IMPORT_FROM_AUTO:
            importFrom = self.getFormat()
        return importFrom

    def getFormat(self):
        for coordFile, _ in self.iterFiles():
            if coordFile.endswith('.pos'):
                return self.IMPORT_FROM_XMIPP
            if coordFile.endswith('.star'):
                return self.IMPORT_FROM_RELION
            if coordFile.endswith('.json') or coordFile.endswith('.box'):
                return self.IMPORT_FROM_EMAN
        return -1

    def getInputMicrographs(self):
        return self.inputMicrographs.get()

    def getMatchingMic(self, coordFile, fileId):
        """ Given a coordinates file check if there is a micrograph
        that this files matches.
        """
        micSet = self.inputMicrographs.get()

        if fileId is None:
            coordBase = removeBaseExt(coordFile)
            for mic in micSet:
                micBase = removeBaseExt(mic.getFileName())
                # temporal use of in
                if coordBase in micBase or micBase in coordBase:
                    return mic
            return None
        else:
            return micSet[fileId]

    def correctCoordinatePosition(self, coord):
        mic = coord.getMicrograph()
        scaleFactor = self.scale.get()
        x = coord.getX()
        y = coord.getY()
        if scaleFactor != 1.:
            x = coord.getX() * scaleFactor
            y = coord.getY() * scaleFactor
        if self.invertX:
            width = mic.getDim()[0]
            x = width - x
        if self.invertY:
            height = mic.getDim()[1]
            y = height - y
        coord.setX(x)
        coord.setY(y)

    def getDefaultBoxSize(self):
        """ This function is used by the wizard to estimate the box size. """
        boxSize = None
        ci = self.getImportClass()
        if hasattr(ci, 'getBoxSize'):
            boxSize = ci.getBoxSize([f for f, _ in self.iterFiles()][0])
            if boxSize is not None and self.scale != 1.:
                boxSize = int(boxSize * self.scale.get())

        return boxSize


class ProtImportCoordinatesPairs(ProtImportFiles):
    """ Protocol to import a set of tilt pair coordinates """
    _label = 'import coordinate pairs'

    IMPORT_FROM_XMIPP = 0
    IMPORT_FROM_EMAN = 1

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2..etc.)
        """
        return ['xmipp', 'eman']

    def _getDefaultChoice(self):
        return self.IMPORT_FROM_XMIPP

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        importChoices = self._getImportChoices()

        form.addSection(label='Import')
        form.addParam('importFrom', EnumParam,
                      choices=importChoices, default=self._getDefaultChoice(),
                      label='Import from',
                      help='Select the type of import.\n'
                           '_Auto_ - detects coordinate file type by extension.\n'
                           '_Xmipp_ - provide *.pos files\n'
                           '_Eman_ - provide info/*.json files')
        form.addParam('copyFiles', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Copy files?",
                      help="By default the files are not copied into the "
                           "project to avoid data duplication and to save "
                           "disk space. Instead of copying, symbolic links are "
                           "created pointing to original files. This approach "
                           "has the drawback that if the project is moved to "
                           "another computer, the links need to be restored.")
        form.addParam('xmippMdFn', PathParam, condition='importFrom==0',
                      important=True,
                      label='Micrograph .xmd file',
                      help='Provide input_micrographs.xmd file that contains '
                           'tilt angles information. This file is usually '
                           'created alongside the coordinates.')
        form.addParam('patternUntilted', PathParam, label=Message.LABEL_PATTERNU,
                      help=Message.TEXT_PATTERN)
        form.addParam('patternTilted', PathParam, label=Message.LABEL_PATTERNT,
                      help=Message.TEXT_PATTERN)
        form.addParam('inputMicrographsTiltedPair', PointerParam,
                      pointerClass='MicrographsTiltPair',
                      label='Input tilt pair micrographs',
                      help='Select the tilt pair micrographs for which you want '
                           'to import coordinates.')
        form.addParam('boxSize', IntParam, label='Box size')
        form.addParam('scale', FloatParam,
                      label='Scale', default=1,
                      help='Factor to scale coordinates')
        form.addParam('invertX', BooleanParam, default=False,
                      label='Invert X')
        form.addParam('invertY', BooleanParam, default=False,
                      label='Invert Y')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.micsFn = self._getExtraPath('input_micrographs.xmd')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        micTiltPairs = self.getInputMicrographs()
        # Get the converted input micrographs in Xmipp format
        writeSetOfMicrographsPairs = importFromPlugin('xmipp3.convert',
                                                      'writeSetOfMicrographsPairs')
        writeSetOfMicrographsPairs(micTiltPairs.getUntilted(),
                                   micTiltPairs.getTilted(),
                                   self.micsFn)
        uCoordSet, tCoordSet, anglesSet = self.readCoordinates()
        # Create CoordinatesTiltPair object
        coordsSet = self._createCoordinatesTiltPair(micTiltPairs,
                                                    uCoordSet, tCoordSet,
                                                    anglesSet, suffix='')
        self._defineOutputs(outputCoordinatesTiltPair=coordsSet)
        self._defineSourceRelation(self.inputMicrographsTiltedPair, coordsSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if not self.patternUntilted.get() or not self.patternTilted.get():
            errors.append(Message.ERROR_PATTERN_EMPTY)
        else:
            filePathsUntilted = glob(expandPattern(self.patternUntilted.get()))
            filePathsTilted = glob(expandPattern(self.patternTilted.get()))

            if len(filePathsUntilted) == 0 or len(filePathsTilted) == 0:
                errors.append(Message.ERROR_PATTERN_FILES)

        return errors

    def _summary(self):
        summary = []

        if not hasattr(self, 'outputCoordinatesTiltPair'):
            msg = 'Output tilt pair coordinates not ready yet'
        else:
            msg = "%s  coordinate pairs from micrographs %s were imported using %s format." % (
                self.outputCoordinatesTiltPair.getSize(),
                self.getObjectTag('inputMicrographsTiltedPair'),
                self._getImportChoices()[self.importFrom.get()])
            if self.scale.get() != 1.:
                msg += " Scale factor %0.2f was applied." % self.scale
            if self.invertX.get():
                msg += " X coordinate was inverted."
            if self.invertY.get():
                msg += " Y coordinate was inverted."

            summary.append(msg)
            summary.append("Output coordinates: %s." % self.getObjectTag('outputCoordinatesTiltPair'))
        return summary

    def _methods(self):
        return self._summary()

    # --------------------------- UTILS functions -----------------------------
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        importFrom = self.importFrom.get()

        if importFrom == self.IMPORT_FROM_XMIPP:
            XmippImport = importFromPlugin('xmipp3.convert', 'XmippImport',
                                           'Xmipp is needed to import .xmd files',
                                           doRaise=True)
            return XmippImport(self, None)
        else:  # import from EMAN
            EmanImport = importFromPlugin('eman2.convert', 'EmanImport',
                                          errorMsg='Eman is needed to import .json or '
                                                   '.box files',
                                          doRaise=True)
            return EmanImport(self, None)

    def getInputMicrographs(self):
        return self.inputMicrographsTiltedPair.get()

    def getMatchingCoord(self, micU, micT):
        """ Given a pair of mics check if there are coordinates
        that match.
        """
        mics = [micU, micT]
        patterns = [self.patternUntilted.get(), self.patternTilted.get()]
        coords = []

        for mic, pattern in izip(mics, patterns):
            micFnBase = removeBaseExt(mic.getFileName())
            for coordFn in self._iterFiles(pattern):
                if coordFn not in coords:
                    coorBase = removeBaseExt(coordFn).replace('_info', '')
                    if coorBase in micFnBase or micFnBase in coorBase:
                        coords.append(coordFn)

        if len(coords) == 2:
            return coords[0], coords[1]
        else:
            return False, False

    def correctCoordinatePosition(self, coord):
        mic = coord.getMicrograph()
        scaleFactor = self.scale.get()
        x = coord.getX()
        y = coord.getY()
        if scaleFactor != 1.:
            x = coord.getX() * scaleFactor
            y = coord.getY() * scaleFactor
        if self.invertX:
            width = mic.getDim()[0]
            x = width - x
        if self.invertY:
            height = mic.getDim()[1]
            y = height - y
        coord.setX(x)
        coord.setY(y)

    def getDefaultBoxSize(self):
        """ This function is used by the wizard to estimate the box size. """
        ci = self.getImportClass()
        pattern = self.patternUntilted.get()
        if hasattr(ci, 'getBoxSize'):
            boxSize = ci.getBoxSize([f for f in self._iterFiles(pattern)][0])
            if boxSize is not None and self.scale != 1.:
                boxSize = int(boxSize * self.scale.get())
            return boxSize

        boxSize = 100
        scale = self.scale.get()
        boxSize = int(boxSize * scale)

        return boxSize

    def readCoordinates(self):
        micTiltPairs = self.getInputMicrographs()
        ci = self.getImportClass()
        uSet = micTiltPairs.getUntilted()
        tSet = micTiltPairs.getTilted()
        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = self._createSetOfCoordinates(uSet, suffix='Untilted')
        tCoordSet = self._createSetOfCoordinates(tSet, suffix='Tilted')
        anglesSet = self._createSetOfAngles()

        def _importCoords(uCoordSet, tCoordSet):
            for micU, micT in izip(uSet, tSet):
                coordFnU, coordFnT = self.getMatchingCoord(micU, micT)
                if coordFnU and coordFnT:

                    def addCoordinateU(coord):
                        coord.setMicrograph(micU)
                        self.correctCoordinatePosition(coord)
                        uCoordSet.append(coord)

                    ci.importCoordinates(coordFnU, addCoordinate=addCoordinateU)

                    def addCoordinateT(coord):
                        coord.setMicrograph(micT)
                        self.correctCoordinatePosition(coord)
                        tCoordSet.append(coord)

                    ci.importCoordinates(coordFnT, addCoordinate=addCoordinateT)

                    def addAngles(ang):
                        anglesSet.append(ang)

                    if self.importFrom.get() == self.IMPORT_FROM_EMAN:
                        ci.importAngles(coordFnT, addAngles=addAngles)

        _importCoords(uCoordSet, tCoordSet)
        boxSize = self.boxSize.get()
        uCoordSet.setBoxSize(boxSize)
        tCoordSet.setBoxSize(boxSize)

        uCoordSet.write()
        tCoordSet.write()

        if self.importFrom.get() == self.IMPORT_FROM_XMIPP:
            readAnglesFromMicrographs = importFromPlugin('xmipp3.convert',
                                                         'readAnglesFromMicrographs')
            if exists(self.xmippMdFn.get()):
                checkAngles = self._compareMicPairs(self.xmippMdFn, uSet, tSet)
                if checkAngles:
                    readAnglesFromMicrographs(self.xmippMdFn, anglesSet)
                else:
                    self.warning('Angles for some micrographs were not found, '
                                 'skipping angles import')
            else:
                self.warning('Micrograph xmd file not provided, so tilt angles will '
                             'not be imported')
        anglesSet.write()

        return uCoordSet, tCoordSet, anglesSet

    def _iterFiles(self, pattern):
        filePaths = glob(expandPattern(pattern))
        for fn in filePaths:
            yield fn

    def _compareMicPairs(self, micFn, uSet, tSet):
        # compare micFn input file and input micsSet
        micMd = md.MetaData(micFn)
        micFnDict = {}
        inputMicsDict = {}
        XmippMdRow = importFromPlugin('xmipp3', 'XmippMdRow')
        for objId in micMd:
            row = XmippMdRow()
            row.readFromMd(micMd, objId)
            micUFn = removeBaseExt(row.getValue(md.MDL_MICROGRAPH))
            micTFn = removeBaseExt(row.getValue(md.MDL_MICROGRAPH_TILTED))
            micFnDict[micUFn] = micTFn

        for micU, micT in izip(uSet, tSet):
            inputMicsDict[removeBaseExt(micU.getFileName())] = removeBaseExt(micT.getFileName())

        for micKey in inputMicsDict:
            if micKey not in micFnDict:
                return False
        return True

    # all funcs below are required for tilt picker to work properly

    def __getOutputSuffix(self):
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(CoordinatesTiltPair):
            suffix = attrName.replace('outputCoordinatesTiltPair', '')
            try:
                counter = int(suffix)
            except:
                counter = 1  # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter + 1) if maxCounter > 0 else ''  # empty if not outputs

    def getSummary(self, coordSet):
        summary = []
        summary.append("Number of particle pairs picked: %s" % coordSet.getSize())
        summary.append("Particle size: %s" % coordSet.getBoxSize())
        return "\n".join(summary)

    def _getBoxSize(self):
        """ Redefine this function to set a specific box size to the output
        coordinates untilted and tilted.
        """
        return None

    def _readAngles(self, micsFn, suffix=''):
        # Read Angles from input micrographs
        anglesSet = self._createSetOfAngles(suffix=suffix)
        readAnglesFromMicrographs = importFromPlugin('xmipp3.convert',
                                                     'readAnglesFromMicrographs')
        readAnglesFromMicrographs(micsFn, anglesSet)
        anglesSet.write()
        return anglesSet

    def _readCoordinates(self, coordsDir, suffix=''):
        micTiltPairs = self.inputMicrographsTiltedPair.get()
        uSuffix = 'Untilted' + suffix
        tSuffix = 'Tilted' + suffix
        uSet = micTiltPairs.getUntilted()
        tSet = micTiltPairs.getTilted()
        # Create Untilted and Tilted SetOfCoordinates
        readSetOfCoordinates = importFromPlugin('xmipp3.convert',
                                                'readSetOfCoordinates')
        uCoordSet = self._createSetOfCoordinates(uSet, suffix=uSuffix)
        readSetOfCoordinates(coordsDir, uSet, uCoordSet)
        uCoordSet.write()
        tCoordSet = self._createSetOfCoordinates(tSet, suffix=tSuffix)
        readSetOfCoordinates(coordsDir, tSet, tCoordSet)
        tCoordSet.write()
        boxSize = self._getBoxSize()
        if boxSize:
            uCoordSet.setBoxSize(boxSize)
            tCoordSet.setBoxSize(boxSize)

        return uCoordSet, tCoordSet

    def registerCoords(self, coordsDir, store=True, readFromExtra=True):
        micTiltPairs = self.inputMicrographsTiltedPair.get()
        suffix = self.__getOutputSuffix()

        uCoordSet, tCoordSet = self._readCoordinates(coordsDir, suffix)

        if readFromExtra:
            micsFn = self._getExtraPath('input_micrographs.xmd')
        else:
            micsFn = self._getPath('input_micrographs.xmd')

        anglesSet = self._readAngles(micsFn, suffix)
        # Create CoordinatesTiltPair object
        outputset = self._createCoordinatesTiltPair(micTiltPairs,
                                                    uCoordSet, tCoordSet,
                                                    anglesSet, suffix)
        summary = self.getSummary(outputset)
        outputset.setObjComment(summary)
        outputName = 'outputCoordinatesTiltPair' + suffix
        outputs = {outputName: outputset}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputMicrographsTiltedPair, outputset)
        if store:
            self._store()
