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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
In this module are protocol base classes related to EM imports of Micrographs, Particles, Volumes...
"""

from os.path import join, exists

from pyworkflow.protocol.params import IntParam, PointerParam, FloatParam, BooleanParam
from pyworkflow.utils.path import removeBaseExt
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking
import xmipp
        
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
        return ['auto', 'xmipp','relion', 'eman', 'dogpicker']

    def _getDefaultChoice(self):
        return  self.IMPORT_FROM_AUTO
    
    def _getFilesCondition(self):
        """ Return an string representing the condition
        when to display the files path and pattern to grab
        files.
        """
        return True
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):


        ProtImportFiles._defineParams(self, form)

    def _defineImportParams(self, form):
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                          label='Input micrographs',
                          help='Select the micrographs for which you want to import coordinates.')

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

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        self._insertFunctionStep('createOutputStep', importFrom,
                                     self.filesPath.get())

    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self, importFrom, *args):
        inputMics = self.inputMicrographs.get()
        coordsSet = self._createSetOfCoordinates(inputMics)
        coordsSet.setBoxSize(self.boxSize.get())
        ci = self.getImportClass()
        for coordFile, fileId in self.iterFiles():
            mic = self.getMatchingMic(coordFile,fileId)
            if mic is not None:
                def addCoordinate(coord):
                    coord.setMicrograph(mic)
                    self.correctCoordinatePosition(coord)
                    coordsSet.append(coord)
                # Parse the coordinates in the given format for this micrograph
                ci.importCoordinates(coordFile, addCoordinate=addCoordinate)
        
        self._defineOutputs(outputCoordinates=coordsSet)
        self._defineSourceRelation(self.inputMicrographs, coordsSet)

    #--------------------------- INFO functions ---------------------------------------------------
    def _summary(self):
        summary = []

        if not hasattr(self, 'outputCoordinates'):
            msg = 'Output coordinates not ready yet'
        else:
            msg = "%s  coordinates from micrographs %s were imported using %s format."%(self.outputCoordinates.getSize(), self.getObjectTag('inputMicrographs'), self._getImportChoices()[self.getImportFrom()])
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
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        filesPath = self.filesPath.get()
        importFrom = self.getImportFrom()
        
        if importFrom == self.IMPORT_FROM_XMIPP:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport
            return XmippImport(self, filesPath)

        elif importFrom == self.IMPORT_FROM_RELION:
            from pyworkflow.em.packages.relion.dataimport import RelionImport
            return RelionImport(self, filesPath)
        
        elif importFrom == self.IMPORT_FROM_EMAN:
            from pyworkflow.em.packages.eman2.dataimport import EmanImport
            return EmanImport(self)
        
        elif importFrom == self.IMPORT_FROM_DOGPICKER:
            from pyworkflow.em.packages.appion.dataimport import DogpickerImport
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
                if coordBase in micBase or micBase in coordBase: #temporal use of in
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
    
        boxSize = 100
        importFrom = self.getImportFrom()
        scale = self.scale.get()

        if importFrom == ProtImportCoordinates.IMPORT_FROM_XMIPP:
            configfile = join(self.filesPath.get(), 'config.xmd')
            existsConfig = exists(configfile)
            if existsConfig:
                md = xmipp.MetaData('properties@' + configfile)
                configobj = md.firstObject()
                boxSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, configobj)
        if importFrom == ProtImportCoordinates.IMPORT_FROM_EMAN:
            # Read the boxSize from the e2boxercache/base.json
            jsonFnbase = join(self.filesPath.get(), 'e2boxercache', 'base.json')
            from pyworkflow.em.packages.eman2 import loadJson
            jsonBoxDict = loadJson(jsonFnbase)
            boxSize = int(jsonBoxDict["box_size"])
        boxSize = (int)(boxSize * scale)
        return boxSize
