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

from pyworkflow.object import Float, Integer
from pyworkflow.protocol.params import PathParam, IntParam, PointerParam, FloatParam, BooleanParam
from pyworkflow.em.data import Coordinate
from pyworkflow.utils.path import removeBaseExt
from pyworkflow.utils.path import join
from os.path import exists




from base import ProtImportFiles



class ProtImportCoordinates(ProtImportFiles):
    """ Protocol to import a set of coordinates """
    _label = 'import coordinates'

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_XMIPP = 1
    IMPORT_FROM_RELION = 2
    IMPORT_FROM_EMAN = 3



    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImportFiles._defineParams(self, form)

        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs', 
                          label='Input micrographs',
                          help='Select the particles that you want to import coordinates.')

        form.addParam('boxSize', IntParam, label='Box size')
        form.addParam('scale', FloatParam,
                      label='Scale', default=1,
                      help='factor to scale coordinates')
        form.addParam('invertX', BooleanParam,
                      label='Invert X')
        form.addParam('invertY', BooleanParam,
                      label='Invert Y')


    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['auto', 'xmipp','relion', 'eman']

    def _getDefaultChoice(self):
        return  self.IMPORT_FROM_AUTO


    def _insertAllSteps(self):
        importFrom = self.importFrom.get()

        ci = self.getImportClass()
        self._insertFunctionStep('importCoordinatesStep', importFrom,
                                     self.importFilePath)

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
        #elif self.importFrom == self.IMPORT_FROM_EMAN:
        #    self.importFilePath = self.sqliteFile.get('').strip()
        #    return EmanImport(self, self.importFilePath)
        else:
            self.importFilePath = ''
            return None
    #--------------------------- INSERT steps functions --------------------------------------------
    def getImportFrom(self):
        importFrom = self.importFrom.get()
        if importFrom == self.IMPORT_FROM_AUTO:
            importFrom = self.getFormat()
        return importFrom


    def importCoordinatesStep(self, importFrom, *args):
         inputMics = self.inputMicrographs.get()
         coordsSet = self._createSetOfCoordinates(inputMics)
         coordsSet.setBoxSize(self.boxSize.get())
         scaleFactor = self.scale.get()
         invertX = self.invertX.get()
         invertY = self.invertY.get()
         ci = self.getImportClass()
         for i, (fileName, fileId) in enumerate(self.iterFiles()):
            for mic in inputMics:
                if removeBaseExt(mic.getFileName()) in removeBaseExt(fileName):#temporal use of in
                    def addCoordinate(coord):
                        coord.setMicrograph(mic)
                        x = coord.getX()
                        y = coord.getY()
                        if scaleFactor != 1.:
                            x = coord.getX() * scaleFactor
                            y = coord.getY() * scaleFactor
                        if invertX:
                            width = mic.getDim()[0]
                            x = width - x
                        if invertY:
                            height = mic.getDim()[1]
                            y = height - y
                        coord.setX(x)
                        coord.setY(y)
                        coordsSet.append(coord)
                    ci.importCoordinates(fileName, addCoordinate=addCoordinate)
                    break





         self._defineOutputs(outputCoordinates=coordsSet)
         self._defineSourceRelation(inputMics, coordsSet)

    #--------------------------- STEPS functions ---------------------------------------------------

    def getFormat(self):
        for i, (fileName, fileId) in enumerate(self.iterFiles()):
            if fileName.endswith('.pos'):
                return self.IMPORT_FROM_XMIPP
            if fileName.endswith('.star'):
                return self.IMPORT_FROM_RELION
            if fileName.endswith('.json'):
                return self.IMPORT_FROM_EMAN
        return -1

    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []

    def _getFilesCondition(self):
        """ Return an string representing the condition
        when to display the files path and pattern to grab
        files.
        """
        return True

    def getDefaultBoxSize(self):
        import xmipp
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


    def _summary(self):
        summary = []

        if not hasattr(self, 'outputCoordinates'):
            msg = 'Output coordinates not ready yet'
        else:
            msg = "%s  coordinates from micrographs %s were imported using %s format."%(self.outputCoordinates.getSize(), self.getObjectTag(self.inputMicrographs.get()), self._getImportChoices()[self.getImportFrom()])
            if self.scale.get() != 1.:
                msg += " Scale factor %d was applied."%self.scale.get()
            if self.invertX.get():
                msg += " X coordinate was inverted."
            if self.invertY.get():
                msg += " Y coordinate was inverted."

            summary.append(msg)
            summary.append("Output coordinates: %s."%self.getObjectTag(self.outputCoordinates))
        return summary

    def _methods(self):
        return self._summary()


    # def importCoordinatesStep(self, micsId, pattern):
    #     """ Copy movies matching the filename pattern
    #     Register other parameters.
    #     """
    #     inputMics = self.inputMicrographs.get()
    #     coordSet = self._createSetOfCoordinates(inputMics)
    #     coordSet.setBoxSize(self.boxSize.get())
    #
    #     coordFiles = self._getFilePaths(pattern)
    #     coordDict = {}
    #     for fn in coordFiles:
    #         # Try to match the micrograph id from filename
    #         # this is set by the user by using #### format in the pattern
    #         match = self._idRegex.match(fn)
    #         if match is None:
    #             raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
    #         ctfId = int(match.group(1))
    #         coordDict[ctfId] = fn
    #
    #     coordinate = Coordinate()
    #     coordinate._araId = Integer()
    #     coordinate._araPeak = Float()
    #
    #     for mic in inputMics:
    #         micId = mic.getObjId()
    #         if micId in coordDict:
    #             araFile = coordDict[micId]
    #             for araId, araPeak, araX, araY in self._parseCoordinateFile(araFile):
    #                 coordinate.setX(araX)
    #                 coordinate.setY(araY)
    #                 coordinate.setMicId(micId)
    #                 coordinate._araId.set(araId)
    #                 coordinate._araPeak.set(araPeak)
    #                 coordinate.cleanObjId()
    #                 coordSet.append(coordinate)
    #         else:
    #             self.warning("Coordinates for micrograph with id %d were not found." % mic.getObjId())
    #
    #     self._defineOutputs(outputCoordinates=coordSet)
    #     self._defineSourceRelation(inputMics, coordSet)
    #
    # def _parseCoordinateFile(self, filename):
    #     """ Parse filename and return iterator over rows data """
    #     f = open(filename)
    #     for line in f:
    #         l = line.strip()
    #         if not l.startswith(';'):
    #             parts = l.split()
    #             yield parts[-4:]
    #     f.close()
