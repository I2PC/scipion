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
from pyworkflow.protocol.params import PathParam, IntParam, PointerParam
from pyworkflow.em.data import Coordinate



from base import ProtImportFiles



class ProtImportCoordinates(ProtImportFiles):
    """ Protocol to import a set of coordinates """
    _label = 'import coordinates'

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_XMIPP3 = 1
    IMPORT_FROM_RELION = 2
    IMPORT_FROM_EMAN = 3



    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtImportFiles._defineParams(self, form)

        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs', 
                          label='Input micrographs',
                          help='Select the particles that you want to import coordinates.')

        form.addParam('boxSize', IntParam, 
                      label='Box size',
                      help='')


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
        self.filesPath = self.filesPath.get()
        if self.importFrom == self.IMPORT_FROM_XMIPP3:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport

            return XmippImport(self, self.filesPath)
        elif self.importFrom == self.IMPORT_FROM_RELION:
            from pyworkflow.em.packages.relion.dataimport import RelionImport
            return RelionImport(self, self.filesPath)
        #elif self.importFrom == self.IMPORT_FROM_EMAN:
        #    self.importFilePath = self.sqliteFile.get('').strip()
        #    return EmanImport(self, self.importFilePath)
        else:
            self.importFilePath = ''
            return None
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importCoordinatesStep',
                                 self.inputMicrographs.get().getObjId(),
                                 self.getPattern())



    def importCoordinatesStep(self, importFrom, *args):
        ci = self.getImportClass()
        ci.importCoordinates()

    #--------------------------- STEPS functions ---------------------------------------------------

        
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
