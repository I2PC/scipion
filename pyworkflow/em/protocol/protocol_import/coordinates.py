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
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking

from base import ProtImport



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
