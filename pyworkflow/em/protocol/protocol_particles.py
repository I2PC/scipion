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
In this module are protocol base classes related to EM Particles
"""

from pyworkflow.protocol.params import PointerParam
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import EMObject
from pyworkflow.utils.properties import Message



class ProtParticles(EMProtocol):
    pass


class ProtProcessParticles(ProtParticles):
    """ This class will serve as a base for all protocol
    that performs some operation on Particles (i.e. filters, mask, resize, etc)
    It is mainly defined by an inputParticles and outputParticles.
    """
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputParticles', PointerParam, important=True,
                      label=Message.LABEL_INPUT_PART, pointerClass='SetOfParticles')
        # Hook that should be implemented in subclasses
        self._defineProcessParams(form)
        form.addParallelSection(threads=2, mpi=1)
        
    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass  


class ProtFilterParticles(ProtProcessParticles):
    """Base class for filters on particles of type ProtPreprocessParticles"""
    pass

class ProtOperateParticles(ProtProcessParticles):
    """Base class for operations on particles of type ProtPreprocessParticles"""
    pass

class ProtMaskParticles(ProtProcessParticles):
    """ This is the base for the branch of mask, 
    between the ProtPreprocessParticles """
    pass


class ProtParticlePicking(ProtParticles):
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCoordinates'):
            summary.append(Message.TEXT_NO_OUTPUT_CO) 
        else:
            #TODO: MOVE following line to manual picking
            summary.append("Number of input micrographs: %d" % self.inputMicrographs.get().getSize())
            summary.append("Number of particles picked: ")
            for _, output in self.iterOutputAttributes(EMObject):
                summary.append('    %d on one set' % output.getSize())
        return summary
    
    def _methods(self):
        methodsMsgs = []
        if not hasattr(self, 'outputCoordinates'):
            return methodsMsgs

        methodsMsgs.append("User picked ")
        for _, output in self.iterOutputAttributes(EMObject):
            methodsMsgs.append('%d particles from %d micrographs with a particle size of %d.' % (output.getSize(), self.inputMicrographs.get().getSize(), output.getBoxSize()))

        return methodsMsgs

    def getCoordsSuffix(self):
        count = len(list(self.iterOutputAttributes(EMObject)))
        suffix = str(count) if count > 1 else ''
        
        return suffix
    
    def getCoords(self):
        suffix = self.getCoordsSuffix()
        outputName = 'outputCoordinates' + suffix
        return getattr(self, outputName)


class ProtExtractParticles(ProtParticles):
    pass
