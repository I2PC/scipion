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
from itertools import izip




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
        
        __threads, __mpi = self._getDefaultParallel()
        
        form.addParallelSection(threads=__threads, mpi=__mpi)
        
    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass
    
    def _getDefaultParallel(self):
        """ Return the default value for thread and MPI
        for the parallel section definition.
        """
        return (2, 1)


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


class ProtExtractParticles(ProtParticles):
    pass


class ProtParticlePicking(ProtParticles):

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                      label=Message.LABEL_INPUT_MIC, important=True,                      
                      help='Select the SetOfMicrographs to be used during picking.')

    #--------------------------- INFO functions ----------------------------------------------------
    def getSummary(self, coordSet):
        summary = []
        summary.append("Number of particles picked: %d" % coordSet.getSize())
        summary.append("Particle size: %d" % coordSet.getBoxSize())
        return "\n".join(summary)

    def getMethods(self, output):
        msg = 'User picked %d particles with a particle size of %s.' % (output.getSize(), output.getBoxSize())
        return msg
    
    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s of size %d." % (self.getObjectTag(self.getInputMicrographs()), self.getInputMicrographs().getSize()))
        
        if self.getOutputsSize() > 1:
            for key, output in self.iterOutputAttributes(EMObject):
                msg = self.getMethods(output)
                methodsMsgs.append("*%s:*\n%s"%(self.getObjectTag(output), msg))
        elif(self.getOutputsSize() == 1):
            output = self.getCoords()
            msg = self.getMethods(output)
            methodsMsgs.append("*%s:*\n%s"%(self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def getInputMicrographs(self):
        return self.inputMicrographs.get()
    
    def getCoords(self):
        count = self.getOutputsSize()
        if count == 0:
            return None
        suffix = str(count) if count > 1 else ''
        outputName = 'outputCoordinates' + suffix
        return getattr(self, outputName)

    def _createOutput(self, outputDir):
        micSet = self.getInputMicrographs()
        count = self.getOutputsSize()
        suffix = str(count + 1) if count > 0 else ''
        outputName = 'outputCoordinates' + suffix
        coordSet = self._createSetOfCoordinates(micSet, suffix)
        self.readSetOfCoordinates(outputDir, coordSet)
        coordSet.setObjComment(self.getSummary(coordSet))
        outputs = {outputName: coordSet}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    def readSetOfCoordinates(self, workingDir, coordSet):
        pass

    def _summary(self):
        summary = []
        if self.getInputMicrographs() is  not None:
            summary.append("Number of input micrographs: %d" % self.getInputMicrographs().getSize())

        if(self.getOutputsSize() > 1):
            for key, output in self.iterOutputAttributes(EMObject):
                label = output.getObjLabel() if output.getObjLabel() != "" else key
                summary.append("*%s:*\n%s"%(key, output.getObjComment()))
        elif(self.getOutputsSize() == 1):
            summary.append(self.getCoords().getObjComment())
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def getCoordsDir(self):
        pass

    def registerCoords(self, coordsDir):
        count = self.getOutputsSize()
        suffix = str(count + 1) if count > 0 else ''
        outputName = 'outputCoordinates' + suffix
        from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates
        inputset = self.getInputMicrographs()
        outputset = self._createSetOfCoordinates(inputset, suffix=suffix)#micrographs are the input set if protocol is not finished
        readSetOfCoordinates(coordsDir, outputset.getMicrographs(), outputset)
        summary = self.getSummary(outputset)
        outputset.setObjComment(summary)
        outputs = {outputName: outputset}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(inputset, outputset)
        self._store()





