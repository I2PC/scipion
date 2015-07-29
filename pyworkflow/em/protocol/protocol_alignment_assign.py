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


from pyworkflow.em import ALIGN_2D
from pyworkflow.protocol.params import PointerParam
from protocol_2d import ProtAlign2D


class ProtAlignmentAssign(ProtAlign2D):
    """ Assign a the alignment calculated for a set of particles to another set.
    This protocol will take into account the differences of pixel size (A/pix)
    between the two sets and multiply by the right factor the shifts.
    The particles with the alignment can also be a subset of the other images
    """
    _label = 'alignment assign'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Select the particles that you want to update the new alignment.')
        form.addParam('inputAlignment', PointerParam, pointerClass='SetOfParticles',
                      label="Input alignments",
                      help='Select the particles with alignment to be apply to the other particles.')

        form.addParallelSection(threads=0, mpi=0)

#--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('createOutputStep')

    def _updateItem(self, item, row):
        """ Implement this function to do some
        update actions over each single item
        that will be stored in the output Set.
        """
        # Add alignment info from corresponding item on inputAlignment
        inputAlignment = self.inputAlignment.get()
        scale = inputAlignment.getSamplingRate()/self.inputParticles.get().getSamplingRate()

        alignedParticle = inputAlignment[item.getObjId()]
        # If alignment is found for this particle set the alignment info on the output particle, if not do not write that item
        if alignedParticle is not None:
            alignment = alignedParticle.getTransform()
            alignment.scaleShifts2D(scale)
            item.setTransform(alignment)
        else:
            item._appendItem = False

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()
        inputAlignment = self.inputAlignment.get()
        outputParticles = self._createSetOfParticles()
        outputParticles.copyInfo(inputParticles)

        outputParticles.copyItems(inputParticles,
                            updateItemCallback=self._updateItem)

        self._defineOutputs(outputParticles=outputParticles)
        self._defineSourceRelation(self.inputParticles, outputParticles)
        self._defineSourceRelation(self.inputAlignment, outputParticles)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            scale = self.inputAlignment.get().getSamplingRate()/self.inputParticles.get().getSamplingRate()
            summary.append("Assigned alignment to %s particles from a total of %s." % (self.outputParticles.getSize(), self.inputParticles.get().getSize()))
            if scale != 1:
                summary.append("Applied scale of %s." % scale)
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputParticles'):
            methods.append("Output particles not ready yet.")
        else:
            scale = self.inputAlignment.get().getSamplingRate()/self.inputParticles.get().getSamplingRate()
            methods.append("We assigned alignment to %s particles from %s and produced %s."
                           % (self.outputParticles.getSize(), self.getObjectTag('inputParticles'), self.getObjectTag('outputParticles')))
            if scale != 1:
                methods.append("Applied scale factor of %s." % scale)
        return methods

    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        #check that input set of aligned particles do have 2D alignment
        errors = [ ]
        inputAlignmentSet = self.inputAlignment.get()
        if not inputAlignmentSet.hasAlignment():
            errors.append("Input alignment set should contains some kind of alignment (2D, 3D or Projection).")
        else:
            # Just for consistency, check that the particles really contains Transform object
            first = inputAlignmentSet.getFirstItem()
            alignment = first.getTransform()
            if alignment is None:
                errors.append('Inconsistency detected in *Input alignment* !!!')
                errors.append('It has alignment: _%s_, but the alignment is missing!!!' % inputAlignmentSet.getAlignment())
            
        # Add some errors if input is not valid
        return errors
