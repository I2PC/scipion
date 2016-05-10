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

import pyworkflow.em as em
import pyworkflow.em.metadata as md

import pyworkflow.protocol.params as params

from pyworkflow.em.convert import ImageHandler
from pyworkflow.utils.properties import Message
from convert import (xmippToLocation, writeSetOfParticles)

        
class XmippProtApplyAlignment(em.ProtAlign2D):
    """ Apply alignment parameters and produce a new set of images. """
    _label = 'apply alignment 2d'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam, important=True,  
                      pointerCondition='hasAlignment2D',
                      label=Message.LABEL_INPUT_PART, pointerClass='SetOfParticles',
                      help='Select the particles that you want to apply the'
                           'alignment parameters.')
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for call cl2d align program"""
        
        # Create a metadata with the geometrical information 
        # as expected by Xmipp
        imgsFn = self._getPath('input_particles.xmd')
        self._insertFunctionStep('convertInputStep', imgsFn)
        self._insertFunctionStep('applyAlignmentStep', imgsFn)
        self._insertFunctionStep('createOutputStep')        

    #--------------------------- STEPS functions --------------------------------------------        
    
    def convertInputStep(self, outputFn):
        """ Create a metadata with the images and geometrical information. """
        writeSetOfParticles(self.inputParticles.get(), outputFn)
        
        return [outputFn]
    
    def applyAlignmentStep(self, inputFn):
        """ Create a metadata with the images and geometrical information. """
        outputStk = self._getPath('aligned_particles.stk')
        args = '-i %(inputFn)s -o %(outputStk)s --apply_transform ' % locals()
        self.runJob('xmipp_transform_geometry', args)
        
        return [outputStk]
     
    def createOutputStep(self):
        particles = self.inputParticles.get()

        # Generate the SetOfAlignmet
        alignedSet = self._createSetOfParticles()
        alignedSet.copyInfo(particles)

        inputMd = self._getPath('aligned_particles.xmd')
        alignedSet.copyItems(particles,
                             updateItemCallback=self._updateItem,
                             itemDataIterator=md.iterRows(inputMd, sortByLabel=md.MDL_ITEM_ID))
        # Remove alignment 2D
        alignedSet.setAlignment(em.ALIGN_NONE)

        # Define the output average

        avgFile = self._getExtraPath("average.xmp")

        imgh = ImageHandler()
        avgImage = imgh.computeAverage(alignedSet)

        avgImage.write(avgFile)

        avg = em.Particle()
        avg.setLocation(1, avgFile)
        avg.copyInfo(alignedSet)

        self._defineOutputs(outputAverage=avg)
        self._defineSourceRelation(self.inputParticles, avg)

        self._defineOutputs(outputParticles=alignedSet)
        self._defineSourceRelation(self.inputParticles, alignedSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors
        
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Applied alignment to %s particles." % self.inputParticles.get().getSize())
        return summary

    def _methods(self):
        if not hasattr(self, 'outputParticles'):
            return ["Output particles not ready yet."]
        else:
            return ["We applied alignment to %s particles from %s and produced %s."
                    % (self.inputParticles.get().getSize(), self.getObjectTag('inputParticles'), self.getObjectTag('outputParticles'))]
    
    #--------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        """ Implement this function to do some
        update actions over each single item
        that will be stored in the output Set.
        """
        # By default update the item location (index, filename) with the new binary data location
        newFn = row.getValue(md.MDL_IMAGE)
        newLoc = xmippToLocation(newFn)
        item.setLocation(newLoc)
        # Also remove alignment info
        item.setTransform(None)

