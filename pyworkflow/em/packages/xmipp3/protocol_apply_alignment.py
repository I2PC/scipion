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
This sub-package contains wrapper around align2d Xmipp program
"""

from pyworkflow.em import *  
from convert import (readSetOfParticles, locationToXmipp, 
                     writeSetOfParticles)

       
        
class XmippProtApplyAlignment(ProtAlign2D):
    """ Apply alignment parameters and produce a new set of images. """
    _label = 'apply alignment'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineAlignParams(self, form):
#         form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
#                       label='Input particles', 
#                       help='Select the particles that you want to apply the'
#                       'alignment parameters.')
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
        
        readSetOfParticles(self._getPath('aligned_particles.xmd'), alignedSet)
        
        self._defineOutputs(outputParticles=alignedSet)
        self._defineSourceRelation(particles, alignedSet)
                

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors
        
    def _summary(self):
        summary = []
        return summary
    
