# **************************************************************************
# *
# * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
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




from pyworkflow.em.protocol import ProtProcessParticles
import pyworkflow.protocol.params as params
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles

#import pyworkflow.em as em
#import pyworkflow.em.metadata as md



#from pyworkflow.em.convert import ImageHandler
#from pyworkflow.utils.properties import Message
#from convert import (xmippToLocation, writeSetOfParticles)

        
class XmippProtApplyTransformationMatrix(ProtProcessParticles):
    """ 
    Apply transformation matrix  of an aligned volume on 
    a set of particles to modify their angular assignment.
    Note:
    These particles are practically related to the 
    aligned volume (but before alignment).
    """
    
    _label = 'apply transformation matrix'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help="Aligned particles that their  "
                           "angular assignment needs to be modified.")
        
      
        
        
        
        form.addParallelSection(threads=1, mpi=4)    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):        
        fnInputParticles = self._getExtraPath('input_particles.xmd')
        self._insertFunctionStep('convertInputStep', fnInputParticles)
        #self._insertFunctionStep('applyAlignmentStep', fnInputParticles)
        #self._insertFunctionStep('createOutputStep')
    #--------------------------- STEPS functions --------------------------------------------        
    
    def convertInputStep(self, fnInputParticles):
        """ Create a metadata with the images and geometrical information. """
        writeSetOfParticles(self.inputParticles.get(), fnInputParticles)
        for part in self.inputParticles.get().iterItems():        
            Mpart = part.getTransform().getMatrix()
        print "Mpart = ", Mpart
       
    
    
    #--------------------------- INFO functions --------------------------------------------
    #def _validate(self):############  check shavad file alignment e volume hast ya na...
    #    errors = []
    #    return errors
        
    #def _summary(self):
    #    summary = []
    #    if not hasattr(self, 'outputParticles'):
    #        summary.append("Output particles not ready yet.")
    #    else:
    #        summary.append("Applied alignment to %s particles." % self.inputParticles.get().getSize())
    #    return summary

    #def _methods(self):
    #    if not hasattr(self, 'outputParticles'):
    #        return ["Output particles not ready yet."]
    #    else:
    #        return ["We applied alignment to %s particles from %s and produced %s."
    #                % (self.inputParticles.get().getSize(), self.getObjectTag('inputParticles'), self.getObjectTag('outputParticles'))]
    
    #--------------------------- UTILS functions --------------------------------------------
    