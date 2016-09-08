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
import numpy as np

        
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
        form.addParam('inputTransformMat', params.PathParam, 
                      label="Tranfsormation matrix",  
                      help="This is a txt file named 'transformation-matrix'.\n"
                           "Note:\n "
                           "Use param --copyGeo in xmipp_volume_align to create "
                           "this file.")       
        
        form.addParallelSection(threads=1, mpi=4)    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):        
        fnOutputParts = self._getExtraPath('output_particles.xmd')
        inputSet = self.inputParticles.get()
        self._insertFunctionStep('applyTransformMatStep', fnOutputParts, inputSet)        
    #--------------------------- STEPS functions --------------------------------------------        
    
    def applyTransformMatStep(self, fnOutputParts, inputSet):        
        data = []
        fhTransformMat = open(self.inputTransformMat.get(), "r")        
        for line in fhTransformMat:
            fields = line.split()
            rowdata = map(float, fields)
            data.extend(rowdata)
            
        count = 0
        TransformMat = [[0 for i in range(4)] for i in range(4)]
        for i in range(4):
            for j in range(4):
                TransformMat[i][j] = data[count]
                count += 1
        TransformMat = np.array (TransformMat)
        TransformMatrix = np.matrix(TransformMat)
        print "Transformation matrix (volume_align output) =\n", TransformMatrix
        #inverseTransform = np.linalg.inv(TransformMatrix)
        #print "Inverse transformation matrix (volume_align output) =\n", inverseTransform
               
       
        outputSet = self._createSetOfParticles()
        for part in inputSet.iterItems():        
            partTransformMat = part.getTransform().getMatrix()            
            partTransformMatrix = np.matrix(partTransformMat)            
            applyMatrix = TransformMatrix * partTransformMatrix
            #applyMatrix = inverseTransform * partTransformMatrix
                      
            part.getTransform().setMatrix(applyMatrix)
            outputSet.append(part)       
        outputSet.copyInfo(inputSet)        
        self._defineOutputs(outputParticles=outputSet)        
        writeSetOfParticles(outputSet, fnOutputParts)
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if not self.inputTransformMat.get():
            errors.append("We need transformation matrix "
                          "to process input particles!!!")
        return errors
        
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Applied alignment to %s particles." % 
                           self.inputParticles.get().getSize())
        return summary

    def _methods(self):
        if not hasattr(self, 'outputParticles'):
            return ["Output particles not ready yet."]
        else:
            return ["We applied alignment to %s particles from %s and produced %s."
                    %(self.inputParticles.get().getSize(), 
                      self.getObjectTag('inputParticles'), 
                      self.getObjectTag('outputParticles'))]
    