# **************************************************************************
# *
# * Authors:  J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Nov 2014
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


from os.path import basename

from pyworkflow.protocol.params import PointerParam, FileParam
from pyworkflow.em.protocol import BatchProtocol
from pyworkflow.em.data import SetOfParticles
import xmipp 



class BatchProtNMACluster(BatchProtocol):
    """ Protocol executed when a cluster is created
    from NMA images and theirs deformations.
    """
    _label = 'nma cluster'
    
    def _defineParams(self, form):
        form.addHidden('inputNmaDimred', PointerParam, pointerClass='EMObject')
        form.addHidden('sqliteFile', FileParam)
        
    #--------------------------- INSERT steps functions --------------------------------------------
        
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('reconstructStep')
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------   
        
    def convertInputStep(self):
        # It is unusual to create the output in the convertInputStep,
        # but just to avoid reading twice the sqlite with particles
        inputSet = self.inputNmaDimred.get().getInputParticles()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputSet)
        
        tmpSet = SetOfParticles(filename=self.sqliteFile.get())        
        partSet.appendFromImages(tmpSet)
        # Register outputs
        self._defineOutputs(outputParticles=partSet)
        self._defineTransformRelation(inputSet, partSet)

    def reconstructStep(self):
        pass
    
    def createOutputStep(self):
        pass
        
        
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        return []
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions --------------------------------------------

    def getOutputMatrixFile(self):
        return self._getExtraPath('output_matrix.txt')
