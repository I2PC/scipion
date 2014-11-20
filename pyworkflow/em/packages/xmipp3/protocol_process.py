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
This sub-package contains classes to use in common processing operations of SetOfParticles, Volume or SetOfVolumes
"""

from pyworkflow.em.protocol import EMProtocol  
from pyworkflow.em.data import Volume
from convert import (writeSetOfParticles, readSetOfParticles, 
                     writeSetOfVolumes, readSetOfVolumes,
                     getImageLocation)



class XmippProcess(EMProtocol):
    """ Class to create a base template for Xmipp protocol that share
    a common structure: build a commmand line and call a program. """
    def __init__(self):
        self._args = "-i %(inputFn)s"
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self):
        """ convert if necessary"""
        pass
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._defineFilenames()
        self._insertFunctionStep("convertInputStep")
        self._insertProcessStep()
        self._insertFunctionStep('createOutputStep')
        
    def _preprocessOutput(self, output):
        """ This function should be implemented
        if some additional modifications needs to be done
        before population of the output elements.
        """
        pass  
      
    def _postprocessOutput(self, output):
        """ This function should be implemented
        if some additional modifications needs to be done
        after the population of output elements.
        """
        pass


class XmippProcessParticles(XmippProcess):
    """ Class to create a base template for Xmipp protocols 
    that process SetOfParticles
    """
    def __init__(self):
        XmippProcess.__init__(self)
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self):
        """ convert if necessary"""
        writeSetOfParticles(self.inputParticles.get(), self.inputFn)

    def createOutputStep(self):
        inputSet = self.inputParticles.get()
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(inputSet)
        
        self._preprocessOutput(imgSet)
        readSetOfParticles(self.outputMd, imgSet)
        self._postprocessOutput(imgSet)
        
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(inputSet, imgSet)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _defineFilenames(self):
        self.inputFn = self._getTmpPath('input_particles.xmd')
        self.outputMd = self._getExtraPath('output_images.xmd')
        self.outputStk = self._getExtraPath('output_images.stk')

class XmippProcessVolumes(XmippProcess):
    """ Class to create a base template for Xmipp protocols that process 
    both volume or a SetOfVolumes objects 
    """
    def __init__(self):
        XmippProcess.__init__(self)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self):
        """ convert if necessary"""
        if not self._isSingleInput():
            writeSetOfVolumes(self.inputVolumes.get(), self.inputFn)
            
    def createOutputStep(self):
        volInput = self.inputVolumes.get()
        
        if self._isSingleInput():
            vol = Volume()
            vol.copyInfo(volInput)
            vol.setLocation(1, self.outputStk)
            self._postprocessOutput(vol)
            self._defineOutputs(outputVol=vol)
        else:
            # ToDo: createSetOfVolumes not work properly when the protocol is resumed.
            volumes = self._createSetOfVolumes()
            volumes.copyInfo(volInput)
            
            self._preprocessOutput(volumes)
            readSetOfVolumes(self.outputMd, volumes)
            self._postprocessOutput(volumes)
            
            self._defineOutputs(outputVol=volumes)
        
        self._defineTransformRelation(volInput, self.outputVol)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _isSingleInput(self):
        return isinstance(self.inputVolumes.get(), Volume)
        
    def _defineFilenames(self):
        """ Prepare the files to process """
        if self._isSingleInput():
            self.inputFn = getImageLocation(self.inputVolumes.get())
            self.outputStk = self._getExtraPath("output_volume.vol")
        else:
            self.inputFn = self._getTmpPath('input_volumes.xmd')
            self.outputStk = self._getExtraPath("output_volumes.stk")
            self.outputMd = self._getExtraPath('output_volumes.xmd')
   
