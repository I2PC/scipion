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

from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
import xmipp3
from convert import createXmippInputImages, readSetOfParticles

from pyworkflow.em.constants import *
from constants import *


class XmippProcess(EMProtocol):
    """ Class to create a base template for Xmipp protocol that share
    a common structure: build a commmand line and call a program. """
    def __init__(self):
        self._args = "-i %(inputFn)s"
        
    def _insertAllSteps(self):
        self._defineFilenames()
        self._insertFunctionStep("convertStep")
        self._insertProcessStep(self.inputFn, self.outputStk, self.outputMd)
        self._insertFunctionStep('createOutput')
        
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        args = self._getCommand(inputFn)
        
        if outputFn != inputFn:
            args += " -o " + outputFn
        
        if outputMd is not None:
            args += (" --save_metadata_stack %(outputMd)s --keep_input_columns --track_origin") % locals()
        
        self._insertRunJobStep(self._program, args)


class XmippProcessParticles(ProtProcessParticles, XmippProcess):
    """ Class to create a base template for Xmipp protocols that process SetOfParticles """
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
        XmippProcess.__init__(self)
    
    def _defineFilenames(self):
        self.inputFn = createXmippInputImages(self, self.inputParticles.get())
        self.outputMd = self._getPath('output_images.xmd')
        self.outputStk = self._getPath('output_images.stk')
    
    def convertStep(self):
        """ convert if necessary"""
        pass
    
    def createOutput(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        readSetOfParticles(self.outputMd, imgSet, imgSet.hasCTF())
        self._processOutput(imgSet)
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(self.inputParticles.get(), imgSet)
        
    def _processOutput(self, outputPcts):
        """ This function should be implemented
        if some additional modifications needs to be done
        on output particles.
        """
        pass


class XmippProcessVolumes(ProtPreprocessVolumes, XmippProcess):
    """ Class to create a base template for Xmipp protocols that process both volume or a SetOfVolumes objects """
    def __init__(self, **args):
        ProtPreprocessVolumes.__init__(self, **args)
        XmippProcess.__init__(self)
    
    def _defineFilenames(self):
        """ Prepare the files to process """
        self.outputStk = self._getPath("output_volumes.stk")
        self.inputFn = self.outputStk
        self.outputMd = None
    
    def convertStep(self):
        """ convert if necessary"""
        volSet = self.inputVolumes.get()
        
        # Check volSet is a volume or a stack
        if isinstance(volSet, Volume):
            self.iniModel  = volSet.getFileName()
            self.singleVolume = True
            ImageHandler().convert(self.iniModel, (1, self.outputStk))
        else:
            volSet.writeStack(self.outputStk)
            self.singleVolume = False
    
    def createOutput(self):
        volSet = self.inputVolumes.get()
        if self.singleVolume:
            vol = Volume()
            vol.copyInfo(volSet)
            vol.setFileName(self.outputStk)
            self._defineOutputs(outputVol=vol)
        else:
            volumes = self._createSetOfVolumes()
            volumes.copyInfo(volSet)
            inVolSet = self.inputVolumes.get()
            for i, vol in enumerate(inVolSet):
                j = i + 1 
                vol.setLocation(j, self.outputStk)
                volumes.append(vol)
            self._defineOutputs(outputVol=volumes)

        self._defineTransformRelation(volSet, self.outputVol)
    
    def _processOutput(self, outputPcts):
        """ This function should be implemented
        if some additional modifications needs to be done
        on output particles.
        """
        pass
    
