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
This sub-package contains protocol for particles filters operations
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from data import *
import xmipp3


        
class XmippDefFilterParticles(DefProcessParticles):
        
    def _addProcessParam(self):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass
      
      
class XmippProtFilter(xmipp3.XmippProtocol, ProtFilterParticles):
    """ Protocol base for Xmipp filters. """
    
    _definition = XmippDefFilterParticles()
    
    def __init__(self):
        ProtFilterParticles.__init__(self)
        self._program = "xmipp_transform_filter"
        self._args = "-i %(inputFn)s --save_metadata_stack %(outputMd)s --keep_input_columns --track_origin "
    
    def _defineFilenames(self):
        self.inputFn = self._getPath('input_images.xmd')    
        self.outputMd = self._getPath('output_images.xmd')
        self.outputStk = self._getPath('output_images.stk')
        
    def _defineSteps(self):
        """ Mainly prepare the command line for call ml(f)2d program"""
        self._defineFilenames()
        self.inputPcts = self.inputParticles.get()       
        inputPctsFn = self._insertConvertStep('inputPcts', XmippSetOfImages, self.inputFn)
        self._insertFilterStep(inputPctsFn, self.outputStk)        
        self._insertFunctionStep('createOutput')
        
    def createOutput(self):
        outputPcts = XmippSetOfParticles(self.outputMd)
        self._processOutput(outputPcts)
        self._defineOutputs(outputParticles=outputPcts)
        
    def _insertFilterStep(self, inputFn, outputFn):
        """ This is the function that will insert the filter step. """
        pass
        
    def _processOutput(self, outputPcts):
        """ This function should be implemented
        if some additional modifications needs to be done
        on output particles.
        """
        pass
        
    def _summary(self):
        summary = []
        return summary
    

class XmippDefFourierFilter(XmippDefFilterParticles):
    """ Definition for XmippProtoFourierFilter """

    def _addProcessParam(self):
        self.addParam('filterType', EnumParam, choices=['low pass', 'high pass', 'band pass'],
                      label="Filter type", default=xmipp3.BAND_PASS,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Select what type of Fourier filter do you want to apply.\n'
                           '<low pass>: all frequency components below <High frequency> are preserved.\n'
                           '<high pass>: all frequency components above <Low frequency> are preserved.\n'
                           '<band pass>: all frequency components between <Low frequency> and <High frequency> are preserved.\n')
        self.addParam('lowFreq', DigFreqParam, default=0.02, condition='filterType != 0',
                      label='Low Frequency (0 < f < 0.5)',
                      help='Low frequency cuttoff to apply the filter.\n')          
        self.addParam('highFreq', DigFreqParam, default=0.35, 
                      label='High Frequency (0 < f < 0.5)', condition='filterType != 1',
                      help='High frequency cuttoff to apply the filter.\n'
                           'Set to 0.5 for a <high pass> filter.')          
        self.addParam('freqDecay', FloatParam, default=0.02, expertLevel=LEVEL_ADVANCED,
                      label='Frequency decay',
                      help='It describes the length of the amplitude decay in a raised cosine') 
     
class XmippProtFourierFilter(XmippProtFilter):
    """ This class implement a protocol for Fourier filter.
    You may do a lowpass filter by setting FreqLow to 0. 
    You may do a high pass filter by setting FreqHigh to 0.5.
    """
    _definition = XmippDefFourierFilter()
    
    def _insertFilterStep(self, inputFn, outputFn):
        lowFreq = self.lowFreq.get()
        highFreq = self.highFreq.get()
        filterType = self.filterType.get()
        freqDecay = self.freqDecay.get()
        outputMd = self.outputMd
        
        if filterType == xmipp3.LOW_PASS:
            filterStr = "low_pass %(highFreq)f"
        elif filterType == xmipp3.HIGH_PASS:
            filterStr = "high_pass %(lowFreq)f"
        else:
            filterStr = "band_pass %(lowFreq)f %(highFreq)f"
        args = self._args + " --fourier " + filterStr + " %(freqDecay)f "
        if outputFn != inputFn:
            args += " -o %(outputFn)s"
                        
        self._insertRunJobStep(self._program, args % locals())
    

class XmippDefGaussianFilter(XmippDefFilterParticles):
    """ Definition for XmippProtoFourierFilter """

    def _addProcessParam(self):
        self.addParam('freqSigma', FloatParam, default=0.04, expertLevel=LEVEL_ADVANCED,
                      label='Frequency sigma',
                      help='Remind that the Fourier frequency is normalized between 0 and 0.5') 
        
class XmippProtGaussianFilter(XmippProtFilter):
    
    _definition = XmippDefGaussianFilter()
    
    def _insertFilterStep(self, inputFn, outputFn):
        freqSigma = self.freqSigma.get()
        outputMd = self.outputMd
        
        args = self._args + " --fourier gaussian %(freqSigma)f "
        if outputFn != inputFn:
            args += " -o %(outputFn)s"
                        
        self._insertRunJobStep(self._program, args % locals())
    