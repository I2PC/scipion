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
import xmipp3
from convert import createXmippInputImages, readSetOfParticles


class XmippProcess(ProtProcessParticles):
    """ Class to create a base template for Xmipp protocol that share
    a common structure: build a commmand line and call a program. """
    def __init__(self):
        ProtProcessParticles.__init__(self)
        self._args = "-i %(inputFn)s --save_metadata_stack %(outputMd)s --keep_input_columns --track_origin "
        
    def _defineFilenames(self):
        self.inputFn = createXmippInputImages(self, self.inputParticles.get())
        self.outputMd = self._getPath('output_images.xmd')
        self.outputStk = self._getPath('output_images.stk')
        
    def _insertAllSteps(self):
        self._defineFilenames()
        self._insertProcessStep(self.inputFn, self.outputStk, self.outputMd)        
        self._insertFunctionStep('createOutput')
        
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        args = self._getCommand(inputFn, outputFn, outputMd)
            
        if outputFn != inputFn:
            args += " -o " + outputFn
                        
        self._insertRunJobStep(self._program, args)
                
    def createOutput(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        readSetOfParticles(self.outputMd, imgSet, imgSet.hasCTF())
        self._processOutput(imgSet)
        self._defineOutputs(outputParticles=imgSet)
        
    def _processOutput(self, outputPcts):
        """ This function should be implemented
        if some additional modifications needs to be done
        on output particles.
        """
        pass
        
            
class XmippProtFilter(ProtFilterParticles, XmippProcess):
    """ Some filters operations such as: Fourier or Gaussian. """
    _label = 'filter particles'
    
    def __init__(self):
        ProtFilterParticles.__init__(self)
        XmippProcess.__init__(self)
        self._program = "xmipp_transform_filter"
        
    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam, choices=['fourier', 'gaussian'], 
                      default=xmipp3.FILTER_FOURIER, display=EnumParam.DISPLAY_LIST,
                      label="Filter type")
        form.addParam('fourierMode', EnumParam, choices=['low pass', 'high pass', 'band pass'],
                      default=xmipp3.FILTER_BAND_PASS,
                      condition='filterType == %d' % xmipp3.FILTER_FOURIER, 
                      label="Fourier mode", 
                      help='Depending on the filter mode some frequency (freq.) components \n'
                           'are keept and some are removed. \n '
                           '_low pass_: components below *High freq.* are preserved. \n '
                           '_high pass_: components above *Low freq.* are preserved. \n '
                           '_band pass_: components between *Low freq.* and *High freq.* are preserved. \n ')
        form.addParam('lowFreq', DigFreqParam, default=0.02, 
                      condition='filterType == %d and fourierMode != %d' % (xmipp3.FILTER_FOURIER, xmipp3.FILTER_LOW_PASS),
                      label='Low frequency (0 < f < 0.5)',
                      help='Low frequency cuttoff to apply the filter. \n ')          
        form.addParam('highFreq', DigFreqParam, default=0.35, 
                      condition='filterType == %d and fourierMode != %d' % (xmipp3.FILTER_FOURIER, xmipp3.FILTER_HIGH_PASS),
                      label='High frequency (0 < f < 0.5)', 
                      help='High frequency cuttoff to apply the filter.')          
        form.addParam('freqDecay', FloatParam, default=0.02, 
                      condition='filterType == %d' % xmipp3.FILTER_FOURIER, 
                      label='Frequency decay',
                      help='It describes the length of the amplitude decay in a raised cosine') 
        form.addParam('freqSigma', FloatParam, default=0.04, 
                      condition='filterType == %d' % xmipp3.FILTER_GAUSSIAN,
                      label='Frequency sigma',
                      help='Remind that the Fourier frequency is normalized between 0 and 0.5')   
    
    def _getCommand(self, inputFn, outputFn, outputMd):
        if self.filterType == xmipp3.FILTER_FOURIER:
            lowFreq = self.lowFreq.get()
            highFreq = self.highFreq.get()
            fourierMode = self.fourierMode.get()
            freqDecay = self.freqDecay.get()
            
            if fourierMode == xmipp3.FILTER_LOW_PASS:
                filterStr = "low_pass %(highFreq)f"
            elif fourierMode == xmipp3.FILTER_HIGH_PASS:
                filterStr = "high_pass %(lowFreq)f"
            else:
                filterStr = "band_pass %(lowFreq)f %(highFreq)f"
            args = (self._args + " --fourier " + filterStr + " %(freqDecay)f ") % locals()
        
        elif self.filterType == xmipp3.FILTER_GAUSSIAN:
            freqSigma = self.freqSigma.get()
            args = (self._args + " --fourier gaussian %(freqSigma)f ") % locals()
            
        else:
            raise Exception("Unknown filter type: %d" % self.filterType.get())
            
        return args


class XmippProtMask(XmippProcess):
    """ This class implement a protocol for applying a mask with Xmipp.
    """
    _label = 'mask particles'
    
    def __init__(self):
        XmippProcess.__init__(self)
        self._program = "xmipp_transform_mask"

    def _defineProcessParams(self, form):
        """ Add common mask parameters that can be used
        in several protocols definitions.
        Params:
            form: the Definition instance.
        """
        form.addParam('maskType', EnumParam, 
                      choices=['raised_cosine', 'circular', 'file'], 
                      default=xmipp3.MASK_RAISED_COSINE, 
                      label="Mask type", display=EnumParam.DISPLAY_COMBO,
                      help='Select which type of mask do you want to apply. \n ')
        
        form.addParam('maskRadius', IntParam, default=-1, 
                      condition='maskType != %d' % xmipp3.MASK_FILE,
                      label='Mask radius (px)',
                      help='This is the radius (in pixels) of the spherical mask ')       
    
        form.addParam('maskOuterRadius', IntParam, default=2, 
                      condition='maskType == %d' % xmipp3.MASK_RAISED_COSINE,
                      label='Mask outer radius (px)',
                      help='Outer radius in pixels for the raised cosine mask ')
 
        form.addParam('maskImage', PointerParam, pointerClass='Mask', 
                      condition='maskType == %d' % xmipp3.MASK_FILE,
                      label="Mask image", 
                      help='The mask file should have the same dimensions as your input particles. \n '
                           'The protein region should be 1 and the solvent should be 0.') 
                    
        form.addParam('fillType', EnumParam, 
                      choices=['value', 'min', 'max', 'avg'], 
                      default=xmipp3.MASK_FILL_VALUE,
                      label="Fill with ", display=EnumParam.DISPLAY_COMBO,
                      help='Select how are you going to fill the pixel values outside the mask. ')
        
        form.addParam('fillValue', IntParam, default=0, 
                      condition='fillType == %d' % xmipp3.MASK_FILL_VALUE,
                      label='Fill value',
                      help='Value to fill the pixel values outside the mask. ')   

    def _getCommand(self, inputFn, outputFn, outputMd):
        fillStr = self.getEnumText('fillType')
        maskType = self.getEnumText('maskType')
        maskRadius = self.maskRadius.get()
        maskBand = maskRadius + self.maskOuterRadius.get()
        
        # Weird logic of Xmipp, minus to define the mask
        maskRadius *= -1
        maskBand *= -1
        
        if self.fillType == xmipp3.MASK_FILL_VALUE:
            fillStr = str(self.fillValue.get())
        
        self._args += " --substitute %(fillStr)s --mask %(maskType)s "
        if self.maskType == xmipp3.MASK_RAISED_COSINE:
            self._args += " %(maskRadius)d %(maskBand)d" 
        elif self.maskType == xmipp3.MASK_CIRCULAR:
            self._args += " %(maskRadius)d"
        elif self.maskType == xmipp3.MASK_FILE:
            self._args += self.maskFn
        else:
            raise Exception("Unrecognized mask type: %d" % self.maskType.get())
            
        return self._args % locals()
    
    def _defineFilenames(self):
        XmippProcess._defineFilenames(self)
        self.maskFn = self._getPath('mask.spi')
        
    def _insertProcessStep(self, inputFn, outputFn, outputMd):
        if self.maskType == xmipp3.MASK_FILE:
            self._insertFunctionStep('convertMask', self.maskImage.get().getLocation())
        XmippProcess._insertProcessStep(self, inputFn, outputFn, outputMd)
    
    def convertMask(self, *args):
        """ Convert the input mask to file. """
        ImageHandler().convert(self.maskImage.get().getLocation(), (None, self.maskFn))

        
