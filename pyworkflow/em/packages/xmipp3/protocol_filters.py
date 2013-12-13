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


      
class XmippProtFilter(ProtFilterParticles):
    """ Protocol base for Xmipp filters. """
    def __init__(self):
        ProtFilterParticles.__init__(self)
        self._program = "xmipp_transform_filter"
        self._args = "-i %(inputFn)s --save_metadata_stack %(outputMd)s --keep_input_columns --track_origin "
    
    def _defineFilenames(self):
        self.inputFn = createXmippInputImages(self, self.inputParticles.get())
#        self.inputFn = self._getPath('input_images.xmd')    
        self.outputMd = self._getPath('output_images.xmd')
        self.outputStk = self._getPath('output_images.stk')
        
    def _defineSteps(self):
        self._defineFilenames()
        inputPctsFn = createXmippInputImages(self, self.inputParticles.get())
#        self.inputPcts = self.inputParticles.get()       
#        inputPctsFn = self._insertConvertStep('inputPcts', XmippSetOfImages, self.inputFn)
        self._insertFilterStep(inputPctsFn, self.outputStk, self.outputMd)        
        self._insertFunctionStep('createOutput')
        
    def createOutput(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        readSetOfParticles(self.outputMd, imgSet, imgSet.hasCTF())
        imgSet.write()
#        outputPcts = XmippSetOfParticles(self.outputMd)
        self._processOutput(imgSet)
        self._defineOutputs(outputParticles=imgSet)
        
    def _insertFilterStep(self, inputFn, outputFn, outputMd):
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
    
     
class XmippProtFourierFilter(XmippProtFilter):
    """ This class implement a protocol for Fourier filter.
    You may do a lowpass filter by setting FreqLow to 0. 
    You may do a high pass filter by setting FreqHigh to 0.5.
    """
    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam, choices=['low pass', 'high pass', 'band pass'],
                      label="Filter type", default=xmipp3.FILTER_BAND_PASS,
                      help='Select what type of Fourier filter do you want to apply. \n '
                           '*low pass*: all frequency components below <High frequency> are preserved. \n '
                           '*high pass*: all frequency components above <Low frequency> are preserved. \n '
                           '*band pass*: all frequency components between <Low frequency> and <High frequency> are preserved. \n ')
        form.addParam('lowFreq', DigFreqParam, default=0.02, 
                      condition='filterType != %d' % xmipp3.FILTER_LOW_PASS,
                      label='Low Frequency (0 < f < 0.5)',
                      help='Low frequency cuttoff to apply the filter. \n ')          
        form.addParam('highFreq', DigFreqParam, default=0.35, 
                      label='High Frequency (0 < f < 0.5)', 
                      condition='filterType != %d' % xmipp3.FILTER_HIGH_PASS,
                      help='High frequency cuttoff to apply the filter. \n '
                           'Set to 0.5 for a <high pass> filter.')          
        form.addParam('freqDecay', FloatParam, default=0.02, 
                      label='Frequency decay',
                      help='It describes the length of the amplitude decay in a raised cosine') 
    
    def _insertFilterStep(self, inputFn, outputFn, outputMd):
        lowFreq = self.lowFreq.get()
        highFreq = self.highFreq.get()
        filterType = self.filterType.get()
        freqDecay = self.freqDecay.get()
        
        if filterType == xmipp3.FILTER_LOW_PASS:
            filterStr = "low_pass %(highFreq)f"
        elif filterType == xmipp3.FILTER_HIGH_PASS:
            filterStr = "high_pass %(lowFreq)f"
        else:
            filterStr = "band_pass %(lowFreq)f %(highFreq)f"
        args = self._args + " --fourier " + filterStr + " %(freqDecay)f "
        if outputFn != inputFn:
            args += " -o %(outputFn)s"
                        
        self._insertRunJobStep(self._program, args % locals())
    

class XmippProtGaussianFilter(XmippProtFilter):
    def _defineProcessParams(self, form):
        form.addParam('freqSigma', FloatParam, default=0.04, 
                      label='Frequency sigma',
                      help='Remind that the Fourier frequency is normalized between 0 and 0.5') 
    
    def _insertFilterStep(self, inputFn, outputFn, outputMd):
        freqSigma = self.freqSigma.get()
        
        args = self._args + " --fourier gaussian %(freqSigma)f "
        if outputFn != inputFn:
            args += " -o %(outputFn)s"
                        
        self._insertRunJobStep(self._program, args % locals())
        

class XmippProtMask(XmippProtFilter):
    """ This class implement a protocol for applying a mask with Xmipp.
    """
    def __init__(self):
        XmippProtFilter.__init__(self)
        self._program = "xmipp_transform_mask"

    def _defineProcessParams(self, form):
        """ Add common mask parameters that can be used
        in several protocols definitions.
        Params:
            form: the Definition instance.
        """
        form.addParam('maskType', EnumParam, 
                      choices=['raised_cosine', 'circular', 'binary_file'], 
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
            
        form.addParam('maskFile', StringParam, default='', 
                      label='Binary mask file', 
                      condition='maskType == %d' % xmipp3.MASK_FILE,
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

    def _insertFilterStep(self, inputFn, outputFn, outputMd):
        args = self._args
        fillStr = self.getEnumText('fillType')
        maskType = self.getEnumText('maskType')
        maskRadius = self.maskRadius.get()
        maskBand = maskRadius + self.maskOuterRadius.get()
        
        if self.fillType == xmipp3.MASK_FILL_VALUE:
            fillStr = str(self.fillValue.get())
        
        args += " --substitute %(fillStr)s --mask %(maskType)s "
        if self.maskType == xmipp3.MASK_RAISED_COSINE:
            args += " -%(maskRadius)d -%(maskBand)d" 
        elif self.maskType == xmipp3.MASK_CIRCULAR:
            args += " -%(maskRadius)d"
        else:
            args += self.maskFile.get()
            
        if outputFn != inputFn:
            args += " -o %(outputFn)s"
                        
        self._insertRunJobStep(self._program, args % locals())

        
class XmippProtResize(XmippProtFilter):
    """ This class implement a protocol to change dimesions of the particles with Xmipp.
    """
    def __init__(self):
        XmippProtFilter.__init__(self)
        self._program = "xmipp_transform_mask"
        
    def _defineProcessParams(self, form):        
        form.addParam('resizeOperation', EnumParam, 
                      choices=['resize', 'crop'], 
                      default=0,
                      label="Resize operation", display=EnumParam.DISPLAY_COMBO,
                      help='Select how do you want to change the size of the particles. \n '
                      '<resize>: you will provide the new size (in pixels) for your particles. \n '
                      '<crop>: you choose how many pixels you want to crop from each border. \n ')
        
        form.addParam('newSize', IntParam, default=0,
                      condition='resizeOperation == 0',
                      label='New image size (px)',
                      help='This is the size in pixels of the particle images.')       
    
        form.addParam('cropSize', IntParam, default=0,
                      condition='resizeOperation == 1',
                      label='Crop size (px)',
                      help='This is the desired output size(in pixels) after cropping.') 

    def _insertFilterStep(self, inputFn, outputFn, outputMd):
        #TODO: Check if we want to separate crop from resize
        raise Exception('Not yet implemented.')
