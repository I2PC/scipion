# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco   (jgomez@cnb.csic.es)
# *              Joaquin Oton   (oton@cnb.csic.es)
# *              Airen Zaldivar (azaldivar@cnb.csic.es)
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

from pyworkflow.em import *  
from pyworkflow.utils import *  
from protocol_process import XmippProcess, XmippProcessParticles, XmippProcessVolumes
from pyworkflow.em.constants import *
from constants import *

RESIZE_SAMPLINGRATE = 0
RESIZE_DIMENSIONS = 1
RESIZE_FACTOR = 2
RESIZE_PYRAMID = 3



class XmippProtResize():
    """ This class implement the common features to change dimensions of either SetOfParticles, Volume or SetOfVolumes objects.
    """
    _inputLabel = None # This should be 'particles' or 'volumes'
    
    def __init__(self, **args):
        self._programWindow = "xmipp_transform_window"
        self._programResize = "xmipp_image_resize"
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        # Resize operation
        form.addParam('doResize', BooleanParam, default=False,
                      label='Resize %s?' % self._inputLabel,
                      help='If you set to *Yes*, you should provide a resize option.')
        form.addParam('resizeOption', EnumParam,
                      choices=['Sampling Rate', 'Dimensions', 'Factor', 'Pyramid'],
                      condition='doResize',
                      default=RESIZE_SAMPLINGRATE,
                      label="Resize option", display=EnumParam.DISPLAY_COMBO,
                      help='Select an option to resize the images: \n '
                      '_Sampling Rate_: Set the desire sampling rate to resize. \n'
                      '_Dimensions_: Set the output dimensions. Resize operation can be done in Fourier space.\n'
                      '_Factor_: Set a resize factor to resize. \n '
                      '_Pyramid_: Use positive level value to expand and negative to reduce. \n')
        form.addParam('resizeSamplingRate', FloatParam, default=1.0,
                      condition='doResize and resizeOption==%d' % (RESIZE_SAMPLINGRATE),
                      label='Resize sampling rate (A/px)',
                      help='Set the new output sampling rate.')
        form.addParam('doFourier', BooleanParam, default=False,
                      condition='doResize and resizeOption==%d' % (RESIZE_DIMENSIONS),
                      label='Use fourier method to resize?',
                      help='If you set to *True*, the final dimensions must be lower than the original ones.')
        form.addParam('resizeDim', IntParam, default=0,
                      condition='doResize and resizeOption==%d' % (RESIZE_DIMENSIONS),
                      label='New image size (px)',
                      help='Size in pixels of the particle images <x> <y=x> <z=x>.')
        form.addParam('resizeFactor', FloatParam, default=0.5,
                      condition='doResize and resizeOption==%d' % (RESIZE_FACTOR),
                      label='Resize factor',
                      help='New size is the old one x resize factor.')
        form.addParam('resizeLevel', IntParam, default=0,
                      condition='doResize and resizeOption==%d' % (RESIZE_PYRAMID),
                      label='Pyramid level',
                      help='Use positive value to expand and negative to reduce.')
        # Window operation
        form.addParam('doWindow', BooleanParam, default=False,
                      label='Apply a window operation?',
                      help='If you set to *Yes*, you should provide a window option.')
        form.addParam('windowOperation', EnumParam,
                      choices=['crop', 'window'],
                      condition='doWindow',
                      default=1,
                      label="Window operation", display=EnumParam.DISPLAY_COMBO,
                      help='Select how do you want to change the size of the particles. \n '
                      '_resize_: you will provide the new size (in pixels) for your particles. \n '
                      '_crop_: you choose how many pixels you want to crop from each border. \n ')
        form.addParam('cropSize', IntParam, default=0,
                      condition='doWindow and windowOperation == 0',
                      label='Crop size (px)',
                      help='This is the amount of pixels cropped in each border. \n '
                           'e.g: if you set 10 pixels, the dimensions of the\n'
                           'object (SetOfParticles, Volume or SetOfVolumes) will be\n'
                           'reduce in 20 pixels (2 borders * 10 pixels)')
        form.addParam('windowSize', IntParam, default=0,
                      condition='doWindow and windowOperation == 1',
                      label='Window size (px)',
                      help='This is the size in pixels of the particle images.')
                
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProcessStep(self):
        isFirstStep = True
        
        if self.doResize:
            isFirstStep = False
            args = self._resizeArgs()
            self._insertFunctionStep("resizeStep", args)
            
        if self.doWindow:
            args = self._windowArgs(isFirstStep)
            self._insertFunctionStep("windowStep", args)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def resizeStep(self, args):
        self.runJob(self._programResize, args)
    
    def windowStep(self, args):
        self.runJob(self._programWindow, args)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        
        if self.doResize and self.resizeOption.get() == RESIZE_SAMPLINGRATE and self.doFourier:
#             imgSet = self.inputParticles.get()
            size = self._getSetSize()
            if self.resizeDim.get() > size:
                errors.append('Fourier resize method cannot be used to increase the dimensions')
                
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeCommonArgs(self):
        samplingRate = self._getSetSampling()
        inputFn = self.inputFn
        if self.resizeOption == RESIZE_SAMPLINGRATE:
            newSamplingRate = self.resizeSamplingRate.get()
            factor = samplingRate / newSamplingRate
            self.samplingRate = newSamplingRate
            args = self._args + " --factor %(factor)f"
        
        elif self.resizeOption == RESIZE_DIMENSIONS:
            size = self.resizeDim.get()
            dim = self._getSetSize()
            self.samplingRate = (samplingRate * float(dim) / float(size))
            
            if self.doFourier:
                args = self._args + " --fourier %(size)d"
            else:
                args = self._args + " --dim %(size)d"
            
        elif self.resizeOption == RESIZE_FACTOR:
            factor = self.resizeFactor.get()                                               
            self.samplingRate = samplingRate / factor
            args = self._args + " --factor %(factor)f"
        
        else:
            level = self.resizeLevel.get()
            factor = pow(2, level)
            self.samplingRate = samplingRate / factor
            args = self._args + " --pyramid %(level)d"
            
        return args % locals()
    
    def _windowCommonArgs(self):
        dim = self._getSetSize()
        if self.getEnumText('windowOperation') == "crop":
            cropSize = self.cropSize.get() * 2
            windowSize = dim - cropSize
            args = " --crop %(cropSize)s "
        else:
            windowSize = self.windowSize.get()
            args = " --size %(windowSize)s"
        
        self.newWindowSize = windowSize
        
        return args % locals()
    
    def _getSize(self, imgSet):
        """ get the size of an object"""
        if isinstance(imgSet, Volume):
            Xdim = imgSet.getDim()[0]
        else:
            Xdim = imgSet.getDimensions()[0]
        return Xdim
    
    def _getSampling(self, imgSet):
        """ get the sampling rate of an object"""
        samplingRate = imgSet.getSamplingRate()
        return samplingRate


class XmippProtCropResizeParticles(ProtProcessParticles, XmippProtResize, XmippProcessParticles):
    """ Crop or resize a set of particles """
    _label = 'crop/resize particles'
    _inputLabel = 'particles'
    
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
        XmippProtResize.__init__(self)
        XmippProcessParticles.__init__(self)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtResize._defineProcessParams(self, form)
     
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self):
        """ convert if necessary"""
        pass
#         imgSet = self.inputParticles.get()
#         imgSet.writeStack(self.outputStk)
    
#     def createOutputStep(self):
#         inImgSet = self.inputParticles.get()
#         outImgSet = self._createSetOfParticles()
#         outImgSet.copyInfo(inImgSet)
#         if self.doResize:
#             outImgSet.setSamplingRate(self.samplingRate)
#         
#         for i, img in enumerate(inImgSet):
#             j = i + 1
#             img.setLocation(j, self.outputStk)
#             if self.doResize:
#                 img.setSamplingRate(self.samplingRate)
#             outImgSet.append(img)
#         
#         self._defineOutputs(outputParticles=outImgSet)
#         self._defineTransformRelation(inImgSet, self.outputParticles)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            sampling = self._getSampling(self.outputParticles)
            size = self._getSize(self.outputParticles)
            if self.doResize:
                summary.append("The sampling rate of the output particles are: %0.3f" % sampling)
            if self.doWindow.get():
                if self.getEnumText('windowOperation') == "crop":
                    summary.append("*Crop operation*: New size %s" % size)
                else:
                    summary.append("*Window operation*: New size %s" % size)
        return summary
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeArgs(self):
        args = self._resizeCommonArgs()
        args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        return args
    
    def _windowArgs(self, isFirstStep):
        if isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += self._windowCommonArgs()
        return args
    
    def _getSetSize(self):
        """ get the size of SetOfParticles object"""
        imgSet = self.inputParticles.get()
        size = self._getSize(imgSet)
        return size
    
    def _getSetSampling(self):
        """ get the sampling rate of SetOfParticles object"""
        imgSet = self.inputParticles.get()
        samplingRate = self._getSampling(imgSet)
        return samplingRate


class XmippProtCropResizeVolumes(ProtPreprocessVolumes, XmippProtResize, XmippProcessVolumes):
    """ Crop or resize a set of volumes """
    _label = 'crop/resize volumes'
    _inputLabel = 'volumes'
    
    def __init__(self, **args):
        ProtPreprocessVolumes.__init__(self, **args)
        XmippProtResize.__init__(self, **args)
        XmippProcessVolumes.__init__(self)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtResize._defineProcessParams(self, form)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):
        volSet = self.inputVolumes.get()
        if self._isSingleInput():
            vol = Volume()
            vol.copyInfo(volSet)
            if self.doResize:
                vol.setSamplingRate(self.samplingRate)
            vol.setFileName(self.outputStk)
            self._defineOutputs(outputVol=vol)
        else:
            volumes = self._createSetOfVolumes()
            volumes.copyInfo(volSet)
            if self.doResize:
                volumes.setSamplingRate(self.samplingRate)
            for i, vol in enumerate(volSet):
                j = i + 1
                vol.setSamplingRate(self.samplingRate)
                vol.setLocation(j, self.outputStk)
                volumes.append(vol)
            self._defineOutputs(outputVol=volumes)

        self._defineTransformRelation(volSet, self.outputVol)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        
        if not hasattr(self, 'outputVol'):
            summary.append("Output volume(s) not ready yet.") 
        else:
            sampling = self._getSampling(self.outputVol)
            size = self._getSize(self.outputVol)
            if self.doResize:
                summary.append("The sampling rate of the output voume(s) are: %0.3f" % sampling)
            if self.doWindow.get():
                if self.getEnumText('windowOperation') == "crop":
                    summary.append("*Crop operation*: New size %d" % size)
                else:
                    summary.append("*Window operation*: New size %d" % size)
        return summary
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeArgs(self):
        args = self._resizeCommonArgs()
        if self._isSingleInput():
            args += " -o %s" % self.outputStk
        else:
            args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        return args
    
    def _windowArgs(self, isFirstStep):
        if isFirstStep:
            if self._isSingleInput():
                args = "-i %s -o %s" % (self.inputFn, self.outputStk)
            else:
                args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += self._windowCommonArgs()
        return args
    
    def _getSetSize(self):
        """ get the size of either Volume or SetOfVolumes objects"""
        imgSet = self.inputVolumes.get()
        size = self._getSize(imgSet)
        return size
    
    def _getSetSampling(self):
        """ get the sampling rate of either Volume or SetOfVolumes objects"""
        imgSet = self.inputVolumes.get()
        samplingRate = self._getSampling(imgSet)
        return samplingRate
