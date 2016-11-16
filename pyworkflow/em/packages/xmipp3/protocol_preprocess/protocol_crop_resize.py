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

from pyworkflow.protocol.params import BooleanParam, EnumParam, FloatParam, IntParam
from pyworkflow.em.data import Volume

from protocol_process import XmippProcessParticles, XmippProcessVolumes



class XmippResizeHelper():
    """ Common features to change dimensions of either SetOfParticles, Volume or SetOfVolumes objects.
    """

    RESIZE_SAMPLINGRATE = 0
    RESIZE_DIMENSIONS = 1
    RESIZE_FACTOR = 2
    RESIZE_PYRAMID = 3

    WINDOW_OP_CROP = 0
    WINDOW_OP_WINDOW = 1

    #--------------------------- DEFINE param functions --------------------------------------------
    @classmethod
    def _defineProcessParams(cls, protocol, form):
        # Resize operation
        form.addParam('doResize', BooleanParam, default=False,
                      label='Resize %s?' % protocol._inputLabel,
                      help='If you set to *Yes*, you should provide a resize option.')
        form.addParam('resizeOption', EnumParam,
                      choices=['Sampling Rate', 'Dimensions', 'Factor', 'Pyramid'],
                      condition='doResize',
                      default=cls.RESIZE_SAMPLINGRATE,
                      label="Resize option", display=EnumParam.DISPLAY_COMBO,
                      help='Select an option to resize the images: \n '
                      '_Sampling Rate_: Set the desire sampling rate to resize. \n'
                      '_Dimensions_: Set the output dimensions. Resize operation can be done in Fourier space.\n'
                      '_Factor_: Set a resize factor to resize. \n '
                      '_Pyramid_: Use positive level value to expand and negative to reduce. \n')
        form.addParam('resizeSamplingRate', FloatParam, default=1.0,
                      condition='doResize and resizeOption==%d' % cls.RESIZE_SAMPLINGRATE,
                      label='Resize sampling rate (A/px)',
                      help='Set the new output sampling rate.')
        form.addParam('doFourier', BooleanParam, default=False,
                      condition='doResize and resizeOption==%d' % cls.RESIZE_DIMENSIONS,
                      label='Use fourier method to resize?',
                      help='If you set to *True*, the final dimensions must be lower than the original ones.')
        form.addParam('resizeDim', IntParam, default=0,
                      condition='doResize and resizeOption==%d' % cls.RESIZE_DIMENSIONS,
                      label='New image size (px)',
                      help='Size in pixels of the particle images <x> <y=x> <z=x>.')
        form.addParam('resizeFactor', FloatParam, default=0.5,
                      condition='doResize and resizeOption==%d' % cls.RESIZE_FACTOR,
                      label='Resize factor',
                      help='New size is the old one x resize factor.')
        form.addParam('resizeLevel', IntParam, default=0,
                      condition='doResize and resizeOption==%d' % cls.RESIZE_PYRAMID,
                      label='Pyramid level',
                      help='Use positive value to expand and negative to reduce.')
        # Window operation
        form.addParam('doWindow', BooleanParam, default=False,
                      label='Apply a window operation?',
                      help='If you set to *Yes*, you should provide a window option.')
        form.addParam('windowOperation', EnumParam,
                      choices=['crop', 'window'],
                      condition='doWindow',
                      default=cls.WINDOW_OP_WINDOW,
                      label="Window operation", display=EnumParam.DISPLAY_COMBO,
                      help='Select how to change the size of the particles.\n'
                      '_cls.RESIZE_: provide the new size (in pixels) for your particles.\n'
                      '_crop_: choose how many pixels to crop from each border.\n')
        form.addParam('cropSize', IntParam, default=0,
                      condition='doWindow and windowOperation == %d' % cls.WINDOW_OP_CROP,
                      label='Crop size (px)',
                      help='Amount of pixels cropped from each border.\n'
                           'e.g: if you set 10 pixels, the dimensions of the\n'
                           'object (SetOfParticles, Volume or SetOfVolumes) will be\n'
                           'reduced in 20 pixels (2 borders * 10 pixels)')
        form.addParam('windowSize', IntParam, default=0,
                      condition='doWindow and windowOperation == %d' % cls.WINDOW_OP_WINDOW,
                      label='Window size (px)',
                      help='Size in pixels of the output object. It will be '
                           'expanded or cutted in all directions such that the '
                           'origin remains the same.')

    #--------------------------- INSERT steps functions --------------------------------------------
    @classmethod
    def _insertProcessStep(cls, protocol):
        isFirstStep = True
        
        if protocol.doResize:
            isFirstStep = False
            args = protocol._resizeArgs()
            protocol._insertFunctionStep("resizeStep", args)
            
        if protocol.doWindow:
            args = protocol._windowArgs(isFirstStep)
            protocol._insertFunctionStep("windowStep", args)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def resizeStep(self, args):
        self.runJob("xmipp_image_resize", args)
    
    def windowStep(self, args):
        self.runJob("xmipp_transform_window", args, numberOfMpi=1)
        
    #--------------------------- INFO functions ----------------------------------------------------
    @classmethod
    def _validate(cls, protocol):
        errors = []
        
        if protocol.doResize and protocol.resizeOption == cls.RESIZE_SAMPLINGRATE and protocol.doFourier:
#             imgSet = self.inputParticles.get()
            size = protocol._getSetSize()
            if protocol.resizeDim > size:
                errors.append('Fourier resize method cannot be used to increase the dimensions')
                
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    @classmethod
    def _resizeCommonArgs(cls, protocol):
        samplingRate = protocol._getSetSampling()
        inputFn = protocol.inputFn
        
        if protocol.resizeOption == cls.RESIZE_SAMPLINGRATE:
            newSamplingRate = protocol.resizeSamplingRate.get()
            factor = samplingRate / newSamplingRate
            args = protocol._args + " --factor %(factor)f"
        elif protocol.resizeOption == cls.RESIZE_DIMENSIONS:
            size = protocol.resizeDim.get()
            dim = protocol._getSetSize()
            factor = float(size) / float(dim)
            newSamplingRate = samplingRate / factor
            
            if protocol.doFourier:
                args = protocol._args + " --fourier %(size)d"
            else:
                args = protocol._args + " --dim %(size)d"
        elif protocol.resizeOption == cls.RESIZE_FACTOR:
            factor = protocol.resizeFactor.get()
            newSamplingRate = samplingRate / factor
            args = protocol._args + " --factor %(factor)f"
        elif protocol.resizeOption == cls.RESIZE_PYRAMID:
            level = protocol.resizeLevel.get()
            factor = 2**level
            newSamplingRate = samplingRate / factor
            args = protocol._args + " --pyramid %(level)d"
            
        protocol.samplingRate = newSamplingRate
        protocol.factor = factor
        
        return args % locals()
    
    @classmethod
    def _windowCommonArgs(cls, protocol):
        op = protocol.getEnumText('windowOperation')
        if op == "crop":
            cropSize2 = protocol.cropSize.get() * 2
            protocol.newWindowSize = protocol._getSetSize() - cropSize2
            return " --crop %d " % cropSize2
        elif op == "window":
            windowSize = protocol.windowSize.get()
            protocol.newWindowSize = windowSize
            return " --size %d " % windowSize


def _getSize(imgSet):
    """ get the size of an object"""
    if isinstance(imgSet, Volume):
        Xdim = imgSet.getDim()[0]
    else:
        Xdim = imgSet.getDimensions()[0]
    return Xdim

def _getSampling(imgSet):
    """ get the sampling rate of an object"""
    samplingRate = imgSet.getSamplingRate()
    return samplingRate


class XmippProtCropResizeParticles(XmippProcessParticles):
    """ Crop or resize a set of particles """
    _label = 'crop/resize particles'
    _inputLabel = 'particles'
    
    def __init__(self, **kwargs):
        XmippProcessParticles.__init__(self, **kwargs)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippResizeHelper._defineProcessParams(self, form)
        form.addParallelSection(threads=0, mpi=8)
        
    def _insertProcessStep(self):
        XmippResizeHelper._insertProcessStep(self)
     
    #--------------------------- STEPS functions ---------------------------------------------------
    def resizeStep(self, args):
        self.runJob("xmipp_image_resize", args)
    
    def windowStep(self, args):
        self.runJob("xmipp_transform_window", args, numberOfMpi=1)
        
    def _preprocessOutput(self, output):
        """ We need to update the sampling rate of the 
        particles if the Resize option was used.
        """
        self.inputHasAlign = self.inputParticles.get().hasAlignment()
        
        if self.doResize:
            output.setSamplingRate(self.samplingRate)
            
    def _updateItem(self, item, row):
        """ Update also the sampling rate and 
        the alignment if needed.
        """
        XmippProcessParticles._updateItem(self, item, row)
        if self.doResize:
            if item.hasCoordinate():
                item.scaleCoordinate(self.factor)
            item.setSamplingRate(self.samplingRate)
            if self.inputHasAlign:
                item.getTransform().scaleShifts(self.factor)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output images not ready yet.") 
        else:
            sampling = _getSampling(self.outputParticles)
            size = _getSize(self.outputParticles)
            if self.doResize:
                summary.append("Output particles have a different sampling "
                               "rate (pixel size): *%0.3f* A/px" % sampling)
                summary.append("Resizing method: *%s*" %
                               self.getEnumText('resizeOption'))
            if self.doWindow:
                if self.getEnumText('windowOperation') == "crop":
                    summary.append("The particles were cropped.")
                else:
                    summary.append("The particles were windowed.")
                summary.append("New size: *%s* px" % size)
        return summary

    def _methods(self):
        methods = ["We took input particles %s of size %d " % (self.getObjectTag('inputParticles'), len(self.inputParticles.get()))]
        if self.doWindow:
            if self.getEnumText('windowOperation') == "crop":
                methods += ["cropped them"]
            else:
                methods += ["windowed them"]
        if self.doResize:
            outputParticles = getattr(self, 'outputParticles', None)
            if outputParticles is None or outputParticles.getDim() is None:
                methods += ["Output particles not ready yet."]
            else:
                methods += ['resized them to %d px using the "%s" method%s' %
                            (outputParticles.getDim()[0],
                             self.getEnumText('resizeOption'),
                             " in Fourier space" if self.doFourier else "")]
        if not self.doResize and not self.doWindow:
            methods += ["did nothing to them"]
        str = "%s and %s. Output particles: %s" % (", ".join(methods[:-1]), methods[-1], self.getObjectTag('outputParticles'))
        return [str]

    def _validate(self):
        return XmippResizeHelper._validate(self)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeArgs(self):
        args = XmippResizeHelper._resizeCommonArgs(self)
        args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        return args
    
    def _windowArgs(self, isFirstStep):
        if isFirstStep:
            args = "-i %s -o %s --save_metadata_stack %s --keep_input_columns" % (self.inputFn, self.outputStk, self.outputMd)
        else:
            args = "-i %s" % self.outputStk
        args += XmippResizeHelper._windowCommonArgs(self)
        return args
    
    def _getSetSize(self):
        """ get the size of SetOfParticles object"""
        imgSet = self.inputParticles.get()
        size = _getSize(imgSet)
        return size
    
    def _getSetSampling(self):
        """ get the sampling rate of SetOfParticles object"""
        imgSet = self.inputParticles.get()
        samplingRate = _getSampling(imgSet)
        return samplingRate
    
    def _getDefaultParallel(self):
        """ Return the default value for thread and MPI
        for the parallel section definition.
        """
        return (0, 1)



class XmippProtCropResizeVolumes(XmippProcessVolumes):
    """ Crop or resize a set of volumes """
    _label = 'crop/resize volumes'
    _inputLabel = 'volumes'
    
    def __init__(self, **kwargs):
        XmippProcessVolumes.__init__(self, **kwargs)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippResizeHelper._defineProcessParams(self, form)
        
    def _insertProcessStep(self):
        XmippResizeHelper._insertProcessStep(self)
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def resizeStep(self, args):
        self.runJob("xmipp_image_resize", args)
    
    def windowStep(self, args):
        self.runJob("xmipp_transform_window", args, numberOfMpi=1)
        
    def _preprocessOutput(self, volumes):
        # We use the preprocess only whne input is a set
        # we do not use postprocess to setup correctly
        # the samplingRate before each volume is added
        if not self._isSingleInput():
            if self.doResize:
                volumes.setSamplingRate(self.samplingRate)

    def _postprocessOutput(self, volumes):
        # We use the postprocess only when input is a volume
        if self._isSingleInput():
            if self.doResize:
                volumes.setSamplingRate(self.samplingRate)            
     
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        
        if not hasattr(self, 'outputVol'):
            summary.append("Output volume(s) not ready yet.") 
        else:
            sampling = _getSampling(self.outputVol)
            size = _getSize(self.outputVol)
            if self.doResize:
                summary.append("Output volume(s) have a different sampling "
                               "rate (pixel size): *%0.3f* A/px" % sampling)
                summary.append("Resizing method: *%s*" %
                               self.getEnumText('resizeOption'))
            if self.doWindow.get():
                if self.getEnumText('windowOperation') == "crop":
                    summary.append("The volume(s) were cropped.")
                else:
                    summary.append("The volume(s) were windowed.")
                summary.append("New size: *%s* px" % size)
        return summary

    def _methods(self):
        if self._isSingleInput():
            methods = ["We took one volume"]
            pronoun = "it"
        else:
            methods = ["We took %d volumes" % self.inputVolumes.get().getSize()]
            pronoun = "them"
        if self.doWindow.get():
            if self.getEnumText('windowOperation') == "crop":
                methods += ["cropped %s" % pronoun]
            else:
                methods += ["windowed %s" % pronoun]
        if self.doResize:
            outputVol = getattr(self, 'outputVol', None)
            if outputVol is None or self.outputVol.getDim() is None:
                methods += ["Output volume not ready yet."]
            else:
                methods += ['resized %s to %d px using the "%s" method%s' %
                            (pronoun, self.outputVol.getDim()[0],
                             self.getEnumText('resizeOption'),
                             " in Fourier space" if self.doFourier else "")]
        if not self.doResize and not self.doWindow:
            methods += ["did nothing to %s" % pronoun]
            # TODO: does this case even work in the protocol?
        return ["%s and %s." % (", ".join(methods[:-1]), methods[-1])]

    def _validate(self):
        return XmippResizeHelper._validate(self)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeArgs(self):
        args = XmippResizeHelper._resizeCommonArgs(self)
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
        args += XmippResizeHelper._windowCommonArgs(self)
        return args
    
    def _getSetSize(self):
        """ get the size of either Volume or SetOfVolumes objects"""
        imgSet = self.inputVolumes.get()
        size = _getSize(imgSet)
        return size
    
    def _getSetSampling(self):
        """ get the sampling rate of either Volume or SetOfVolumes objects"""
        imgSet = self.inputVolumes.get()
        samplingRate = _getSampling(imgSet)
        return samplingRate
