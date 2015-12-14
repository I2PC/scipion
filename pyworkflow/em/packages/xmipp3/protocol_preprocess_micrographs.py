# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This sub-package contains the XmippPreprocessMicrographs protocol
"""

from os.path import basename
from pyworkflow.utils import getExt, replaceExt
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, FloatParam, LabelParam
from pyworkflow.em.protocol import ProtPreprocessMicrographs



class XmippProtPreprocessMicrographs(ProtPreprocessMicrographs):
    """Protocol to preprocess a set of micrographs in the project. You can crop borders, remove bad pixels, etc. """
    _label = 'preprocess micrographs'


    def __init__(self, **args):        
        ProtPreprocessMicrographs.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE params functions --------------------------------------------
    
    def _defineParams(self, form):
        form.addSection(label='Preprocess')
        form.addParam('inputMicrographs', PointerParam, 
                      pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True,
                      help='Select the SetOfMicrograph to be preprocessed.')
        
        form.addParam('orderComment', LabelParam, label="Operations are performed in the order shown below", important=True)
        form.addParam('doCrop', BooleanParam, default=False,
                      label='Crop borders?', 
                      help='Crop a given amount of pixels from each border.')
        form.addParam('cropPixels', IntParam, default=10, condition='doCrop',
                      label='Pixels to crop',
                      help='Amount of pixels you want to crop from borders.')
        form.addParam('doLog', BooleanParam, default=False,
                      label='Take logarithm?', 
                      help='Depending on your acquisition system you may need to take the logarithm '
                           'of the pixel values in order to have a linear relationship between '
                           'the gray values in the image and those in the volume. a - b ln(x+c) '
                           'by default 4.431-0.4018*LN((P1+336.6)) is applied (right one for nikon coolscan 9000)')
        line = form.addLine('Log', condition='doLog', 
                            help='Parameters in a - b ln(x+c).')
        line.addParam('logA', FloatParam, default=4.431, label='a')
        line.addParam('logB', FloatParam, default=0.402, label='b')
        line.addParam('logC', FloatParam, default=336.6, label='c')
        
        form.addParam('doRemoveBadPix', BooleanParam, default=False,
                      label='Remove bad pixels?',
                      help='Values will be thresholded to this multiple of standard deviations. '
                           'Typical values are about 5, i.e., pixel values beyond 5 times the '
                           'standard deviation will be substituted by the local median. '
                           'Set this option to -1 for not applying it.')
        form.addParam('mulStddev', IntParam, default=5, condition='doRemoveBadPix',
                      label='Multiple of Stddev',
                      help='Multiple of standard deviation.')    
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast?',
                      help='Multiply by -1')
        form.addParam('doDownsample', BooleanParam, default=False,
                      label='Downsample micrographs?',
                      help='Downsample micrographs by a given factor.')
        form.addParam('downFactor', FloatParam, default=2., condition='doDownsample',
                      label='Downsampling factor',
                      help='Non-integer downsample factors are possible. Must be larger than 1.')
        form.addParam('doSmooth', BooleanParam, default=False,
                      label='Gaussian filter',
                      help="Apply a Gaussian filter in real space")
        form.addParam('sigmaConvolution', FloatParam, default=2, condition="doSmooth",
                      label='Gaussian sigma (px)',
                      help="The larger this value, the more the effect will be noticed")
        form.addParam('doNormalize', BooleanParam, default=False,
                      label='Normalize micrograph?',
                      help='Normalize micrographs to be zero mean and standard deviation one')
    
        form.addParallelSection(threads=2, mpi=1)

    def _defineInputs(self):
        """ Store some of the input parameter in a dictionary for
        an easy replacement in the programs command line. 
        """
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get()
        
        # Parameters needed to preprocess the micrographs
        self.params = {'downFactor': self.downFactor.get(),
                       'cropPixels': 2 * self.cropPixels.get(),
                       'logA': self.logA.get(),
                       'logB': self.logB.get(),
                       'logC': self.logC.get(),
                       'stddev': self.mulStddev.get(),
                       'sigmaConvolution': self.sigmaConvolution.get()}
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Insert one or many processing steps per micrograph. """
        self._defineInputs()
        # For each micrograph insert the steps to preprocess it
        pre = []
        for mic in self.inputMics:
            fnOut = self._getOutputMicrograph(mic)
            stepId = self._insertStepsForMicrograph(mic.getFileName(), fnOut)
            pre.append(stepId)
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep', prerequisites=pre)
    
    def _insertStepsForMicrograph(self, inputMic, outputMic):
        self.params['inputMic'] = inputMic
        self.params['outputMic'] = outputMic
        self.lastStepId = None
        self.prerequisites = [] # First operation should not depend on nothing before
        # Crop
        self.__insertOneStep(self.doCrop, "xmipp_transform_window",
                            " -i %(inputMic)s --crop %(cropPixels)d -v 0")
        # Take logarithm
        self.__insertOneStep(self.doLog, "xmipp_transform_filter",
                            " -i %(inputMic)s --log --fa %(logA)f --fb %(logB)f --fc %(logC)f")
        # Remove bad pixels
        self.__insertOneStep(self.doRemoveBadPix, "xmipp_transform_filter",
                            " -i %(inputMic)s --bad_pixels outliers %(stddev)f -v 0")
        # Invert
        self.__insertOneStep(self.doInvert, "xmipp_image_operate",
                            "-i %(inputMic)s --mult -1")
        # Downsample
        self.__insertOneStep(self.doDownsample, "xmipp_transform_downsample",
                            "-i %(inputMic)s --step %(downFactor)f --method fourier")
        # Smooth
        self.__insertOneStep(self.doSmooth, "xmipp_transform_filter",
                            "-i %(inputMic)s --fourier real_gaussian %(sigmaConvolution)f")
        # Normalize
        self.__insertOneStep(self.doNormalize, "xmipp_transform_normalize",
                            "-i %(inputMic)s --method OldXmipp")
        return self.lastStepId
    
    def __insertOneStep(self, condition, program, arguments):
        """Insert operation if the condition is met.
        Possible conditions are: doDownsample, doCrop...etc"""
        if condition.get():
            # If the input micrograph and output micrograph differss, 
            # add the -o option
            if self.params['inputMic'] != self.params['outputMic']:
                arguments += " -o %(outputMic)s"
            # Insert the command with the formatted parameters
            self.lastStepId = self._insertRunJobStep(program, arguments % self.params,
                                                     prerequisites=self.prerequisites)
            self.prerequisites = [self.lastStepId] # next should depend on this step
            # Update inputMic for next step as outputMic
            self.params['inputMic'] = self.params['outputMic']
    
    #--------------------------- STEPS functions ---------------------------------------------------
    
    def createOutputStep(self):        
        outputMics = self._createSetOfMicrographs()
        outputMics.copyInfo(self.inputMics)
        
        if self.doDownsample.get():
            outputMics.setDownsample(self.downFactor.get())

        for mic in self.inputMics:
            # Update micrograph name and append to the new Set
            mic.setFileName(self._getOutputMicrograph(mic))
            outputMics.append(mic)

        self._defineOutputs(outputMicrographs=outputMics)
        self._defineTransformRelation(self.inputMics, outputMics)
        
    #--------------------------- INFO functions ----------------------------------------------------
    
    def _validate(self):
        validateMsgs = []
        # Some prepocessing option need to be marked
        if not(self.doCrop or self.doDownsample or self.doLog or self.doRemoveBadPix or self.doInvert or self.doNormalize or self.doSmooth):
            validateMsgs.append('Some preprocessing option need to be selected.')
        return validateMsgs
    
    def _citations(self):
        return ["Sorzano2009d"]
    
    def _hasOutput(self):
        return (getattr(self, 'outputMicrographs', False) and
                self.outputMicrographs.hasValue())
                
    def _summary(self):
        if not self._hasOutput():
            return ["*Output Micrographs* not ready yet."]

        summary = []
        summary.append("Micrographs preprocessed: *%d*" % self.inputMicrographs.get().getSize())
        if self.doCrop:
            summary.append("Cropped *%d* pixels per each border." % self.cropPixels)
        if self.doLog:
            summary.append("Formula applied: %f - %f ln(x + %f)" % (self.logA, self.logB, self.logC,))
        if self.doRemoveBadPix:
            summary.append("Multiple of standard deviation to remove pixels: %d" % self.mulStddev)
        if self.doInvert:
            summary.append("Contrast inverted")
        if self.doDownsample:
            summary.append("Downsampling factor: %0.2f" % self.downFactor)
        if self.doSmooth:
            summary.append("Gaussian filtered with sigma=%f (px)"%self.sigmaConvolution.get())
        if self.doNormalize:
            summary.append("Normalized to mean 0 and variance 1")
        return summary
    
    def _methods(self):
        if not self._hasOutput():
            return ['*Output micrographs* not ready yet.']

        txt = "The micrographs in set %s have " % self.getObjectTag('inputMicrographs')

        if self.doCrop:
            txt += "been cropped by %d pixels " % self.cropPixels
        if self.doLog:
            txt += ("changed from transmisivity to density with the formula: "
                    "%f - %f * ln(x + %f) " % (self.logA, self.logB, self.logC))
        if self.doRemoveBadPix:
            txt += "had pixels removed, the ones with standard deviation beyond %d " % self.mulStddev.get()
        if self.doRemoveBadPix:
            txt += "contrast inverted "
        if self.doDownsample:
            txt += "been downsampled with a factor of %0.2f " % self.downFactor.get()
        if self.doSmooth:
            txt += "been Gaussian filtered with a sigma of %0.2f pixels "%self.sigmaConvolution.get()
        if self.doNormalize:
            txt += "been normalized to mean 0 and variance 1"

        return [txt, "The resulting set of micrographs is %s" %
                self.getObjectTag('outputMicrographs')]

    #--------------------------- UTILS functions --------------------------------------------
    def _getOutputMicrograph(self, mic):
        """ Return the name of the output micrograph, given
        the input Micrograph object.
        """
        fn = mic.getFileName()
        extFn = getExt(fn)
        if extFn != ".mrc":
            fn = replaceExt(fn, "mrc")
        
        fnOut = self._getExtraPath(basename(fn))
        
        return fnOut
