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


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from data import *


class XmippDefPreprocessMicrograph(Form):
    """Create the definition of parameters for
    the XmippPreprocessMicrographs protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Preprocess')
        self.addParam('inputMicrographs', PointerParam, label="Input micrographs", isImportant=True,
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrograph to be preprocessed.')
        
        self.addParam('doCrop', BooleanParam, default=False,
                      label='Crop borders?', 
                      help='Crop a given amount of pixels from each border.')
        self.addParam('cropPixels', IntParam, default=10, condition='doCrop',
                      label='Pixels to crop',
                      help='Amount of pixels you want to crop from borders.')
        self.addParam('doLog', BooleanParam, default=False,
                      label='Take logarithm?', 
                      help='Depending on your acquisition system you may need to take the logarithm '
                           'of the pixel values in order to have a linear relationship between '
                           'the gray values in the image and those in the volume. a - b ln(x+c) '
                           'by default 4.431-0.4018*LN((P1+336.6)) is applied (right one for nikon coolscan 9000)')
        self.addParam('logA', FloatParam, default=4.431, condition='doLog',
                      label='a',
                      help='Parameter a in a - b ln(x+c).')
        self.addParam('logB', FloatParam, default=0.4018, condition='doLog',
                      label='b',
                      help='Parameter b in a - b ln(x+c).')
        self.addParam('logC', FloatParam, default=336.6, condition='doLog',
                      label='c',
                      help='Parameter c in a - b ln(x+c).')
        self.addParam('doRemoveBadPix', BooleanParam, default=False,
                      label='Remove bad pixels?',
                      help='Values will be thresholded to this multiple of standard deviations. '
                           'Typical values are about 5, i.e., pixel values beyond 5 times the '
                           'standard deviation will be substituted by the local median. '
                           'Set this option to -1 for not applying it.')
        self.addParam('mulStddev', IntParam, default=5, condition='doRemoveBadPix',
                      label='Multiple of Stddev',
                      help='Multiple of standard deviation.')    
        self.addParam('doDownsample', BooleanParam, default=False,
                      label='Downsample micrographs?',
                      help='Downsample micrographs by a given factor.')
        self.addParam('downFactor', FloatParam, default=2., condition='doDownsample',
                      label='Downsampling factor',
                      help='Non-integer downsample factors are possible. Must be larger than 1.')
    
        self.addParallelSection(threads=2, mpi=1)
        
class XmippProtPreprocessMicrographs(ProtPreprocessMicrographs):
    """Protocol to preprocess a set of micrographs in the project"""
    _definition = XmippDefPreprocessMicrograph()
    _label = 'Xmipp Micrographs Preprocessing'
    

    def __init__(self, **args):        
        ProtPreprocessMicrographs.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    def _defineSteps(self):
        '''for each micrograph insert the steps to preprocess it
        '''
        print "En defineSteps"
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        
        IOTable = {}
        
        # Parameters needed to preprocess the micrographs
        self.params = {'downFactor': self.downFactor.get(),
                       'cropPixels': self.cropPixels.get(),
                       'logA': self.logA.get(),
                       'logB': self.logB.get(),
                       'logC': self.logC.get(),
                       'stddev': self.mulStddev.get()}
        
        # For each micrograph insert the steps to preprocess it
        pre = []
        for mic in self.inputMics:
            fn = mic.getFileName()
            fnOut = self._getPath(os.path.basename(fn))
            stepId = self.__insertStepsForMicrograph(fn, fnOut)
            pre.append(stepId)
            IOTable[fn] = fnOut
        
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput', IOTable, prerequisites=pre)
        
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
            
    def __insertStepsForMicrograph(self, inputMic, outputMic):
        self.params['inputMic'] = inputMic
        self.params['outputMic'] = outputMic
        self.lastStepId = None
        self.prerequisites = [] # First operation should not depend on nothing before
        
        # Downsample
        self.__insertOneStep(self.doDownsample, "xmipp_transform_downsample",
                            "-i %(inputMic)s --step %(downFactor)f --method fourier")
        # Crop
        self.__insertOneStep(self.doCrop, "xmipp_transform_window",
                            " -i %(inputMic)s --crop %(cropPixels)d -v 0")
        # Take logarithm
        self.__insertOneStep(self.doLog, "xmipp_transform_filter",
                            " -i %(inputMic)s --log --fa %(logA)f --fb %(logB)f --fc %(logC)f")
        # Remove bad pixels
        self.__insertOneStep(self.doRemoveBadPix, "xmipp_transform_filter",
                            " -i %(inputMic)s --bad_pixels outliers %(stddev)f -v 0")
                
        return self.lastStepId
    
    def createOutput(self, IOTable):
        
        path = self._getPath('micrographs.sqlite')
        micSet = SetOfMicrographs()
        micSet.setFileName(path)

        mapsId = {}
            
        for mic in self.inputMics:
            outMicFn = Micrograph()
            outMicFn.setFileName(IOTable[mic.getFileName()])
            # Updating micrograph name
            outMicFn.setId(mic.getId())
            micSet.append(outMicFn)
            #mapsId[mic.getId()] = xmicFn.getId()
            
        micSet.copyInfo(self.inputMics)
        
        if self.doDownsample.get():
            micSet.setDownsample(self.downFactor.get())

        micSet.write()
        self._defineOutputs(outputMicrographs=micSet)

    def _summary(self):
        summary = []
        if not self.inputMicrographs.hasValue():
            summary.append("No <Input Micrographs> selected.")
        else:
            summary.append("Input micrographs: " + self.inputMicrographs.get().getNameId())
            if self.doDownsample.get():
                summary.append("Downsampling factor: %d" % self.downFactor.get())
            if self.doCrop.get():
                summary.append("Number of pixels cropped: %d" % self.cropPixels.get())
            if self.doLog.get():
                summary.append("Formula applied: %f - %f ln(x + %f)" % (self.logA.get(), self.logB.get(), self.logC.get(),))
            if self.doRemoveBadPix.get():
                summary.append("Number of pixels removed: %d" % self.mulStddev.get())
        return summary
    
    def _validate(self):
        validateMsgs = []
        # Some prepocessing option need to be marked
        if not(self.doCrop or self.doDownsample or self.doLog or self.doRemoveBadPix):
            validateMsgs.append('Some preprocessing option need to be selected.')
        return validateMsgs
