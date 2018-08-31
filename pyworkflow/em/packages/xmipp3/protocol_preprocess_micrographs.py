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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import basename
from pyworkflow.utils import getExt, replaceExt
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED
import pyworkflow.protocol.constants as cons
from pyworkflow.protocol.params import (PointerParam, BooleanParam, IntParam,
                                        FloatParam, LabelParam)
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtPreprocessMicrographs
from pyworkflow.em.data import SetOfMicrographs
from pyworkflow.object import Set


class XmippProtPreprocessMicrographs(ProtPreprocessMicrographs):
    """Protocol to preprocess a set of micrographs in the project.
    You can crop borders, remove bad pixels, etc. """
    _label = 'preprocess micrographs'


    def __init__(self, **args):        
        ProtPreprocessMicrographs.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE params functions -----------------------
    
    def _defineParams(self, form):
        form.addSection(label='Preprocess')
        form.addParam('inputMicrographs', PointerParam, 
                      pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True,
                      help='Select the SetOfMicrograph to be preprocessed.')
        
        form.addParam('orderComment', LabelParam,
                      label="Operations are performed in the order shown below",
                      important=True)
        form.addParam('doCrop', BooleanParam, default=False,
                      label='Crop borders?', 
                      help='Crop a given amount of pixels from each border.')
        form.addParam('cropPixels', IntParam, default=10, condition='doCrop',
                      label='Pixels to crop',
                      help='Amount of pixels you want to crop from borders.')
        form.addParam('doLog', BooleanParam, default=False,
                      label='Take logarithm?', 
                      help='Depending on your acquisition system you may need '
                           'to take the logarithm of the pixel values in '
                           'order to have a linear relationship betweenthe '
                           'gray values in the image and those in the volume. '
                           'a - b ln(x+c) by default 4.431-0.4018*'
                           'LN((P1+336.6)) is applied (right one for nikon '
                           'coolscan 9000)')
        line = form.addLine('Log', condition='doLog', 
                            help='Parameters in a - b ln(x+c).')
        line.addParam('logA', FloatParam, default=4.431, label='a')
        line.addParam('logB', FloatParam, default=0.402, label='b')
        line.addParam('logC', FloatParam, default=336.6, label='c')
        
        form.addParam('doRemoveBadPix', BooleanParam, default=False,
                      label='Remove bad pixels?',
                      help='Values will be thresholded to this multiple of '
                           'standard deviations. Typical values are about 5, '
                           'i.e., pixel values beyond 5 times the standard '
                           'deviation will be substituted by the local median. '
                           'Set this option to -1 for not applying it.')
        form.addParam('mulStddev', IntParam, default=5,
                      condition='doRemoveBadPix',
                      label='Multiple of Stddev',
                      help='Multiple of standard deviation.')    
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast?',
                      help='Multiply by -1')
        form.addParam('doDownsample', BooleanParam, default=False,
                      label='Downsample micrographs?',
                      help='Downsample micrographs by a given factor.')
        form.addParam('downFactor', FloatParam, default=2.,
                      condition='doDownsample',
                      label='Downsampling factor',
                      help='Non-integer downsample factors are possible. '
                           'Must be larger than 1.')
        form.addParam('doDenoise', BooleanParam, default=False,
                      label='Denoising',
                      help="Apply a denoising method")
        form.addParam('maxIteration', IntParam, default=50, condition='doDenoise',
                      label='Max. number of iterations',
                      help='Max. number of iterations. Higher number = better '
                           'output but slower calculation. Must be larger '
                           'than 1.')
        form.addParam('doSmooth', BooleanParam, default=False,
                      label='Gaussian filter',
                      help="Apply a Gaussian filter in real space")
        form.addParam('sigmaConvolution', FloatParam, default=2,
                      condition="doSmooth",
                      label='Gaussian sigma (px)',
                      help="The larger this value, the more the effect will "
                           "be noticed")
        form.addParam('doHighPass', BooleanParam, default=False,
                      label='Highpass filter',
                      help="Apply a highpass filter in real space")
        form.addParam('highCutoff', FloatParam, default=0.002,
                      condition="doHighPass",
                      label='Cutoff frequency',
                      help="In normalized frequencies (<0.5). For example, "
                           "if you want to remove patterns larger than "
                           "500 pixels, use 1/500=0.002")
        form.addParam('highRaised', FloatParam, default=0.001,
                      condition="doHighPass", expertLevel=LEVEL_ADVANCED,
                      label='Transition bandwidth',
                      help="In normalized frequencies (<0.5). For example, "
                           "if you want to remove patterns larger than "
                           "1000 pixels, use 1/1000=0.001")
        form.addParam('doLowPass', BooleanParam, default=False,
                      label='Lowpass filter',
                      help="Apply a lowpass filter in real space")
        form.addParam('lowCutoff', FloatParam, default=0.4,
                      condition="doLowPass",
                      label='Cutoff frequency',
                      help="In normalized frequencies (<0.5). For example, "
                           "if you want to remove the crystalline ice at a frequency of 4A "
                           "and the pixel size is 0.5A, then the cutoff should be 0.5/4=0.125")
        form.addParam('lowRaised', FloatParam, default=0.001,
                      condition="doLowPass", expertLevel=LEVEL_ADVANCED,
                      label='Transition bandwidth',
                      help="In normalized frequencies (<0.5). The number of pixels in Fourier "
                           "will be approximately lowRaised*Xdim")
        form.addParam('doNormalize', BooleanParam, default=False,
                      label='Normalize micrograph?',
                      help='Normalize micrographs to be zero mean and '
                           'standard deviation one')
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
                       'sigmaConvolution': self.sigmaConvolution.get(),
                       'highCutoff': self.highCutoff.get(),
                       'highRaised': self.highRaised.get(),
                       'lowCutoff': self.lowCutoff.get(),
                       'lowRaised': self.lowRaised.get(),
                       'maxIterTV': self.maxIteration.get()}

    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):
        self._defineInputs()
        self.insertedDict = {}
        preprocessSteps = self._insertNewMicsSteps(self.insertedDict,
                                                   self.inputMicrographs.get())
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=preprocessSteps, wait=True)

    def createOutputStep(self):
        pass

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _insertNewMicsSteps(self, insertedDict, inputMics):
        deps = []
        for mic in inputMics:
            if mic.getObjId() not in insertedDict:
                fnOut = self._getOutputMicrograph(mic)
                stepId = self._insertStepsForMicrograph(mic.getFileName(), fnOut)
                deps.append(stepId)
                insertedDict[mic.getObjId()] = stepId
        return deps

    def _stepsCheck(self):
        # Input micrograph set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new micrographs to process from the input set
        micsFile = self.inputMicrographs.get().getFileName()
        micsSet = SetOfMicrographs(filename=micsFile)
        micsSet.loadAllProperties()
        self.SetOfMicrographs = [m.clone() for m in micsSet]
        self.streamClosed = micsSet.isStreamClosed()
        micsSet.close()
        newMics = any(m.getObjId() not in self.insertedDict
                      for m in self.inputMics)
        outputStep = self._getFirstJoinStep()
        if newMics:
            fDeps = self._insertNewMicsSteps(self.insertedDict,
                                             self.inputMics)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()


    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m.clone() for m in self.SetOfMicrographs
                   if int(m.getObjId()) not in doneList and self._isMicDone(m)]

        # We have finished when there is not more input micrographs (stream closed)
        # and the number of processed micrographs is equal to the number of inputs
        self.finished = self.streamClosed and (len(doneList) + len(newDone)) == len(self.SetOfMicrographs)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        outSet = self._loadOutputSet(SetOfMicrographs, 'micrographs.sqlite')

        for mic in newDone:
            micOut = em.data.Micrograph()
            if self.doDownsample:
                micOut.setSamplingRate(self.inputMicrographs.get().getSamplingRate() * self.downFactor.get())
            micOut.setObjId(mic.getObjId())
            micOut.setFileName(self._getOutputMicrograph(mic))
            micOut.setMicName(mic.getMicName())
            outSet.append(micOut)

        self._updateOutputSet('outputMicrographs', outSet, streamMode)
        
        if len(doneList)==0: #firstTime
            self._defineTransformRelation(self.inputMicrographs, outSet) 
            
        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)


    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputMicrographs.get()
        outputSet.copyInfo(inputs)
        if self.doDownsample:
            outputSet.setSamplingRate(self.inputMicrographs.get().getSamplingRate() * self.downFactor.get())
        return outputSet


    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
        outputSet.setStreamState(state)
        if self.hasAttribute(outputName):
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protcol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            # Here the defineOutputs function will call the write() method
            self._defineOutputs(**{outputName: outputSet})
            self._store(outputSet)
        # Close set databaset to avoid locking it
        outputSet.close()

    
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
        # Denoise
        self.__insertOneStep(self.doDenoise, "xmipp_transform_filter",
                            "-i %(inputMic)s --denoiseTV --maxIterTV %(maxIterTV)d")
        # Smooth
        self.__insertOneStep(self.doSmooth, "xmipp_transform_filter",
                            "-i %(inputMic)s --fourier real_gaussian %(sigmaConvolution)f")
        # Highpass
        self.__insertOneStep(self.doHighPass, "xmipp_transform_filter",
                            "-i %(inputMic)s --fourier high_pass %(highCutoff)f %(highRaised)f")
        # Lowpass
        self.__insertOneStep(self.doLowPass, "xmipp_transform_filter",
                            "-i %(inputMic)s --fourier low_pass %(lowCutoff)f %(lowRaised)f")
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

        
    #--------------------------- INFO functions --------------------------------
    
    def _validate(self):
        validateMsgs = []
        # Some prepocessing option need to be marked
        if not(self.doCrop or self.doDownsample or self.doLog or self.doRemoveBadPix or self.doInvert
               or self.doNormalize or self.doDenoise or self.doSmooth or self.doHighPass or self.doLowPass):
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
        if self.doDenoise:
            summary.append("Denoising applied with %d iterations" % self.maxIteration.get())
        if self.doSmooth:
            summary.append("Gaussian filtered with sigma=%f (px)"%self.sigmaConvolution.get())
        if self.doHighPass:
            summary.append("Highpass filtered with cutoff=%f (1/px)"%self.highCutoff.get())
        if self.doLowPass:
            summary.append("Lowpass filtered with cutoff=%f (1/px)"%self.lowCutoff.get())
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
        if self.doDenoise:
            txt += "been Denoised with %d iterations " % self.maxIteration.get()
        if self.doSmooth:
            txt += "been Gaussian filtered with a sigma of %0.2f pixels "%self.sigmaConvolution.get()
        if self.doHighPass:
            txt += "been highpass filtered with a cutoff of %f (1/px) "%self.highCutoff.get()
        if self.doLowPass:
            txt += "been lowpass filtered with a cutoff of %f (1/px) "%self.lowCutoff.get()
        if self.doNormalize:
            txt += "been normalized to mean 0 and variance 1"

        return [txt, "The resulting set of micrographs is %s" %
                self.getObjectTag('outputMicrographs')]

    #--------------------------- UTILS functions -------------------------------
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

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]
        return doneList

    def _getAllDone(self):
        return self._getExtraPath('DONE_all.TXT')

    def _writeDoneList(self, micList):
        """ Write to a text file the items that have been done. """
        with open(self._getAllDone(), 'a') as f:
            for mic in micList:
                f.write('%d\n' % mic.getObjId())

    def _isMicDone(self, mic):
        """ A movie is done if the marker file exists. """
        return os.path.exists(self._getMicDone(mic))

    def _getMicDone(self, mic):
        fn = mic.getFileName()
        extFn = getExt(fn)
        if extFn != ".mrc":
            fn = replaceExt(fn, "mrc")
        return self._getExtraPath('%s' % basename(fn))