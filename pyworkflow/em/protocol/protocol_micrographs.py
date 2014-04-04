# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
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
In this module are protocol base classes related to EM Micrographs

"""
from pyworkflow.em.protocol import *

class ProtMicrographs(EMProtocol):
    pass

class ProtCTFMicrographs(ProtMicrographs):
    """ Base class for all protocols that estimates the CTF"""
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.methodsInfo = String()
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        
        form.addParam('inputMicrographs', PointerParam, important=True,
                      label=Message.LABEL_INPUT_MIC, pointerClass='SetOfMicrographs')
#        form.addParam('ampContrast', FloatParam, default=0.1,
#                      label='Amplitude Contrast',
#                      help='It should be a positive number, typically between 0.05 and 0.3.')
        form.addParam('lowRes', FloatParam, default=0.05,
                      label=Message.LABEL_LOW_RES,
                      help=Message.TEXT_LOW_RES)
        form.addParam('highRes', FloatParam, default=0.35,
                      label=Message.LABEL_HIGH_RES, 
                      help=Message.TEXT_HIGH_RES)
        form.addParam('minDefocus', FloatParam, default=0.5,
                      label=Message.LABEL_MIN_FOCUS,
                      help=Message.TEXT_MIN_FOCUS,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('maxDefocus', FloatParam, default=10.,
                      label=Message.LABEL_MAX_FOCUS,
                      help=Message.TEXT_MAX_FOCUS,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('windowSize', IntParam, default=256,
                      label=Message.LABEL_WINDOW_SIZE,
                      help=Message.TEXT_WINDOW_SIZE,
                      expertLevel=LEVEL_ADVANCED)
        
        form.addParallelSection(threads=2, mpi=0)       
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform ctf estimation on a set of micrographs.
        """
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        acquisition = self.inputMics.getAcquisition()
         
        self._params = {'voltage': acquisition.getVoltage(),
                        'sphericalAberration': acquisition.getSphericalAberration(),
                        'magnification': acquisition.getMagnification(),
                        'ampContrast': acquisition.getAmplitudeContrast(),
                        'samplingRate': self.inputMics.getSamplingRate(),
                        'scannedPixelSize': self.inputMics.getScannedPixelSize(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        # Convert from microns to Amstrongs
                        'minDefocus': self.minDefocus.get() * 1e+4, 
                        'maxDefocus': self.maxDefocus.get() * 1e+4
                       }
        
        self._prepareCommand()
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each micrograph insert the steps to process it
        for micFn, micDir, _ in self._iterMicrographs():
            # CTF estimation with Xmipp
            stepId = self._insertFunctionStep('_estimateCTF', micFn, micDir,
                                              prerequisites=[]) # Make estimation steps indepent between them
            deps.append(stepId)
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(Message.ERROR_NO_EST_CTF)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCTF'):
            summary.append(Message.TEXT_NO_CTF_READY)
        else:
            summary.append("CTF estimation of %d micrographs." % self.inputMicrographs.get().getSize())
        return summary
    
    def _methods(self):
        methods = []
        
        if not hasattr(self, 'outputCTF'):
            methods.append(Message.TEXT_NO_CTF_READY)
        else:
            methods.append(self.methodsInfo.get())
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getExtraPath(removeBaseExt(mic.getFileName()))        
        
    def _iterMicrographs(self):
        """ Iterate over micrographs and yield
        micrograph name and a directory to process.
        """
        for mic in self.inputMics:
            micFn = mic.getFileName()
            micDir = self._getExtraPath(removeBaseExt(micFn)) 
            yield (micFn, micDir, mic)  
                
    def _prepareCommand(self):
        """ This function should be implemented to prepare the
        arguments template if doesn't change for each micrograph
        After this method self._program and self._args should be set. 
        """
        pass
    
    def _defocusMaxMin(self, list):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        minimum = float(min(list))/10000
        maximum = float(max(list))/10000
        msg = "The range of micrograph's experimental defocus are %(minimum)0.3f - %(maximum)0.3f microns" % locals()
        self.methodsInfo.set(msg)


class ProtPreprocessMicrographs(ProtMicrographs):
    pass


class ProtProcessMovies(ProtPreprocessMicrographs):
    """Protocol base for protocols to process movies from direct detectors cameras"""
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputMovies', PointerParam, important=True,
                      label=Message.LABEL_INPUT_MOVS, pointerClass='SetOfMovies')
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        movSet = self.inputMovies.get()
        self._micList = []
        for mov in movSet:
            movFn = mov.getFirstItem().getFileName()
            self._insertFunctionStep('processMoviesStep', movFn)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def processMoviesStep(self, movFn):
        movName = removeBaseExt(movFn)
        
        self._createMovWorkingDir(movName)
        movDir = self._movWorkingDir(movName)
        
        self._enterDir(movDir)
        self._defineProgram()
        movRelFn = os.path.relpath(movFn, movDir)
        args = "%s" % movRelFn
        self.runJob(self._program, args)
        self._leaveDir()
        
        micJob = join(movDir, "justtest.mrc")
        micFn = self._getExtraPath(movName + ".mrc")
        moveFile(micJob, micFn)
        self._micList.append(micFn)
    
    def createOutputStep(self):
        micSet = self._createSetOfMicrographs()
        movSet = self.inputMovies.get()
        micSet.setAcquisition(movSet.getAcquisition())
        micSet.setSamplingRate(movSet.getSamplingRate())
        
        for m in self._micList:
            mic = Micrograph()
            mic.setFileName(m)
            micSet.append(mic)
        self._defineOutputs(outputMicrographs=micSet)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _createMovWorkingDir(self, movFn):
        """create a new directory for the movie and change to this directory.
        """
        workDir = self._movWorkingDir(movFn)
        makePath(workDir)   # Create a directory for a current iteration
    
    def _movWorkingDir(self, movFn):
        """ Define which is the directory for the current movie"""
        movDir = '%s' % movFn
        workDir = self._getTmpPath(movFn)
        return workDir
    
    
    
class ProtOpticalAlignment(ProtProcessMovies):
    """ Protocol to align movies, from direct detectors cameras, into micrographs.
    """
    _label = 'optical alignment'
    
    def _defineProgram(self):
        XMP_OPT_ALIGN = 'xmipp_optical_alignment'
        self._program = join(os.environ['OPT_ALIGN_HOME'], XMP_OPT_ALIGN)