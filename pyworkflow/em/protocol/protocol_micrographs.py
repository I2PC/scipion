# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

import os
from os.path import join, basename, exists, dirname, relpath

from pyworkflow.object import String
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, LEVEL_EXPERT
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, TextParam, BooleanParam, FileParam
from pyworkflow.utils.path import copyTree, copyFile, removeBaseExt, makePath, moveFile
from pyworkflow.utils.properties import Message
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Micrograph, SetOfImages, SetOfCTF



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

        line = form.addLine('Resolution', 
                            help='Give a value in digital frequency (i.e. between 0.0 and 0.5). '
                                 'These cut-offs prevent the typical peak at the center of the PSD and high-resolution'
                                 'terms where only noise exists, to interfere with CTF estimation. The default lowest '
                                 'value is 0.05 but for micrographs with a very fine sampling this may be lowered towards 0.'
                                 'The default highest value is 0.35, but it should '+'be increased for micrographs with '
                                 'signals extending beyond this value. However, if your micrographs extend further than '
                                 '0.35, you should consider sampling them at a finer rate.')
        line.addParam('lowRes', FloatParam, default=0.05,
                      label='Lowest' )
        line.addParam('highRes', FloatParam, default=0.35,
                      label='Highest')
        
        line = form.addLine('Defocus search range (microns)', expertLevel=LEVEL_ADVANCED,
                            help='Select _minimum_ and _maximum_ values for defocus search range (in microns).'
                                 'Underfocus is represented by a positive number.')
        line.addParam('minDefocus', FloatParam, default=0.5, 
                      label='Min')
        line.addParam('maxDefocus', FloatParam, default=10.,
                      label='Max')
        
        form.addParam('windowSize', IntParam, default=256, expertLevel=LEVEL_ADVANCED,
                      label='Window size', 
                      help='The PSD is estimated from small patches of this size. Bigger patches '
                           'allow identifying more details. However, since there are fewer windows, '
                           'estimations are noisier.')
        
        form.addParallelSection(threads=2, mpi=1)       
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform ctf estimation on a set of micrographs.
        """
        
        self._defineValues()
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
    def _defineValues(self):
        """ This function get some parameters of the micrographs"""
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
    
    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getExtraPath(removeBaseExt(mic.getFileName()))        
    
    def _iterMicrographs(self):
        """ Iterate over micrographs and yield
        micrograph name and a directory to process.
        """
        for mic in self.inputMics:
            micFn = mic.getFileName()
            micDir = self._getMicrographDir(mic) 
            yield (micFn, micDir, mic)  
    
    def _prepareCommand(self):
        """ This function should be implemented to prepare the
        arguments template if doesn't change for each micrograph
        After this method self._program and self._args should be set. 
        """
        pass
    
    def _defocusMaxMin(self, defocusList):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        minimum = float(min(defocusList))/10000
        maximum = float(max(defocusList))/10000
        msg = "The range of micrograph's experimental defocus are %(minimum)0.3f - %(maximum)0.3f microns" % locals()
        self.methodsInfo.set(msg)


class ProtRecalculateCTF(ProtMicrographs):
    """ Base class for all protocols that re-calculate the CTF"""
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.methodsInfo = String()
        self.summaryInfo = String()
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        
        form.addParam('inputCtf', PointerParam, important=True,
              label="input the SetOfCTF to recalculate", pointerClass='SetOfCTF')
        form.addParam('inputValues', TextParam)
        form.addHidden('sqliteFile', FileParam)
        
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform ctf re-estimation on a set of CTFs.
        """
        self._insertFunctionStep('createSubsetOfCTF')
        inputValsId = self._insertFunctionStep('copyInputValues')
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each psd insert the steps to process it
        self.values = self._splitFile(self.inputValues.get())
        for line in self.values:
            # CTF Re-estimation with Xmipp
            copyId = self._insertFunctionStep('copyFiles', line, prerequisites=[inputValsId])
            stepId = self._insertFunctionStep('_estimateCTF',line, prerequisites=[copyId]) # Make estimation steps independent between them
            deps.append(stepId)
        # Insert step to create output objects
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def createSubsetOfCTF(self):
        """ Create a subset of CTF and Micrographs analyzing the CTFs. """
        
        self._loadDbNamePrefix() # load self._dbName and self._dbPrefix
        
        self.setOfCtf = self._createSetOfCTF("_subset")
        modifiedSet = SetOfCTF(filename=self._dbName, prefix=self._dbPrefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                mic = ctf.getMicrograph()
                self.setOfCtf.append(ctf)
#         self.setOfCtf.write()
           
    def copyInputValues(self):
        """ Copy a parameter file that contain the info of the
        micrographs to recalculate its CTF to a current directory"""
        srcFile = self.inputValues.get()
        baseFn = basename(srcFile)
        dstFile = self._getTmpPath(baseFn)
        copyFile(srcFile, dstFile)
    
    def copyFiles(self, line):
        """Copy micrograph's directory tree"""
        objId = self._getObjId(line)
        ctfModel = self.setOfCtf.__getitem__(objId)
        mic = ctfModel.getMicrograph()
        
        prevDir = self._getPrevMicDir(mic)
        micDir = self._getMicrographDir(mic)
        # Create micrograph dir under extra directory
        print "creating path micDir=", micDir
        makePath(micDir)
        if not exists(micDir):
            raise Exception("No created dir: %s " % micDir)
        copyTree(prevDir, micDir)
    
    def _estimateCTF(self, line):
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
        if not (hasattr(self, 'outputCTF') and hasattr(self, 'outputMicrographs')):
            summary.append(Message.TEXT_NO_CTF_READY)
        else:
            summary.append(self.summaryInfo.get())
        return summary
    
    def _methods(self):
        methods = []
        
        if not (hasattr(self, 'outputCTF') and hasattr(self, 'outputMicrographs')):
            methods.append(Message.TEXT_NO_CTF_READY)
        else:
            methods.append(self.methodsInfo.get())
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _defineValues(self, line):
        """ This function get the acquisition info of the micrographs"""
        
        objId = self._getObjId(line)
        ctfModel = self.setOfCtf.__getitem__(objId)
        mic = ctfModel.getMicrograph()
        
        acquisition = mic.getAcquisition()
        scannedPixelSize = mic.getSamplingRate() * acquisition.getMagnification() / 10000
        self._params = {'voltage': acquisition.getVoltage(),
                        'sphericalAberration': acquisition.getSphericalAberration(),
                        'magnification': acquisition.getMagnification(),
                        'ampContrast': acquisition.getAmplitudeContrast(),
                        'scannedPixelSize': scannedPixelSize,
                        'samplingRate': mic.getSamplingRate()
                       }
    
    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getExtraPath(removeBaseExt(mic.getFileName()))        
    
    def _prepareCommand(self):
        """ This function should be implemented to prepare the
        arguments template if doesn't change for each micrograph
        After this method self._program and self._args should be set. 
        """
        pass
    
    def _getPrevMicDir(self, mic):
        
        objFn = self.inputCtf.get().getFileName()
        directory = dirname(objFn)
        return join(directory, "extra", removeBaseExt(mic.getFileName()))
    
    def _getObjId(self, values):
        return int(values[0])
    
    def _splitFile(self, filename):
        """ This method split the parameter filename into lines"""
        values = []
        f1 = open(filename)
        for l in f1:
            split = l.split()
            values.append(split)
        f1.close()
        return values
    
    def _defocusMaxMin(self, values):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        minimum = float(min(values))/10000
        maximum = float(max(values))/10000
        msg = "The range of micrograph's experimental defocus are %(minimum)0.3f - %(maximum)0.3f microns" % locals()
        self.methodsInfo.set(msg)
    
    def _ctfCounter(self, values):
        """ This function return the number of CTFs that was recalculated.
        """
        numberOfCTF = len(values)/2
        msg = "CTF Re-estimation of %(numberOfCTF)d micrographs" % locals()
        self.summaryInfo.set(msg)
    
    def _loadDbNamePrefix(self):
        """ Setup filename and prefix for db connection. """
        self._dbName = self.sqliteFile.get()
        self._dbPrefix = ""
        if self._dbPrefix.endswith('_'):
            self._dbPrefix = self._dbPrefix[:-1] 

class ProtPreprocessMicrographs(ProtMicrographs):
    pass


class ProtProcessMovies(ProtPreprocessMicrographs):
    """Protocol base for protocols to process movies from direct detectors cameras"""
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputMovies', PointerParam, important=True,
                      label=Message.LABEL_INPUT_MOVS, pointerClass='SetOfMovies')
        form.addParam('doGPU', BooleanParam, default=False,
                      label="Use GPU (vs CPU)", help="Set to true if you want the GPU implementation")
        form.addParam('GPUCore', IntParam, default=0,
                      label="Choose GPU core",
                      condition="doGPU",
                      help="GPU may have several cores. Set it to zero if you do not know what we are talking about")
        form.addParam('winSize', IntParam, default=150,
                      label="Window size", expertLevel=LEVEL_EXPERT,
                      help="Window size (shifts are assumed to be constant within this window).")
        line = form.addLine('Drop Frames (NOT IMPLEMENTED):',
                      help='Drop first and last frames. set to 0 in orser to keep all\n\n'
                           'NOT IMPLEMENTED YET.')
        line.addParam('firstFrames', IntParam, default='0',
                      label='First')
        line.addParam('lastFrames',  IntParam, default='0',
                      label='Last')
        form.addParallelSection(threads=1, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        movSet = self.inputMovies.get()
        for mov in movSet:
            #TODO What milks is this?
            mov.load()
            movFn = mov.getFirstItem().getFileName()
            self._insertFunctionStep('processMoviesStep', movFn)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def processMoviesStep(self, movFn):
        movName = removeBaseExt(movFn)
        self._createMovWorkingDir(movName)
        self._defineProgram()
        self._micList = []
        micFn = self._getExtraPath( movName+ "_aligned.spi")
        args = '-i %s -o %s --winSize %d'%(movFn, micFn, self.winSize.get())
        if False:
            args += ' --dropFirst %d --dropLast %d' % (self.firstFrames.get(),self.lastFrames.get())
        if self.doGPU:
            args += ' --gpu %d' % self.GPUCore.get()
        self.runJob(self._program, args)

        # micJob = join(movDir, "justtest.mrc")
        # micFn = self._getExtraPath(movName + ".mrc")
        # moveFile(micJob, micFn)
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
        workDir = self._getTmpPath(movFn)
        return workDir
    
    
class ProtOpticalAlignment(ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'movies optical alignment'
    
    def _defineProgram(self):
        if self.doGPU:
            self._program = 'xmipp_optical_alignment_gpu'
        else:
            self._program = 'xmipp_optical_alignment_cpu'
