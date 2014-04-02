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
In this module are protocol base classes related to 2D processing

"""
from pyworkflow.em.protocol import *

class ProtPreprocessVolumes(EMProtocol):
    """ This class will serve as a base for all protocol
    that performs some operation on Volumes (i.e. filters, mask, resize, etc)
    It is mainly defined by an inputVolumes and outputVolumes.
    """
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputVolumes', PointerParam, important=True,
                      label=Message.LABEL_INPUT_VOLS, pointerClass='Volume, SetOfVolumes',
                      help='Can be a density volume or a SetOfVolumes.')
        # Hook that should be implemented in subclasses
        self._defineProcessParams(form)
        form.addParallelSection(threads=2, mpi=1)
        
    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass


class ProtFilterVolumes(ProtPreprocessVolumes):
    """ This is the base for the branch of filters, 
    between the ProtPreprocessVolumes """
    pass

class ProtMaskVolumes(ProtPreprocessVolumes):
    """ This is the base for the branch of mask, 
    between the ProtPreprocessVolumes """
    pass

#class ProtInitialVolume(EMProtocol):
#    pass

class ProtRefine3D(EMProtocol):
    pass

class ProtClassify3D(EMProtocol):
    pass

class ProtValidate3D(EMProtocol):
    pass

class ProtCreateMask3D(EMProtocol):
    pass

class ProtProcessMovies(EMProtocol):
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


class ProtInitialVolume(EMProtocol):
    """Protocol base for Initial volumes protocols"""
    pass

class ProtAlignVolume(EMProtocol):
    """Protocol base for Align volumes protocols"""
    pass