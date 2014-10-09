# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from os.path import join

from pyworkflow.protocol.params import IntParam, BooleanParam, LEVEL_EXPERT
from pyworkflow.utils.path import moveFile
from pyworkflow.em.data import Micrograph
from pyworkflow.em.protocol import ProtProcessMovies
            
    
class ProtOpticalAlignment(ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'movie optical alignment'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('doGPU', BooleanParam, default=False,
                      label="Use GPU (vs CPU)", help="Set to true if you want the GPU implementation")
        form.addParam('GPUCore', IntParam, default=0,
                      label="Choose GPU core",
                      condition="doGPU",
                      help="GPU may have several cores. Set it to zero if you do not know what we are talking about")
        form.addParam('winSize', IntParam, default=150,
                      label="Window size", expertLevel=LEVEL_EXPERT,
                      help="Window size (shifts are assumed to be constant within this window).")
        line = form.addLine('Skip Frames:',
                      help='Drop first and last frames. set to 0 in order to keep all\n'
                           'First frame is 1\n')
        line.addParam('firstFrame', IntParam, default='0',
                      label='First')
        line.addParam('lastFrame',  IntParam, default='0',
                      label='Last')
        form.addParallelSection(threads=1, mpi=1)

    
    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):
        movSet = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(movSet)

        for movie in self.inputMovies.get():
            micName = self._getMicName(movie.getObjId())
            mic = Micrograph()
            # All micrograph are copied to the 'extra' folder after each step
            mic.setFileName(self._getExtraPath(micName))
            micSet.append(mic)

        self._defineOutputs(outputMicrographs=micSet)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    
    def _processMovie(self, movieId, movieName, movieFolder):
        """ Process the movie actions, remember to:
        1) Generate all output files inside movieFolder (usually with cwd in runJob)
        2) Copy the important result files after processing (movieFolder will be deleted!!!)
        """
        # Get the program (either with gpu or not)
        program = self._getProgram()
        # Prepare the arguments
        micName = self._getMicName(movieId)
        winSize = self.winSize.get()
        args = '-i %(movieName)s -o %(micName)s --winSize %(winSize)d' % locals()
        firstFrame = self.firstFrame.get()
        lastFrame = self.lastFrame.get()
        if firstFrame or lastFrame:
            args += ' --nst %d --ned %d' % (firstFrame, lastFrame)
        if self.doGPU:
            args += ' --gpu %d' % self.GPUCore.get()

        self.runJob(program, args, cwd=movieFolder)
        # Move output micrograph to 'extra' folder
        moveFile(join(movieFolder, micName), self._getExtraPath()) 
             
    def _getProgram(self):
        if self.doGPU:
            return 'xmipp_optical_alignment_gpu'
        else:
            return 'xmipp_optical_alignment_cpu'
    
  #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        numThreads = self.numberOfThreads;
        if numThreads>1:
            if self.doGPU:
                errors.append("GPU and Parallelization can not be used together")
        return errors

    def _citations(self):

        return ['Abrishami2014a']

    def _methods(self):
        """ METHODS TO DO"""
        pass

    def _summary(self):
        summary = []
        summary.append('Number of input movies: *%d*' % self.inputMovies.get().getSize())
        summary.append('Using a window size of: *%d*' % self.winSize.get())

        if self.doGPU:
            summary.append('- Used GPU for processing')

        return summary
