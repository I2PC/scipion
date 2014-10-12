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

from pyworkflow.protocol.params import IntParam, StringParam, BooleanParam, LEVEL_EXPERT, LEVEL_ADVANCED, EnumParam
from pyworkflow.utils.path import moveFile
from pyworkflow.em.data import Micrograph
from pyworkflow.em.protocol import ProtProcessMovies

# Alignment methods enum
AL_OPTICAL = 0
AL_DOSEFGPU = 1
AL_DOSEFGPUOPTICAL = 2
AL_AVERAGE = 3

class ProtOpticalAlignment(ProtProcessMovies):
    """ Aligns movies, from direct detectors cameras, into micrographs.
    """
    _label = 'movie optical alignment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        form.addParam('alignMethod', EnumParam, choices=['optical flow', 'dosefgpu',
                                                         'dosefgpu + optical flow', 'average'],
                      label="Alignment method", default=AL_OPTICAL,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Method to use for alignment of the movies')
        group = form.addGroup('Common parameters')
        group.addParam('GPUCore', IntParam, default=0,
                      label="Choose GPU core",
                      help="GPU may have several cores. Set it to zero if you do not know what we are talking about")
        line = group.addLine('Used in alignment',
                            help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('alignFrame0', IntParam, default=0, label='Fisrt')
        line.addParam('alignFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        group = form.addGroup('Optical Flow parameters',condition="alignMethod==%d or alignMethod==%d " % (AL_OPTICAL, AL_DOSEFGPUOPTICAL))
        group.addParam('doGPU', BooleanParam, default=False,
                      label="Use GPU (vs CPU)", help="Set to true if you want the GPU implementation")
        group.addParam('winSize', IntParam, default=150,
                      label="Window size", expertLevel=LEVEL_EXPERT,
                      help="Window size (shifts are assumed to be constant within this window).")
        #---------------------------------- DosefGPU Params--------------------------------
        group = form.addGroup('DosefGPU parameters',condition="alignMethod==%d or alignMethod==%d " % (AL_DOSEFGPU, AL_DOSEFGPUOPTICAL))
        line = group.addLine('Used in final sum',
                             help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line.addParam('sumFrame0', IntParam, default=0, label='First')
        line.addParam('sumFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', IntParam, default=0, label='X')
        line.addParam('cropOffsetY', IntParam, default=0, label='Y')
        line = group.addLine('Crop dimensions (px)',
                      help='How many pixels to crop from offset\n'
                           'If equal to 0, use maximum size.')
        line.addParam('cropDimX', IntParam, default=0, label='X')
        line.addParam('cropDimY', IntParam, default=0, label='Y')
        group.addParam('binFactor', IntParam, default=1,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')
        group.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""
-bft       150               BFactor in pix^2.
-pbx       96                Box dimension for searching CC peak.
-fod       2                 Number of frame offset for frame comparision.
-nps       0                 Radius of noise peak.
-sub       0                 1: Save as sub-area corrected sum. 0: Not.
-srs       0                 1: Save uncorrected sum. 0: Not.
-ssc       0                 1: Save aligned stack. 0: Not.
-scc       0                 1: Save CC Map. 0: Not.
-slg       1                 1: Save Log. 0: Not.
-atm       1                 1: Align to middle frame. 0: Not.
-dsp       1                 1: Save quick results. 0: Not.
-fsc       0                 1: Calculate and log FSC. 0: Not.
                      """)
        form.addParallelSection(threads=1, mpi=1)


    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        alMethod = self.alignMethod.get()
        #if alMethod == AL_DOSEFGPU:
            # Also create a Set of Movies with the alignment parameters
            #movieSet = self._createSetOfMovies()
            #movieSet.copyInfo(inputMovies)

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
        # Read common parameters
        firstFrame = self.alignFrame0.get()
        lastFrame = self.alignFrameN.get()
        gpuId = self.GPUCore.get()
        alMethod = self.alignMethod.get()
        # For simple average execution
        if alMethod == AL_AVERAGE:
            args = '-i %(movieName)s -o %(micName)s' % locals()
            args += ' --nst %d --ned %d --simpleAverage' % (firstFrame, lastFrame)
            self.runJob(program, args, cwd=movieFolder)
        # For Optical Flow execution
        if alMethod == AL_OPTICAL:
            winSize = self.winSize.get()
            args = '-i %(movieName)s -o %(micName)s --winSize %(winSize)d' % locals()
            args += ' --nst %d --ned %d' % (firstFrame, lastFrame)
            if self.doGPU:
                args += ' --gpu %d' % gpuId
            self.runJob(program, args, cwd=movieFolder)
        # For DosefGPU Execution
        if alMethod == AL_DOSEFGPU:
            #logFile = self._getLogFile(movieId)
            #gainFile = self.inputMovies.get().getGain()
            args = {'-crx': self.cropOffsetX.get(),
                    '-cry': self.cropOffsetY.get(),
                    '-cdx': self.cropDimX.get(),
                    '-cdy': self.cropDimY.get(),
                    '-bin': self.binFactor.get(),
                    '-nst': self.alignFrame0.get(),
                    '-ned': self.alignFrameN.get(),
                    '-nss': self.sumFrame0.get(),
                    '-nes': self.sumFrameN.get(),
                    '-gpu': gpuId,
                    #'-flg': logFile,
                    }
            command = '%(movieName)s -fcs %(micName)s ' % locals()
            command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
            command += ' ' + self.extraParams.get()
            import pyworkflow.em.packages.dosefgpu as dosefgpu
            self.runJob(program, command, cwd=movieFolder,
                        env=dosefgpu.getEnviron())


        # Move output micrograph to 'extra' folder
        moveFile(join(movieFolder, micName), self._getExtraPath())

    def _getProgram(self):
        alMethod = self.alignMethod.get()
        if alMethod == AL_AVERAGE:
            return 'xmipp_optical_alignment_cpu'
        if alMethod == AL_OPTICAL:
            if self.doGPU:
                return 'xmipp_optical_alignment_gpu'
            else:
                return 'xmipp_optical_alignment_cpu'
        elif alMethod == AL_DOSEFGPU:
            return 'dosefgpu_driftcorr'

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
