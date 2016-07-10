# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
Protocol wrapper around the motioncor2 tool for movie alignment
"""

import os, sys
from os.path import join, basename, exists

from pyworkflow.utils import makePath, moveFile, removeExt, getExt, Environ
from pyworkflow.protocol.params import StringParam, IntParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtProcessMovies


class ProtMotionCor(ProtProcessMovies):
    """
    Wrapper protocol for MotionCor2: Drift correction and Dose weighting program
    Written by Shawn Q. Zheng @ David Agard Lab, UCSF.

    If dose weighting is enabled, the program produces both weighted and
    unweighted (should be used for CTF estimation) aligned sums.
    """
    _label = 'align movies'
             
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        
        group = form.addGroup('Frame range')
        line2 = group.addLine('Used in alignment and final sum',
                             help='First and last frames used in alignment and total sum.\n'
                                  'The first frame in the stack is *0*.' )
        line2.addParam('frame0', IntParam, default=0, label='First')
        line2.addParam('frameN', IntParam, default=0, label='Last',
                       help='If *0*, use maximum value')

        form.addParam('patch', StringParam, default='5 5',
                      label='Number of patches',
                      help='Number of patches to be used for patch based alignment, '
                            'Set to *0 0* to do only global motion correction.')
        
        form.addParam('frameDose', FloatParam, default='0.0',
                      label='Frame dose (e/A^2)',
                      help='Frame dose in e/A^2. If set to 0.0, dose weighting will be skipped.')

        form.addParam('group', IntParam, default='1',
                      label='Group N frames',
                      help='Group every specified number of frames by adding them together.'
                           'The alignment is then performed on the summed frames.'
                           'By default, no grouping is performed.')

        form.addParam('tol', FloatParam, default='0.5',
                      label='Tolerance (px)',
                      help='Tolerance for iterative alignment, default 0.5px.')

        group2 = form.addGroup('Crop and binning')
        line = group2.addLine('Crop dimensions (px)',
                              expertLevel=LEVEL_ADVANCED,
                              help='How many pixels to crop from offset\n'
                                   'If equal to 0, use maximum size.')
        line.addParam('cropDimX', IntParam, default=0, label='X')
        line.addParam('cropDimY', IntParam, default=0, label='Y')

        group2.addParam('binFactor', FloatParam, default=1.0,
                        expertLevel=LEVEL_ADVANCED,
                        label='Binning factor',
                        help='Binning before processing, performed in Fourier space, default 1.0.')

        form.addParam('gpuId', StringParam, default='0',
                      expertLevel=LEVEL_ADVANCED,
                      label='GPU IDs',
                      help='GPU devices IDs. For multiple GPUs set to *0 1 2*.')

        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="""
        -Bft       100              BFactor for alignment, in px^2.
        -Iter      5                Maximum iterations for iterative alignment.
        -MaskCent  0 0              Center of subarea that will be used for alignment,
                                    default 0 0 corresponding to the frame center.
        -MaskSize  1.0 1.0          The size of subarea that will be used for alignment,
                                    default 1.0 1.0 corresponding full size.
        -Align     1                Generate aligned sum (1) or simple sum (0).
        -FmRef     0                Specify which frame to be the reference to which
                                    all other frames are aligned, default 0 is aligned to the first frame,
                                    other value aligns to the central frame.
        -Tilt      0 0              Tilt angle range for a dose fractionated tomographic tilt series stack,
                                    e.g. -60 60
                            """)

        form.addParallelSection(threads=0, mpi=0)
              
    def _processMovie(self, movieId, movieName, movieFolder):
        inputName = movieName
        rootName = movieName.replace('frames', '')
        micName = self._getNameExt(rootName, '_aligned_sum', 'mrc')
        logName = removeExt(rootName)
        gainFile = self.inputMovies.get().getGain()
        gpuId = self.gpuId.get()

        args = {'-Patch': self.patch.get(),
                '-Tol': self.tol.get(),
                '-FtBin': self.binFactor.get(),
                '-FmDose': self.frameDose.get(),
                '-PixSize': self.inputMovies.get().getSamplingRate(),
                '-kV': self.inputMovies.get().getAcquisition().getVoltage(),
                '-Throw': self.frame0.get(),
                '-Trunc': self.frameN.get(),
                '-Group': self.group.get(),
                '-Crop': '%d %d' % (self.cropDimX.get(), self.cropDimY.get()),
                '-Gpu': gpuId
                }
        
        if gainFile is not None:
            args['-Gain'] = gainFile

        if getExt(inputName) == 'tif' or getExt(inputName) == 'tiff':
            command = '-InTiff %(inputName)s -OutMrc %(micName)s -LogFile %(logName)s_ ' % locals()
        else:
            command = '-InMrc %(inputName)s -OutMrc %(micName)s -LogFile %(logName)s_ ' % locals()
        command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
        command += ' ' + self.extraParams.get()

        try:
            self.runJob(self._getProgram(), command, cwd=movieFolder, env=self.getEnviron())
        except:
            print >> sys.stderr, self._getProgram(), " failed for movie %(inputName)s" % locals()

        # Move the micrograph and alignment text files
        # before clean of movie folder
        makePath(self._getExtraPath('logs'))
        moveFile(join(movieFolder, micName), self._getExtraPath())
        if self.frameDose.get() != 0.0:
            wtSum = self._getNameExt(rootName, '_aligned_sum_DW', 'mrc')
            moveFile(join(movieFolder, wtSum), self._getExtraPath())

        if self.patch.get() == '0 0':
            logFile = self._getNameExt(rootName, '_0-Full', 'log')
            moveFile(join(movieFolder, logFile), self._getExtraPath('logs'))
        else:
            logFileNames = ['_0-Patch-Frame', '_0-Patch-Full', '_0-Patch-Patch']
            for fn in logFileNames:
                fullName = self._getNameExt(rootName, fn, 'log')
                moveFile(join(movieFolder, fullName), self._getExtraPath('logs'))

    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs(suffix='_unweighted')
        micSet.copyInfo(inputMovies)
        micSet.setObjLabel('unweighted aligned sums')

        if self.frameDose.get() == 0.0:
            for movie in inputMovies:
                movieId = movie.getObjId()
                rootName = movie.getFileName().replace('frames', '')
                micName = self._getNameExt(rootName, '_aligned_sum', 'mrc')

                mic = micSet.ITEM_TYPE()
                mic.setObjId(movieId)
                mic.setFileName(self._getExtraPath(micName))
                micSet.append(mic)

            self._defineOutputs(outputMicrographs=micSet)
            self._defineTransformRelation(inputMovies, micSet)
        else:
            micSet2 = self._createSetOfMicrographs(suffix='_weighted')
            micSet2.copyInfo(inputMovies)
            micSet2.setObjLabel('dose-weighted aligned sums')

            for movie in inputMovies:
                movieId = movie.getObjId()
                rootName = movie.getFileName().replace('frames', '')
                micName = self._getNameExt(rootName, '_aligned_sum', 'mrc')
                micName2 = self._getNameExt(rootName, '_aligned_sum_DW', 'mrc')

                mic = micSet.ITEM_TYPE()
                mic.setObjId(movieId)
                mic.setFileName(self._getExtraPath(micName))
                micSet.append(mic)

                mic2 = micSet2.ITEM_TYPE()
                mic2.setObjId(movieId)
                mic2.setFileName(self._getExtraPath(micName2))
                micSet2.append(mic2)

            self._defineOutputs(outputMicrographs=micSet, ouputMicrographs=micSet2)
            self._defineTransformRelation(inputMovies, micSet)
            self._defineTransformRelation(inputMovies, micSet2)
        
    #--------------------------- UTILS functions ---------------------------------------------------

    def _summary(self):
        summary = []
        firstFrame = self.frame0.get()
        lastFrame = self.frameN.get()
        if self.inputMovies.get():
            summary.append('Number of input movies: *%d*' % self.inputMovies.get().getSize())
        if lastFrame == 0:
            summary.append('Frames used in alignment: *%d* to *%s* (first frame is 0)' % (firstFrame, 'Last Frame'))
        else:
            summary.append('Frames used in alignment: *%d* to *%d* (first frame is 0)' % (firstFrame, lastFrame))

        return summary

    def _citations(self):
        return ['Zheng2016']

    def _methods(self):
        methods = []
        if self.patch.get() == '0 0':
            methods.append('Alignment method: global motion correction')
        else:
            methods.append('Alignment method: global + local motion correction for *%s* patches' % self.patch.get())
        if self.frameDose.get() != 0.0:
            methods.append('- Dose weighting enabled. Dose per frame: *%0.2f e/A^2*' % self.frameDose.get())
        if self.group.get() != 1:
            methods.append('- Group *%d* consecutive frames' % self.group.get())
        methods.append('- Tolerance of alignment accuracy: *%0.2f px*' % self.tol.get())

        return methods
    
    def _validate(self):
        errors = []
        # Check that the program exists
        if not exists(self._getProgram()):
            errors.append("Binary '%s' does not exits. \n"
                          "Check configuration file: ~/.config/scipion/scipion.conf\n"
                          "and set MOTIONCOR2_home variable properly." % self._getProgram())
            print "os.environ['MOTIONCOR2_HOME']", os.environ['MOTIONCOR2_HOME']
        return errors

    def _getProgram(self):
        """ Return the program binary that will be used. """
        binary = 'motioncor2'
        program = join(os.environ['MOTIONCOR2_HOME'], basename(binary))

        return program

    def getEnviron(self):
        """ Return the envirion settings to run motioncor2 program. """
        """ Setup the environment variables needed to launch Motioncor2. """
        environ = Environ(os.environ)
        environ.update({
            'LD_LIBRARY_PATH': join(os.environ.get('MOTIONCOR2_CUDA_LIB', ''))
            }, position=Environ.BEGIN)
        return environ
