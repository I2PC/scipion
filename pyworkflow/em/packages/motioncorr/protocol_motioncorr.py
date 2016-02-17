# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
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
Protocol wrapper around the MotionCorr for movie alignment
"""

import os, sys

import pyworkflow.em as em
import pyworkflow.utils.path as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
from pyworkflow.em.protocol import ProtAlignMovies

from convert import parseMovieAlignment


class ProtMotionCorr(ProtAlignMovies):
    """
    Wrapper protocol to Dose Fractionation Tool: Flat fielding and Drift correction
    Wrote by Xueming Li @ Yifan Cheng Lab, UCSF   
    """

    _label = 'motioncorr alignment'
    CONVERT_TO_MRC = 'mrc'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineAlignmentParams(self, form):
        form.addParam('gpuMsg', params.LabelParam, default=True,
                      label='WARNING! You need to have installed CUDA'
                            ' libraries and a nvidia GPU')

        ProtAlignMovies._defineAlignmentParams(self, form)

        group = form.addGroup('GPU', expertLevel=cons.LEVEL_ADVANCED)
        group.addParam('GPUCore', params.IntParam, default=0,
                       label="Choose GPU core",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on.")
        group.addParam('extraParams', params.StringParam, default='',
                       label='Additional parameters',
                       help="""
-bft       150               BFactor in pix^2.
-pbx       96                Box dimension for searching CC peak.
-fod       2                 Number of frame offset for frame comparison.
-nps       0                 Radius of noise peak.
-sub       0                 1: Save as sub-area corrected sum. 0: Not.
-srs       0                 1: Save uncorrected sum. 0: Not.
-scc       0                 1: Save CC Map. 0: Not.
-slg       1                 1: Save Log. 0: Not.
-atm       1                 1: Align to middle frame. 0: Not.
-dsp       1                 1: Save quick results. 0: Not.
-fsc       0                 1: Calculate and log FSC. 0: Not.
                            """)

        form.addParallelSection(threads=0, mpi=0)

    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movie):
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outputMicFn = self._getOutputMicName(movie)
        outputMovieFn = self._getOutputMovieName(movie)

        # Get the number of frames and the range to be used for alignment and sum
        numberOfFrames = movie.getNumberOfFrames()
        a0, aN = self._getFrameRange(numberOfFrames, 'align')
        s0, sN = self._getFrameRange(numberOfFrames, 'sum')

        logFile = self._getLogFile(movie)
        args = {'-crx': self.cropOffsetX.get(),
                '-cry': self.cropOffsetY.get(),
                '-cdx': self.cropDimX.get(),
                '-cdy': self.cropDimY.get(),
                '-bin': self.binFactor.get(),
                '-nst': a0,
                '-ned': aN,
                '-nss': s0,
                '-nes': sN,
                '-gpu': self.GPUCore.get(),
                '-flg': logFile,
                }

        # FIXME: Always produce the average micrograph?
        command = '%s -fcs %s ' % (movie.getBaseName(), outputMicFn)
        command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])

        if inputMovies.getGain():
            command += " -fgr " + inputMovies.getGain()

        if inputMovies.getDark():
            command += " -fdr " + inputMovies.getDark()

        if self.doSaveMovie:
            command += " -fct %s -ssc 1" % outputMovieFn

        command += ' ' + self.extraParams.get()
        program = 'dosefgpu_driftcorr'

        try:
            #self.info("Running: %s %s" % (program, command))
            self.runJob(program, command, cwd=movieFolder)
        except:
            print >> sys.stderr, program, " failed for movie %(movieName)s" % locals()

        # FIXME: Why the following lines?
        #putils.cleanPattern(os.path.join(movieFolder, movieName))
        #putils.moveTree(self._getTmpPath(), self._getExtraPath())

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        bin = int(self.binFactor.get())

        if not (bin == 1 or bin == 2):
            errors.append("Binning factor can only be 1 or 2")

        return errors

    #--------------------------- UTILS functions ---------------------------------------------------
    def _getMovieLogFile(self, movie):
        return 'micrograph_%06d_Log.txt' % movie.getObjId()

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
         The shifts should refer to the original micrograph without any binning.
         In case of a bining greater than 1, the shifts should be scaled.
        """
        return parseMovieAlignment(self._getMovieLogFile(movie))



