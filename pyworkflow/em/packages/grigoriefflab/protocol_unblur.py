# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
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
This module contains the protocol for CTF estimation with Unblur3
"""

#import os
#import sys
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.protocol.params import  (IntParam, PointerParam,
                                         BooleanParam, FloatParam,
                                         LEVEL_ADVANCED)
from pyworkflow.utils.properties import Message



class ProtUnblur(ProtProcessMovies):
    """Unblur is used to align the frames of movies recorded
    on an electron microscope to reduce image blurring due
    to beam-induced motion. It reads stacks of movies that
    are stored in MRC/CCP4 format. Unblur generates frame
    sums that can be used in subsequent image processing
    steps and optionally applies an exposure-dependent
    filter to maximize the signal at all resolutions
    in the frame averages."""
    _label = 'Unblur'
    
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        form.addParam('alignFrameRange', IntParam,
                      default=-1,
                      label='Number of frames per movie ',
                      help='Number of frames per movie. -1 -> all frames ')
        form.addParam('doApplyDoseFilter', BooleanParam, default=True,
              label='Apply Dose filter',
              help='Apply Dose filter')
        form.addParam('exposurePerFrame', FloatParam,
                      label='Exposure per frame (e/A^2)',
                      help='Exposure per frame (e/A^2)')

        #group = form.addGroup('Expert Options')
        form.addParam('minShiftInitSearch', FloatParam,
                      default=2.,
                      label='Min. Shift Initial search (A)',
                      help='Minimum shift for initial search (Angstroms)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('OutRadShiftLimit', FloatParam,
                      default=200.,
                      label='Outer radius shift limit (A)',
                      help='Outer radius shift limit (Angstroms)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('HWVertFourMask', FloatParam,
                      default=1.,
                      label='Half-width vertical Fourier mask',
                      help='Half-width of central vertical line of Fourier mask',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('HWHoriFourMask', FloatParam,
                      default=1.,
                      label='Half-width horizontal Fourier mask',
                      help='Half-width of central horizontal line of Fourier mask',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('terminationShiftThreshold', FloatParam,
                      default=0.1,
                      label='Termination shift threshold',
                      help='Termination shift threshold',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('maximumNumberIterations', IntParam,
                      default=10,
                      label='Maximum number of iterations',
                      help='Maximum number of iterations',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('doRestoreNoisePower', BooleanParam,
                      default=True,
                      label='Restore Noise Power? ',
                      help='Restore Noise Power? ',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('doVerboseOutput', BooleanParam,
                      default=False,
                      label='Verbose Output?',
                      help='Verbose Output?',
                      expertLevel=LEVEL_ADVANCED)
        form.addParallelSection(threads=4, mpi=1)

#
#Input stack filename                [my_movie.mrc] : kk.mrc
#Number of frames per movie                    [34] :
#Output aligned sum file       [my_aligned_sum.mrc] :
#Output shifts file                 [my_shifts.txt] :
#Pixel size of images (A)                       [1] :
#Apply Dose filter?                            [NO] : YES
#Exposure per frame (e/A^2)                   [1.0] :
#Acceleration voltage (kV)                  [300.0] :
#Set Expert Options?                           [NO] : YES
#Output FRC file                       [my_frc.txt] :
#Minimum shift for initial search (Angstroms)
#[2.0]                                              :
#Outer radius shift limit (Angstroms)       [200.0] :
#B-factor to apply to images (A^2)           [1500] :
#Half-width of central vertical line of Fourier mask
#[1]                                                :
#Half-width of central horizontal line of Fourier mask
#[1]                                                :
#Termination shift threshold                  [0.1] :
#Maximum number of iterations                  [10] :
#Restore Noise Power?                         [YES] :
#Verbose Output?                               [NO] : YES

    #--------------------------- STEPS functions ---------------------------------------------------
    def _processMovie(self, movieId, movieName, movieFolder):
        HERHEREHERE
        """call program here"""
        inputName = movieName
        micName = self._getMicName(movieId)
        logFile = self._getLogFile(movieId)
        gainFile = self.inputMovies.get().getGain()

    def createOutputStep(self):
        """save whatever"""
        pass

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []