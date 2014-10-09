# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
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
Protocol wrapper around the ResMap tool for local resolution
"""

from os.path import join, exists, basename
from glob import glob

from pyworkflow.utils import makePath, moveFile, copyFile
from pyworkflow.protocol.params import StringParam, IntParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtProcessMovies

from convert import parseMovieAlignment



class ProtDosefGpu(ProtProcessMovies):
    """
    Wrapper protocol to 
    
    Dose Fractionation Tool: Flat fielding and Drift correction
    Wrote by Xueming Li @ Yifan Cheng Lab, UCSF   
    """
    _label = 'align movies'
             
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        
        group = form.addGroup('Frame range')
        line1 = group.addLine('Used in alignment',
                             help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line1.addParam('alignFrame0', IntParam, default=0, label='Fisrt')
        line1.addParam('alignFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')
        
        line2 = group.addLine('Used in final sum',
                             help='First and last frames used in alignment.\n'
                                  'The first frame in the stack is *0*.' )
        line2.addParam('sumFrame0', IntParam, default=0, label='First')
        line2.addParam('sumFrameN', IntParam, default=0, label='Last',
                      help='If *0*, use maximum value')        
        
        form.addParam('gpuId', StringParam, default='0',
                      expertLevel=LEVEL_ADVANCED,
                      label='GPU id',
                      help='GPU device ID')
        
        form.addParam('extraParams', StringParam, default='',
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
        
        form.addSection('Crop and binning')

        line = form.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', IntParam, default=0, label='X')
        line.addParam('cropOffsetY', IntParam, default=0, label='Y')
    
    line = form.addLine('Crop dimensions (px)',
                      help='How many pixels to crop from offset\n'
                           'If equal to 0, use maximum size.')
        line.addParam('cropDimX', IntParam, default=0, label='X')
        line.addParam('cropDimY', IntParam, default=0, label='Y')


        form.addParam('binFactor', IntParam, default=1, 
                      label='Binning factor',
                      help='1x or 2x. Bin stack before processing.')
              
    form.addParallelSection(threads=1, mpi=1)
                     
    def _processMovie(self, movieId, movieName, movieFolder):
        inputName = movieName
        micName = self._getMicName(movieId)
        logFile = self._getLogFile(movieId)
        gainFile = self.inputMovies.get().getGain()
        gpuId = self.gpuId.get()

# TODO Check better way to handle gain correction
#         if gainFile is not None:
#             # Apply the gain correction to flat the raw movie
#             correctedName = movieName.replace('.mrc', '_corrected.mrc')
#             
#             self.runJob('dosefgpu_flat', 
#                         '%(inputName)s %(correctedName)s %(gainFile)s %(gpuId)s' % locals(),
#                         cwd=movieFolder)
#            
#            inputName = correctedName
        
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
                '-flg': logFile,
                }
        
        #TODO: check the gain can be handle in dosefgpu_driftcoor program
        #if gainFile is not None:
        #    args['-fgr'] = gainFile
        
        command = '%(inputName)s -fcs %(micName)s ' % locals()
        command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
        command += ' ' + self.extraParams.get()

        self.runJob('dosefgpu_driftcorr', command, cwd=movieFolder)
        # Move the micrograph and alignment text file
        # before clean of movie folder
        moveFile(join(movieFolder, micName), self._getExtraPath())
        moveFile(join(movieFolder, logFile), self._getExtraPath())        
        
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        # Also create a Set of Movies with the alignment parameters
        movieSet = self._createSetOfMovies()
        movieSet.copyInfo(inputMovies)
        
        for movie in inputMovies:
            movieId = movie.getObjId()
            micName = self._getMicName(movieId)
            movieFolder = self._getMovieFolder(movieId)
          
            mic = micSet.ITEM_TYPE()
            mic.setObjId(movieId)
            mic.setFileName(self._getExtraPath(micName))
            micSet.append(mic)
            
            # Parse the alignment parameters and store the log files
            alignedMovie = movie.clone()
            logFile = self._getExtraPath(self._getLogFile(movieId))
            alignment = parseMovieAlignment(logFile)
            alignedMovie.setAlignment(alignment)
            movieSet.append(alignedMovie)
            
        self._defineOutputs(outputMicrographs=micSet)
        self._defineTransformRelation(inputMovies, micSet)
        
        self._defineOutputs(outputMovies=movieSet)
        self._defineTransformRelation(inputMovies, movieSet)
        
    #--------------------------- UTILS functions ---------------------------------------------------
        
    def _getLogFile(self, movieId):
        return 'micrograph_%06d_Log.txt' % movieId
            
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
                
        return errors
    
