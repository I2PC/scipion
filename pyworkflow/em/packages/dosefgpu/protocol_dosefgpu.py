# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from os.path import join

from pyworkflow.utils import makePath, moveFile
from pyworkflow.protocol.params import StringParam, IntParam, LEVEL_ADVANCED
from pyworkflow.em.protocol import ProtProcessMovies



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
-frs       FileName.mrc      Uncorrected sum
-fcs       FileName.mrc      Corrected sum
-fct       FileName.mrc      Corrected stack
-fcm       FileName.mrc      CC map
-flg       FileName.txt      Log file
                      """)
        
        form.addSection('Crop and binning')
        line = form.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', IntParam, default=0, label='X')
        line.addParam('cropOffsetY', IntParam, default=0, label='Y')
        form.addParam('cropDimension', IntParam, default=0, 
                      label='Crop dimension',
                      help='How many pixel to crop from offset\n'
                           'If equal to 0, use maximum size.')
        form.addParam('binFactor', IntParam, default=1, 
                      label='Binning factor',
                      help='1x or 2x. Bin stack before processing.')
        
                     
    def _processMovie(self, movieId, movieName):
        inputName = movieName
        micrographName = self._getMicName(movieId)
        gainFile = self.inputMovies.get().getGain()
        gpuId = self.gpuId.get()

        if gainFile is not None:
            # Apply the gain correction to flat the raw movie
            correctedName = movieName.replace('.mrc', '_corrected.mrc')
            
            self.runJob('dosefgpu_flat', 
                        '%(inputName)s %(correctedName)s %(gainFile)s %(gpuId)s' % locals())
            
            inputName = correctedName
        
        args = {'-crx': self.cropOffsetX.get(),
                '-cry': self.cropOffsetY.get(),
                '-crd': self.cropDimension.get(),
                '-bin': self.binFactor.get(),
                '-nst': self.alignFrame0.get(),
                '-ned': self.alignFrameN.get(),
                '-nss': self.sumFrame0.get(),
                '-nes': self.sumFrameN.get(),
                '-gpu': gpuId,
                }
        
        command = '%(inputName)s -fcs %(micrographName)s ' % locals()
        command += ' '.join(['%s %s' % (k, v) for k, v in args.iteritems()])
        command += ' ' + self.extraParams.get()

        self.runJob('dosefgpu_driftcorr', command)
        
    def printJob(self, program, args):
        self.info(">>> Running: %s %s" % (program, args))
        
    def createOutputStep(self):
        inputMovies = self.inputMovies.get()
        micSet = self._createSetOfMicrographs()
        micSet.copyInfo(inputMovies)
        # Also create a Set of Movies with the alignment parameters
        movieSet = self._createSetOfMovies()
        movieSet.copyInfo(inputMovies)
        
        # Create a folder to store the resulting micrographs
        micsFolder = self._getPath('micrographs')
        makePath(micsFolder)
        
        for movie in inputMovies():
            movieId = movie.getObjId()
            micName = self._getMicName(movieId)
            micNameSrc = join(self._getMovieFolder(movieId), micName)
            micNameDst = join(micsFolder, micName)
            # Move the resulting micrograph before delete of movies folder
            moveFile(micNameSrc, micNameDst)            
            mic = micSet.ITEM_TYPE()
            mic.setFileName(micNameDst)
            
            # Parse the alignment parameters and store the log files
            
            micSet.append(mic)
            
            
        self._defineOutputs(outputMicrographs=micSet)
        self._defineTransformRelation(inputMovies, micSet)
        
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
                
        return errors
    
