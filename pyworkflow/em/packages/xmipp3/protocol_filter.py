# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
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
This sub-package contains protocol for particles filters operations
"""

from pyworkflow.protocol.params import FloatParam, EnumParam, DigFreqParam
from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from protocol_process import XmippProcess, XmippProcessParticles, XmippProcessVolumes

from pyworkflow.em.constants import *
from constants import *


# Some Filter Modes constants to be used locally
# the special cases of low pass, high pass and band pass 
# should preserve the em.constants values 0, 1 and 2 respectively 
# for properly working of the wizards
#Fourier filters
FM_LOW_PASS  = FILTER_LOW_PASS   #0
FM_HIGH_PASS = FILTER_HIGH_PASS  #1
FM_BAND_PASS = FILTER_BAND_PASS  #2
FM_CTF       = 3
# Real Space Filters
FM_MEDIAN = 0
#WaveLets decomposition base
FM_DAUB4   = 0
FM_DAUB12  = 1
FM_DAUB20  = 2
#Wavelet filters
FM_REMOVE_SCALE      = 0
FM_BAYESIAN          = 1
FM_SOFT_THRESHOLDING = 2
FM_ADAPTIVE_SOFT     = 3
FM_CENTRAL           = 4


class XmippProtFilter():
    """ Some filters operations such as: Fourier or Gaussian. """
    tmpCTF = "ctf.xmd"

    def __init__(self, **args):
        self._program = "xmipp_transform_filter"

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        form.addParam('filterSpace', EnumParam, choices=['fourier', 'real', 'wavelets'],
                      default=FILTER_SPACE_FOURIER,
                      label="Filter space")
        form.addParam('filterModeFourier', EnumParam, choices=['low pass', 'high pass', 'band pass', 'ctf'],
                      default=FM_BAND_PASS,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label="Filter mode",
                      help='Depending on the filter mode some frequency (freq.) components \n'
                           'are keept and some are removed. \n '
                           '_low pass_: components below *High freq.* are preserved. \n '
                           '_high pass_: components above *Low freq.* are preserved. \n '
                           '_band pass_: components between *Low freq.* and *High freq.* '
                           'are preserved. \n '
                           'ctf:  apply first CTF in CTFset to all the particles'
        )

        form.addParam('filterModeReal', EnumParam, choices=['median'],
                      default=FM_MEDIAN,
                      condition='filterSpace == %d' % FILTER_SPACE_REAL,
                      label="Filter mode",
                      help='median: replacing each pixel with the median of neighboring pixels. \n'
        )
        form.addParam('filterModeWaveLets', EnumParam, choices=['daub4','daub12','daub20'],
                      default=FM_DAUB4,
                      condition='filterSpace == %d' % FILTER_SPACE_WAVELET,
                      label="Filter mode",
                      help='DAUB4: filter using the DAUB4 wavelet transform. \n '
        )

        #String that identifies filter in Fourier Space
        fourierCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_FOURIER,
                                                           self.getModesCondition('filterModeFourier',
                                                                                  FM_LOW_PASS,
                                                                                  FM_HIGH_PASS,
                                                                                  FM_BAND_PASS))
        #String that identifies filter in Real Space
        realCondition    = 'filterSpace == %d and (%s)' % (FILTER_SPACE_REAL,
                                                           self.getModesCondition('filterModeReal',
                                                                                  FM_MEDIAN))
        #String that identifies filter in Real Space
        waveLetCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_WAVELET,
                                                           self.getModesCondition('filterModeWaveLets',
                                                                                  FM_DAUB4,
                                                                                  FM_DAUB12,
                                                                                  FM_DAUB20))
        #fourier

        form.addParam('freqInAngstrom', BooleanParam, default=True,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label='Provide frequencies in Angstroms?',
                      help='If *Yes*, the frequencies values for filter\n'
                           'should be provided in Angstroms. If *No*, the\n'
                           'values should be in digital frequencies (between 0 and 0.5).')
        # Frequencies in Angstroms
        line = form.addLine('Frequency (A)',
                            condition=fourierCondition + ' and freqInAngstrom',
                            help='Range to apply the filter')
        line.addParam('lowFreqA', FloatParam, default=60,
                      condition=self.getModesCondition('filterModeFourier',FM_BAND_PASS
                                                       , FM_HIGH_PASS),
                      label='Lowest')
        line.addParam('highFreqA', FloatParam, default=10,
                      condition=self.getModesCondition('filterModeFourier',FM_BAND_PASS
                                                       , FM_LOW_PASS),
                      label='Highest')
        
        # Digital frequencies
        line = form.addLine('Frequency (dig)',
                            condition=fourierCondition + ' and (not freqInAngstrom)',
                            help='Range to apply the filter')
        line.addParam('lowFreqDig', DigFreqParam, default=0.02,
                      condition=self.getModesCondition('filterModeFourier',FM_BAND_PASS
                                                       , FM_HIGH_PASS),
                      label='Lowest')
        line.addParam('highFreqDig', DigFreqParam, default=0.35,
                      condition=self.getModesCondition('filterModeFourier',FM_BAND_PASS
                                                       , FM_LOW_PASS),
                      label='Highest')        
        
        form.addParam('freqDecay', FloatParam, default=0.02,
                      condition=fourierCondition,
                      label='Frequency decay',
                      help='It describes the length of the amplitude decay in a raised cosine')

        form.addParam('inputCTF', PointerParam,
                      condition='filterModeFourier == %d' % FM_CTF,
                      label='CTF Object',
                      pointerClass='CTFModel',
                      help='Object with CTF information')

        #wavelets
        form.addParam('waveLetMode',  EnumParam, choices=['remove_scale',
                                                          'bayesian(not implemented)',
                                                          'soft_thresholding',
                                                          'adaptive_soft',
                                                          'central'],
                      default=FM_REMOVE_SCALE,
                      condition='filterSpace == %d' % FILTER_SPACE_WAVELET,
                      label='mode',
                      help='filter mode to be applied in wavelet space')

    def getModesCondition(self,filterMode, *filterModes):
        command = ' or '.join(['%s==%d' % (filterMode,fm) for fm in filterModes])
        #print("command", command)
        return command
    #        return ' or '.join(['%s==%d' % (fmKey,fmValue) for fmKey,fmValue in zip(filterMode,filterModes)])

    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertProcessStep(self):

        inputFn = self.inputFn
        ctfModel = self._getTmpPath(self.tmpCTF)
        
        if self.filterSpace == FILTER_SPACE_FOURIER:
            if self.freqInAngstrom:
                lowFreq = self.lowFreqA.get()
                highFreq = self.highFreqA.get()
                samplingStr = ' --sampling %f' % self.getInputSampling()
            else:
                lowFreq = self.lowFreqDig.get()
                highFreq = self.highFreqDig.get()
                samplingStr = ''
                
            filterMode = self.filterModeFourier.get()
                        
            freqDecay = self.freqDecay.get()

            if filterMode == FM_LOW_PASS:
                filterStr = "low_pass %(highFreq)f %(freqDecay)f "
            elif filterMode == FM_HIGH_PASS:
                filterStr = "high_pass %(lowFreq)f %(freqDecay)f "
            elif filterMode == FM_BAND_PASS:
                filterStr = "band_pass %(lowFreq)f %(highFreq)f %(freqDecay)f "
            elif filterMode == FM_CTF:
                filterStr = "ctf %(ctfModel)s "

            else:
                raise Exception("Unknown fourier filter mode: %d" % filterMode)

            args = (self._args + " --fourier " + filterStr + samplingStr ) % locals()

        elif self.filterSpace == FILTER_SPACE_REAL:
            filterMode = self.filterModeReal.get()
            if filterMode == FM_MEDIAN:
                filterStr = "--median "
            else:
                raise Exception("Unknown real filter mode: %d" % filterMode)

        elif self.filterSpace == FILTER_SPACE_WAVELET:
            filterMode = self.filterModeWavelets.get()
            filterStr = "--wavelets "
            if filterMode == FM_DAUB4:
                filterStr += "DAUB4 "
            elif filterMode == FM_DAUB12:
                filterStr += "DAUB12 "
            elif filterMode == FM_DAUB20:
                filterStr += "DAUB20 "
            else:
                raise Exception("Unknown wavelets filter mode: %d" % filterMode)
            waveLetMode = self.waveLetMode.get()
            if waveLetMode == FM_REMOVE_SCALE:
                filterStr += "remove_scale "
            elif waveLetMode == FM_BAYESIAN:
                raise Exception("Bayesian filter not implemented")
            elif waveLetMode == FM_SOFT_THRESHOLDING:
                filterStr += "soft_thresholding "
            elif waveLetMode == FM_ADAPTIVE_SOFT:
                filterStr += "adaptive_soft "
            elif waveLetMode == FM_CENTRAL:
                filterStr += "central "
        else:
            raise Exception("Unknown filter space: %d" % self.filterSpace.get())

        if filterMode == FM_CTF:
            self._insertFunctionStep("convertCTFXmippStep", ctfModel)# save CTF model
        self._insertFunctionStep("filterStep", args)
        
    def getInputSampling(self):
        """ Function to return the sampling rate of input objects. 
        Should be implemented for filter volumes and particles.
        """
        pass


class XmippProtFilterParticles(ProtFilterParticles,
                               XmippProcessParticles,
                               XmippProtFilter):
    """ Apply Fourier filters to a set of particles  """
    _label = 'filter particles'

    def __init__(self, **args):
        ProtFilterParticles.__init__(self, **args)
        XmippProcessParticles.__init__(self)
        XmippProtFilter.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtFilter._defineProcessParams(self, form)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertCTFXmippStep(self, ctfModel):
        from convert import writeCTFModel
        #defocus comes from here
        ctf = self.inputCTF.get()
        inputSet = self.inputParticles.get()
        #is there any voltage, sampling and amplitud contrast available?
        acquisition = inputSet.getAcquisition()
        sampling = inputSet.getSamplingRate()
        # Spherical aberration in mm
        ctf._xmipp_ctfVoltage = Float(acquisition.getVoltage())
        ctf._xmipp_ctfSphericalAberration = Float(acquisition.getSphericalAberration())
        ctf._xmipp_ctfQ0 = Float(acquisition.getAmplitudeContrast())
        ctf._xmipp_ctfSamplingRate = Float(sampling)
        writeCTFModel(self.inputCTF.get(), ctfModel)

    #--------------------------- STEPS functions ---------------------------------------------------
    def filterStep(self, args):
        args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        self.runJob("xmipp_transform_filter", args)
        
    def getInputSampling(self):
        return self.inputParticles.get().getSamplingRate()
    
    
class XmippProtFilterVolumes(ProtFilterVolumes, XmippProcessVolumes, XmippProtFilter):
    """ Apply Fourier filters to a set of volumes """
    _label = 'filter volumes'    #--------------------------- UTILS functions ---------------------------------------------------
    
    def __init__(self, **args):
        ProtFilterVolumes.__init__(self, **args)
        XmippProcessVolumes.__init__(self)
        XmippProtFilter.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtFilter._defineProcessParams(self, form)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def filterStep(self, args):
        if self._isSingleInput():
            args += " -o %s" % self.outputStk
        else:
            args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        
        self.runJob("xmipp_transform_filter", args)

    def getInputSampling(self):
        return self.inputVolumes.get().getSamplingRate()
