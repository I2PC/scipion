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
Protocols for particles filter operations.
"""

from pyworkflow.object import Float
from pyworkflow.protocol.params import (
    FloatParam, EnumParam, DigFreqParam, BooleanParam, PointerParam)
from pyworkflow.em.constants import FILTER_LOW_PASS, FILTER_HIGH_PASS, FILTER_BAND_PASS
from pyworkflow.em.protocol import ProtFilterParticles, ProtFilterVolumes
from pyworkflow.em.packages.xmipp3.constants import (
    FILTER_SPACE_FOURIER, FILTER_SPACE_REAL, FILTER_SPACE_WAVELET)
from protocol_process import XmippProcessParticles, XmippProcessVolumes



class XmippFilterHelper():
    """ Filter operations such as: Fourier or Gaussian. """
    tmpCTF = "ctf.xmd"

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
    #Wavelets decomposition base
    FM_DAUB4   = 0
    FM_DAUB12  = 1
    FM_DAUB20  = 2
    #Wavelet filters
    FM_REMOVE_SCALE      = 0
    FM_BAYESIAN          = 1
    FM_SOFT_THRESHOLDING = 2
    FM_ADAPTIVE_SOFT     = 3
    FM_CENTRAL           = 4

    #--------------------------- DEFINE param functions --------------------------------------------
    @classmethod
    def _defineProcessParams(cls, form):
        form.addParam('filterSpace', EnumParam, choices=['fourier', 'real', 'wavelets'],
                      default=FILTER_SPACE_FOURIER,
                      label="Filter space")
        form.addParam('filterModeFourier', EnumParam, choices=['low pass', 'high pass', 'band pass', 'ctf'],
                      default=cls.FM_BAND_PASS,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label="Filter mode",
                      help='Depending on the filter mode some frequency (freq.) components\n'
                           'are kept and some are removed.\n '
                           '_low pass_: components below *High freq.* are preserved.\n '
                           '_high pass_: components above *Low freq.* are preserved.\n '
                           '_band pass_: components between *Low freq.* and *High freq.* '
                           'are preserved. \n'
                           'ctf: apply first CTF in CTFset to all the particles'
        )

        form.addParam('filterModeReal', EnumParam, choices=['median'],
                      default=cls.FM_MEDIAN,
                      condition='filterSpace == %d' % FILTER_SPACE_REAL,
                      label="Filter mode",
                      help='median: replace each pixel with the median of neighboring pixels.\n'
        )
        form.addParam('filterModeWavelets', EnumParam, choices=['daub4','daub12','daub20'],
                      default=cls.FM_DAUB4,
                      condition='filterSpace == %d' % FILTER_SPACE_WAVELET,
                      label="Filter mode",
                      help='DAUB4: filter using the DAUB4 wavelet transform.\n '
        )

        #String that identifies filter in Fourier Space
        fourierCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_FOURIER,
                                                           cls.getModesCondition('filterModeFourier',
                                                                                  cls.FM_LOW_PASS,
                                                                                  cls.FM_HIGH_PASS,
                                                                                  cls.FM_BAND_PASS))
        #String that identifies filter in Real Space
        realCondition    = 'filterSpace == %d and (%s)' % (FILTER_SPACE_REAL,
                                                           cls.getModesCondition('filterModeReal',
                                                                                  cls.FM_MEDIAN))
        #String that identifies filter in Real Space
        waveletCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_WAVELET,
                                                           cls.getModesCondition('filterModeWavelets',
                                                                                  cls.FM_DAUB4,
                                                                                  cls.FM_DAUB12,
                                                                                  cls.FM_DAUB20))
        #fourier

        form.addParam('freqInAngstrom', BooleanParam, default=True,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label='Provide frequencies in Angstroms?',
                      help='If *Yes*, the frequency values for the filter\n'
                           'should be provided in Angstroms. If *No*, the\n'
                           'values should be in digital frequencies (between 0 and 0.5).')
        # Frequencies in Angstroms
        line = form.addLine('Frequency (A)',
                            condition=fourierCondition + ' and freqInAngstrom',
                            help='Range of frequencies to use in the filter')
        line.addParam('lowFreqA', FloatParam, default=60,
                      condition=cls.getModesCondition('filterModeFourier',
                                                       cls.FM_BAND_PASS, cls.FM_HIGH_PASS),
                      label='Lowest')
        line.addParam('highFreqA', FloatParam, default=10,
                      condition=cls.getModesCondition('filterModeFourier',
                                                       cls.FM_BAND_PASS, cls.FM_LOW_PASS),
                      label='Highest')

        # Digital frequencies
        line = form.addLine('Frequency (dig)',
                            condition=fourierCondition + ' and (not freqInAngstrom)',
                            help='Range of frequencies to use in the filter')
        line.addParam('lowFreqDig', DigFreqParam, default=0.02,
                      condition=cls.getModesCondition('filterModeFourier',
                                                       cls.FM_BAND_PASS, cls.FM_HIGH_PASS),
                      label='Lowest')
        line.addParam('highFreqDig', DigFreqParam, default=0.35,
                      condition=cls.getModesCondition('filterModeFourier',
                                                       cls.FM_BAND_PASS, cls.FM_LOW_PASS),
                      label='Highest')

        form.addParam('freqDecay', FloatParam, default=0.02,
                      condition=fourierCondition,
                      label='Frequency decay',
                      help='Length of the amplitude decay in a raised cosine')

        form.addParam('inputCTF', PointerParam,
                      condition='filterModeFourier == %d' % cls.FM_CTF,
                      label='CTF Object',
                      pointerClass='CTFModel',
                      help='Object with CTF information')

        #wavelets
        form.addParam('waveletMode',  EnumParam, choices=['remove_scale',
                                                          'bayesian(not implemented)',
                                                          'soft_thresholding',
                                                          'adaptive_soft',
                                                          'central'],
                      default=cls.FM_REMOVE_SCALE,
                      condition='filterSpace == %d' % FILTER_SPACE_WAVELET,
                      label='mode',
                      help='filter mode to be applied in wavelet space')

    @classmethod
    def getModesCondition(cls, filterMode, *filterModes):
        command = ' or '.join('%s==%d' % (filterMode,fm) for fm in filterModes)
        #print("command", command)
        return command
    #        return ' or '.join(['%s==%d' % (fmKey,fmValue) for fmKey,fmValue in zip(filterMode,filterModes)])

    #--------------------------- INSERT steps functions --------------------------------------------
    @classmethod
    def _insertProcessStep(cls, protocol):

        if protocol.filterSpace == FILTER_SPACE_FOURIER:
            if protocol.freqInAngstrom:
                lowFreq = protocol.lowFreqA.get()
                highFreq = protocol.highFreqA.get()
                samplingStr = ' --sampling %f ' % protocol.getInputSampling()
            else:
                lowFreq = protocol.lowFreqDig.get()
                highFreq = protocol.highFreqDig.get()
                samplingStr = ''

            mode = protocol.filterModeFourier.get()

            freqDecay = protocol.freqDecay.get()

            if mode == cls.FM_LOW_PASS:
                filterStr = " low_pass %f %f " % (highFreq, freqDecay)
            elif mode == cls.FM_HIGH_PASS:
                filterStr = " high_pass %f %f " % (lowFreq, freqDecay)
            elif mode == cls.FM_BAND_PASS:
                filterStr = " band_pass %f %f %f " % (lowFreq, highFreq, freqDecay)
            elif mode == cls.FM_CTF:
                ctfModel = protocol._getTmpPath(protocol.tmpCTF)
                filterStr = " ctf %s " % ctfModel
                # Save CTF model too
                protocol._insertFunctionStep("convertCTFXmippStep", ctfModel)
            else:
                raise Exception("Unknown fourier filter mode: %d" % mode)

            args = " --fourier " + filterStr + samplingStr
        elif protocol.filterSpace == FILTER_SPACE_REAL:
            mode = protocol.filterModeReal.get()
            if mode == cls.FM_MEDIAN:
                filterStr = " --median "
            else:
                raise Exception("Unknown real filter mode: %d" % mode)

            args = filterStr
        elif protocol.filterSpace == FILTER_SPACE_WAVELET:
            filterMode = protocol.filterModeWavelets.get()
            filterStr = " --wavelet "
            if filterMode == cls.FM_DAUB4:
                filterStr += " DAUB4 "
            elif filterMode == cls.FM_DAUB12:
                filterStr += " DAUB12 "
            elif filterMode == cls.FM_DAUB20:
                filterStr += " DAUB20 "
            else:
                raise Exception("Unknown wavelets filter mode: %d" % filterMode)

            waveletMode = protocol.waveletMode.get()
            if waveletMode == cls.FM_REMOVE_SCALE:
                filterStr += " remove_scale "
            elif waveletMode == cls.FM_BAYESIAN:
                raise Exception("Bayesian filter not implemented")
            elif waveletMode == cls.FM_SOFT_THRESHOLDING:
                filterStr += " soft_thresholding "
            elif waveletMode == cls.FM_ADAPTIVE_SOFT:
                filterStr += " adaptive_soft "
            elif waveletMode == cls.FM_CENTRAL:
                filterStr += " central "

            args = filterStr
        else:
            raise Exception("Unknown filter space: %d" % protocol.filterSpace.get())

        protocol._insertFunctionStep("filterStep",
                                 (protocol._args % {'inputFn': protocol.inputFn}) + args)

    def getInputSampling(self):
        """ Function to return the sampling rate of input objects.
        Should be implemented for filter volumes and particles.
        """
        pass


class XmippProtFilterParticles(ProtFilterParticles, XmippProcessParticles):
    """ Apply Fourier filters to a set of particles  """
    _label = 'filter particles'

    def __init__(self, **kwargs):
        ProtFilterParticles.__init__(self, **kwargs)
        XmippProcessParticles.__init__(self, **kwargs)
        self._program = "xmipp_transform_filter"

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippFilterHelper._defineProcessParams(form)
        
    def _insertProcessStep(self):
        XmippFilterHelper._insertProcessStep(self)

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


class XmippProtFilterVolumes(ProtFilterVolumes, XmippProcessVolumes):
    """ Apply Fourier filters to a set of volumes """
    _label = 'filter volumes'

    #--------------------------- UTILS functions ---------------------------------------------------

    def __init__(self, **kwargs):
        ProtFilterVolumes.__init__(self, **kwargs)
        XmippProcessVolumes.__init__(self, **kwargs)
        self._program = "xmipp_transform_filter"

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippFilterHelper._defineProcessParams(form)
        
    def _insertProcessStep(self):
        XmippFilterHelper._insertProcessStep(self)
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def filterStep(self, args):
        if self._isSingleInput():
            args += " -o %s" % self.outputStk
        else:
            args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)

        self.runJob("xmipp_transform_filter", args)

    def getInputSampling(self):
        return self.inputVolumes.get().getSamplingRate()
