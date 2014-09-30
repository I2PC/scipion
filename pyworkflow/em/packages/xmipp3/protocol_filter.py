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
FM_LOW_PASS = FILTER_LOW_PASS
FM_HIGH_PASS = FILTER_HIGH_PASS
FM_BAND_PASS = FILTER_BAND_PASS



class XmippProtFilter():
    """ Some filters operations such as: Fourier or Gaussian. """
    
    def __init__(self, **args):
        self._program = "xmipp_transform_filter"
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        form.addParam('filterSpace', EnumParam, choices=['fourier', 'real'], 
                      default=FILTER_SPACE_FOURIER,
                      label="Filter space")
        form.addParam('fourierMode', EnumParam, choices=['low pass', 'high pass', 'band pass'],
                      default=FM_BAND_PASS,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label="Fourier mode", 
                      help='Depending on the filter mode some frequency (freq.) components \n'
                           'are keept and some are removed. \n '
                           '_low pass_: components below *High freq.* are preserved. \n '
                           '_high pass_: components above *Low freq.* are preserved. \n '
                           '_band pass_: components between *Low freq.* and *High freq.* are preserved. \n ')
        
        freqCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_FOURIER, 
                                                       self.getModesCondition(FM_HIGH_PASS, FM_LOW_PASS, FM_BAND_PASS))
        line = form.addLine('Frequency',
                            condition=freqCondition, 
                            help='Range to apply the filter')        
        line.addParam('lowFreq', DigFreqParam, default=0.02, 
                      condition=self.getModesCondition(FM_BAND_PASS, FM_HIGH_PASS),
                      label='Lowest')
        line.addParam('highFreq', DigFreqParam, default=0.35, 
                      condition=self.getModesCondition(FM_BAND_PASS, FM_LOW_PASS),
                      label='Highest')
        
        form.addParam('freqDecay', FloatParam, default=0.02, 
                      condition=freqCondition,
                      label='Frequency decay',
                      help='It describes the length of the amplitude decay in a raised cosine')
         
        form.addParam('freqSigma', FloatParam, default=0.04, 
                      condition='filterSpace == %d' % FILTER_SPACE_REAL,
                      label='Frequency sigma',
                      help='Remind that the Fourier frequency is normalized between 0 and 0.5')

    def getModesCondition(self, *filterModes):
        return ' or '.join(['fourierMode==%d' % fm for fm in filterModes])
    
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertProcessStep(self):
        
        inputFn = self.inputFn
        
        if self.filterSpace == FILTER_SPACE_FOURIER:
            lowFreq = self.lowFreq.get()
            highFreq = self.highFreq.get()
            fourierMode = self.fourierMode.get()
            freqDecay = self.freqDecay.get()
            
            if fourierMode == FM_LOW_PASS:
                filterStr = "low_pass %(highFreq)f"
            elif fourierMode == FM_HIGH_PASS:
                filterStr = "high_pass %(lowFreq)f"
            else:
                filterStr = "band_pass %(lowFreq)f %(highFreq)f"
            args = (self._args + " --fourier " + filterStr + " %(freqDecay)f ") % locals()
        
        elif self.filterSpace == FILTER_SPACE_REAL:
            freqSigma = self.freqSigma.get()
            args = (self._args + " --fourier gaussian %(freqSigma)f ") % locals()
            
        else:
            raise Exception("Unknown filter space: %d" % self.filterSpace.get())
        
        self._insertFunctionStep("filterStep", args)


class XmippProtFilterParticles(ProtFilterParticles, XmippProcessParticles, XmippProtFilter):
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
    def filterStep(self, args):
        args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        self.runJob("xmipp_transform_filter", args)


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

