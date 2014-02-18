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

from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from protocol_process import XmippProcess, XmippProcessParticles, XmippProcessVolumes

from pyworkflow.em.constants import *
from constants import *



class XmippProtFilter():
    """ Some filters operations such as: Fourier or Gaussian. """
    _label = 'filter particles'
    
    def __init__(self, **args):
        self._program = "xmipp_transform_filter"
        
    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam, choices=['fourier', 'gaussian'], 
                      default=FILTER_FOURIER, display=EnumParam.DISPLAY_LIST,
                      label="Filter type")
        form.addParam('fourierMode', EnumParam, choices=['low pass', 'high pass', 'band pass'],
                      default=FILTER_BAND_PASS,
                      condition='filterType == %d' % FILTER_FOURIER,
                      label="Fourier mode", 
                      help='Depending on the filter mode some frequency (freq.) components \n'
                           'are keept and some are removed. \n '
                           '_low pass_: components below *High freq.* are preserved. \n '
                           '_high pass_: components above *Low freq.* are preserved. \n '
                           '_band pass_: components between *Low freq.* and *High freq.* are preserved. \n ')
        form.addParam('lowFreq', DigFreqParam, default=0.02, 
                      condition='filterType == %d and fourierMode != %d' % (FILTER_FOURIER, FILTER_LOW_PASS),
                      label='Low frequency (0 < f < 0.5)',
                      help='Low frequency cuttoff to apply the filter. \n ')          
        form.addParam('highFreq', DigFreqParam, default=0.35, 
                      condition='filterType == %d and fourierMode != %d' % (FILTER_FOURIER, FILTER_HIGH_PASS),
                      label='High frequency (0 < f < 0.5)', 
                      help='High frequency cuttoff to apply the filter.')          
        form.addParam('freqDecay', FloatParam, default=0.02, 
                      condition='filterType == %d' % FILTER_FOURIER, 
                      label='Frequency decay',
                      help='It describes the length of the amplitude decay in a raised cosine') 
        form.addParam('freqSigma', FloatParam, default=0.04, 
                      condition='filterType == %d' % FILTER_GAUSSIAN,
                      label='Frequency sigma',
                      help='Remind that the Fourier frequency is normalized between 0 and 0.5')   
    
    def _getCommand(self, inputFn):
        if self.filterType == FILTER_FOURIER:
            lowFreq = self.lowFreq.get()
            highFreq = self.highFreq.get()
            fourierMode = self.fourierMode.get()
            freqDecay = self.freqDecay.get()
            
            if fourierMode == FILTER_LOW_PASS:
                filterStr = "low_pass %(highFreq)f"
            elif fourierMode == FILTER_HIGH_PASS:
                filterStr = "high_pass %(lowFreq)f"
            else:
                filterStr = "band_pass %(lowFreq)f %(highFreq)f"
            args = (self._args + " --fourier " + filterStr + " %(freqDecay)f ") % locals()
        
        elif self.filterType == FILTER_GAUSSIAN:
            freqSigma = self.freqSigma.get()
            args = (self._args + " --fourier gaussian %(freqSigma)f ") % locals()
            
        else:
            raise Exception("Unknown filter type: %d" % self.filterType.get())
            
        return args


class XmippProtFilterParticles(ProtFilterParticles, XmippProtFilter, XmippProcessParticles):
    """ Apply some filter to SetOfParticles """
    _label = 'filter particles'
    
    def __init__(self, **args):
        ProtFilterParticles.__init__(self)
        XmippProtFilter.__init__(self, **args)
        XmippProcessParticles.__init__(self, **args)


class XmippProtFilterVolumes(ProtFilterVolumes, XmippProtFilter, XmippProcessVolumes):
    """ Apply some filter to SetOfParticles """
    _label = 'filter volumes'
    
    def __init__(self, **args):
        ProtFilterVolumes.__init__(self)
        XmippProtFilter.__init__(self, **args)
        XmippProcessVolumes.__init__(self, **args)

