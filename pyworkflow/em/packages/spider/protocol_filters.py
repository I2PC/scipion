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
from pyworkflow.utils import removeBaseExt
from constants import *
from spider import SpiderShell, SpiderProtocol
        

      
class SpiderProtFilter(ProtFilterParticles, SpiderProtocol):
    """ Protocol for Spider filters. """
    _label = 'filters'
    
    def __init__(self):
        ProtFilterParticles.__init__(self)
        SpiderProtocol.__init__(self)
        self._op = "FQ"
        self._params = {'ext': 'stk', 
                        'particles': 'particles_filtered',
                        'particlesSel': 'particles_filtered_sel'}

    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam, choices=['Top-hat', 'Gaussian', 'Fermi', 'Butterworth', 'Raised cosine'],
                      label="Filter type", default=3,
                      help="""Select what type of filter do you want to apply.
                      
<Top-hat>   <low-pass>. truncation. Filter is a "top-hat" function that truncates
                the Fourier transform at spatial frequency: SPF.
            <high-pass>. truncation. Filter is inverse "top-hat" function that 
                passes the Fourier transform beyond spatial frequency: SPF.
            
<Gaussian>  <low-pass>. 
                Filter is the Gaussian function: EXP(-F**2 / (2 * SPF**2)), 
                where F is the frequency.
            <high-pass>. Filter is complement of the Gaussian function: 
                1 - EXP(-F**2 / (2 * SPF**2)), where F is the frequency.
              
<Fermi>     <low-pass>. Filter is: 1 / (1 + EXP((F - SPF) / T)) which negotiates 
                between "Top-hat" and Gaussian characteristics, depending on 
                the value of the temperature: T (see below).
            <high-pass>. Filter is: 1 / (1 + EXP((F - SPF) / -T)) 
                (Same as in low-pass, but with T replaced by -T).

<Butterworth> <low-pass>. Filter is: 1 / (SQRT(1 + F / RAD)**(2 * ORDER)) 
                where 
                ORDER = (2 * log(eps/SQRT(**2-1)) ) / (log(Flow/Fup)) 
                RAD = Flow / ((eps)**(2 / ORDER)) 
                In the Butterworth filter the ORDER determines the filter 
                fall off and RAD corresponds to the cut-off radius. 
                Frequencies below the lower frequency are preserved, 
                frequencies above the upper frequency are removed, 
                with a smooth transition in between lower and upper 
                limiting frequencies.
              <high-pass>. Filter is: 1 - (1 / (SQRT(1 + F / RAD)**(2 * ORDER))) 
                Frequencies below the lower frequency are removed, 
                frequencies above upper frequency are preserved, 
                with a smooth transition in between lower and upper 
                limiting frequencies.

<Raised cosine> <low-pass>. Filter is: 
                    0.5 * (COS(PI * (F - Flow) / (Flow - Fup)) + 1) 
                       if Flow < F < Fup, 
                    1  if F < Flow, and 
                    0  if F > Fup. 
                <high-pass>. Filter is: 
                    0.5 * (-COS(PI*(F - Flow) / (Flow - Fup)) + 1) 
                       if Flow < F < Fup 
                    1  if F < Flow, and 
                    1  if F > Fup. 
                           """)
        form.addParam('filterMode', EnumParam, choices=['low-pass', 'high-pass'],
                      label='Filter mode', default=0,
                      )
        form.addParam('usePadding', BooleanParam, default=True, 
                      label='Use padding?',
                      help="If <No> padding is applied, this may lead to artifacts near boundary of image,\n"
                           "we suggest use of padding to avoid this.\n")  
        form.addParam('filterRadius', DigFreqParam, default=0.12, 
                      label='Filter radius (0 < f < 0.5)',
                      condition='filterType <= %d' % FILTER_GAUSSIAN,
                      help='Low frequency cuttoff to apply the filter.\n')  
        form.addParam('lowFreq', DigFreqParam, default=0.1, 
                 
         label='Low Frequency (0 < f < 0.5)',
                      condition='filterType > %d' % FILTER_GAUSSIAN,
                      help='Low frequency cuttoff to apply the filter.\n')          
        form.addParam('highFreq', DigFreqParam, default=0.2, 
                      label='High Frequency (0 < f < 0.5)', 
                      condition='filterType > %d' % FILTER_GAUSSIAN,
                      help='High frequency cuttoff to apply the filter.\n'
                           'Set to 0.5 for a <high pass> filter.')          
        form.addParam('temperature', FloatParam, default=0.3, 
                      label='Temperature T:',
                      condition='filterType == %d' % FILTER_FERMI,
                      help='Enter a temperature parameter T The filter falls off roughly within \n'
                           'this reciprocal distance (in terms of frequency units).')     
        
    def _insertAllSteps(self):
        # Define some names
        self.particlesStk = self._getPath('%(particles)s.%(ext)s' % self._params)
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('filterParticles', self.filterType.get())

    def filterParticles(self, filterType):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        OP = self._op
        args = []

        if not self.usePadding:
            OP += ' NP'
        
        if filterType <= FILTER_GAUSSIAN:
            args.append(self.filterRadius.get())
        else:
            args.append('%f %f' % (self.lowFreq.get(), self.highFreq.get()))
            
        if filterType == FILTER_FERMI:
            args.append(self.temperature.get())
        
        # Map to spected filter number in Spider for operation FQ    
        filterNumber = filterType * 2 + 1
        # Consider low-pass or high-pass
        filterNumber += self.filterMode.get()
        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())

        self._enterWorkingDir() # Do operations inside the run working dir
        
        spi = SpiderShell(ext=self._params['ext']) # Create the Spider process to send commands        
        particlesStk = removeBaseExt(self.particlesStk)
        
        # Run a loop for filtering
        locStr = particlesStk + '@******[part]'
        cmds = ['do lb5 [part] = 1,%d' % n,
                OP, locStr, locStr, filterNumber] + args + ['lb5']
        
        for c in cmds:
            spi.runCmd(c)
            
        for i in range(1, n+1):
            img = Image()
            img.setLocation(i, self.particlesStk)
            imgSet.append(img)
            
        spi.close()
            
        self._leaveWorkingDir() # Go back to project dir
            
        self._defineOutputs(outputParticles=imgSet)
        
    def _summary(self):
        summary = []
        return summary
    
