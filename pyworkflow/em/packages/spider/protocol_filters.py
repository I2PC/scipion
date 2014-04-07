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
from spider import *
        
      
class SpiderProtFilter(ProtFilterParticles, SpiderProtocol):
    """ Apply Fourier filters to an image or volume 
    using Spider FQ or FQ NP. 
    
    To improve boundary quality the image is padded 
    with the average value to twice the original size 
    during filtration if padding is selected.  
    
    See more documentation in: 
    [[http://spider.wadsworth.org/spider_doc/spider/docs/man/fq.html][FQ Spider online manual]]
    """
    _label = 'filter particles'
    
    def __init__(self, **args):
        ProtFilterParticles.__init__(self)
        SpiderProtocol.__init__(self)
        self._op = "FQ"
        self._params = {'ext': 'stk', 
                        'particles': 'particles_filtered',
                        'particlesSel': 'particles_filtered_sel'}

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam, choices=['Top-hat', 'Gaussian', 'Fermi', 'Butterworth', 'Raised cosine'],
                      label="Filter type", default=3,
                      help="""Select what type of filter do you want to apply.
                      
*Top-hat*: Filter is a "top-hat" function 
that truncates the Fourier transform at spatial frequency.
            
*Gaussian*: Filter is the Gaussian function: EXP(-F**2 / (2 * SPF**2)), 
where F is the frequency.
              
*Fermi*: Filter is: 1 / (1 + EXP((F - SPF) / T)) which negotiates 
between "Top-hat" and Gaussian characteristics, depending on 
the value of the temperature.

*Butterworth* Filter is: 1 / (SQRT(1 + F / RAD)**(2 * ORDER)) 
The ORDER determines the filter fall off and RAD corresponds 
to the cut-off radius. 

*Raised cosine* Filter is: 

See detailed description of the filter in [[http://spider.wadsworth.org/spider_doc/spider/docs/man/fq.html][FQ Spider online manual]]
                           """)
        form.addParam('filterMode', EnumParam, choices=['low-pass', 'high-pass'],
                      label='Filter mode', default=0,
                      )
        form.addParam('usePadding', BooleanParam, default=True, 
                      label='Use padding?',
                      help='If set to *Yes*, to improve boundary quality\n'
                           'the image is padded with the average value to\n'
                           'twice the original size during filtration.\n\n'
                           'If *No* padding is applied, this may lead to\n'
                           'artifacts near boundary of image.')
        form.addParam('filterRadius', DigFreqParam, default=0.12, 
                      label='Filter radius (0 < f < 0.5)',
                      condition='filterType <= %d' % FILTER_GAUSSIAN,
                      help='Low frequency cuttoff to apply the filter.\n')  
        line = form.addLine('Frequency', help='Range to apply the filter. Expected values between 0 and 0.5.')
        line.addParam('lowFreq', DigFreqParam, default=0.1, 
                      condition='filterType > %d' % FILTER_GAUSSIAN,
                      label='Lowest')
        line.addParam('highFreq', 
                      DigFreqParam, default=0.2, 
                      condition='filterType > %d' % FILTER_GAUSSIAN,
                      label='Highest')
#        form.addParam('lowFreq', DigFreqParam, default=0.1, 
#                 
#         label='Low Frequency (0 < f < 0.5)',
#                      condition='filterType > %d' % FILTER_GAUSSIAN,
#                      help='Low frequency cuttoff to apply the filter.\n')          
#        form.addParam('highFreq', DigFreqParam, default=0.2, 
#                      label='High Frequency (0 < f < 0.5)', 
#                      condition='filterType > %d' % FILTER_GAUSSIAN,
#                      help='High frequency cuttoff to apply the filter.\n'
#                           'Set to 0.5 for a <high pass> filter.')          
        form.addParam('temperature', FloatParam, default=0.3, 
                      label='Temperature T:',
                      condition='filterType == %d' % FILTER_FERMI,
                      help='Enter a temperature parameter T The filter falls off roughly within \n'
                           'this reciprocal distance (in terms of frequency units).')     
        
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        # Define some names
        self.particlesStk = self._getPath('%(particles)s.%(ext)s' % self._params)
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('filterStep', self.filterType.get())

    #--------------------------- STEPS functions --------------------------------------------       
    def filterStep(self, filterType):
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
            
        spi.close()
            
        self._leaveWorkingDir() # Go back to project dir
            
        for i, img in enumerate(particles):
            img.setLocation(i+1, self.particlesStk)
            imgSet.append(img)            
            
        self._defineOutputs(outputParticles=imgSet)
        
#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
    
    
    
