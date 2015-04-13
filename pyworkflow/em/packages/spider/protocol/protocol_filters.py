# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Tapu Shaikh            (shaikh@ceitec.muni.cz)
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
Particle filter operations.
"""

from pyworkflow.em.protocol import ProtFilterParticles  
from pyworkflow.protocol.params import EnumParam, BooleanParam, DigFreqParam, FloatParam
from pyworkflow.utils.path import removeBaseExt

from ..constants import FILTER_SPACE_REAL, FILTER_FERMI, FILTER_BUTTERWORTH, FILTER_LOWPASS, FILTER_HIGHPASS
from ..spider import SpiderShell
from protocol_base import SpiderProtocol
        
      
class SpiderProtFilter(ProtFilterParticles, SpiderProtocol):
    """ Apply Fourier filters to an image or volume 
    using Spider FQ or FQ NP. 
    
    To improve boundary quality the image is padded 
    with the average value to twice the original size 
    during filtration if padding is selected.  
    
    See more documentation at: 
    [[http://spider.wadsworth.org/spider_doc/spider/docs/man/fq.html][SPIDER's FQ online manual]]
    """
    _label = 'filter particles'
    
    def __init__(self, **kwargs):
        ProtFilterParticles.__init__(self, **kwargs)
        SpiderProtocol.__init__(self, **kwargs)
        self._op = "FQ"
        self._params = {'ext': 'stk', 
                        'particles': 'particles_filtered',
                        'particlesSel': 'particles_filtered_sel'}

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineProcessParams(self, form):
        form.addParam('filterType', EnumParam, choices=['Top-hat', 'Gaussian', 'Fermi', 'Butterworth', 'Raised cosine'],
                      label="Filter type", default=FILTER_BUTTERWORTH,
                      help="""Select what type of filter do you want to apply.
                      
*Top-hat*: Filter is a "top-hat" function 
that truncates the Fourier transform at spatial frequency.
            
*Gaussian*: Filter is the Gaussian function: EXP(-F**2 / (2 * SPF**2)), 
where F is the frequency.
              
*Fermi*: Filter is: 1 / (1 + EXP((F - SPF) / T)) which negotiates 
between "Top-hat" and Gaussian characteristics, depending on 
the value of the temperature T.

*Butterworth* Filter is: 1 / (SQRT(1 + F / RAD)**(2 * ORDER)) 
The ORDER determines the filter fall off and RAD corresponds 
to the cut-off radius. 

*Raised cosine* Filter is: 0.5 * (COS(PI * (F - Flow) / (Flow - Fup)) + 1) 
if Flow < F < Fup, 1 if F < Flow, and 0 if F > Fup

See detailed description of the filter at [[http://spider.wadsworth.org/spider_doc/spider/docs/man/fq.html][SPIDER's FQ online manual]]
                           """)
        form.addParam('filterMode', EnumParam, choices=['low-pass', 'high-pass'],
                      label='Filter mode', default=FILTER_LOWPASS,
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
                      condition='filterType <= %d or filterType == %d' % (FILTER_SPACE_REAL,FILTER_FERMI),
                      help='Low frequency cutoff to apply the filter.\n')  
        
        line = form.addLine('Frequency', 
                            condition='filterType > %d and filterType != %d' % (FILTER_SPACE_REAL,FILTER_FERMI),
                            help='Range to apply the filter. Expected values between 0 and 0.5.')
        line.addParam('lowFreq', DigFreqParam, condition='filterMode==%d' % FILTER_HIGHPASS, default=0.1,
                    label='Lowest')
        line.addParam('highFreq', DigFreqParam, condition='filterMode==%d' % FILTER_LOWPASS,
                    default=0.2, label='Highest')
         
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
        
        if filterType <= FILTER_SPACE_REAL:
            args.append(self.filterRadius.get())
        else:
            args.append('%f %f' % (self.lowFreq, self.highFreq))
            
        if filterType == FILTER_FERMI:
            args.append(self.temperature.get())
        
        # Map to spected filter number in Spider for operation FQ    
        filterNumber = filterType * 2 + 1
        # Consider low-pass or high-pass
        filterNumber += self.filterMode.get()
        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())

        self._enterWorkingDir() # Do operations inside the run working dir
        
        spi = SpiderShell(ext=self.getExt()) # Create the Spider process to send commands        
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
        self._defineTransformRelation(particles, imgSet)
        
#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        pixelSize = self.inputParticles.get().getSamplingRate()
        
        summary = []
        summary.append('Used filter: *%s %s*' % (self.getEnumText('filterType'), self.getEnumText('filterMode')))
 
        if self.filterType <= FILTER_SPACE_REAL or self.filterType == FILTER_FERMI: 
            summary.append('Filter radius: *%s px^-1*' % self.filterRadius)
            radiusAngstroms = pixelSize/self.filterRadius
            summary.append('Filter radius: *%s Angstroms*' % radiusAngstroms )
        else:
            summary.append('Frequency range: *%s - %s*' % (self.lowFreq, self.highFreq))

        if self.filterType == FILTER_FERMI:
            summary.append('Temperature factor: *%s*' % self.temperature)

        summary.append('Padding set to: *%s*' % self.usePadding)
        return summary
    
    def _methods(self):
        methods = []
        msg = '\nIput particles %s were %s filtered using a %s filter' % (self.getObjectTag('inputParticles') ,
                                                                          self.getEnumText('filterMode'),
                                                                          self.getEnumText('filterType'))

        if self.filterType <= FILTER_SPACE_REAL or self.filterType == FILTER_FERMI: 
            msg += ', using a radius of %s px^-1' % self.filterRadius
        else:
            msg += ', using a frequency range of %s to %s px^-1' % (self.lowFreq, self.highFreq)

        if self.filterType == FILTER_FERMI: 
            msg += ' and a temperature factor of of %s px^-1' % self.temperature
            
        if self.usePadding:
            msg += ', padding the images by a factor of two.'
        else:
            msg += ' with no padding.'

        methods.append(msg)
        methods.append('Output particles: %s' % self.getObjectTag('outputParticles'))

        return methods
