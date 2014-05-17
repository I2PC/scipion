# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
#                Tapu Shaikh            (shaikh@ceitec.muni.cz)
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
from spider import SpiderProtocol

      
class SpiderProtAlignPairwise(ProtAlign2D, SpiderProtocol):
    """ Reference-free alignment shift and rotational alignment of an image series. """
    _label = 'align pairwise'
    
    def __init__(self, **args):
        ProtAlign2D.__init__(self, **args)
        SpiderProtocol.__init__(self, **args)
        self.alignDir = 'pairwise'
         
        self._params = {'ext': 'stk',
                        'particles': 'particles_aligned',
                        'particlesSel': 'particles_aligned_sel',
                        'average': join(self.alignDir, 'rfreeavg001') # If change the name, change it in template batch
                        }    
    
    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineAlignParams(self, form):
        form.addParam('innerRadius', IntParam, default=5,
                      label='Inner radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <diam>/2 (in pixel units) will be analyzed.')
        form.addParam('diameter', IntParam, default=88, 
                      label='Outer diameter(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <diam>/2 (in pixel units) will be analyzed.')
        form.addParam('searchRange', IntParam, default=8, 
                      label='Search range(px):',
                      help='In the translational alignment, shifts of up to\n'
                           '<searchRange> will be allowed.')
        form.addParam('stepSize', IntParam, default=2, 
                      label='Step size(px):',
                      help='In the translational alignment, shifts will be analyzed\n'
                           'in units of <stepSize> (in pixel units).')        
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('alignParticlesStep', 
                                 self.innerRadius.get(), self.diameter.get())
        #self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------       
    
    def alignParticlesStep(self, innerRadius, diameter):
        """ Execute the pairwise.msa script to alignm the particles.
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        xdim = particles.getDimensions()[0]
        
        self._enterWorkingDir() # Do operations inside the run working dir
        
        # Align images
        makePath(self.alignDir)

        self._params.update({
                             '[idim-header]': xdim,
                             '[inner-rad]': self.innerRadius.get(),
                             '[obj-diam]': self.diameter.get(),
                             '[search-range]': self.searchRange.get(),
                             '[step-size]': self.stepSize.get(),
                             '[selection_list]': 'particles_aligned_sel',
                             '[unaligned_image]': 'particles_aligned@*****',
                            })
        
        self.runScript('mda/pairwise.msa', self._params['ext'], self._params)
                    
        self._leaveWorkingDir() # Go back to project dir
            
    def createOutputStep(self):
        """ Register the output (the alignment and the aligned particles.)
        """
        particles = self.inputParticles.get()
        n = particles.getSize()

        # Create the output average image
        avg = Particle()
        avg.copyInfo(particles)
        avg.setLocation(NO_INDEX, self._getFileName('average'))
        self._defineOutputs(outputAverage=avg)
        
        # Create the output aligned images        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(particles)
        
        outputStk = self._getFileName('particles')
        for i in range(1, n+1):
            img = Image()
            img.setLocation(i, outputStk)
            imgSet.append(img)
        self._defineOutputs(outputParticles=imgSet)
        
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        particles = self.inputParticles.get()
        xdim = particles.getDimensions()[0]
        r = xdim / 2
        innerRadius = self.innerRadius.get()
        if innerRadius > r or innerRadius < 0:
            errors.append("<innerRadius> should be between 0 and %d" % r)
        
        return errors
    
    def _citations(self):
        return ['Marco1996']
    
    def _summary(self):
        summary = []
        return summary
    

