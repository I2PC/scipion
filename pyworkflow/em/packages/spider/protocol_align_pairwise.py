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
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.utils import removeExt, removeBaseExt, makePath, getLastFile
from constants import *
from spider import SpiderShell, SpiderProtocol, runSpiderTemplate
from convert import locationToSpider
from glob import glob
        

      
class SpiderProtAlignPairwise(ProtAlign, SpiderProtocol):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    Described in Marco S, Chagoyen M, de la Fraga LG, Carazo JM, Carrascosa JL (1996)
    "A variant to the Random Approximation of the reference-free algorithm."
    Ultramicroscopy. Vol 66: pg. 5-10.
    """
    _label = 'align pairwise'
    
    def __init__(self, **args):
        ProtAlign.__init__(self, **args)
        SpiderProtocol.__init__(self, **args)
        self.alignDir = 'pairwise'
         
        self._params = {'ext': 'stk',
                        'particles': 'particles_aligned',
                        'particlesSel': 'particles_aligned_sel',
                        'average': join(self.alignDir, 'rfreeavg001') # If change the name, change it in template batch
                        }    
    
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
    def _defineSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('alignParticles', 
                                 self.innerRadius.get(), self.diameter.get())
        self._insertFunctionStep('createOutput')

    def alignParticles(self, innerRadius, diameter):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        xdim = particles.getDimensions()[0]
        
        self._enterWorkingDir() # Do operations inside the run working dir
        
        # Align images
        makePath(self.alignDir)

        self._params.update({
                             'dim': xdim,
                             'diameter': self.diameter.get(),
                             'innerRadius': self.innerRadius.get(),
                             'searchRange': self.searchRange.get(),
                             'stepSize': self.stepSize.get()
                            })
        runSpiderTemplate('pairwise.txt', self._params['ext'], self._params)
                    
        
#        spi = SpiderShell(ext='stk', log='script2.log') # Create the Spider process to send commands
#        
#        lastAvg = getLastFile(avgPattern)
#        avg = Particle()
#        avg.copyInfo(particles)
#        avg.setLocation(NO_INDEX, self._getPath(lastAvg))
#        self._defineOutputs(outputAverage=avg)
#        
#        lastDoc = getLastFile(docPattern)
#        f = open(lastDoc) # Open last doc
#        particlesPrefix = self._params['particles']
#        
#        for i, line in enumerate(f):
#            line = line.strip()
#            if len(line) and not line.startswith(';'):
#                angle, shiftX, shiftY = [float(s) for s in line.split()[2:]]
#                inLoc = locationToSpider(i, particlesPrefix)
#                incore = "_1"
#                spi.runFunction('RT SQ', inLoc, incore, angle, (shiftX, shiftY))
#                spi.runFunction('CP', incore, inLoc)
#            
#        spi.close()
            
        self._leaveWorkingDir() # Go back to project dir
            
    def createOutput(self):
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
        imgSet.write()
        self._defineOutputs(outputParticles=imgSet)
        
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        particles = self.inputParticles.get()
        xdim = particles.getDimensions()[0]
        r = xdim / 2
        innerRadius = self.innerRadius.get()
        if innerRadius > r or innerRadius < 0:
            errors.append("<innerRadius> should be between 0 and %d" % r)
        
        return errors
    

