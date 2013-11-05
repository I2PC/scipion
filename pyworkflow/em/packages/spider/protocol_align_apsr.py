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
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.utils import removeExt, removeBaseExt, makePath, getLastFile
from constants import *
from spider import SpiderShell, SpiderProtocol, runSpiderTemplate
from convert import locationToSpider
from glob import glob
        

      
class SpiderProtAlignAPSR(ProtAlign, SpiderProtocol):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    Uses Spider AP SR command.
    """
    def __init__(self, **args):
        ProtAlign.__init__(self, **args)
        SpiderProtocol.__init__(self, **args)
        
        self._params = {'ext': 'stk',
                        'particles': 'particles_aligned',
                        'particlesSel': 'particles_aligned_sel',
                        }    
    
    def _defineAlignParams(self, form):
        form.addParam('innerRadius', IntParam, default=5,
                      label='Inner radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <outerRadius> (in pixel units) will be analyzed.')
        form.addParam('outerRadius', IntParam, default=50, 
                      label='Outer radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <outerRadius> (in pixel units) will be analyzed.')
        
    def _defineSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('alignParticles', 
                                 self.innerRadius.get(), self.outerRadius.get())
        self._insertFunctionStep('createOutput')

    def alignParticles(self, innerRadius, outerRadius):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        xdim = particles.getDimensions()[0]
        
        self._enterWorkingDir() # Do operations inside the run working dir
        
        # Align images
        alignDir = 'align'
        makePath(alignDir)
        avgPattern = join(alignDir, 'avg_iter***')
        docPattern = join(alignDir, 'doc_iter***')

        self._params.update({
                             'innerRadius': self.innerRadius.get(),
                             'outerRadius': self.outerRadius.get(),
                             'numberOfParticles': n,
                             'dim': xdim,
                             'avgPattern': avgPattern,
                             'docPattern': docPattern,
                            })
        runSpiderTemplate('ap_sr.txt', 'stk', self._params)
                    
        
        spi = SpiderShell(ext='stk', log='script2.log') # Create the Spider process to send commands
        
        lastAvg = getLastFile(avgPattern)
        avg = Particle()
        avg.copyInfo(particles)
        avg.setLocation(NO_INDEX, self._getPath(lastAvg))
        self._defineOutputs(outputAverage=avg)
        
        lastDoc = getLastFile(docPattern)
        f = open(lastDoc) # Open last doc
        particlesPrefix = self._params['particles']
        
        for i, line in enumerate(f):
            line = line.strip()
            if len(line) and not line.startswith(';'):
                angle, shiftX, shiftY = [float(s) for s in line.split()[2:]]
                inLoc = locationToSpider(i, particlesPrefix)
                incore = "_1"
                spi.runFunction('RT SQ', inLoc, incore, angle, (shiftX, shiftY))
                spi.runFunction('CP', incore, inLoc)
            
        spi.close()
            
        self._leaveWorkingDir() # Go back to project dir
            
    def createOutput(self):
        """ Register the output (the alignment and the aligned particles.)
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        
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
        outerRadius = self.outerRadius.get()
        if innerRadius > r or innerRadius < 0:
            errors.append("<innerRadius> should be between 0 and %d" % r)        
        if outerRadius > r or outerRadius < 0:
            errors.append("<outerRadius> should be between 0 and %d" % r)
        
        return errors
    

