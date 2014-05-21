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
This sub-package contains protocol for particles filters operations
"""
from os.path import join

from pyworkflow.em import ProtAlign2D, IntParam, Image, Particle, NO_INDEX
from pyworkflow.utils import getLastFile, makePath

from ..constants import *
from ..spider import SpiderShell, runSpiderTemplate
from ..convert import locationToSpider
from protocol_base import SpiderProtocol


      
class SpiderProtAlign(ProtAlign2D, SpiderProtocol):
    """ Base class for Spider alignment protocols. """
    
    def __init__(self, script, alignDir, **args):
        ProtAlign2D.__init__(self, **args)
        SpiderProtocol.__init__(self, **args)
        self._script = script
        self._alignDir = alignDir
        
        self._params = {'ext': 'stk',
                        'particles': 'input_particles',
                        'particlesSel': 'input_particles_sel',
                        'average': join(self._alignDir, 'rfreeavg001'), # If change the name, change it in template batch
                        'particlesAligned': join(self._alignDir, 'stkaligned')
                        }  
        
    def getAlignDir(self):
        return self._alignDir 

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineAlignParams(self, form):
        line = form.addLine('Radius (px):', 
                            help='In the rotational alignment, only rings between\n'
                                 '_innerRadius_ and _outerRadius_ (in pixel units)\n'
                                 'will be analyzed.')
        line.addParam('innerRadius', IntParam, default=5, label='Inner')
        line.addParam('outerRadius', IntParam, default=44, label='Outer')

        
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles', 
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        self._insertFunctionStep('alignParticlesStep', 
                                 self.innerRadius.get(), self.outerRadius.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------       

    def createOutputStep(self):
        """ Register the output (the alignment and the aligned particles.)
        """
        # Create the output average image
        avg = Particle()
        avg.copyInfo(self.inputParticles.get())
        avg.setLocation(NO_INDEX, self.getAverage())
        self._defineOutputs(outputAverage=avg)
        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        outputStk = self._getFileName('particlesAligned')
        imgSet.readStack(outputStk)
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
        if innerRadius > r or innerRadius < 1:
            errors.append("*innerRadius* should be between 1 and %d." % r)        
        if outerRadius > r or outerRadius < 1:
            errors.append("*outerRadius* should be between 1 and %d." % r)
        if innerRadius >= outerRadius:
            errors.append("*innerRadius* should be less than *outerRadius*.")
        
        return errors
    

