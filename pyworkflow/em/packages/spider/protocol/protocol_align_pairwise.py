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

import pyworkflow.utils as pwutils
from pyworkflow.em import IntParam

from ..spider import getScript
from protocol_align_base import SpiderProtAlign


      
class SpiderProtAlignPairwise(SpiderProtAlign):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    
    This alignment scheme aligns a pair of images at a time and then averages 
    them. Then, the averages of each of those pairs is aligned and averages, 
    and then pairs of those pairs, etc. Compared to [[http://spider.wadsworth.org/spider_doc/spider/docs/man/apsr.html][AP SR]], this alignment 
    scheme appears to be less random, which chooses seed images as alignment 
    references.
    
    For more information, see Step 2b at [[http://spider.wadsworth.org/spider_doc/spider/docs/techs/MSA/index.html#pairwise][SPIDER's MDA online manual]]
    """
    _label = 'align pairwise'
    
    def __init__(self, **args):
        SpiderProtAlign.__init__(self, 'mda/pairwise.msa', 'pairwise', **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineAlignParams(self, form):
        SpiderProtAlign._defineAlignParams(self, form)
        
        form.addParam('searchRange', IntParam, default=8, 
                      label='Search range (px):',
                      help='In the translational alignment, shifts of up to\n'
                           '_searchRange_ (in pixel units) will be allowed.')
        form.addParam('stepSize', IntParam, default=2, 
                      label='Step size (px):',
                      help='Alignments will be evaluated in units of _stepSize_ \n'
                           '(in pixel units) up to a maximum of +/- _searchRange_.')        
        form.addParallelSection(threads=2, mpi=0)    
    
    #--------------------------- STEPS functions --------------------------------------------       
    
    def alignParticlesStep(self, innerRadius, outerRadius):
        """ Execute the pairwise.msa script to alignm the particles.
        """
        particles = self.inputParticles.get()
        xdim = particles.getDimensions()[0]
        
        self._params.update({
                             '[idim-header]': xdim,
                             '[cg-option]': self.cgOption.get(),
                             '[inner-rad]': innerRadius,
                             '[outer-rad]': outerRadius, # convert radius to diameter
                             '[search-range]': self.searchRange.get(),
                             '[step-size]': self.stepSize.get(),
                             '[selection_list]': self._params['particlesSel'],
                             '[unaligned_image]': self._params['particles'] + '@******',
                             '[nummps]': self.numberOfThreads.get(),
                            })
        
        copy1Script = getScript('mda/center1.msa')
        newScript = pwutils.replaceBaseExt(copy1Script, self.getExt())
        pwutils.copyFile(copy1Script, self._getPath(newScript))
        self.runTemplate(self.getScript(), self.getExt(), self._params)
                
    def getAverage(self):
        return self._getPath(self.getAlignDir(), 'rfreeavg001.%s' % self.getExt()) 
       
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _citations(self):
        return ['Marco1996']
    
    def _summary(self):
        summary = []
        summary.append('Radius range (px): *%s - %s*' % (self.innerRadius, self.outerRadius))
        summary.append('Search range (px): *%s*' % self.searchRange)
        summary.append('Step size (px): *%s*' % self.stepSize)
        
        return summary
    
    def _methods(self):
        msg  = "Input particles %s " % self.getObjectTag('inputParticles')
        msg += "were subjected to a pairwise reference-free alignment using the "
        msg += "'pyramidal system for prealignment construction' ([Marco1996]), "
        msg += "using radii %s to %s pixels. " % (self.innerRadius, self.outerRadius)
        msg += "Particles were then aligned to this initial reference-free average "
        msg += "using SPIDER command _AP SH_ using a "
        msg += "search range of %s pixels and a step size of %s pixels. " % (self.searchRange, 
                                                                            self.stepSize)
        if self.hasAttribute('outputParticles'):
            msg += "Output particles: %s" % self.getObjectTag('outputParticles')
        
        return [msg]
