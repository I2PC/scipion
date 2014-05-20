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

from pyworkflow.em import IntParam

from protocol_align_base import SpiderProtAlign


      
class SpiderProtAlignPairwise(SpiderProtAlign):
    """ Reference-free alignment shift and rotational alignment of an image series. """
    _label = 'align pairwise'
    
    def __init__(self, **args):
        SpiderProtAlign.__init__(self, 'mda/pairwise.msa', 'pairwise', **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineAlignParams(self, form):
        SpiderProtAlign._defineAlignParams(self, form)
        
        form.addParam('searchRange', IntParam, default=8, 
                      label='Search range (px):',
                      help='In the translational alignment, shifts of up to\n'
                           '_searchRange_ will be allowed.')
        form.addParam('stepSize', IntParam, default=2, 
                      label='Step size(px):',
                      help='In the translational alignment, shifts will be analyzed\n'
                           'in units of _stepSize_ (in pixel units).')        
    
    #--------------------------- STEPS functions --------------------------------------------       
    
    def alignParticlesStep(self, innerRadius, outerRadius):
        """ Execute the pairwise.msa script to alignm the particles.
        """
        particles = self.inputParticles.get()
        xdim = particles.getDimensions()[0]
        
        self._params.update({
                             '[idim-header]': xdim,
                             '[inner-rad]': innerRadius,
                             '[obj-diam]': outerRadius * 2, # convert radius to diameter
                             '[search-range]': self.searchRange.get(),
                             '[step-size]': self.stepSize.get(),
                             '[selection_list]': self._params['particlesSel'],
                             '[unaligned_image]': self._params['particles'] + '@******',
                            })
        
        self.runScript(self.getScript(), self.getExt(), self._params)
                
    def getAverage(self):
        return self._getPath(self.getAlignDir(), 'rfreeavg001.%s' % self.getExt()) 
       
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _citations(self):
        return ['Marco1996']
    
    def _summary(self):
        summary = []
        return summary
    

