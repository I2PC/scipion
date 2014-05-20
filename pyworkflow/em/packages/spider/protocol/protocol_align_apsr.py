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

from pyworkflow.utils.path import getLastFile

from protocol_align_base import SpiderProtAlign


      
class SpiderProtAlignAPSR(SpiderProtAlign):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    Uses Spider AP SR command.
    """
    _label = 'align apsr'
    
    def __init__(self, **args):
        SpiderProtAlign.__init__(self, 'mda/apsr4class.msa', 'apsr', **args)
        
    def alignParticlesStep(self, innerRadius, outerRadius):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        self._params.update({
                             '[x55]': innerRadius,
                             '[x60]': outerRadius * 2, # convert radius to diameter
                             '[group_particles]': self._params['particlesSel'],
                             '[unaligned]': self._params['particles'] + '@******',
                             '[aligned_stack]': self._params['particlesAligned'],
                            })
        
        self.runScript(self.getScript(), self.getExt(), self._params)
                
    def getAverage(self):
        pattern = self._getPath(self.getAlignDir(), 'iteravg*.%s' % self.getExt())
        return getLastFile(pattern)
    
    def _summary(self):
        summary = []
        return summary
    
