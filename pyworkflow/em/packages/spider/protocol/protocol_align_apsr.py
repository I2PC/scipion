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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.utils.path import getLastFile
from protocol_align_base import SpiderProtAlign


      
class SpiderProtAlignAPSR(SpiderProtAlign):
    """ 
    Reference-free alignment (both translational and rotational)
    of an image series. Uses Spider AP SR command.
    
    See detailed description at:
    [[http://spider.wadsworth.org/spider_doc/spider/docs/man/apsr.html][SPIDER's AP SR online manual]]
    """
    _label = 'align apsr'
    
    def __init__(self, **args):
        SpiderProtAlign.__init__(self, 'mda/apsr4class.msa', 'apsr', **args)

    def _defineAlignParams(self, form):
        SpiderProtAlign._defineAlignParams(self, form)

        # Hide the center of gravity option from the GUI since it is not
        # used in this particular alignment method
        cgOption = form.getParam('cgOption')
        cgOption.config(condition='False')

        form.addParallelSection(threads=1, mpi=0)
        
    def alignParticlesStep(self, innerRadius, outerRadius):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        self._params.update({
                             '[inner-rad]': innerRadius,
                             '[outer-rad]': outerRadius,
                             '[group_particles]': self._params['particlesSel'],
                             '[unaligned]': self._params['particles'] + '@******',
                             '[aligned_stack]': self._params['particlesAligned'],
                             '[nummps]': self.numberOfThreads.get()
                            })
        
        self.runTemplate(self.getScript(), self.getExt(), self._params)
                
    def getAverage(self):
        pattern = self._getPath(self.getAlignDir(), 'iteravg*.%s' % self.getExt())
        return getLastFile(pattern)
    
    def _summary(self):
        summary = []
        summary.append('Radius range: *%s - %s*' % (self.innerRadius, self.outerRadius))
        
        return summary
    
    def _methods(self):
        if hasattr(self, 'outputParticles'):
            msg  = "Input particles %s were " % self.getObjectTag('inputParticles')
            msg += "initially subjected to reference-free alignment using SPIDER's "
            msg += "_AP SR_ command, using radii %s to %s pixels. " % (self.innerRadius, 
                                                                       self.outerRadius)
            msg += "Output particles: %s" % self.getObjectTag('outputParticles') 
        else:
            msg = "Output not ready yet."
        
        return [msg]

    def _citations(self):
        return ['Penczek1992']
