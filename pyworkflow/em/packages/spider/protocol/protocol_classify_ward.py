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

from protocol_classify_base import SpiderProtClassifyCluster
      
      

class SpiderProtClassifyWard(SpiderProtClassifyCluster):
    """ Finds clusters of images/elements in factor space
    (or a selected subspace) by using Diday’s method of moving centers,
    and applies hierarchical ascendant classification (HAC) (using
    Ward's method) to the resulting cluster centers.
    Uses the Spider CL HC program
    """
    _label = 'classify ward'
    
    def __init__(self, **kwargs):
        SpiderProtClassifyCluster.__init__(self, 'mda/hierarchical.msa', 'HC', **kwargs)
        
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        summary.append('Number of factors: *%s*' % self.numberOfFactors)
        return summary
    
    def _methods(self):
        msg  = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were subjected to Ward's method  "
        msg += "(SPIDER command [[http://spider.wadsworth.org/spider_doc/spider/docs/man/clhc.html][CL HC]]) "
        msg += "using %s factors. " % self.numberOfFactors
        return [msg]
#        return self._summary()  # summary is quite explicit and serve as methods
    
