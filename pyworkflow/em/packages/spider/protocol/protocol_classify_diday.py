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

      

class SpiderProtClassifyDiday(SpiderProtClassifyCluster):
    """ Performs automatic clustering using Diday's method and
    Hierarchical Ascendant Classification (HAC) using Ward's criterion
    on factors produced by CA or PCA. Uses the Spider CL CLA program
    """
    _label = 'classify diday'
    
    def __init__(self, **kwargs):
        SpiderProtClassifyCluster.__init__(self, 'mda/cluster.msa', 'CLA',  **kwargs)

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
        msg += "were subjected to Diday's method of moving centers "
        msg += "(SPIDER command [[http://spider.wadsworth.org/spider_doc/spider/docs/man/clcla.html][CL CLA]]) "
        msg += "using %s factors. " % self.numberOfFactors
        return [msg]
#        # I don't know what's he's talking about.  The summary is totally empty as of now.
#        return self._summary()  # summary is quite explicit and serve as methods
    
