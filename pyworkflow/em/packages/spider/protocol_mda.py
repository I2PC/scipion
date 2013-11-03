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
This sub-package contains Spider protocol for PCA.
"""


from pyworkflow.em import *  
from pyworkflow.utils import removeExt, removeBaseExt, makePath, moveFile, copyFile, basename
from constants import *
from spider import SpiderShell, SpiderDocFile
from convert import locationToSpider
from glob import glob

      
# TODO: Remove from ProtAlign, and put in other category     
class SpiderWfMDA(ProtClassify):
    """ Ward's method, using 'CL HC' 
    """
    def __init__(self):
        ProtClassify.__init__(self)
        self._params = {'ext': 'stk',
                        'particles': 'particles_input',
                        'dendroPs': 'dendrogram',
                        'dendroDoc': 'docdendro',
                        'averages': 'averages'
                        }
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Input images to perform MDA')
        
        #form.addParam('maskType', )
              
        form.addParam('filterProtocolClass', ProtocolClassParam, 
                      protocolClassName='ProtFilterParticles',
                      label="Filter protocol", 
                      help='Select which Filter Protocol do you want to use')        
        
        
    def _getFileName(self, key):
        #TODO: Move to a base Spider protocol
        template = '%(' + key + ')s.%(ext)s'
        return self._getPath(template % self._params)
    
    def _defineSteps(self):
        self._insertFunctionStep('convertInput', self.inputParticles.get().getFileName())
    
    def convertInput(self, inputFilename):
        """ Convert the input particles to a Spider stack. """
        particles = self.inputParticles.get()
        ih = ImageHandler()
        particlesStk = self._getFileName('particles')
        
        for i, p in enumerate(particles):
            ih.convert(p.getLocation(), (i+1, particlesStk))
            
    def _summary(self):
        summary = []
        return summary
    

