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
from spider import SpiderShell
from convert import locationToSpider
from glob import glob

      
# TODO: Remove from ProtAlign, and put in other category     
class SpiderProtClassifyWard(ProtClassify):
    """ Ward's method, using 'CL HC' 
    """
    def __init__(self):
        ProtClassify.__init__(self)
        self._params = {'ext': 'stk',
                        'inputImage': 'input_image',
                        'outputMask': 'output_mask'
                        }
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Input images to perform PCA')
        
        #form.addParam('maskType', )
              
        form.addParam('imcFile', StringParam, default="",
                      label='IMC file generated in CA-PCA')        
        form.addParam('numberOfFactors', IntParam, default=10,
                      label='Number of factors',
                      help='After running, examine the eigenimages and decide which ones to use.\n'
                           'Typically all but the first few are noisy.')
        
        
    def _defineSteps(self):
        self._insertFunctionStep('classifyWard', self.imcFile.get(), self.numberOfFactors.get())
        #self._insertFunctionStep('createOutput')
#; classification, hierarchical
#cl hc
#[cas_prefix]_IMC  ; INPUT
#(1-x27)  ; factors to use
#(0)      ; no factor weighting
#(5)      ; clustering criterion (5==Ward's method)
#Y        ; dendrogram PostScript file?
#[ps_dendrogram]   ; OUTPUT
#Y        ; dendrogram document file?
#[dendrogram_doc]  ; OUTPUT
    
    def classifyWard(self, imcFile, numberOfFactors):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        self._params.update(locals()) # Store input params in dict
        
        
        # Copy file to working directory, it could be also a link
        imcLocalFile = basename(imcFile)
        copyFile(imcFile, self._getPath(imcLocalFile))
        imcLocalFile = removeExt(imcLocalFile)

        self._enterWorkingDir() # Do operations inside the run working dir


        spi = SpiderShell(ext=self._params['ext'], log='script.stk') # Create the Spider process to send commands 
        spi.runFunction('CL HC', imcLocalFile, '1-%d' % numberOfFactors, 0, 5, 
                        'Y', 'dendogram', 'Y', 'docdendro')
        spi.close()
        
        self._leaveWorkingDir() # Go back to project dir

            
    def _summary(self):
        summary = []
        return summary
    

