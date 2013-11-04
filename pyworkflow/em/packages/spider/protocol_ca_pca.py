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
from pyworkflow.em.packages.spider.spider import PcaFile
"""
This sub-package contains protocol for 
Correspondence Analysis or Principal Component Analysis
"""

from pyworkflow.em import *  
from pyworkflow.utils import removeExt, removeBaseExt, makePath, moveFile
from constants import *
from spider import SpiderShell, SpiderDocFile, SpiderProtocol, copyTemplate, runSpiderTemplate
from convert import locationToSpider
from glob import glob



class SpiderProtCAPCA(SpiderProtocol):
    """ Correspondence Analysis or Principal Component Analysis. """
    def __init__(self):
        SpiderProtocol.__init__(self)
        self._params = {'ext': 'stk',
                        'spiderParticles': 'particles_spider',
                        'spiderSel': 'particles_selfile',
                        'outputParticles': 'particles_output',
                        'spiderMask': 'mask',
                        'imcFile': join('CA', 'cas_IMC'),
                        'seqFile': join('CA', 'cas_SEQ'),
                        'eigFile': join('CA', 'cas_EIG'),
                        'eigenimages': join('CA', 'eigenimg'),
                        'reconstituted': join('CA', 'reconstituted')
                        }
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input particles to perform CA or PCA.')        
        form.addParam('analysisType', EnumParam, default=CA, choices=['CA', 'PCA', 'IPCA'],
                      label='Analysis type',
                      help='Select which type of analysis you want to perform')
        form.addParam('addConstant', FloatParam, default=0,
                      condition="analysisType==%d" % CA, 
                      label='Additive constant',
                      help='Additive constant, 0 means automatic.')       
        form.addParam('numberOfFactors', IntParam, default=25,
                      label='Number of eigenfactors',
                      help='Number of eigenfactors to calculate.')
        form.addParam('maskType', EnumParam, choices=['circular', 'from file'], default=0,
                      label='Mask type', help='Select which type of mask do you want to apply.')
        form.addParam('maskRadius', IntParam, default=-1,
                      label='Mask radius (pix)', condition='maskType==0',
                      help='If -1, the entire radius (in pixels) will be considered.')
        form.addParam('maskImage', PointerParam, label="Mask image", condition='maskType==1',
                      pointerClass='Mask', help="Select a mask file")       
        
    def _defineSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInput', 'inputParticles',
                                 self._getFileName('spiderParticles'), self._getFileName('spiderSel'))
        self._insertFunctionStep('convertMask')
        self._insertFunctionStep('runCAPCA', self.analysisType.get(), 
                                 self.numberOfFactors.get(), self.maskType.get())
        
    def convertMask(self):
        """ Convert the input mask if needed and
        copy some spider needed scripts. 
        """
        # Copy mask if selected
        if self.maskType > 0: # mask from file
            maskFn = self._getFileName('spiderMask')
            ih.convert(self.maskImage.get().getLocation(), 
                       (1, maskFn))
        # Copy template scripts
        copyTemplate('ploteigen.gnu', self._getPath())
        copyTemplate('eigendoc.py', self._getPath())
        
    def runCAPCA(self, analysisType, numberOfFactors, maskType):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        dim = self.inputParticles.get().getDimensions()[0]
        
        self._params.update({'dim': dim,
                             'radius': self.maskRadius.get(),
                             'analysisType': analysisType + 1, # Index starts at 0
                             'addConstant': self.addConstant.get(),
                             'numberOfFactors': numberOfFactors
                             })
                                
        self._enterWorkingDir() # Do operations inside the run working dir
        
#        spi = SpiderShell(ext=self._params['ext'], log='script.stk') # Create the Spider process to send commands 
#        
#        spi.runScript('ca_pca.txt', self._params)
#        
#        spi.close(end=False)   

        runSpiderTemplate('ca_pca.txt', self._params['ext'], self._params)
        
        self._leaveWorkingDir() # Go back to project dir
        
        # Generate outputs
        imc = PcaFile()
        imc.filename.set(self._getFileName('imcFile'))
        seq = PcaFile()
        seq.filename.set(self._getFileName('seqFile'))
        
        self._defineOutputs(imcFile=imc, seqFile=seq)
        
            
    def _summary(self):
        summary = []
        return summary
    

