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
from spider import SpiderShell, SpiderDocFile, SpiderProtocol
from convert import locationToSpider
from glob import glob

      
class SpiderWfMDA(ProtClassify, SpiderProtocol):
    """ Ward's method, using 'CL HC' 
    """
    def __init__(self, **args):
        ProtClassify.__init__(self, **args)
        EMProtocol.__init__(self, **args)
        
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
        #form.addSection(label='1.Filter')      
        form.addParam('filter', ProtocolClassParam, 
                      protocolClassName='SpiderProtFilter',
                      label="1.Filter protocol", 
                      help='Select which Filter Protocol do you want to use')      
        
        #form.addSection(label='2.Align')      
        form.addParam('align', ProtocolClassParam, allowSubclasses=True,
                      protocolClassName='ProtAlign',
                      label="2.Align protocol", 
                      help='Select which Filter Protocol do you want to use')  
           
        #form.addSection(label='3.Dimension reduction')      
        form.addParam('dimred', ProtocolClassParam, 
                      protocolClassName='SpiderProtCAPCA',
                      label="3. Dimension reduction", 
                      help='Select which Filter Protocol do you want to use')      
        
        #form.addSection(label='4.Classification')      
        form.addParam('classify', ProtocolClassParam, 
                      protocolClassName='SpiderProtClassifyWard',
                      label="4. Classification protocol", 
                      help='Select which Filter Protocol do you want to use')         

    
    def _defineSteps(self):
        self._insertFunctionStep('workflowStep', self.inputParticles.get().getFileName())

            
    def workflowStep(self, inputFileName):
        self.filterInstance.inputParticles.set(self.inputParticles.get())
        self.runProtocol(self.filterInstance)
        
        self.alignInstance.inputParticles.set(self.filterInstance.outputParticles)
        self.runProtocol(self.alignInstance)
         
#        protMask = getattr(self, '')
#        protMask.inputImage.set(protAPSR.outputAverage)
#        self.runProtocol(protMask)       
              
#        protCAPCA = SpiderProtCAPCA()
#        protCAPCA.maskType.set(1)
#        protCAPCA.maskImage.set(protMask.outputMask)
        self.dimredInstance.inputParticles.set(self.alignInstance.outputParticles)
        self.runProtocol(self.dimredInstance)
        
        self.classifyInstance.inputParticles.set(self.alignInstance.outputParticles)
        self.classifyInstance.pcaFilePointer.set(self.dimredInstance.imcFile)
        self.runProtocol(self.classifyInstance)
        
        # TODO: Fix this, since is a deeper problem with pointers
        #self.alignInstance.outputParticles.setStore(True)
        #self.dimred.imcFile.setStore(True)
        
    
            
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        return []
    

