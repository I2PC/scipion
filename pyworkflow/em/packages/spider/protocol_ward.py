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
                        'particles': 'particles_input',
                        'dendroPs': 'dendrogram',
                        'dendroDoc': 'docdendro',
                        'averages': 'averages'
                        }
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Input images to perform PCA')
        
        #form.addParam('maskType', )
              
        form.addParam('imcFile', StringParam, default="",
                      label='IMC file generated in CA-PCA')        
        form.addParam('numberOfFactors', IntParam, default=10,
                      label='Number of factors',
                      help='After running, examine the eigenimages and decide which ones to use.\n'
                           'Typically all but the first few are noisy.')
        
        
    def _getFileName(self, key):
        #TODO: Move to a base Spider protocol
        template = '%(' + key + ')s.%(ext)s'
        return self._getPath(template % self._params)
    
    def _defineSteps(self):
        self._insertFunctionStep('convertInput', self.inputParticles.getFileName())
        self._insertFunctionStep('classifyWard', self.imcFile.get(), self.numberOfFactors.get())
        self._insertFunctionStep('buildDendrogram', True)
    
    def convertInput(self, inputFilename):
        """ Convert the input particles to a Spider stack. """
        particles = self.inputParticles.get()
        ih = ImageHandler()
        
        for i, p in enumerate(particles):
            ih.convert(p.getLocation(), (i+1, self.particlesStk))
            
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
                        'Y', self._params['dendroPs'], 'Y', self._params['dendroDoc'])
        spi.close()
        
        self._leaveWorkingDir() # Go back to project dir
        
        
    def buildDendrogram(self, writeAverages=False):
        """ Parse Spider docfile with the information to build the dendogram.
        Params:
            dendroFile: docfile with a row per image. 
                 Each row contains the image id and the height.
        """ 
        dendroFile = self._getFileName('dendroDoc')
        f = open(dendroFile)
        values = []
        for line in f:
            line = line.strip()
            if not line.startswith(';'):
                values.append(float(line.split()[3]))
        f.close()
        self.dendroValues = values
        #TODO: COPY the input particles
        #self.dendroImages = self.inputParticles.
        
        return self._buildDendrogram(0, len(values)-1, 1, writeAverages)
    
    def _buildDendrogram(self, leftIndex, rightIndex, index, writeAverages=False):
        """ This function is recursively called to create the dendogram graph(binary tree)
        and also to write the average image files.
        Params:
            leftIndex, rightIndex: the indinxes within the list where to search.
            index: the index of the class average.
            writeImages: flag to select when to write averages.
        From self:
            self.dendroValues: the list with the heights of each node
            self.dendroImages: image stack filename to read particles
            self.dendroAverages: stack name where to write averages
        It will search for the max in values list (between minIndex and maxIndex).
        Nodes to the left of the max are left childs and the other right childs.
        """
        maxValue = self.dendroValues[leftIndex]
        maxIndex = 0
        for i, v in enumerate(self.dendroValues[leftIndex+1:rightIndex]):
            if v > maxValue:
                maxValue = v
                maxIndex = i+1
        
        m = maxIndex + leftIndex
        node = {'height': maxValue, 'childs': [], 
                'length': 1, 'index': index}#len(self.dendroValues[leftIndex:rightIndex])}
        
        ih = ImageHandler()

        if writeAverages:
            node['image'] = ih.read((m+1, self.dendroImages))
            
        def addChildNode(left, right, index):
            if right > left:
                child = self._buildDendrogram(left, right, index, writeAverages)
                node['childs'].append(child)
                node['length'] += child['length'] 
                if writeAverages:
                    node['image'] += child['image']
                    del child['image']
                
        if rightIndex > leftIndex + 1:
            addChildNode(leftIndex, m, 2*index)
            addChildNode(m+1, rightIndex, 2*index+1)
            if writeAverages:
                #TODO: node['image'] /= float(node['length'])
                ih.write(node['image'], (index, self.dendroAverages))
        return node

            
    def _summary(self):
        summary = []
        return summary
    

