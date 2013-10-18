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
This sub-package contains protocol for particles filters operations
"""


from pyworkflow.em import *  
from pyworkflow.utils import removeExt, removeBaseExt, makePath
from constants import *
from spider import SpiderShell
from convert import locationToSpider
from glob import glob
        
class SpiderDefAlignAPSR(Form):        
    """Create the definition of parameters for
    the XmippProtAlign2d protocol"""
    def __init__(self):
        Form.__init__(self)

        self.addSection(label='Input')
        self.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select the input particles to be aligned.')        
        self.addParam('innerRadius', IntParam, default=5,
                      label='Inner radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <outerRadius> (in pixel units) will be analyzed.')
        self.addParam('outerRadius', IntParam, default=50, 
                      label='Outer radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <outerRadius> (in pixel units) will be analyzed.')

      
class SpiderProtAlignAPSR(ProtAlign):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    Uses Spider AP SR command.
    """
    _definition = SpiderDefAlignAPSR()
    
    def __init__(self):
        ProtAlign.__init__(self)
        self._op = "FQ"
    
    def _defineSteps(self):
        # Define some names
        self.inputStk = self._getPath('input_particles.stk')
        self.outputStk = self._getPath('output_particles.stk')
        # Insert processing steps
        self._insertFunctionStep('convertInput')
        self._insertFunctionStep('alignParticles', 
                                 self.innerRadius.get(), self.outerRadius.get())
        #self._insertFunctionStep('createOutput')
        
    def convertInput(self):
        """ Convert the input particles to a Spider stack. """
        particles = self.inputParticles.get()
        ih = ImageHandler()
        
        for i, p in enumerate(particles):
            ih.convert(p.getLocation(), (i+1, self.inputStk))

    def alignParticles(self, innerRadius, outerRadius):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        imgSet.load() # Force to create mapper before entering working dir

        self._enterWorkingDir() # Do operations inside the run working dir
        
        spi = SpiderShell(ext='stk', log='script.log') # Create the Spider process to send commands        
        #inputStk = removeBaseExt(self.inputStk)
        outputStk = removeBaseExt(self.outputStk)
        inputStk = outputStk
        
        # Create the images docfile
        selfile = 'input_particles_sel'
        spi.runFunction('DOC CREATE', selfile, 0, '1-%d' % n)
        # Generate a filtered disc
        spi.runFunction('PT', '_1', (100,100), 'C', (50,50), 45, 'N')
        spi.runFunction('FQ', '_1', '_2', 3, 0.02)
        # Ensure the stack has Big endian as expected by SPIDER
        for i in range(1, n+1):
            spi.runFunction('CP', 'input_particles@%d' % i, 'output_particles@%d' % i)
        # Align images
        alignDir = 'align'
        makePath(alignDir)
        avgPattern = join(alignDir, 'avg_iter***')
        docPattern = join(alignDir, 'doc_iter***')
        spi.runFunction('AP SR', 'output_particles@******', '1-%d' % n, 
                        90, (1,45), '_2', avgPattern, docPattern)
        
        spi.close()
        
        
        import xmipp
        from protlib_import import readMdFromSpider
        
        avgFiles = glob(avgPattern)
        md = xmipp.MetaData()
        for avg in sorted(avgFiles):
            objId = md.addObject()
            md.setValue(xmipp.MDL_IMAGE, self._getPath(avg), objId)
        md.write('averages.xmd')
        
        md.read('output_particles.stk')
        docFiles = glob(docPattern)
        for doc in docFiles:
            mdDoc = readMdFromSpider(doc, 'anglePsi shiftX shiftY')
            md.merge(mdDoc)
            md.write(doc.replace('.stk', '.xmd'))
            
            
        
            
        
        self._leaveWorkingDir() # Go back to project dir
            
#        for i in range(1, n+1):
#            img = Image()
#            img.setLocation(i, self.outputStk)
#            imgSet.append(img)
#        imgSet.write()
#        self._defineOutputs(outputParticles=imgSet)
        
    def _summary(self):
        summary = []
        return summary
    

