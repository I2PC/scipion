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
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.utils import removeExt, removeBaseExt, makePath, getLastFile
from constants import *
from spider import SpiderShell, runSpiderTemplate
from convert import locationToSpider
from glob import glob
        

PART_INPUT = 'particles_input'
PART_BIGEN = 'particles_big' 
PART_OUTPUT = 'particles_output'

      
class SpiderProtAlignAPSR(ProtAlign):
    """ Reference-free alignment shift and rotational alignment of an image series. 
    Uses Spider AP SR command.
    """
    def _defineAlignParams(self, form):
        form.addParam('innerRadius', IntParam, default=5,
                      label='Inner radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <outerRadius> (in pixel units) will be analyzed.')
        form.addParam('outerRadius', IntParam, default=50, 
                      label='Outer radius(px):',
                      help='In the rotational alignment, only rings between\n'
                           '<innerRadius> and <outerRadius> (in pixel units) will be analyzed.')
        
    
    def _defineSteps(self):
        self._op = "FQ"
        # Define some names
        self.inputStk = self._getPath('%s.stk' % PART_INPUT)
        self.bigenStk = self._getPath('%s.stk' % PART_BIGEN)
        self.outputStk = self._getPath('%s.stk' % PART_OUTPUT)
        # Insert processing steps
        self._insertFunctionStep('convertInput')
        self._insertFunctionStep('alignParticles', 
                                 self.innerRadius.get(), self.outerRadius.get())
        self._insertFunctionStep('createOutput')
        #self._insertFunctionStep('createOutput')
        
    def convertInput(self):
        """ Convert the input particles to a Spider stack. """
        particles = self.inputParticles.get()
        ih = ImageHandler()
        inputStk = self._getPath('%s.stk' % PART_INPUT)
        
        for i, p in enumerate(particles):
            ih.convert(p.getLocation(), (i+1, self.inputStk))

    def alignParticles(self, innerRadius, outerRadius):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        
        self._enterWorkingDir() # Do operations inside the run working dir
        
#        spi = SpiderShell(ext='stk', log='script.log') # Create the Spider process to send commands        
#        
        # Create the images docfile
        #selfile = 'input_particles_sel'
        #spi.runFunction('DOC CREATE', selfile, 0, '1-%d' % n)
        particles = self.inputParticles.get()
        n = particles.getSize()
        
#        # Generate a filtered disc
#        spi.runFunction('PT', '_1', (100,100), 'C', (50,50), 45, 'N')
#        spi.runFunction('FQ', '_1', '_2', 3, 0.02)
#        # Ensure the stack has Big endian as expected by SPIDER
#        for i in range(1, n+1):
#            spi.runFunction('CP', '%s@%d' % (PART_INPUT, i), '%s@%d' % (PART_BIGEN, i))
#        spi.close() # This is done to ensure that AP SR is done
                    # we need a better way to wait on completion
        # Align images
        alignDir = 'align'
        makePath(alignDir)
        avgPattern = join(alignDir, 'avg_iter***')
        docPattern = join(alignDir, 'doc_iter***')

#        spi.runFunction('AP SR', 
#                        '%s@******' % PART_BIGEN, '1-%d' % n, 
#                        90, (1,45), '_2', avgPattern, docPattern)
        params = {'innerRadius': self.innerRadius.get(),
                  'outerRadius': self.outerRadius.get(),
                  'numberOfParticles': n,
                  'dim': 100,
                  'particles': PART_INPUT,
                  'particlesBig': PART_BIGEN,
                  'avgPattern': avgPattern,
                  'docPattern': docPattern,
                  }
        runSpiderTemplate('ap_sr.txt', 'stk', params)
        
                    
        
        spi = SpiderShell(ext='stk', log='script2.log') # Create the Spider process to send commands
        # TODO: Remove this, only testing for creating Xmipp metadata
#        import xmipp
#        from protlib_import import readMdFromSpider
#        
        lastAvg = getLastFile(avgPattern)
        avg = Particle()
        avg.copyInfo(particles)
        avg.setLocation(NO_INDEX, self._getPath(lastAvg))
        self._defineOutputs(outputAverage=avg)
        
#        md = xmipp.MetaData()
#        for avg in sorted(avgFiles):
#            objId = md.addObject()
#            md.setValue(xmipp.MDL_IMAGE, self._getPath(avg), objId)
#        md.write('averages.xmd')
        
        #md.read('%s.stk' % PART_OUTPUT)
        lastDoc = getLastFile(docPattern)
        f = open(lastDoc) # Open last doc
        
        for i, line in enumerate(f):
            line = line.strip()
            if len(line) and not line.startswith(';'):
                angle, shiftX, shiftY = [float(s) for s in line.split()[2:]]
                inLoc = locationToSpider(i, PART_BIGEN)
                outLoc = locationToSpider(i, PART_OUTPUT)
                spi.runFunction('RT SQ', inLoc, outLoc, (angle, 1), (shiftX, shiftY))
            
        spi.close()
#        for doc in docFiles:
#            mdDoc = readMdFromSpider(doc, 'anglePsi shiftX shiftY')
#            md.merge(mdDoc)
#            md.write(doc.replace('.stk', '.xmd'))
            
        self._leaveWorkingDir() # Go back to project dir
            
    def createOutput(self):
        """ Register the output (the alignment and the aligned particles.)
        """
        particles = self.inputParticles.get()
        n = particles.getSize()
        
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(particles)
        
        outputStk = self._getPath('%s.stk' % PART_OUTPUT)
        for i in range(1, n+1):
            img = Image()
            img.setLocation(i, outputStk)
            imgSet.append(img)
        imgSet.write()
        self._defineOutputs(outputParticles=imgSet)
        
    def _summary(self):
        summary = []
        return summary
    

