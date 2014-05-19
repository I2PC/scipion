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

from pyworkflow.em import *  
from pyworkflow.utils import removeExt
import pyworkflow.utils.graph as graph

from ..spider import SpiderDocFile
from ..constants import *
from protocol_base import SpiderProtClassify

      

class SpiderProtClassifyKmeans(SpiderProtClassify):
    """ Diday's method, using 'CL CLA' 
    """
    _label = 'classify kmeans'
    
    def __init__(self, **kwargs):
        SpiderProtClassify.__init__(self, **kwargs)
        
        self._params = {'ext': 'stk',
                        'particles': 'input_particles',
                        'particlesSel': 'input_particles_sel',
                        'classDir': 'KM',
                        }        

    #--------------------------- DEFINE param functions --------------------------------------------  
     
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles',
                      help='Input images to perform PCA')
        form.addParam('pcaFile', PointerParam, pointerClass='PcaFile',
                      label="PCA file", 
                      help='IMC or SEQ file generated in CA-PCA')        
        form.addParam('numberOfFactors', IntParam, default=10,
                      label='Number of factors',
                      help='After running, examine the eigenimages and decide which ones to use.\n'
                           'Typically all but the first few are noisy.')
        form.addParam('numberOfClasses', IntParam, default=4, 
                      label='Number of classes',
                      help='Desired number of classes.')
        
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):    
        
        pcaFile = self.pcaFile.get().filename.get()
        
        self._insertFunctionStep('convertInput', 'inputParticles',
                                 self._getFileName('particles'), self._getFileName('particlesSel'))
        
        self._insertFunctionStep('classifyKmeansStep', pcaFile, 
                                 self.numberOfFactors.get(), self.numberOfClasses.get())
        
        self._insertFunctionStep('createOutputStep')
            
    #--------------------------- STEPS functions --------------------------------------------    
       
    def classifyKmeansStep(self, imcFile, numberOfFactors, numberOfClasses):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        # Copy file to working directory, it could be also a link
        imcLocalFile = basename(imcFile)
        copyFile(imcFile, self._getPath(imcLocalFile))
        self.info("Copied file '%s' to '%s' " % (imcFile, imcLocalFile))
        # Spider automatically add _IMC to the ca-pca result file
        imcBase = removeExt(imcLocalFile).replace('_IMC', '')
        
        self._params.update({'x20': numberOfClasses,
                             'x27': numberOfFactors,
                             '[cas_prefix]': imcBase,
                             '[particles]': self._params['particles'] + '@******',
                             '[class_dir]': self._params['classDir'],
                             })

        self.runScript('mda/kmeans.msa', self._params['ext'], self._params)

    def createOutputStep(self):
        """ Create the SetOfClass from the docfiles with the images-class
        assigment, the averages for each class.
        """
        particles = self.inputParticles.get()
        sampling = particles.getSamplingRate()
        classes2D = self._createSetOfClasses2D(particles)
            
        for classId in range(1, self.numberOfClasses.get()+1):
            class2D = Class2D()
            class2D.setObjId(classId)
            
            avgImg = Particle()
            avgImg.setSamplingRate(sampling)
            avgFn = self._getPath(self._params['classDir'], 'classavg%03d.stk' % classId)
            avgImg.setLocation(1, avgFn)
            
            class2D.setRepresentative(avgImg)
            classes2D.append(class2D)
            
            docClass = self._getPath(self._params['classDir'], 
                                     'docclass%03d.stk' % classId)
            doc = SpiderDocFile(docClass)
            
            for values in doc.iterValues():
                imgId = int(values[0])
                img = particles[imgId]
                class2D.append(img)
                
        self._defineOutputs(outputClasses=classes2D)
         
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
    
    #--------------------------- UTILS functions --------------------------------------------
    
    
