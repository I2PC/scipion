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

from pyworkflow.em import Particle, Class2D
from pyworkflow.protocol.params import IntParam

from ..spider import SpiderDocFile
from protocol_classify_base import SpiderProtClassify

      

class SpiderProtClassifyKmeans(SpiderProtClassify):
    """ Diday's method, using 'CL CLA' 
    """
    _label = 'classify kmeans'
    
    def __init__(self, **kwargs):
        SpiderProtClassify.__init__(self, 'mda/kmeans.msa', 'KM', **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------  
     
    def _defineParams(self, form):
        SpiderProtClassify._defineParams(self, form)

        form.addParam('numberOfClasses', IntParam, default=4, 
                      label='Number of classes',
                      help='Desired number of classes.')
        
    def getNumberOfClasses(self):
        return self.numberOfClasses.get()
            
    #--------------------------- STEPS functions --------------------------------------------    
       
    def _updateParams(self):
        self._params.update({'x20': self.getNumberOfClasses(),
                             '[particles]': self._params['particles'] + '@******',
                             })

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
            avgFn = self._getPath(self.getClassDir(), 'classavg%03d.stk' % classId)
            avgImg.setLocation(1, avgFn)
            
            class2D.setRepresentative(avgImg)
            classes2D.append(class2D)
            
            docClass = self._getPath(self.getClassDir(), 
                                     'docclass%03d.stk' % classId)
            doc = SpiderDocFile(docClass)
            
            for values in doc.iterValues():
                imgId = int(values[0])
                img = particles[imgId]
                class2D.append(img)
                
        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(particles, classes2D)
         
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        summary.append('Number of classes: *%s*' % self.getNumberOfClasses())
        summary.append('Number of factors: *%s*' % self.numberOfFactors)
        return summary
    
    def _methods(self):
        msg  = "\nInput particles %s " % self.getObjectTag('inputParticles')
        msg += "were divided into % classes using K-means classification " % self.getNumberOfClasses()
        msg += "(SPIDER command [[http://spider.wadsworth.org/spider_doc/spider/docs/man/clkm.html][CL KM]]) "
        msg += "using %s factors. " % self.numberOfFactors
        return [msg]
    
    #--------------------------- UTILS functions --------------------------------------------
    
    
