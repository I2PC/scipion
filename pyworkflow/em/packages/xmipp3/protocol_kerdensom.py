# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, dirname, exists
from pyworkflow.em import *  
import xmipp

import xmipp3
from convert import writeSetOfParticles, readSetOfClasses2D, xmippToLocation
from glob import glob


class KendersomBaseClassify(ProtClassify2D):
    """ Class to create a base template for Kendersom and rotational spectra protocols that share
    a common structure. 
    """
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles', pointerCondition='hasAlignment',
                      help='Select the input images from the project.'
                           'It should be a SetOfParticles class')
        form.addParam('useMask', BooleanParam, default=False,
                      label='Use a Mask ?', 
                      help='If you set to *Yes*, you should provide a mask')
        form.addParam('Mask', PointerParam , condition='useMask',
                      label="Mask", pointerClass='Mask',
                      help='Mask image will serve to enhance the classification')
        
        line = form.addLine('Dimension of the map', 
                            help='Josue tiene que meter el help')
        line.addParam('SomXdim', IntParam, default=7,
                      label='X')
        line.addParam('SomYdim', IntParam, default=7,
                      label='Y')
               
        form.addParam('SomReg0', IntParam, default=1000, expertLevel=LEVEL_ADVANCED,
                      label='Initial regularization factor', 
                      help='The kerdenSOM algorithm anneals from an initial high regularization factor'
                      'to a final lower one, in a user-defined number of steps.'
                      'If the output map is too smooth, lower the regularization factors'
                      'If the output map is not organized, higher the regularization factors'
                      'See [[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM][KerDenSOM]]')
        form.addParam('SomReg1', IntParam, default=200, expertLevel=LEVEL_ADVANCED,
                      label='Final regularization factor:')
        form.addParam('SomSteps', IntParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label='Regularization steps:',
                      help='Number of steps to lower the regularization factor')
        form.addParam('extraParams', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label="Additional parameters:", 
                      help='Additional parameters for kerdensom program. \n For a complete description'
                      'See [[http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM][KerDenSOM]]')
        self._addParams(form)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _prepareParams(self):
        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('images.xmd') 
        self._insertFunctionStep('convertInputStep')  
        
        if self.useMask:
            mask = self.Mask.get().getFileName()
        else:
            mask = None
            
        self._params = {'oroot': self._getExtraPath("kerdensom"),
                        'imgsFn': self.imgsFn,
                        'mask': mask,
                        'SomXdim': self.SomXdim.get(),
                        'SomYdim': self.SomYdim.get(),
                        'SomReg0': self.SomReg0.get(),
                        'SomReg1': self.SomReg1.get(),
                        'SomSteps': self.SomSteps.get(),
                        'extraParams': self.extraParams.get(),
                        'vectors': self._getExtraPath("vectors.xmd"),
                        'classes': self._getExtraPath("classes.stk"),
                        'averages': self._getExtraPath("averages.stk"),
                        'kvectors': self._getExtraPath("kerdensom_vectors.xmd"),
                        'kclasses': self._getExtraPath("kerdensom_classes.xmd")
                       }
    
    def _insertAllSteps(self):
        self._prepareParams()
        self._insertProccesStep()
        self._insertFunctionStep('rewriteClassBlockStep')
        self._insertFunctionStep('createOutputStep')
    
    def _insertKerdensomStep(self):
        args = '-i %(vectors)s --oroot %(oroot)s --xdim %(SomXdim)d --ydim %(SomYdim)d' + \
               ' --deterministic_annealing %(SomSteps)f %(SomReg0)f %(SomReg1)f %(extraParams)s'
        self._insertRunJobStep("xmipp_classify_kerdensom", args % self._params)
#        deleteFiles([self._getExtraPath("vectors.xmd"),self._getExtraPath("vectors.vec")], True)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn) 

    #--------------------------- INFO functions ----------------------------------------------------
    def rewriteClassBlockStep(self):
        firstImage = self.inputParticles.get().getFirstItem()
        fnClasses = self._params['kclasses']
        mdClasses = "classes@%s" % fnClasses
        fnClassStack = self._params['classes']
        fnAverageStack = self._params['averages']      
        
        md = xmipp.MetaData(mdClasses)
        image = ImageHandler().createImage()
        
        counter = 1
        
        for objId in md:
            imageName =  "%06d@%s" % (counter, fnClassStack)
            averageName = "%06d@%s" % (counter, fnAverageStack)
            
            if md.getValue(xmipp.MDL_CLASS_COUNT, objId) > 0:
                # compute the average of images assigned to this class
                classPrefix = 'class%06d' % counter
                classMd = '%s_images@%s' % (classPrefix, fnClasses)
                classRoot = self._getTmpPath(classPrefix)
                self.runJob('xmipp_image_statistics', 
                            '-i %s --save_image_stats %s -v 0' % (classMd, classRoot))
                image.read(classRoot + 'average.xmp')
            else:
                # Create emtpy image as average
                image.read(firstImage.getLocation()) # just to take the right dimensions and datatype
                image.initConstant(0.)
                
            image.write(averageName)
            md.setValue(xmipp.MDL_IMAGE, imageName, objId)
            md.setValue(xmipp.MDL_IMAGE2, averageName, objId)
            
            counter += 1
            
        md.write(mdClasses, xmipp.MD_APPEND)
        
    def _preprocessClass(self, classItem, classRow):
        classItem.average = Particle()
        classItem.average.setLocation(xmippToLocation(classRow.getValue(xmipp.MDL_IMAGE2)))
        
    def createOutputStep(self):
        """ Store the kenserdom object 
        as result of the protocol.
        """
        imgSet = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(imgSet)
        readSetOfClasses2D(classes2DSet, self._params['kclasses'], 
                           preprocessClass=self._preprocessClass)
        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(self.inputParticles, classes2DSet)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        if self.SomReg0 < self.SomReg1:
            errors.append("Regularization must decrease over iterations:")
            errors.append("    Initial regularization must be larger than final")
        if self.useMask:
            mask = self.Mask.get()
            if mask is None:
                errors.append("You have selected to use a mask. Select one.")
        return errors
    
    def _summary(self):
        return self._methods()

    def _methods(self):
        messages = []  
        if not hasattr(self, 'outputClasses'):
            messages.append("Output classification not ready yet.")
        elif self.inputParticles.get() is None:
            messages.append('Input not selected yet.')
        else:    
            messages.append("*Kendersom classification*")
            messages.append('%s particles from %s were classified to obtain %s classes %s.'
                            % (self.inputParticles.get().getSize(), self.getObjectTag('inputParticles'), self.outputClasses.getSize(), self.getObjectTag('outputClasses')))
            if self.useMask:
                messages.append('Mask %s was used in classification.' % self.getObjectTag('Mask'))
        return messages


class XmippProtKerdensom(KendersomBaseClassify):
    """
    Classifies a set of images using  Kohonen's Self-Organizing Feature Maps (SOM) 
    and Fuzzy c-means clustering technique (FCM) .
    
    The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    
    KerdenSOM is an excellent tool for classification, especially when
    using a large number of data and classes and when the transition between
    the classes is almost continuous, with no clear separation between them.
    
    The input images must be previously aligned.
    """
    _label = 'kerdensom'
    
    def __init__(self, **args):
        KendersomBaseClassify.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _addParams(self, form):
        pass
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertProccesStep(self):
        self._insertImgToVectorStep()
        self._insertKerdensomStep()
        self._insertVectorToImgStep()
    
    def _insertImgToVectorStep(self):
        """ Insert runJob for convert into a vector Md """
        args = ' -i %(imgsFn)s -o %(vectors)s '
        if self.useMask:
            args += ' --mask binary_file %(mask)s'
        
        self._insertRunJobStep("xmipp_image_vectorize", args % self._params)
   
    def _insertVectorToImgStep(self):
        args = ' -i %(kvectors)s -o %(classes)s' 
        if self.useMask:
            args += ' --mask binary_file %(mask)s'
        self._insertRunJobStep("xmipp_image_vectorize", args % self._params)
#        deleteFiles([self._getExtraPath("kerdensom_vectors.xmd"),self._getExtraPath("kerdensom_vectors.vec")], True)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        return KendersomBaseClassify._validate(self)
    
    def _summary(self):
        return KendersomBaseClassify._summary(self)
    
    def _methods(self):
        return KendersomBaseClassify._methods(self)
    
    def _citations(self):
        return ['PascualMontano2001', 'PascualMontano2002']
