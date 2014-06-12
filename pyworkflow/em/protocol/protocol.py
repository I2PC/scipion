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
In this module are protocol base classes related to EM.
Them should be sub-classes in the different sub-packages from
each EM-software packages.
"""

import os
import shutil
from glob import glob

from pyworkflow.object import String, Float
from pyworkflow.protocol import *
from pyworkflow.protocol.params import *
from pyworkflow.em.constants import *
from pyworkflow.em.data import *
from pyworkflow.utils.path import removeBaseExt, join, basename, cleanPath
from pyworkflow.utils.properties import Message, Icon 



class EMProtocol(Protocol):
    """ Base class to all EM protocols.
    It will contains some common functionalities. 
    """
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        
    def __createSet(self, SetClass, template, suffix):
        """ Create a set and set the filename using the suffix. 
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        cleanPath(setFn)
        setObj = SetClass(filename=setFn)
        return setObj
    
    def _createSetOfMicrographs(self, suffix=''):
        return self.__createSet(SetOfMicrographs, 'micrographs%s.sqlite', suffix)
    
    def _createSetOfCoordinates(self, micSet, suffix=''):
        coordSet = self.__createSet(SetOfCoordinates, 'coordinates%s.sqlite', suffix)
        coordSet.setMicrographs(micSet)       
        return coordSet
    
    def _createSetOfParticles(self, suffix=''):
        return self.__createSet(SetOfParticles, 'particles%s.sqlite', suffix)
    
    def _createSetOfClasses2D(self, imgSet, suffix=''):
        classes = self.__createSet(SetOfClasses2D, 'classes2D%s.sqlite', suffix)
        classes.setImages(imgSet)
        return classes
    
    def _createSetOfClasses3D(self, imgSet, suffix=''):
        classes =  self.__createSet(SetOfClasses3D, 'classes3D%s.sqlite', suffix)
        classes.setImages(imgSet)
        return classes
    
    def _createSetOfClassesVol(self, suffix=''):
        return self.__createSet(SetOfClassesVol, 'classesVol%s.sqlite', suffix)
    
    def _createSetOfVolumes(self, suffix=''):
        return self.__createSet(SetOfVolumes, 'volumes%s.sqlite', suffix)
    
    def _createSetOfCTF(self, suffix=''):
        return self.__createSet(SetOfCTF, 'ctfs%s.sqlite', suffix)
    
    def _createSetOfMovies(self, suffix=''):
        return self.__createSet(SetOfMovies, 'movies%s.sqlite', suffix)
    
    def _createSetOfAlignment(self, particles, suffix=''):
        alignment = self.__createSet(SetOfAlignment, 'alignment%s.sqlite', suffix)
        alignment.setParticles(particles)
        
        return alignment
    
    def _defineSourceRelation(self, srcObj, dstObj):
        """ Add a DATASOURCE relation between srcObj and dstObj """
        self._defineRelation(RELATION_SOURCE, srcObj, dstObj)
    
    def _defineTransformRelation(self, srcObj, dstObj):
        self._defineRelation(RELATION_TRANSFORM, srcObj, dstObj)
        # A transform relation allways implies a source relation
        self._defineSourceRelation(srcObj, dstObj)
        
    def _defineCtfRelation(self, srcObj, dstObj):
        self._defineRelation(RELATION_CTF, srcObj, dstObj)
        # A ctf relation allways implies a source relation
        self._defineSourceRelation(srcObj, dstObj)
    
    def _insertChild(self, key, child):
        if isinstance(child, Set):
            child.write()
        Protocol._insertChild(self, key, child)
        
    def _validateDim(self, obj1, obj2, errors, label1='Input 1', label2='Input 2'):
        """ Validate that obj1 and obj2 has the same dimensions.
        Params:
            obj1, obj2: input objects that can be Images or SetOfImages subclasses.
            errors: an error list where to append the error if not same dimensions
            label1, label2: labels for input objects used for the error message.
        """
        d1 = obj1.getDim()
        d2 = obj2.getDim()
        if d1 != d2:
            msg = '%s and %s have not the same dimensions, \n' % (label1, label2)
            msg += 'which are %d and %d, respectively' % (d1, d2)
            errors.append(msg)

 
class ProtSets(EMProtocol):
    pass


class ProtUserSubSet(ProtSets):
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self._inputType = String(args.get('inputType', None))
        self._outputType = String(args.get('outputType', None))
        self.setObjLabel(args.get('label', self.getClassName()))
        
    def _defineParams(self, form):
        form.addHidden('inputObject', PointerParam, pointerClass='EMObject')
        form.addHidden('sqliteFile', FileParam)
        form.addHidden('outputClassName', StringParam)
        
    def _insertAllSteps(self):
        self._insertFunctionStep('createSetOfImagesStep')
    
    def _createSubSetFromImages(self, inputImages):
        className = inputImages.getClassName()
        createFunc = getattr(self, '_create' + className)
        fn, prefix = self.sqliteFile.get().split(',')
        if prefix.endswith('_'):
            prefix = prefix[:-1]
        modifiedSet = inputImages.getClass()(filename=fn, prefix=prefix)
        
        output = createFunc()
        output.copyInfo(inputImages)
        for img in modifiedSet:
            if img.isEnabled():                
                output.append(img)
        # Register outputs
        outputDict = {'output' + className.replace('SetOf', ''): output}
        self._defineOutputs(**outputDict)
    
    def _createSubSetFromCTF(self, inputCTFs):
        """ Create a subset of Micrographs analyzing the CTFs. """
        output = self._createSetOfMicrographs()
        firstMic = inputCTFs.getFirstItem().getMicrograph()
        SetOfImages.copyInfo(output, firstMic)
        
        fn, prefix = self.sqliteFile.get().split(',')
        if prefix.endswith('_'):
            prefix = prefix[:-1]
        modifiedSet = SetOfCTF(filename=fn, prefix=prefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                mic = ctf.getMicrograph()
                print "ctfId = ", ctf.getObjId()
                print "micId = ", mic.getObjId()
                output.append(mic)
                
        self._defineOutputs(outputMicrographs=output)
        
    def createSetOfImagesStep(self):
        inputObj = self.inputObject.get()
        
        if isinstance(inputObj, SetOfImages):
            self._createSubSetFromImages(inputObj)
            
        elif isinstance(inputObj, SetOfClasses):
            self._createSubSetFromImages(inputObj.getImages())
           
        elif isinstance(inputObj, SetOfCTF):
            self._createSubSetFromCTF(inputObj) 
            
    def defineOutputSet(self, outputset):
        outputs = {'output' + self.getOutputType(): outputset}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.getInput(), outputset)
    
    def _summary(self):
        summary = []
#         input = self.getInputType()
#         if hasattr(self.getInput(), "getSize"):
#             input = "%d %s"%(self.getInput().getSize(), input)
#         summary.append("From input set of %s created subset of %s %s"%(input, self.getOutputSet().getSize(), self.getOutputType()))
        return summary

    def _methods(self):
        return self._summary()


class ProtJoinSets(ProtSets):
    """ Protocol to join two sets. """
    _label = 'join sets'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputSets', MultiPointerParam, label="Input set of images", important=True, 
                      pointerClass='SetOfImages', minNumObjects=2, maxNumObjects=0,
                      help='Select two or more set of images (micrographs, particles o volumes) to be joined.'
                           'If you select 3 sets with 100, 200, 200 images, the final set should contans a '
                           'total of 500 images. All sets should have the same sampling rate.'
                           )
#         form.addParam('inputSet1', PointerParam, label="Input set 1", 
#                       pointerClass='SetOfImages', 
#                       help='Select the set of images (it can be a set of micrographs, particles o volumes) from the project.'
#                            'They should 2 or more object classes')
#         form.addParam('inputSet2', PointerParam, label="Input set 2", 
#                       pointerClass='SetOfImages', 
#                       help='Select the set of images (it can be a set of micrographs, particles o volumes) from the project.'
#                            'They should 2 or more object classes')
    
    #--------------------------- INSERT steps functions --------------------------------------------   
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutput')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutput(self):
        #Read Classname and generate corresponding SetOfImages (SetOfParticles, SetOfVolumes, SetOfMicrographs)
        self.inputSets.printAll()
        
        self.inputType = str(self.inputSets[0].get().getClassName())
        #self.inputType = str(self.inputSet1.get().getClassName())
        outputSetFunction = getattr(self, "_create%s" % self.inputType)
        outputSet = outputSetFunction()
        
        #Copy info from input (sampling rate, etc)
        outputSet.copyInfo(self.inputSets[0].get())
        #outputSet.copyInfo(self.inputSet1.get())
        #inputSets = [self.inputSet1, self.inputSet2]
       
        for itemSet in self.inputSets:
            for itemObj in itemSet.get():
                itemObj.cleanObjId()
                outputSet.append(itemObj)
        
        self._defineOutputs(outputImages=outputSet)
        
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        classList = []
        #inputSets = [self.inputSet1, self.inputSet2]
        for itemSet in self.inputSets:
            itemClassName = itemSet.get().getClassName()
            if len(classList) == 0 or itemClassName not in classList:
                classList.append(itemClassName)
            
        errors = []
        if len(classList) > 1:
            errors.append("Object should have same type")
            errors.append("Types of objects found: " + ", ".join(classList))
        return errors   

    def _summary(self):
        summary = []
#         if not hasattr(self, 'outputImages'):
#             summary.append("Output set not ready yet.")
#         else:
#             summary.append("Input sets of type %s:" % self.outputImages.getClassName())
#             inputSets = [self.inputSet1, self.inputSet2]
#             for itemSet in inputSets:
#                 summary.append("%s" % itemSet.get().getNameId())
        return summary
        
    def _methods(self):
        methods = []
        if not hasattr(self, 'outputImages'):
            methods.append("Protocol has not finished yet.")
        else:
            m = "We have joint the following sets: "
            #inputSets = [self.inputSet1, self.inputSet2]
            for itemSet in self.inputSets:
                m += "%s, " % itemSet.get().getNameId()
            methods.append(m[:-2])
        
        return methods
