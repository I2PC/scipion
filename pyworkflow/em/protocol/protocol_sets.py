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
This module contains protocols related to Set operations such us:
- subsets
- joins
- split
... etc
"""

from protocol import EMProtocol
from pyworkflow.protocol.params import (PointerParam, FileParam, StringParam,
                                        MultiPointerParam, IntParam)
from pyworkflow.em.data import SetOfImages, SetOfCTF, SetOfClasses



class ProtSets(EMProtocol):
    """ Base class for all protocols related to subsets. """
    pass


class ProtUserSubSet(ProtSets):
    """ Create subsets from the GUI.
    This protocol will be executed mainly calling the script 'pw_create_image_subsets.py'
    from the ShowJ gui. The enabled/disabled changes will be stored in a temporary sqlite
    file that will be read to create the new subset.
    """
     
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addHidden('inputObject', PointerParam, pointerClass='EMObject')
        form.addHidden('otherObject', PointerParam, pointerClass='EMObject', allowsNull=True)
        form.addHidden('sqliteFile', FileParam)
        form.addHidden('outputClassName', StringParam)
        
    def _insertAllSteps(self):
        self._insertFunctionStep('createSetOfImagesStep')
    
    def _createSubSetFromImages(self, inputImages):
        className = inputImages.getClassName()
        createFunc = getattr(self, '_create' + className)
        modifiedSet = inputImages.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        
        output = createFunc()
        output.copyInfo(inputImages)
        output.appendFromImages(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)
        
    def _createSubSetFromClasses(self, inputClasses):
        outputClassName = self.outputClassName.get()
        
        if (outputClassName.startswith('SetOfParticles') or
            outputClassName.startswith('SetOfVolumes')):
            # We need to distinguish two cases:
            # a) when we want to create images by grouping class images
            # b) create a subset from a particular class images
            from pyworkflow.mapper.sqlite import SqliteFlatDb
            db = SqliteFlatDb(dbName=self._dbName, tablePrefix=self._dbPrefix)
            itemClassName = db.getSelfClassName()
            if itemClassName.startswith('Class'):
                if outputClassName.endswith('Representatives'):
                    self._createRepresentativesFromClasses(inputClasses, 
                                                           outputClassName.split(',')[0])
                else:
                    self._createImagesFromClasses(inputClasses)
            else:
                self._createSubSetFromImages(inputClasses.getImages())
        elif outputClassName.startswith('SetOfClasses'):
            self._createClassesFromClasses(inputClasses)
        else:
            raise Exception("Unrecognized output type: '%s'" % outputClassName)  
              
    def _createSubSetFromCTF(self, inputCTFs):
        """ Create a subset of Micrographs analyzing the CTFs. """
        output = self._createSetOfMicrographs()
        setOfMics = inputCTFs.getMicrographs()
        SetOfImages.copyInfo(output, setOfMics)
        
        modifiedSet = SetOfCTF(filename=self._dbName, prefix=self._dbPrefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                mic = ctf.getMicrograph()
                output.append(mic)
                
        self._defineOutputs(outputMicrographs=output)
        
    def _createRepresentativesFromClasses(self, inputClasses, outputClassName):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        inputImages = inputClasses.getImages()
        createFunc = getattr(self, '_create' + outputClassName)
        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating REPRESENTATIVES of images from classes,  sqlite file: %s" % self._dbName)
        
        output = createFunc()
        output.copyInfo(inputImages)
        for cls in modifiedSet:
            if cls.isEnabled():
                img = cls.getRepresentative()
                img.copyObjId(cls)
                output.append(img)
        # Register outputs
        self._defineOutput('Representatives', output)
        
    def _createImagesFromClasses(self, inputClasses):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        inputImages = inputClasses.getImages()
        className = inputImages.getClassName()
        createFunc = getattr(self, '_create' + className)
        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating subset of images from classes,  sqlite file: %s" % self._dbName)
        
        output = createFunc()
        output.copyInfo(inputImages)
        output.appendFromClasses(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)
 
    def _createClassesFromClasses(self, inputClasses):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        #inputImages = inputClasses.getImages()
        className = inputClasses.getClassName()
        createFunc = getattr(self, '_create' + className)
        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating subset of classes from classes,  sqlite file: %s" % self._dbName)
        
        output = createFunc(inputClasses.getImages())
        output.appendFromClasses(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)
               
    def createSetOfImagesStep(self):
        inputObj = self.inputObject.get()
        
        self._loadDbNamePrefix() # load self._dbName and self._dbPrefix
        
        if isinstance(inputObj, SetOfImages):
            self._createSubSetFromImages(inputObj)
            
        elif isinstance(inputObj, SetOfClasses):
            self._createSubSetFromClasses(inputObj)
           
        elif isinstance(inputObj, SetOfCTF):
            self._createSubSetFromCTF(inputObj) 
            
        elif isinstance(inputObj, EMProtocol):
            from pyworkflow.mapper.sqlite import SqliteFlatDb
            db = SqliteFlatDb(dbName=self._dbName, tablePrefix=self._dbPrefix)
            setClassName = db.getProperty('self') # get the set class name
            setObj = globals()[setClassName](filename=self._dbName, prefix=self._dbPrefix)
            otherObj = self.otherObject.get()
            
            if isinstance(setObj, SetOfClasses):
                setObj.setImages(otherObj)
                self._createSubSetFromClasses(setObj)
                
            elif isinstance(setObj, SetOfImages):
                setObj.copyInfo(otherObj) # copy info from original images
                self._createSubSetFromImages(setObj)
                
            # Clean the pointer to other, to avoid dependency in the graph
            self.otherObject.set(None)
                
            
            
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

    def _defineOutput(self, className, output):
        outputDict = {'output' + className.replace('SetOf', ''): output}
        self._defineOutputs(**outputDict) 
        
    def _loadDbNamePrefix(self):
        """ Setup filename and prefix for db connection. """
        self._dbName, self._dbPrefix = self.sqliteFile.get().split(',')
        if self._dbPrefix.endswith('_'):
            self._dbPrefix = self._dbPrefix[:-1] 
            

class ProtJoinSets(ProtSets):
    """ Protocol to join two or more sets of images.
    This protocol allows to select two or more set of images 
    and will produce another set joining all elements of the 
    selected sets. It will validate that all sets are of the
    same type of elements (Micrographs, Particles or Volumes) 
    """
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
    
    #--------------------------- INSERT steps functions --------------------------------------------   
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        #Read Classname and generate corresponding SetOfImages (SetOfParticles, SetOfVolumes, SetOfMicrographs)
        self.inputType = str(self.inputSets[0].get().getClassName())
        outputSetFunction = getattr(self, "_create%s" % self.inputType)
        outputSet = outputSetFunction()
        
        #Copy info from input (sampling rate, etc)
        outputSet.copyInfo(self.inputSets[0].get())
       
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

        if not hasattr(self, 'outputImages'):
            summary.append("Protocol has not finished yet.")
        else:
            m = "We have joint the following sets: "
            #inputSets = [self.inputSet1, self.inputSet2]
            for itemSet in self.inputSets:
                m += "%s, " % itemSet.get().getNameId()
            summary.append(m[:-2])
        
        return summary
        
    def _methods(self):
        return self._summary()
            

class ProtSplitSet(ProtSets):
    """ Protocol to split a set in two or more subsets. 
    """
    _label = 'split sets'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputSet', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfImages', 
                      help='Select the set of images that you want to split. '
                           )
        form.addParam('numberOfSets', IntParam, label="Number of subsets", default=2,
                      help='Select how many subsets do you want to create. '
                           )    
    #--------------------------- INSERT steps functions --------------------------------------------   
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        inputClassName = str(inputSet.getClassName())
        outputSetFunction = getattr(self, "_create%s" % inputClassName)
        n = self.numberOfSets.get()
        
        # Create as many subsets as requested by the user
        subsets = [outputSetFunction(suffix=str(i)) for i in range(1, n+1)]
        # Iterate over the images in the input set and assign
        # diffent subsets
        for i, img in enumerate(self.inputSet.get()):
            subsets[i % n].append(img)
            
        key = 'output' + inputClassName.replace('SetOf', '') + '%02d'
        for i in range(1, n+1):
            subset = subsets[i-1]
            subset.copyInfo(inputSet)
            self._defineOutputs(**{key % i: subset})
            self._defineTransformRelation(inputSet, subset)
        
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.inputSet.get().getSize() < self.numberOfSets:
            errors.append("The number of subsets requested is greater than")
            errors.append("the number of elements in the input set.")
        return errors   
    
    
class ProtIntersectSet(ProtSets):
    """    
    Create a set with the intersection (common elements) 
    from two different sets.
    
    Usually, there is a bigger set with all elements...
    and a smaller one (obtained from classification, cleanning...). 
    In that case, the desired result is the elements from the 
    original set that are present in the smaller set. 
    
    Both set should be of the same type of elements 
    (micrographs, particles, volumes)
    """
    _label = 'intersect sets'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputFullSet', PointerParam, label="Full images set", important=True, 
                      pointerClass='SetOfImages', 
                      help='Even if the intersection can be applied to two subsets,\n'
                           'the most common use-case is to retrieve a subset of  \n'
                           'elements from an original full set.\n' 
                           '*Note*: the images of the result set will be the same \n'
                           'ones of this input set.'
                           )
        form.addParam('inputSubSet', PointerParam, label="Subset of images", important=True, 
                      pointerClass='SetOfImages', 
                      help='The elements that are in this (normally smaller) set and \n'
                           'in the full set will be included in the result set'
                           )
        
    #--------------------------- INSERT steps functions --------------------------------------------   
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        inputFullSet = self.inputFullSet.get()
        inputSubSet = self.inputSubSet.get()
        
        inputClassName = inputFullSet.getClassName()
        outputSetFunction = getattr(self, "_create%s" % inputClassName)

        outputSet = outputSetFunction()
        outputSet.copyInfo(inputFullSet)    
        # Iterate over the images in the smaller set
        # and take the info from the full set
        for img in inputSubSet:
            #TODO: this can be improved if you perform an
            # intersection directly in sqlite
            origImg = inputFullSet[img.getObjId()]
            outputSet.append(origImg)
            
        key = 'output' + inputClassName.replace('SetOf', '') 
        self._defineOutputs(**{key: outputSet})
        self._defineTransformRelation(inputFullSet, outputSet)
        self._defineSourceRelation(inputSubSet, outputSet)
        
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.inputFullSet.get().getClassName() != self.inputSubSet.get().getClassName():
            errors.append("Both the full set and the subset should be of the same")
            errors.append("type of elements (micrographs, particles, volumes).")
        return errors   

    