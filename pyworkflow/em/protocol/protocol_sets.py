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
- unions
- split
... etc
"""
import os
import math
from os.path import basename
from protocol import EMProtocol
from pyworkflow.protocol.params import (PointerParam, FileParam, StringParam, BooleanParam,
                                        MultiPointerParam, IntParam)
from pyworkflow.em.data import (SetOfImages, SetOfCTF, SetOfClasses, SetOfVolumes, Volume,
                                SetOfClasses2D, SetOfClasses3D) #we need to import this to be used dynamically 

from pyworkflow.em.data_tiltpairs import MicrographsTiltPair, TiltPair
from pyworkflow.utils.path import createLink

class ProtSets(EMProtocol):
    """ Base class for all protocols related to subsets. """
    pass


class ProtUserSubSet(ProtSets):
    """ Create subsets from the GUI.
    This protocol will be executed mainly calling the script 'pw_create_image_subsets.py'
    from the ShowJ gui. The enabled/disabled changes will be stored in a temporary sqlite
    file that will be read to create the new subset.
    """
    _label = 'create subset'
     
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addHidden('inputObject', PointerParam, pointerClass='EMObject')
        form.addHidden('other', StringParam, allowsNull=True)
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
        self._defineSourceRelation(inputImages, output)

    def _createVolumeFromSet(self, volfile):
        src = volfile
        if ':' in src:
            src = src.split(':')[0]
        dest = self._getPath(basename(src))
        createLink(src, dest)
        vol = Volume()
        vol.setFileName(dest)
        inputObject = self.inputObject.get()
        samplingRate = inputObject.getSamplingRate()
        vol.setSamplingRate(samplingRate)
        self._defineOutputs(outputVolume=vol)
        self._defineTransformRelation(inputObject, vol)



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
              
    def _createMicsSubSetFromCTF(self, inputCTFs):
        """ Create a subset of Micrographs analyzing the CTFs. """
        outputMics = self._createSetOfMicrographs()
        setOfMics = inputCTFs.getMicrographs()
        SetOfImages.copyInfo(outputMics, setOfMics)
        
        modifiedSet = SetOfCTF(filename=self._dbName, prefix=self._dbPrefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                mic = ctf.getMicrograph()
                outputMics.append(mic)
                
        self._defineOutputs(outputMicrographs=outputMics)
        self._defineTransformRelation(setOfMics, outputMics)
        
    def _createSubSetOfCTF(self, inputCtf):
        """ Create a subset of CTF and Micrographs analyzing the CTFs. """
        
        
        setOfCtf = self._createSetOfCTF("_subset")
        
        modifiedSet = SetOfCTF(filename=self._dbName, prefix=self._dbPrefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                mic = ctf.getMicrograph()
                setOfCtf.append(ctf)
                
         # Register outputs
        self._defineOutput(self.outputClassName.get(), setOfCtf)
        self._defineSourceRelation(inputCtf, setOfCtf)
        
        
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
        
    def _createSubSetFromMicrographsTiltPair(self, micrographsTiltPair):
        """ Create a subset of Micrographs Tilt Pair. """
        output = MicrographsTiltPair(filename=self._getPath('micrographs_pairs.sqlite'))
        print "self._dbName=%s" % self._dbName
        modifiedSet = MicrographsTiltPair(filename=self._dbName, prefix=self._dbPrefix)

        for micPairI in modifiedSet:
            untilted = micPairI.getUntilted()
            tilted = micPairI.getTilted()
            if micPairI.isEnabled():
                print "is enabled"
                micPairO = TiltPair()
                micPairO.setUntilted(untilted)
                micPairO.setTilted(tilted)
                output.append(micPairO)
        # Register outputs
        outputDict = {'outputMicrographsTiltPair': output}
        self._defineOutputs(**outputDict) 
        
    def createSetOfImagesStep(self):
        inputObj = self.inputObject.get()
        
        self._loadDbNamePrefix() # load self._dbName and self._dbPrefix


        if (isinstance(inputObj, SetOfVolumes) or isinstance(inputObj, SetOfClasses3D)) \
                and '.' in self.other.get():#is volume file, not id
                print 'other ' + self.other.get()
                volfile = self.other.get()
                self._createVolumeFromSet(volfile)

        elif isinstance(inputObj, SetOfImages):

                self._createSubSetFromImages(inputObj)
            
        elif isinstance(inputObj, SetOfClasses):
            self._createSubSetFromClasses(inputObj)
           
        elif isinstance(inputObj, SetOfCTF):
            outputClassName = self.outputClassName.get()
            if outputClassName.startswith('SetOfMicrographs'):
                self._createMicsSubSetFromCTF(inputObj)
            else:
                self._createSubSetOfCTF(inputObj)
                
        elif isinstance(inputObj, MicrographsTiltPair):
            self._createSubSetFromMicrographsTiltPair(inputObj)
        
        elif isinstance(inputObj, EMProtocol):
            setObj = self.getSetObject()
            otherObj = self.getProject().mapper.selectById(int(self.other.get()))

            
            if isinstance(setObj, SetOfClasses):
                setObj.setImages(otherObj)
                self._createSubSetFromClasses(setObj)
                
            elif isinstance(setObj, SetOfImages):
                setObj.copyInfo(otherObj) # copy info from original images
                self._createSubSetFromImages(setObj)
                

    
    def getSetObject(self):            
        from pyworkflow.mapper.sqlite import SqliteFlatDb
        db = SqliteFlatDb(dbName=self._dbName, tablePrefix=self._dbPrefix)
        setClassName = db.getProperty('self') # get the set class name
        setObj = globals()[setClassName](filename=self._dbName, prefix=self._dbPrefix)
        return setObj    
            
#     def defineOutputSet(self, outputset):
#         outputs = {'output' + self.getOutputType(): outputset}
#         self._defineOutputs(**outputs)
#         self._defineSourceRelation(self.getInput(), outputset)
    
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
        
        _dbName, self._dbPrefix = self.sqliteFile.get().split(',')
        self._dbName = self._getPath('subset.sqlite')
        os.rename(_dbName, self._dbName)

        if self._dbPrefix.endswith('_'):
            self._dbPrefix = self._dbPrefix[:-1]



            

class ProtUnionSet(ProtSets):
    """ Protocol to join two or more sets of images.
    This protocol allows to select two or more set of images
    and will produce another set joining all elements of the 
    selected sets. It will validate that all sets are of the
    same type of elements (Micrographs, Particles or Volumes) 
    """
    _label = 'union sets'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputSets', MultiPointerParam, label="Input set", important=True,
                      pointerClass='EMSet', minNumObjects=2, maxNumObjects=0,
                      help='Select two or more sets (of micrographs, particles, volumes, etc.) to be united.'
                           'If you select 3 sets with 100, 200, 200 elements, the final set will contain a '
                           'total of 500 elements.')
        # TODO: See what kind of restrictions we add (like "All sets should have the same sampling rate.")

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
                print "itemObj", itemObj
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
            m = "We have unioned the following sets: "
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
        
        form.addParam('inputSet', PointerParam, pointerClass='EMSet',
                      label="Input set", important=True,
                      help='Select the set of elements (images, etc) that you want to split.'
                      )
        form.addParam('numberOfSets', IntParam, default=2,
                      label="Number of subsets",
                      help='Select how many subsets do you want to create.'
                      )
        form.addParam('randomize', BooleanParam, default=False,
                      label="Randomize elements",
                      help='Put the elements at random in the different subsets.'
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

        # Iterate over the elements in the input set and assign
        # to different subsets.
        elements = self.inputSet.get()

        ns = [len(elements) // n + (1 if i < len(elements) % n else 0) for i in range(n)]
        pos, i = 0, 0  # index of current bucket and index of position inside it
        for elem in elements.__iter__(random=self.randomize):
            if i >= ns[pos]:
                pos += 1
                i = 0
            subsets[pos].append(elem)
            i += 1

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
    
    
class ProtSubSet(ProtSets):
    """    
    Create a set with the elements of an original set that are also
    referenced in another set.
    
    Usually there is a bigger set with all the elements, and a smaller
    one obtained from classification, cleaning, etc. The desired result
    is a set with the elements from the original set that are also present
    somehow in the smaller set (in the smaller set they may be downsampled
    or processed in some other way).
    
    Both sets should be of the same kind (micrographs, particles, volumes)
    or related (micrographs and CTFs for example).
    """
    _label = 'subset'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputFullSet', PointerParam, label="Full set of items",
                      important=True, pointerClass='EMSet',
                      help='Even if the operation can be applied to two arbitrary sets,\n'
                           'the most common use-case is to retrieve a subset of\n'
                           'elements from an original full set.\n' 
                           '*Note*: the elements of the resulting set will be the same\n'
                           'ones as this input set.')
        form.addParam('inputSubSet', PointerParam, label="Subset of items",
                      important=True, pointerClass='EMSet',
                      help='The elements that are in this (normally smaller) set and\n'
                           'in the full set will be included in the resulting set.')

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
        # Iterate over the elements in the smaller set
        # and take the info from the full set
        for elem in inputSubSet:
            # TODO: this can be improved if we perform intersection directly in sqlite
            origElem = inputFullSet[elem.getObjId()]
            if origElem is not None:
                outputSet.append(origElem)
            
        key = 'output' + inputClassName.replace('SetOf', '') 
        self._defineOutputs(**{key: outputSet})
        self._defineTransformRelation(inputFullSet, outputSet)
        self._defineSourceRelation(inputSubSet, outputSet)
        
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        """Make sure the input data make sense."""

        # self.inputFullSet and self.inputSubSet .get().getClassName() can be SetOf...
        #   Alignment
        #   Angles
        #   Averages
        #   Classes
        #   ClassesVol
        #   Coordinates
        #   CTF
        #   Micrographs
        #   MovieParticles
        #   Movies
        #   Particles
        #   Volumes

        c1 = self.inputFullSet.get().getClassName()
        c2 = self.inputSubSet.get().getClassName()

        if c1 == c2:
            return []

        # Avoid combinations that make no sense.
        for classA, classesIncompatible in [
            ('SetOfParticles', {'SetOfMicrographs', 'SetOfMovies', 'SetOfVolumes'}),
            ('SetOfCoordinates', {'SetOfMicrographs', 'SetOfMovies', 'SetOfVolumes'}),
            ('SetOfVolumes', {'SetOfMicrographs', 'SetOfMovies', 'SetOfParticles', 'SetOfCoordinates'})]:
            if ((c1 == classA and c2 in classesIncompatible) or
                (c2 == classA and c1 in classesIncompatible)):
                return ["The full set and the subset are of incompatible classes",
                        "%s and %s." % (c1, c2)]

        return []  # no errors
