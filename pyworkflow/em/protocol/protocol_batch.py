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
This module contains protocols that are launched
throught other GUI (such as showj) and that
are called "batch" protocols.
"""

import os

from pyworkflow.protocol.params import PointerParam, FileParam, StringParam, IntParam
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import SetOfImages, SetOfCTF, SetOfClasses, SetOfClasses3D, SetOfVolumes
from pyworkflow.em.data_tiltpairs import TiltPair, MicrographsTiltPair



class BatchProtocol(EMProtocol):
    """ Base class to all protocols that are launched
    throught other GUIs (such as showj) and that
    are called "batch" protocols. They should not be
    executed from normal "form" of other protocols.
    """
    pass


class ProtUserSubSet(BatchProtocol):
    """ Create subsets from the GUI.
    This protocol will be executed mainly calling the script 'pw_create_image_subsets.py'
    from the ShowJ gui. The enabled/disabled changes will be stored in a temporary sqlite
    file that will be read to create the new subset.
    """
    _label = 'create subset'
     
    def __init__(self, **args):
        BatchProtocol.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addHidden('inputObject', PointerParam, pointerClass='XmippProtDimredNMA')
        form.addHidden('otherObj', PointerParam, pointerClass='EMObject', allowsNull=True)
        form.addHidden('sqliteFile', FileParam)
        form.addHidden('outputClassName', StringParam)
        form.addHidden('volId', IntParam, allowsNull=True)
        
    def _insertAllSteps(self):
        self._insertFunctionStep('createSetStep')
    
    def _createSubSetFromImages(self, inputImages):
        className = inputImages.getClassName()
        createFunc = getattr(self, '_create' + className)
        modifiedSet = inputImages.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        
        output = createFunc()
        output.copyInfo(inputImages)
        output.appendFromImages(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)
        return output

    def _createSubSetFromClasses(self, inputClasses):
        outputClassName = self.outputClassName.get()
        
        if (outputClassName.startswith('SetOfAverages') or
            outputClassName.startswith('SetOfVolumes') or
            outputClassName.startswith('SetOfParticles')):
            # We need to distinguish two cases:
            # a) when we want to create images by grouping class images
            # b) create a subset from a particular class images
            from pyworkflow.mapper.sqlite import SqliteFlatDb
            db = SqliteFlatDb(dbName=self._dbName, tablePrefix=self._dbPrefix)
            itemClassName = db.getSelfClassName()
            if itemClassName.startswith('Class'):
                if outputClassName.startswith('SetOfParticles'):
                    return self._createImagesFromClasses(inputClasses)
                else:
                    return self._createRepresentativesFromClasses(inputClasses, outputClassName.split(',')[0])
            else:
                return self._createSubSetFromImages(inputClasses.getImages())
        elif outputClassName.startswith('SetOfClasses'):
            return self._createClassesFromClasses(inputClasses)
        else:
            raise Exception("Unrecognized output type: '%s'" % outputClassName)  
              
    def _createMicsSubSetFromCTF(self, inputCTFs):
        """ Create a subset of Micrographs analyzing the CTFs. """
        outputMics = self._createSetOfMicrographs()
        setOfMics = inputCTFs.getMicrographs()
        outputMics.copyInfo(setOfMics)
        
        modifiedSet = SetOfCTF(filename=self._dbName, prefix=self._dbPrefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                mic = ctf.getMicrograph()
                outputMics.append(mic)
                
        self._defineOutputs(outputMicrographs=outputMics)
        self._defineTransformRelation(setOfMics, outputMics)
        return outputMics
        
    def _createSubSetOfCTF(self, inputCtf):
        """ Create a subset of CTF and Micrographs analyzing the CTFs. """
        
        setOfCtf = self._createSetOfCTF("_subset")
        
        modifiedSet = SetOfCTF(filename=self._dbName, prefix=self._dbPrefix)
        
        for ctf in modifiedSet:
            if ctf.isEnabled():
                setOfCtf.append(ctf)
                
        # Register outputs
        self._defineOutput(self.outputClassName.get(), setOfCtf)
        self._defineSourceRelation(inputCtf, setOfCtf)
        return setOfCtf
        
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
        return output

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
        for sampleClass in inputClasses:
            alignment = sampleClass.getAlignment()
            break
        output.setAlignment(alignment)
        output.appendFromClasses(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)
        return output
 
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
        return output
        
    def _createSubSetFromMicrographsTiltPair(self, micrographsTiltPair):
        """ Create a subset of Micrographs Tilt Pair. """
        output = MicrographsTiltPair(filename=self._getPath('micrographs_pairs.sqlite'))
        print "self._dbName=%s" % self._dbName
        modifiedSet = MicrographsTiltPair(filename=self._dbName, prefix=self._dbPrefix)

        for micPairI in modifiedSet:
            untilted = micPairI.getUntilted()
            tilted = micPairI.getTilted()
            if micPairI.isEnabled():

                micPairO = TiltPair()
                micPairO.setUntilted(untilted)
                micPairO.setTilted(tilted)
                output.append(micPairO)
        # Register outputs
        outputDict = {'outputMicrographsTiltPair': output}
        self._defineOutputs(**outputDict)
        return output

    def createSetStep(self):
        setObj = self.createSetObject()
        inputObj = self.inputObject.get()

        if self.volId.get():

            if isinstance(setObj, SetOfVolumes):
                volSet = SetOfVolumes(filename=self._dbName)
                output = volSet[self.volId.get()]
            else:
                classSet = SetOfClasses3D(filename=self._dbName)
                output = classSet[self.volId.get()].getRepresentative()
            self._defineOutputs(outputVolume=output)

        elif isinstance(inputObj, SetOfImages):
                output = self._createSubSetFromImages(inputObj)

        elif isinstance(inputObj, SetOfClasses):
            output = self._createSubSetFromClasses(inputObj)

        elif isinstance(inputObj, SetOfCTF):
            outputClassName = self.outputClassName.get()
            if outputClassName.startswith('SetOfMicrographs'):
                self._createMicsSubSetFromCTF(inputObj)
            else:
                self._createSubSetOfCTF(inputObj)

        elif isinstance(inputObj, MicrographsTiltPair):
            self._createSubSetFromMicrographsTiltPair(inputObj)

        elif isinstance(inputObj, EMProtocol):
            otherObj = self.otherObj.get()

            if isinstance(setObj, SetOfClasses):
                setObj.setImages(otherObj)
                output = self._createSubSetFromClasses(setObj)

            elif isinstance(setObj, SetOfImages):
                setObj.copyInfo(otherObj) # copy info from original images
                output = self._createSubSetFromImages(setObj)
        else:
            className = inputObj.getClassName()
            createFunc = getattr(self, '_create' + className)
            modifiedSet = inputObj.getClass()(filename=self._dbName, prefix=self._dbPrefix)

            output = createFunc()
            for item in modifiedSet:
                if item.isEnabled():
                    output.append(item)

            if hasattr(modifiedSet, 'copyInfo'):
                modifiedSet.copyInfo(output)
            # Register outputs
            self._defineOutput(className, output)


        if isinstance(inputObj, EMProtocol):
            for key, attr in inputObj.iterInputAttributes():
                print attr
                self._defineSourceRelation(attr.get(), output)
        else:
            if not isinstance(inputObj, SetOfCTF):#otherwise setted before
                self._defineSourceRelation(inputObj, output)
    
    def createSetObject(self):
        _dbName, self._dbPrefix = self.sqliteFile.get().split(',')
        self._dbName = self._getPath('subset.sqlite')
        os.rename(_dbName, self._dbName)

        if self._dbPrefix.endswith('_'):
            self._dbPrefix = self._dbPrefix[:-1]

        from pyworkflow.mapper.sqlite import SqliteFlatDb
        db = SqliteFlatDb(dbName=self._dbName, tablePrefix=self._dbPrefix)
        setClassName = db.getProperty('self') # get the set class name
        from pyworkflow.em import getObjects
        setObj = getObjects()[setClassName](filename=self._dbName, prefix=self._dbPrefix)
        return setObj

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return self._summary()

    def _defineOutput(self, className, output):
        outputDict = {'output' + className.replace('SetOf', ''): output}
        self._defineOutputs(**outputDict) 
        

