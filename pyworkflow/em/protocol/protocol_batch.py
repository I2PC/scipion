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
through other GUI (such as showj) and that
are called "batch" protocols.
"""

import os
from itertools import izip

from pyworkflow.protocol.params import PointerParam, FileParam, StringParam
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import (SetOfImages, SetOfCTF, SetOfClasses,
                                SetOfClasses3D, SetOfVolumes, EMObject, EMSet,
                                SetOfNormalModes, SetOfParticles, FSC,
                                Class2D, Class3D, SetOfMicrographs)
from pyworkflow.em.data_tiltpairs import (TiltPair, MicrographsTiltPair,
                                          ParticlesTiltPair)
from pyworkflow.em.data import Mask
from pyworkflow.utils import moveFile



class BatchProtocol(EMProtocol):
    """ Base class to all protocols that are launched
    through other GUIs (such as showj) and that
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
        form.addHidden('inputObject', PointerParam, pointerClass='Object')
        form.addHidden('other', StringParam, allowsNull=True)
        form.addHidden('sqliteFile', FileParam)
        form.addHidden('outputClassName', StringParam)

    def _insertAllSteps(self):
        self._insertFunctionStep('createSetStep')

    def createSetStep(self):
        setObj = self.createSetObject()
        inputObj = self.inputObject.get()
        other = self.other.get()

        if other and ',Volume' in other:
            volId = int(other.split(',')[0])

            if isinstance(setObj, SetOfVolumes):
                volSet = SetOfVolumes(filename=self._dbName)
                output = volSet[volId]
            else:
                classSet = SetOfClasses3D(filename=self._dbName)
                output = classSet[volId].getRepresentative()
            self._defineOutputs(outputVolume=output)

        elif isinstance(inputObj, SetOfImages):
                output = self._createSubSetFromImages(inputObj)

        elif isinstance(inputObj, SetOfClasses):
            output = self._createSubSetFromClasses(inputObj)

        elif isinstance(inputObj, SetOfCTF):
            outputClassName = self.outputClassName.get()
            if outputClassName.startswith('SetOfMicrographs'):
                output = self._createMicsSubSetFromCTF(inputObj)
            else:
                output = self._createSubSetOfCTF(inputObj)
        
        elif isinstance(inputObj, MicrographsTiltPair):
            output = self._createSubSetFromMicrographsTiltPair(inputObj)

        elif isinstance(inputObj, ParticlesTiltPair):
            output = self._createSubSetFromParticlesTiltPair(inputObj)

        elif isinstance(inputObj, EMProtocol):
            otherid = self.other.get()
            otherObj = self.getProject().mapper.selectById(int(otherid))

            if isinstance(setObj, SetOfClasses):
                setObj.setImages(otherObj)
                self._createSubSetFromClasses(setObj)

            elif isinstance(setObj, SetOfImages):
                setObj.copyInfo(otherObj) # copy info from original images
                self._createSubSetFromImages(setObj)

            elif isinstance(setObj, SetOfNormalModes):
                self._createSimpleSubset(otherObj)
            
        else:
            output = self._createSimpleSubset(inputObj)
    
    def _createSimpleSubset(self, inputObj):
        className = inputObj.getClassName()
        createFunc = getattr(self, '_create' + className)
        modifiedSet = inputObj.getClass()(filename=self._dbName, prefix=self._dbPrefix)

        output = createFunc()
        for item in modifiedSet:
            if item.isEnabled():
                output.append(item)

        if hasattr(modifiedSet, 'copyInfo'):
            output.copyInfo(inputObj)

        # Register outputs
        self._defineOutput(className, output)

        if inputObj.hasObjId():
            self._defineTransformRelation(inputObj, output)

        return output

    def _createSubSetFromImages(self, inputImages,
                                copyInfoCallback=None):
        className = inputImages.getClassName()
        setClass = inputImages.getClass()
        inputClass = False

        createFunc = getattr(self, '_create' + className)
        modifiedSet = setClass(filename=self._dbName, prefix=self._dbPrefix)

        output = createFunc()

        if copyInfoCallback is None:
            modifiedSet.loadAllProperties()
            output.copyInfo(modifiedSet)
        else:
            copyInfoCallback(output)

        output.appendFromImages(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)

        if inputImages.hasObjId():
            self._defineTransformRelation(inputImages, output)

        # Define an informative summary of the subset operation
        sizeIn = inputImages.getSize()
        sizeOut = output.getSize()
        sizeDiff = sizeIn - sizeOut
        msg = 'A subset of _%s_ was created, ' % output.getClassName()
        msg += 'discarding *%d* items (%0.1f %%) from the input set.' % (sizeDiff, sizeDiff*100./sizeIn)
        self.summaryVar.set(msg)

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
                    return self._createRepresentativesFromClasses(inputClasses,
                                                                  outputClassName.split(',')[0])
            else:
                def callback(output):
                    self._copyInfoAndSetAlignment(inputClasses, output)

                return self._createSubSetFromImages(inputClasses.getImages(),
                                                    copyInfoCallback=callback)

        elif outputClassName.startswith('SetOfClasses'):
            return self._createClassesFromClasses(inputClasses)
        else:
            raise Exception("Unrecognized output type: '%s'" % outputClassName)  
    
    def _createMicsSubSetFromCTF(self, inputCTFs):
        """ Create a subset of Micrographs analyzing the CTFs. """
        outputMics = self._createSetOfMicrographs()
        setOfMics = inputCTFs.getMicrographs()
        if setOfMics is None:
            raise Exception('Could not create SetOfMicrographs subset from'
                            'this SetOfCTF, the micrographs were not set.')
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
        self._defineTransformRelation(inputCtf, setOfCtf)
        return setOfCtf
        
    def _createRepresentativesFromClasses(self, inputClasses, outputClassName):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        inputImages = inputClasses.getImages()
        createFunc = getattr(self, '_create' + outputClassName)
        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating REPRESENTATIVES of images from classes,  sqlite file: %s" % self._dbName)

        count = 0
        output = createFunc()
        output.copyInfo(inputImages)
        for cls in modifiedSet:
            if cls.isEnabled():
                img = cls.getRepresentative()
                img.copyObjId(cls)
                output.append(img)
                count += 1
        # Register outputs
        self._defineOutput('Representatives', output)
        if inputClasses.hasObjId():
            self._defineSourceRelation(inputClasses, output)
        else:
            self._defineSourceRelation(inputImages, output)
        
        selectmsg = 'we selected %s items' % count if count > 1 else 'was selected 1 item'
        msg = 'From input %s of size %s %s to create output %s'%(inputClasses.getClassName(), 
                                                                 inputClasses.getSize(), 
                                                                 selectmsg, 
                                                                 output.getClassName())
        self.summaryVar.set(msg)
        return output


    def _copyInfoAndSetAlignment(self, inputClasses, output):
        """ This method is used when creating subset of images from classes.
        We need to copy the information from the original classes images
        and also set the proper alignment contained in the classes.
        """
        inputImages = inputClasses.getImages()
        # Copy all info form the original 'classified' images
        output.copyInfo(inputImages)
        # Take the alignment of the first class
        cls = inputClasses.getFirstItem()
        output.setAlignment(cls.getAlignment())

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
        self._copyInfoAndSetAlignment(inputClasses, output)
        output.appendFromClasses(modifiedSet)
        # Register outputs
        self._defineOutput(className, output)
        if inputClasses.hasObjId():
            self._defineSourceRelation(inputClasses, output)
        self._defineTransformRelation(inputImages, output)
        count = len([cls for cls in modifiedSet if cls.isEnabled()])
        selectmsg = 'we selected %s items' % count if count > 1 else 'was selected 1 item'
        msg = 'From input %s of size %s %s to create output %s of size %s'%(inputClasses.getClassName(), 
                                                                            inputClasses.getSize(),  
                                                                            selectmsg, 
                                                                            output.getClassName(), 
                                                                            output.getSize())
        self.summaryVar.set(msg)
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
        if inputClasses.hasObjId():
            self._defineTransformRelation(inputClasses, output)
        else:
            self._defineSourceRelation(inputClasses.getImages(), output)
        count = len([cls for cls in modifiedSet if cls.isEnabled()])
        selectmsg = 'we selected %s items' % count if count > 1 else 'was selected 1 item'
        msg = 'From input %s of size %s %s to create output %s'%(inputClasses.getClassName(), inputClasses.getSize(),  selectmsg, output.getClassName())
        self.summaryVar.set(msg)
        return output
        
    def _createSubSetFromMicrographsTiltPair(self, micrographsTiltPair):
        """ Create a subset of Micrographs Tilt Pair. """
        output = MicrographsTiltPair(filename=self._getPath('micrographs_pairs.sqlite'))
        modifiedSet = MicrographsTiltPair(filename=self._dbName, prefix=self._dbPrefix)
        inputU = micrographsTiltPair.getUntilted()
        inputT = micrographsTiltPair.getTilted()
        outputU = SetOfMicrographs(filename=self._getPath('mics_untilted.sqlite'))
        outputT = SetOfMicrographs(filename=self._getPath('mics_tilted.sqlite'))
        outputU.copyInfo(inputU)
        outputT.copyInfo(inputT)

        for micPair, u, t in izip(modifiedSet, inputU, inputT):
            if micPair.isEnabled():
                output.append(micPair)
                outputU.append(u)
                outputT.append(t)
        output.setUntilted(outputU)
        output.setTilted(outputT)
        outputU.write()
        outputT.write()
        # Register outputs
        outputDict = {'outputMicrographsTiltPair': output}
        self._defineOutputs(**outputDict)
        self._defineTransformRelation(micrographsTiltPair, output)
        return output

    def _createSubSetFromParticlesTiltPair(self, particlesTiltPair):
        """ Create a subset of Particles Tilt Pair. """
        output = ParticlesTiltPair(filename=self._getPath('particles_pairs.sqlite'))
        
        inputU = particlesTiltPair.getUntilted()
        inputT = particlesTiltPair.getTilted()
        outputU = SetOfParticles(filename=self._getPath('particles_untilted.sqlite'))
        outputT = SetOfParticles(filename=self._getPath('particles_tilted.sqlite'))
        outputU.copyInfo(inputU)
        outputT.copyInfo(inputT)

        modifiedSet = ParticlesTiltPair(filename=self._dbName, prefix=self._dbPrefix)

        for pair, u, t in izip(modifiedSet, inputU, inputT):
            if pair.isEnabled():
                output.append(pair)
                outputU.append(u)
                outputT.append(t)
        # Register outputs
        output.setUntilted(outputU)
        output.setTilted(outputT)
        outputU.write()
        outputT.write()
        
        outputDict = {'outputParticlesTiltPair': output}
        self._defineOutputs(**outputDict)
        self._defineTransformRelation(particlesTiltPair, output)
        return output

    def createSetObject(self):
        _dbName, self._dbPrefix = self.sqliteFile.get().split(',')
        self._dbName = self._getPath('subset.sqlite')
        os.rename(_dbName, self._dbName)

        if self._dbPrefix.endswith('_'):
            self._dbPrefix = self._dbPrefix[:-1]

        from pyworkflow.em import loadSetFromDb

        # Ignoring self._dbPrefix here, since we want to load
        # the top-level set in the sqlite file
        setObj = loadSetFromDb(self._dbName)

        return setObj

    def _summary(self):
        summary = []
        msg = self.summaryVar.get()
        if  msg is None:
            msg = self.getDefaultSummary()
        summary.append(msg)
        return summary

    def getDefaultSummary(self):
        inputStr = ''
        inputObj = self.inputObject.get()
        if inputObj is not None:
            inputStr += inputObj.getClassName()
            if isinstance(inputObj, EMSet):
                inputStr += ' of size %s' % inputObj.getSize()
        output = ''
        for _, attr in self.iterOutputAttributes(EMObject):
            output += attr.getClassName()
            if isinstance(attr, EMSet):
                output += ' of size %s' % attr.getSize()

        msg = 'From input %s created output %s ' % (inputStr, output)

        return msg

    def _methods(self):
        return self._summary()

    def _defineOutput(self, className, output):
        outputDict = {'output' + className.replace('SetOf', ''): output}
        self._defineOutputs(**outputDict) 
        

class ProtCreateMask(BatchProtocol):
    
    _label='create mask'

    def _defineParams(self, form):
        form.addHidden('inputObj', PointerParam, pointerClass='EMObject')
        form.addHidden('maskFile', StringParam)

    def _insertAllSteps(self):
        self._insertFunctionStep('createMaskStep')

    def createMaskStep(self):
        inputObj = self.inputObj.get()
        maskSrc=self.maskFile.get()
        basename = os.path.basename(maskSrc)
        maskDst = self._getPath(basename)
        moveFile(maskSrc, maskDst)
        samplingRate = None
        if(hasattr(inputObj, "getSamplingRate")):
            samplingRate = inputObj.getSamplingRate()
        else:
            for key, attr in inputObj.iterInputAttributes():
                if hasattr(attr.get(), "getSamplingRate"):
                    samplingRate = attr.get().getSamplingRate()
        if  not samplingRate:
            raise Exception("sampling rate required")
        
        mask = Mask()
        mask.setFileName(maskDst)
        mask.setSamplingRate(samplingRate)
        self._defineOutputs(outputMask=mask)
        self._defineSourceRelation(self.inputObj, self.outputMask)

    def _summary(self):
        summary = []
        summary.append('From input %s created mask %s'%(self.getObjectTag("inputObj"), self.getObjectTag("outputMask")))
        return summary
        
    def _methods(self):
        return self._summary()



class ProtCreateFSC(BatchProtocol):

    _label='create fsc'

    def _defineParams(self, form):
        pass
        form.addHidden('inputObj', PointerParam, pointerClass='EMObject')

    def setInputObj(self, obj):
        self.inputObj.set(obj)

    def _summary(self):
        summary = []
        summary.append('From input %s created fsc %s'%(self.getObjectTag("inputObj"), self.getObjectTag("outputMask")))
        return summary

    def _methods(self):
        return self._summary()



