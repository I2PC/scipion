# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This sub-package implement projection matching using xmipp 3.1
"""

import xmipp
from pyworkflow.object import Integer
from pyworkflow.utils.path import makePath, copyFile
from pyworkflow.utils import getFloatListFromValues, getBoolListFromValues, getStringListFromValues
from pyworkflow.em import ProtRefine3D, ProtClassify3D

from projmatch_initialize import *
from projmatch_form import _defineProjectionMatchingParams
from projmatch_steps import *



class XmippProtProjMatch(ProtRefine3D, ProtClassify3D):
    """ 3D reconstruction and classification using multireference projection matching"""

    _label = 'projection matching'
    
    FILENAMENUMBERLENGTH = 6

    def __init__(self, **args):        
        ProtRefine3D.__init__(self, **args)
        ProtClassify3D.__init__(self, **args)
        self.numberOfCtfGroups = Integer(1)
        self._lastIter = Integer(0)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._loadInputInfo()
        # Setup the dictionary with filenames templates to 
        # be used by _getFileName
        createFilenameTemplates(self)
        # Load the values from several params generating a list
        # of values per iteration or references
        initializeLists(self)

    def _loadInputInfo(self):
        from pyworkflow.em.packages.xmipp3 import getImageLocation
        
        reference = self.input3DReferences.get() # Input can be either a single volume or a set of volumes.
        
        if isinstance(reference, Volume): # Treat the case of a single volume
            self.referenceFileNames = [getImageLocation(reference)]
        else:
            self.referenceFileNames = [getImageLocation(vol) for vol in reference]
            
        self.numberOfReferences = len(self.referenceFileNames)
        self.resolSam = reference.getSamplingRate()

        
    #--------------------------- DEFINE param functions --------------------------------------------   
        
    def _defineParams(self, form):
        """ Since the form definition is very very large,
        we have do it in a separated function.
        """
        _defineProjectionMatchingParams(self, form)
         
         
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):
        self._initialize()
        # Insert initial steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('executeCtfGroupsStep')
#         insertExecuteCtfGroupsStep(self)
#         insertInitAngularReferenceFileStep(self)
        self._insertFunctionStep('initAngularReferenceFileStep')
        # Steps per iteration
        self._insertItersSteps()
        # Final steps
        self._insertFunctionStep('createOutputStep')
        
    def _insertItersSteps(self):
        """ Insert several steps needed per iteration. """
        
        for iterN in self.allIters():
            dirsStep = self._insertFunctionStep('createIterDirsStep', iterN)
            # Insert some steps per reference volume
            projMatchSteps = []
            for refN in self.allRefs():
                # Mask the references in the iteration
                insertMaskReferenceStep(self, iterN, refN, prerequisites=[dirsStep])
                # Create the library of projections
                insertAngularProjectLibraryStep(self, iterN, refN)
                # Projection matching steps
                projMatchStep = self._insertProjectionMatchingStep(iterN, refN)
                projMatchSteps.append(projMatchStep)
                
            # Select the reference that best fits each image
            self._insertFunctionStep('assignImagesToReferencesStep', iterN, prerequisites=projMatchSteps)
            
            insertAngularClassAverageStep(self, iterN, refN)
    
            # Reconstruct each reference with new averages
            for refN in self.allRefs():
                # Create new class averages with images assigned
                insertReconstructionStep(self, iterN, refN)
                
                if self.doComputeResolution and self._doSplitReferenceImages[iterN]:
                    # Reconstruct two halves of the data
                    insertReconstructionStep(self, iterN, refN, 'Split1')
                    insertReconstructionStep(self, iterN, refN, 'Split2')
                    # Compute the resolution
                    insertComputeResolutionStep(self, iterN, refN)
                    
                insertFilterVolumeStep(self, iterN, refN)

            # Calculate both angles and shifts devitations for this iteration
            self._insertFunctionStep('calculateDeviationsStep', iterN)

    
    def _insertProjectionMatchingStep(self, iterN, refN):
        args = getProjectionMatchingArgs(self, iterN)
        return self._insertFunctionStep('projectionMatchingStep', iterN, refN, args)
    
    #--------------------------- STEPS functions --------------------------------------------       

    def convertInputStep(self):
        """ Generated the input particles metadata expected 
        by projection matching. And copy the generated file to be
        used as initial docfile for further iterations.
        """
        from pyworkflow.em.packages.xmipp3 import writeSetOfParticles
        writeSetOfParticles(self.inputParticles.get(), self.selFileName, 
                            blockName=self.blockWithAllExpImages)
        #copyFile(self.selFileName, self._getFileName('inputParticlesDoc'))
        
    def createIterDirsStep(self, iterN):
        """ Create the necessary directory for a given iteration. """
        iterDirs = [self._getFileName(k, iter=iterN) for k in ['iterDir', 'projMatchDirs', 'libraryDirs']]
    
        for d in iterDirs:
            makePath(d)
            
        return iterDirs
    
    def volumeConvertStep(self, reconstructedFilteredVolume, maskedFileName):
        runVolumeConvertStep(self, reconstructedFilteredVolume, maskedFileName)
    
    def executeCtfGroupsStep(self, **kwargs):
        runExecuteCtfGroupsStep(self, **kwargs)
    
    def transformMaskStep(self, program, args, **kwargs):
        runTransformMaskStep(self, program, args, **kwargs)
    
    def angularProjectLibraryStep(self, iterN, refN, args, stepParams, **kwargs):
        runAngularProjectLibraryStep(self, iterN, refN, args, stepParams, **kwargs)
        
    def initAngularReferenceFileStep(self):
        runInitAngularReferenceFileStep(self)
        
    def projectionMatchingStep(self, iterN, refN, args):
        runProjectionMatching(self, iterN, refN, args)
    
    def assignImagesToReferencesStep(self, iterN):
        runAssignImagesToReferences(self, iterN)
        
    def cleanVolumeStep(self, vol1, vol2):
        cleanPath(vol1, vol2)
    
    def reconstructionStep(self, iterN, refN, program, method, args, suffix, **kwargs):
        runReconstructionStep(self, iterN, refN, program, method, args, suffix, **kwargs)
    
    def storeResolutionStep(self, resolIterMd, resolIterMaxMd, sampling):
        runStoreResolutionStep(self, resolIterMd, resolIterMaxMd, sampling)
    
    def calculateFscStep(self, iterN, refN, args, constantToAdd, **kwargs):
        runCalculateFscStep(self, iterN, refN, args, constantToAdd, **kwargs)
    
    def filterVolumeStep(self, iterN, refN, constantToAddToFiltration, **kwargs):
        runFilterVolumeStep(self, iterN, refN, constantToAddToFiltration, **kwargs)
    
    def createOutputStep(self):
        runCreateOutputStep(self)

    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        if self.doCTFCorrection and not self.doAutoCTFGroup and not exists(self.setOfDefocus.get()):
            errors.append("Error: for non-automated ctf grouping, please provide a docfile!")
        if self.numberOfMpi<=1:
            errors.append("The number of MPI processes has to be larger than 1")
        
        xDimImg = self.inputParticles.get().getXDim()
        xDimVol, _, _ = self.input3DReferences.get().getDim()
        if xDimImg != xDimVol:
            errors.append("The dimensions of the volume(s) and particles must be equal!!!!")
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
    
    def allIters(self):
        """ Iterate over all iterations. """
        for i in range(1, self.numberOfIterations.get()+1):
            yield i
            
    def allRefs(self):
        """ Iterate over all references. """
        for i in range(1, self.numberOfReferences+1):
            yield i
            
    def allCtfGroups(self):
        """ Iterate over all CTF groups. """
        for i in range(1, self.numberOfCtfGroups.get() + 1):
            yield i
            
    def itersFloatValues(self, attributeName, firstValue=-1):
        """ Take the string of a given attribute and
        create a list of floats that will be used by 
        the iteratioins. An special first value will be
        added to the list for iteration 0.
        """
        valuesStr = self.getAttributeValue(attributeName)
        if valuesStr is None:
            raise Exception('None value for attribute: %s' % attributeName)
        return [firstValue] + getFloatListFromValues(valuesStr, length=self.numberOfIterations.get())
    
    def itersBoolValues(self, attributeName, firstValue=False):
        """ Take the string of a given attribute and
        create a list of booleans that will be used by 
        the iteratioins. An special first value will be
        added to the list for iteration 0.
        """
        valuesStr = self.getAttributeValue(attributeName)
        if valuesStr is None:
            raise Exception('None value for attribute: %s' % attributeName)
        return [firstValue] + getBoolListFromValues(valuesStr, length=self.numberOfIterations.get())
        
    def itersStringValues(self, attributeName, firstValue='c1'):
        """ Take the string of a given attribute and
        create a list of strings that will be used by 
        the iteratioins. An special first value will be
        added to the list for iteration 0.
        """
        valuesStr = self.getAttributeValue(attributeName)
        if valuesStr is None:
            raise Exception('None value for attribute: %s' % attributeName)
        return [firstValue] + getStringListFromValues(valuesStr, length=self.numberOfIterations.get())
        
    def _getBlockFileName(self, blockName, blockNumber, filename, length=None):
        l = length or self.FILENAMENUMBERLENGTH
        
        return blockName + str(blockNumber).zfill(l) + '@' + filename
    
    def _getExpImagesFileName(self, filename):
        return self.blockWithAllExpImages + '@' + filename
    
    def _getRefBlockFileName(self, ctfBlName, ctfBlNumber, refBlName, refBlNumber, filename, length=None):
        l = length or self.FILENAMENUMBERLENGTH
        
        return ctfBlName + str(ctfBlNumber).zfill(l) + '_' + refBlName + str(refBlNumber).zfill(l) + '@' + filename

    def _getFourierMaxFrequencyOfInterest(self, iterN, refN):
        """ Read the corresponding resolution metadata and return the
        desired resolution.
        """
        md = xmipp.MetaData(self._getFileName('resolutionXmdMax', iter=iterN, ref=refN))
        return md.getValue(xmipp.MDL_RESOLUTION_FREQREAL, md.firstObject())
    
    def calculateDeviationsStep(self, it):
        """ Calculate both angles and shifts devitations for all iterations
        """
    
        SL = xmipp.SymList()
        mdIter = xmipp.MetaData()
        #for it in self.allIters():
        mdIter.clear()
        SL.readSymmetryFile(self._symmetry[it])
        md1 = xmipp.MetaData(self.docFileInputAngles[it])
        md2 = xmipp.MetaData(self.docFileInputAngles[it-1])
        #ignore disabled,
        md1.removeDisabled()
        md2.removeDisabled()

        #first metadata file may not have shiftx and shifty
        if not md2.containsLabel(xmipp.MDL_SHIFT_X):
            md2.addLabel(xmipp.MDL_SHIFT_X)
            md2.addLabel(xmipp.MDL_SHIFT_Y)
            md2.fillConstant(xmipp.MDL_SHIFT_X,0.)
            md2.fillConstant(xmipp.MDL_SHIFT_Y,0.)
        oldLabels=[xmipp.MDL_ANGLE_ROT,
                   xmipp.MDL_ANGLE_TILT,
                   xmipp.MDL_ANGLE_PSI,
                   xmipp.MDL_SHIFT_X,
                   xmipp.MDL_SHIFT_Y]
        newLabels=[xmipp.MDL_ANGLE_ROT2,
                   xmipp.MDL_ANGLE_TILT2,
                   xmipp.MDL_ANGLE_PSI2,
                   xmipp.MDL_SHIFT_X2,
                   xmipp.MDL_SHIFT_Y2]
        md2.renameColumn(oldLabels,newLabels)
        md2.addLabel(xmipp.MDL_SHIFT_X_DIFF)
        md2.addLabel(xmipp.MDL_SHIFT_Y_DIFF)
        md2.addLabel(xmipp.MDL_SHIFT_DIFF)
        mdIter.join1(md1, md2, xmipp.MDL_IMAGE, xmipp.INNER_JOIN)
        SL.computeDistance(mdIter,False,False,False)
        xmipp.activateMathExtensions()
        #operate in sqlite
        shiftXLabel     = xmipp.label2Str(xmipp.MDL_SHIFT_X)
        shiftX2Label    = xmipp.label2Str(xmipp.MDL_SHIFT_X2)
        shiftXDiff      = xmipp.label2Str(xmipp.MDL_SHIFT_X_DIFF)
        shiftYLabel     = xmipp.label2Str(xmipp.MDL_SHIFT_Y)
        shiftY2Label    = xmipp.label2Str(xmipp.MDL_SHIFT_Y2)
        shiftYDiff      = xmipp.label2Str(xmipp.MDL_SHIFT_Y_DIFF)
        shiftDiff       = xmipp.label2Str(xmipp.MDL_SHIFT_DIFF)
        #timeStr = str(dtBegin)
        operateString   =       shiftXDiff+"="+shiftXLabel+"-"+shiftX2Label
        operateString  += "," + shiftYDiff+"="+shiftYLabel+"-"+shiftY2Label
        mdIter.operate(operateString)
        operateString  =  shiftDiff+"=sqrt(" \
                          +shiftXDiff+"*"+shiftXDiff+"+" \
                          +shiftYDiff+"*"+shiftYDiff+");"
        mdIter.operate(operateString)
        iterFile = self._mdDevitationsFn(it)
        mdIter.write(iterFile,xmipp.MD_APPEND)

        self._setLastIter(it)

    
    def _mdDevitationsFn(self, it):
        mdFn = self._getPath('deviations.xmd')
        return "iter_%03d@" % it + mdFn

    def _setLastIter(self, iterN):
        self._lastIter.set(iterN)
        self._store(self._lastIter)

    def getLastIter(self):
        return self._lastIter.get()
    
    def _fillParticlesFromIter(self, partSet, iteration):
        print("_fillParticlesFromIter")
        import pyworkflow.em.metadata as md
        
        imgSet = self.inputParticles.get()
        imgFn = "all_exp_images@" + self._getFileName('docfileInputAnglesIters', iter=iteration, ref=1)
        partSet.copyInfo(imgSet)
        partSet.setAlignmentProj()
        
        partSet.copyItems(imgSet,
                            updateItemCallback=self._createItemMatrix,
                            itemDataIterator=md.iterRows(imgFn, sortByLabel=md.MDL_ITEM_ID))
    
    def _createItemMatrix(self, item, row):
        from pyworkflow.em.packages.xmipp3.convert import createItemMatrix
        import pyworkflow.em as em
        
        createItemMatrix(item, row, align=em.ALIGN_PROJ)
    
    def _getIterParticles(self, it, clean=False):
        import pyworkflow.em as em
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by 
        converting from this iteration data.star file.
        """
        
        dataParticles = self._getFileName('particlesScipion', iter=it)
        
        if clean:
            cleanPath(dataParticles)
            
        if not exists(dataParticles):
            partSet = em.SetOfParticles(filename=dataParticles)
            self._fillParticlesFromIter(partSet, it)
            partSet.write()
            partSet.close()
        else:
            partSet = em.SetOfParticles(filename=dataParticles)
            imgSet = self.inputParticles.get()
            partSet.copyInfo(imgSet)
            partSet.setAlignmentProj()

        return partSet

