# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
from pyworkflow.utils import getFloatListFromValues, getBoolListFromValues
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
    
    def executeCtfGroupsStep(self):
        makePath(self.ctfGroupDirectory)
        self._log.info("Created CTF directory: '%s'" % self.ctfGroupDirectory)
    #     printLog("executeCtfGroups01"+ CTFDatName, _log) FIXME: print in log this line
    
        if not self.doCTFCorrection:
            md = xmipp.MetaData(self.selFileName)
            block_name = self._getBlockFileName(ctfBlockName, 1, self._getFileName('imageCTFpairs'))
            md.write(block_name)
            self._log.info("Written a single CTF group to file: '%s'" % block_name)
            self.numberOfCtfGroups.set(1)
        else:
            raise Exception("Ctf groupping not implemented yet")
            # TODO: update self.numberOfCtfGroups with the appropiated value
        self._store(self.numberOfCtfGroups)
    
    def angularProjectLibraryStep(self, iterN, refN, args, stepParams, **kwargs):
        runAngularProjectLibraryStep(self, iterN, refN, args, stepParams, **kwargs)
        
    def initAngularReferenceFileStep(self):
        '''Create Initial angular file. Either fill it with zeros or copy input'''
        #NOTE: if using angles, self.selFileName file should contain angles info
        md = xmipp.MetaData(self.selFileName) 
        
        # Ensure this labels are always 
        md.addLabel(xmipp.MDL_ANGLE_ROT)
        md.addLabel(xmipp.MDL_ANGLE_TILT)
        md.addLabel(xmipp.MDL_ANGLE_PSI)
        
        expImages = self._getFileName('inputParticlesDoc')
        ctfImages = self._getFileName('imageCTFpairs')
        
        md.write(self._getExpImagesFileName(expImages))
        blocklist = xmipp.getBlocksInMetaDataFile(ctfImages)
        
        mdCtf = xmipp.MetaData()
        mdAux = xmipp.MetaData()
        readLabels = [xmipp.MDL_ITEM_ID, xmipp.MDL_IMAGE]
        
        for block in blocklist:
            #read ctf block from ctf file
            mdCtf.read(block + '@' + ctfImages, readLabels)
            #add ctf columns to images file
            mdAux.joinNatural(md, mdCtf)
            # write block in images file with ctf info
            mdCtf.write(block + '@' + expImages, xmipp.MD_APPEND)
            
        return [expImages]
    
    def projectionMatchingStep(self, iterN, refN, args):
        runProjectionMatching(self, iterN, refN, args)
    
    def assignImagesToReferencesStep(self, iterN):
        runAssignImagesToReferences(self, iterN)
        
    def cleanVolumeStep(self, vol1, vol2):
        cleanPath(vol1, vol2)
    
    def reconstructionStep(self, iterN, refN, program, method, args, mpi, threads, **kwargs):
        runReconstructionStep(self, iterN, refN, program, method, args, mpi, threads, **kwargs)
    
    def storeResolutionStep(self, resolIterMd, resolIterMaxMd, sampling):
        runStoreResolutionStep(self, resolIterMd, resolIterMaxMd, sampling)
    
    def calculateFscStep(self, iterN, refN, args, constantToAdd, **kwargs):
        runCalculateFscStep(self, iterN, refN, args, constantToAdd, **kwargs)
    
    def filterVolumeStep(self, iterN, refN, constantToAddToFiltration, **kwargs):
        runFilterVolumeStep(self, iterN, refN, constantToAddToFiltration, **kwargs)
    
    def createOutputStep(self):
        print "output generated..........."


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
    
    