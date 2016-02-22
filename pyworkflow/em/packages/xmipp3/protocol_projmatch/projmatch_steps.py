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
Since the Projection Matching protocol of Xmipp 3 has a very large
form definition, we have separated in this sub-module.
"""

import math
from os.path import exists, join

import xmipp
from pyworkflow.object import Float
from pyworkflow.em.data import Volume, SetOfClasses3D
from pyworkflow.utils import getMemoryAvailable, replaceExt, removeExt, cleanPath, makePath, copyFile

from pyworkflow.em.packages.xmipp3.convert import createClassesFromImages
from pyworkflow.em.packages.xmipp3.utils import isMdEmpty


ctfBlockName = 'ctfGroup'
refBlockName = 'refGroup'


def runExecuteCtfGroupsStep(self, **kwargs):
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
        self._log.info('*********************************************************************')
        self._log.info('* Make CTF groups')
        
    #    remove all entries not present in sel file by
    #    join between selfile and metadatafile
        mdCtfData = xmipp.MetaData()
        mdCtfData.read(self.ctfDatName)
    
        mdSel = xmipp.MetaData();
        mdSel.read(self.selFileName)
        mdCtfData.intersection(mdSel, xmipp.MDL_IMAGE)
        tmpCtfDat = self.ctfDatName
        mdCtfData.write(tmpCtfDat)
        args = ' --ctfdat %(tmpCtfDat)s -o %(ctffile)s --wiener --wc %(wiener)s --pad %(pad)s'
        args += ' --sampling_rate %(sampling)s'
        
        params = {'tmpCtfDat' : tmpCtfDat,
                  'ctffile' : self._getFileName('ctfGroupBase') + ':stk',
                  'wiener' : self.wienerConstant.get(),
                  'pad' : self.paddingFactor.get(),
                  'sampling' : self.resolSam
                  }
        
        if self.inputParticles.get().isPhaseFlipped():
            args += ' --phase_flipped '
    
        if self.doAutoCTFGroup:
            args += ' --error %(ctfGroupMaxDiff)s --resol %(ctfGroupMaxResol)s'
            params['ctfGroupMaxDiff'] = self.ctfGroupMaxDiff.get()
            params['ctfGroupMaxResol'] = self.ctfGroupMaxResol.get()
        else:
            if exists(self.setOfDefocus.get()):
                args += ' --split %(setOfDefocus)s'  
                params['setOfDefocus'] = self.setOfDefocus.get()
        
        self.runJob("xmipp_ctf_group", args % params, numberOfMpi=1, **kwargs)
        
        auxMD = xmipp.MetaData("numberGroups@" + self._getFileName('cTFGroupSummary'))
        self.numberOfCtfGroups.set(auxMD.getValue(xmipp.MDL_COUNT, auxMD.firstObject()))
    
    self._store(self.numberOfCtfGroups)


# # Functions outside th loop loop for xmipp_projection_matching
def runInitAngularReferenceFileStep(self):
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


# Functions in loop for xmipp_projection_matching
def insertMaskReferenceStep(self, iterN, refN, **kwargs):
    maskRadius = self.maskRadius.get()
    if maskRadius < 0:
        maskRadius, _, _ = self.input3DReferences.get().getDim()
        maskRadius = maskRadius / 2
    maskedFileName = self._getFileName('maskedFileNamesIters', iter=iterN, ref=refN)
    reconstructedFilteredVolume = self.reconstructedFilteredFileNamesIters[iterN-1][refN]
    
    if self.getEnumText('maskType') != 'None':
        
        args = ' -i %(reconstructedFilteredVolume)s -o %(maskedFileName)s'
        if self.getEnumText('maskType') == 'circular':
            args += ' --mask circular -%(maskRadius)s'
            maskProgram = "xmipp_transform_mask"
        else:
            maskFn = self.inputMask.get().getFileName()
            args += ' --mult %(maskFn)s'
            maskProgram = "xmipp_image_operate"
        
        # Here is used _insertFunctionStep instead of _insertRunJobStep cause xmipp_transform_mask is not implemented with mpi
        self._insertFunctionStep('transformMaskStep', maskProgram, args % locals(), **kwargs)
    else:
        self._insertFunctionStep('volumeConvertStep', reconstructedFilteredVolume, maskedFileName)


def runTransformMaskStep(self, program, args, **kwargs):
    self.runJob(program, args, numberOfMpi=1, **kwargs)


def runVolumeConvertStep(self, reconstructedFilteredVolume, maskedFileName):
    from pyworkflow.em.convert import ImageHandler
    img = ImageHandler()
    img.convert(reconstructedFilteredVolume, maskedFileName)


def insertAngularProjectLibraryStep(self, iterN, refN, **kwargs):
    args =  ' -i %(maskedFileNamesIter)s --experimental_images %(experimentalImages)s'
    args += ' -o %(projectLibraryRootName)s --sampling_rate %(samplingRate)s --sym %(symmetry)s'
    args += 'h --compute_neighbors' 

    ###need one block per reference
    # Project all references
    
    xDim, yDim, zDim = self.input3DReferences.get().getDim()
    memoryUsed = (xDim * yDim * zDim * 8) / pow(2,20)
    if memoryUsed == 0:
        memoryUsed = 1 # If this value is 0, produce an division error in runAngularProjectLibraryStep
    
    stepParams = {'method' : self.getEnumText('projectionMethod')}
    
    expImages = self._getExpImagesFileName(self.docFileInputAngles[iterN-1])
    projectLibraryRootName = self._getFileName('projectLibraryStk', iter=iterN, ref=refN)
    
    params = {'maskedFileNamesIter' : self._getFileName('maskedFileNamesIters', iter=iterN, ref=refN),
              'experimentalImages' : expImages,
              'projectLibraryRootName' : projectLibraryRootName,
              'samplingRate' : self._angSamplingRateDeg[iterN],
              'symmetry' : self._symmetry[iterN],
              }
        
    if self.maxChangeInAngles < 181:
        args += ' --near_exp_data --angular_distance %(maxChangeInAngles)s'
    else:
        args += ' --angular_distance -1'
    
    if self._perturbProjectionDirections[iterN]:
        args +=' --perturb %(perturb)s'
        params['perturb'] = math.sin(math.radians(self._angSamplingRateDeg[iterN])) / 4.

    if self.doRestricSearchbyTiltAngle:
        args += ' --min_tilt_angle %(tilt0)s --max_tilt_angle %(tiltF)s'
        params['tilt0'] = self.tilt0.get()
        params['tiltF'] = self.tiltF.get()
 
    if self.doCTFCorrection:
        params['ctfGroupSubsetFileName'] = self._getFileName('imageCTFpairs')
        args += ' --groups %(ctfGroupSubsetFileName)s'

    if len(self.symmetryGroupNeighbourhood.get()) > 1:
        params['symmetryGroupNeighbourhood'] = self.symmetryGroupNeighbourhood.get()
        args += ' --sym_neigh %(symmetryGroupNeighbourhood)s'

    if self._onlyWinner[iterN]:
        args += ' --only_winner'
    
    if stepParams['method'] == 'fourier':
        memoryUsed = memoryUsed * 6
        stepParams['paddingAngularProjection'] = self.paddingAngularProjection.get()
        stepParams['kernelAngularProjection'] = self.getEnumText('kernelAngularProjection')
        stepParams['constantToAdd'] = self._constantToAddToFiltration[iterN]
    
    stepParams['memoryUsed'] = memoryUsed
        
    self._insertFunctionStep('angularProjectLibraryStep', iterN, refN, args % params, stepParams, **kwargs)

    if not self.doCTFCorrection:
        src = removeExt(projectLibraryRootName) + '_sampling.xmd'
        dst = removeExt(projectLibraryRootName) + ('_group%06d_sampling.xmd' % 1)
        self._insertCopyFileStep(src, dst)


def runAngularProjectLibraryStep(self, iterN, refN, args, stepParams, **kwargs):
    args += ' --method %(method)s'
    if stepParams['method'] == 'fourier':
        if self._fourierMaxFrequencyOfInterest[iterN] == -1:
            fourierMaxFrequencyOfInterest = self._getFourierMaxFrequencyOfInterest(iterN-1, refN)
            fourierMaxFrequencyOfInterest = self.resolSam / fourierMaxFrequencyOfInterest + stepParams['constantToAdd']
            
            if fourierMaxFrequencyOfInterest > 0.5:
                fourierMaxFrequencyOfInterest = 0.5
            elif fourierMaxFrequencyOfInterest < 0.:
                fourierMaxFrequencyOfInterest = 0.001
        else:
            fourierMaxFrequencyOfInterest = self._fourierMaxFrequencyOfInterest[iterN]
        
        stepParams['fourierMaxFrequencyOfInterest'] = fourierMaxFrequencyOfInterest
        args += ' %(paddingAngularProjection)s %(fourierMaxFrequencyOfInterest)s %(kernelAngularProjection)s'
    
    processorsToUse = self.numberOfMpi.get() * self.numberOfThreads.get()
    if processorsToUse > 1:
        memoryAvailable = getMemoryAvailable()
        processorsToUse = min(processorsToUse, memoryAvailable/stepParams['memoryUsed'])
    
    if self.numberOfMpi > 1 and processorsToUse > 1:
        stepParams['mpiJobSize'] = self.mpiJobSize.get()
        args += ' --mpi_job_size %(mpiJobSize)s'
    
    self._log.info('* Create projection library')
    self.runJob('xmipp_angular_project_library', args % stepParams, numberOfMpi=processorsToUse, **kwargs)


def getProjectionMatchingArgs(self, iterN):
    """ Get the arguments for the projection matching program that
    does not vary with the CTF groups. 
    """
    
    args = ' --Ri %(innerRadius)s'
    args += ' --Ro %(outerRadius)s --max_shift %(maxShift)s --search5d_shift %(search5dShift)s'
    args += ' --search5d_step %(search5dStep)s --append'
    
    params = {
              'innerRadius' : self._innerRadius[iterN],
              'outerRadius' : self._outerRadius[iterN],
              'maxShift' : self._maxChangeOffset[iterN],
              'search5dShift' : self._search5DShift[iterN],
              'search5dStep' : self._search5DStep[iterN],
              }
    
    if self.doScale:
        args += ' --scale %(scaleStep)s %(scaleNumberOfSteps)s'
        params['scaleStep'] = self._scaleStep[iterN]
        params['scaleNumberOfSteps'] = self._scaleNumberOfSteps[iterN]
        
    if self.doCTFCorrection and self._referenceIsCtfCorrected[iterN]:
        args += ' --pad %(pad)s'
        params['pad'] = self.paddingFactor.get()

    if self.numberOfMpi > 1:
        params['mpiJobSize'] = self.mpiJobSize.get()
        args += ' --mpi_job_size %(mpiJobSize)s'
        
    return args % params  
    
def runProjectionMatching(self, iterN, refN, args, **kwargs):
    """ Loop over all CTF groups and launch a projection matching for each one.
    Note: Iterate ctf groups in reverse order to have same order as 
          in add_to docfiles from angular_class_average. #FIXME: check why reverse order is needed
    """
    projMatchRootName = self._getFileName('projMatchRootNames', iter=iterN, ref=refN)
    refname = self._getFileName('projectLibraryStk', iter=iterN, ref=refN)
    
    numberOfCtfGroups = self.numberOfCtfGroups.get()
#     ctfGroupName = self._getPath(self.ctfGroupDirectory, '%(ctfGroupRootName)s')
    #remove output metadata
    cleanPath(projMatchRootName)
    
    for ctfN in reversed(list(self.allCtfGroups())):
        self._log.info('CTF group: %d/%d' % (ctfN, numberOfCtfGroups))
        ctfArgs = ' -i %(inputdocfile)s -o %(outputname)s --ref %(refname)s'        
        
        inputdocfile = self._getBlockFileName(ctfBlockName, ctfN, self.docFileInputAngles[iterN-1])
        outputname = self._getBlockFileName(ctfBlockName, ctfN, projMatchRootName)
        baseTxtFile = removeExt(refname)
        neighbFile = baseTxtFile + '_sampling.xmd'
        cleanPath(neighbFile)
        neighbFileb = baseTxtFile + '_group' + str(ctfN).zfill(self.FILENAMENUMBERLENGTH) + '_sampling.xmd'
        copyFile(neighbFileb, neighbFile)
        print "copied file ", neighbFileb, "to", neighbFile
        
        threads = self.numberOfThreads.get()
        trhArgs = ' --mem %(mem)s --thr %(thr)s'
        thrParams = {
                  'mem' : self.availableMemory.get() * threads,
                  'thr' : threads,
                  }
        
        if self.doCTFCorrection and self._referenceIsCtfCorrected[iterN]:
            ctfArgs += ' --ctf %s' % self._getBlockFileName('', ctfN, self._getFileName('stackCTFs'))
        
        progArgs = ctfArgs % locals() + args + trhArgs % thrParams
        self.runJob('xmipp_angular_projection_matching', progArgs, **kwargs)


def runAssignImagesToReferences(self, iterN, **kwargs):
    ''' assign the images to the different references based on the crosscorrelation coeficient
        #if only one reference it just copy the docfile generated in the previous step
        '''
    numberOfCtfGroups = self.numberOfCtfGroups.get()
    #first we need a list with the references used. That is,
    #read all docfiles and map referecendes to a mdl_order
    mdAux = xmipp.MetaData()
    mdSort = xmipp.MetaData()
    md = xmipp.MetaData()
    md1 = xmipp.MetaData()
    mdout = xmipp.MetaData()
    mdout.setComment("Metadata with images, the winner reference as well as the ctf group")
    
    mycounter = 1L
    for ctfN in self.allCtfGroups():
        ctfFilePrefix = self._getBlockFileName(ctfBlockName, ctfN, '')
        for refN in self.allRefs():
            projMatchRootName = self._getFileName('projMatchRootNames', iter=iterN, ref=refN)
            inputdocfile = ctfFilePrefix + projMatchRootName
            md.read(inputdocfile)
            for id in md:
                t = md.getValue(xmipp.MDL_REF, id)
                i = mdSort.addObject()
                mdSort.setValue(xmipp.MDL_REF, t, i)
    
    mdSort.removeDuplicates()
     
    for id in mdSort:
        mdSort.setValue(xmipp.MDL_ORDER, mycounter, id)
        mycounter += 1
    ####################
    outputdocfile = self.docFileInputAngles[iterN]
    cleanPath(outputdocfile)
         
    mdout2 = xmipp.MetaData()
    for ctfN in self.allCtfGroups():
        mdAux.clear()
        ctfFilePrefix = self._getBlockFileName(ctfBlockName, ctfN, '')
        for refN in self.allRefs():
            projMatchRootName = self._getFileName('projMatchRootNames', iter=iterN, ref=refN)
            inputdocfile = ctfFilePrefix + projMatchRootName
            md.clear()
            md.read(inputdocfile)
            #In practice you should not get duplicates
            md.removeDuplicates()
            md.setValueCol(xmipp.MDL_REF3D, refN)
            md.setValueCol(xmipp.MDL_DEFGROUP, ctfN)
            #MD.setValueCol(xmipp.MDL_CTF_MODEL,ctfFilePrefix[:-1])
            mdAux.unionAll(md)
        mdAux.sort()
        md.aggregate(mdAux, xmipp.AGGR_MAX, xmipp.MDL_IMAGE, xmipp.MDL_MAXCC, xmipp.MDL_MAXCC)
        #if a single image is assigned to two references with the same 
        #CC use it in both reconstruction
        #recover atribbutes after aggregate function
         
        md1.joinNatural(md, mdAux)
        mdout.joinNatural(md1, mdSort)
        mdout.write(ctfFilePrefix + outputdocfile, xmipp.MD_APPEND)
        mdout2.unionAll(mdout)
        
    mdout2.write(self.blockWithAllExpImages + '@' + outputdocfile, xmipp.MD_APPEND)
    #Aggregate for ref3D and print warning if all images has been assigned to a single volume
    #ROB
    md1.clear()
    md1.aggregate(mdout2, xmipp.AGGR_COUNT, xmipp.MDL_REF3D, xmipp.MDL_REF3D, xmipp.MDL_COUNT)
    import sys
    md1_size = md1.size()
    numberOfReferences = self.numberOfReferences
    if md1_size != numberOfReferences:
        print >> sys.stderr,"********************************************"
        print >> sys.stderr, md1
        print >> sys.stderr,"ERROR: Some 3D references do not have assigned any projection assigned to them"
        print >> sys.stderr,"Consider reducing the number of 3D references"
        print >> sys.stderr,"Number of References:" ,numberOfReferences
        print >> sys.stderr,"Number of Empty references", numberOfReferences - md1_size
        raise Exception('runAssignImagesToReferences failed')
 
 
    #!a original_angles too
     
    #we are done but for the future it is convenient to create more blocks
    #with the pairs ctf_group reference    
    for ctfN in self.allCtfGroups():
        ctfFilePrefix = self._getBlockFileName(ctfBlockName, ctfN, '')
        #print 'read file: ', ctfFilePrefix+outputdocfile
        mdAux.read(ctfFilePrefix + outputdocfile)
        for refN in self.allRefs():
            auxOutputdocfile = self._getRefBlockFileName(ctfBlockName, ctfN, refBlockName, refN, '')
            #select images with ref3d=iRef3D
            mdout.importObjects(mdAux, xmipp.MDValueEQ(xmipp.MDL_REF3D, refN))
            mdout.write(auxOutputdocfile + outputdocfile, xmipp.MD_APPEND)


def insertAngularClassAverageStep(self, iterN, refN, **kwargs):

 
    refname = self._getFileName('projectLibraryStk', iter=iterN, ref=refN)
    baseTxtFile = refname[:-len('.stk')] 
     
    docFileInputAngles  = self.docFileInputAngles[iterN]
    projLibraryDoc = self._getFileName('projectLibraryDoc', iter=iterN, ref=refN)
    outClasses = self._getFileName('outClasses', iter=iterN, ref=refN)
    # FIXME: Why is necessary ask if docFileInputAngles is empty. check if is a validation step
#     if xmipp.isMdEmpty(docFileInputAngles):
#         print "Empty metadata file: %s" % docFileInputAngles
#         return
    
    params = {'docFileInputAngles' : docFileInputAngles,
              'projLibraryDoc' : projLibraryDoc,
              'outClasses' : outClasses
              }
    
    args = ' -i ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]\$@'
    args += '%(docFileInputAngles)s --lib %(projLibraryDoc)s -o %(outClasses)s'
    
    # FIXME: This option no exist in the form
#     if self.doSaveImagesAssignedToClasses:
#         args += ' --save_images_assigned_to_classes'  
    
    if self.getEnumText('discardImages') == 'maxCC':
        args += ' --limit0 %(discard)s'
        params['discard'] = self._minimumCrossCorrelation[iterN]
    elif self.getEnumText('discardImages') == 'percentage':
        args += ' --limitRper %(discard)s'
        params['discard'] = self._discardPercentage[iterN]
    elif self.getEnumText('discardImages') == 'classPercentage':
        args += ' --limitRclass %(discard)s'
        params['discard'] = self._discardPercentagePerClass[iterN]
    #else 'none'

    
    # On-the fly apply Wiener-filter correction and add all CTF groups together
    if self.doCTFCorrection:
        args += ' --wien %(wien)s --pad %(pad)s'
        params['wien'] = self._getFileName('stackWienerFilters')
        params['pad'] = self.paddingFactor.get()
                    
    if self._doAlign2D[iterN]:
        args += ' --iter %(alignIter)s --Ri %(innerRadius)s --Ro %(outerRadius)s'
        params['alignIter'] = self._align2DIterNr[iterN]
        params['innerRadius'] = self._innerRadius[iterN]
        params['outerRadius'] = self._outerRadius[iterN]
        
    if self.doComputeResolution and self._doSplitReferenceImages[iterN]:
        args += ' --split'

    processorsToUse = self.numberOfMpi.get() * self.numberOfThreads.get()                   
    if self.numberOfMpi > 1:
        params['mpiJobSize'] = self.mpiJobSize.get()
        args += ' --mpi_job_size %(mpiJobSize)s'
 
    self._insertRunJobStep('xmipp_angular_class_average', args % params, processorsToUse, **kwargs)


def insertReconstructionStep(self, iterN, refN, suffix='', **kwargs):
    
    method = self.getEnumText('reconstructionMethod')
    reconsXmd = 'reconstructionXmd' + suffix
    reconsVol = 'reconstructedFileNamesIters' + suffix
    
    args = ' -i %(reconsXmd)s -o %(reconsVol)s --sym %(symmetry)s'
    params = {'reconsXmd' : self._getFileName(reconsXmd, iter=iterN, ref=refN),
              'reconsVol' : self._getFileName(reconsVol, iter=iterN, ref=refN),
              'symmetry' : self._symmetry[iterN]
              }
    
    if method == 'wbp':
        program = 'xmipp_reconstruct_wbp'
        args += ' --doc %(reconsXmd)s --weight --use_each_image ' + self.wbpReconstructionExtraCommand.get()
        
    elif method == 'art':
        program = 'xmipp_reconstruct_art'
        args += ' --WLS'
        
        if len(self._artLambda) >= 1:
            args += ' -l %(artLambda)s ' + self.artReconstructionExtraCommand.get()
            params['artLambda'] = self._artLambda[iterN]
    
    elif method == 'fourier':
        program = 'xmipp_reconstruct_fourier'
        args += ' --weight --padding %(pad)s %(pad)s'
        params['pad'] = self.paddingFactor.get()
        
    self._insertFunctionStep('reconstructionStep', iterN, refN, program, method, args % params, suffix, **kwargs)


def runReconstructionStep(self, iterN, refN, program, method, args, suffix, **kwargs):
    #if input metadata is empty create a Blanck image
    reconsXmd = 'reconstructionXmd' + suffix
    reconsVol = 'reconstructedFileNamesIters' + suffix
    mdFn = self._getFileName(reconsXmd, iter=iterN, ref=refN)
    volFn = self._getFileName(reconsVol, iter=iterN, ref=refN)
    maskFn = self._getFileName('maskedFileNamesIters', iter=iterN, ref=refN)
    if method=="art":
        mpi = 1
        threads = 1
    else:
        mpi = self.numberOfMpi.get()
        threads = self.numberOfThreads.get()
        args += ' --thr %d' % threads
    
    if isMdEmpty(mdFn):
        img = xmipp.Image()
        img.read(maskFn)
        #(x,y,z,n) = img.getDimensions()
        self._log.warning("Metadata '%s' is empty. \n Creating a random volume file '%s'" % (mdFn, volFn))
        #createEmptyFile(ReconstructedVolume,x,y,z,n)
        img.initRandom()
        img.write(volFn)
    else:
        if method == 'fourier':
            if self._fourierMaxFrequencyOfInterest[iterN] == -1:
                fourierMaxFrequencyOfInterest = self._getFourierMaxFrequencyOfInterest(iterN-1, refN)
                fourierMaxFrequencyOfInterest = self.resolSam / fourierMaxFrequencyOfInterest + self._constantToAddToMaxReconstructionFrequency[iterN]
                
                if fourierMaxFrequencyOfInterest > 0.5:
                    fourierMaxFrequencyOfInterest = 0.5
                elif fourierMaxFrequencyOfInterest < 0.:
                    fourierMaxFrequencyOfInterest = 0.001
            else:
                fourierMaxFrequencyOfInterest = self._fourierMaxFrequencyOfInterest[iterN]
            
            args += ' --max_resolution %s' % fourierMaxFrequencyOfInterest
    
        if mpi > 1:
            args += ' --mpi_job_size %s' % self.mpiJobSize.get()
        
        self._log.info('*********************************************************************')
        self._log.info('* Reconstruct volume using %s' % method)
        self.runJob( program, args, numberOfMpi=mpi, numberOfThreads=threads, **kwargs)


def insertComputeResolutionStep(self, iterN, refN, **kwargs):
    vol1 = self._getFileName('reconstructedFileNamesItersSplit1', iter=iterN, ref=refN)
    vol2 = self._getFileName('reconstructedFileNamesItersSplit2', iter=iterN, ref=refN)
    resolIterMd = self._getFileName('resolutionXmd', iter=iterN, ref=refN)
    resolIterMaxMd = self._getFileName('resolutionXmdMax', iter=iterN, ref=refN)
    samplingRate = self.resolSam
    resolutionXmdCurrIter = self._getFileName('resolutionXmd', iter=iterN, ref=refN)
    # Prevent high-resolution correlation because of discrete mask from wbp
    outRadius = self._outerRadius[iterN]
    if outRadius < 0:
        outRadius, _, _ = self.input3DReferences.get().getDim()
        outRadius = outRadius / 2
    innerRadius = outRadius - 2
    outputVolumes = [vol1, vol2]
    for vol in outputVolumes:
        args = ' -i %(vol)s --mask  raised_cosine -%(innerRadius)s -%(outRadius)s'
        self._insertFunctionStep('transformMaskStep', "xmipp_transform_mask", args % locals(), **kwargs)
    
    args = ' --ref %(vol1)s -i %(vol2)s --sampling_rate %(samplingRate)s -o %(resolutionXmdCurrIter)s' % locals()
    
    self._insertFunctionStep('calculateFscStep', iterN, refN, args, self._constantToAddToMaxReconstructionFrequency[iterN], **kwargs)
    self._insertFunctionStep('storeResolutionStep', resolIterMd, resolIterMaxMd, samplingRate)
    #if cleanup=true delete split volumes 
    if self.cleanUpFiles:
        self._insertFunctionStep('cleanVolumeStep', vol1, vol2)


def runCalculateFscStep(self, iterN, refN, args, constantToAdd, **kwargs):
    if self.getEnumText('reconstructionMethod') == 'fourier':
        if self._fourierMaxFrequencyOfInterest[iterN] == -1:
            fourierMaxFrequencyOfInterest = self._getFourierMaxFrequencyOfInterest(iterN-1, refN)
            normalizedFreq = self.resolSam / fourierMaxFrequencyOfInterest + constantToAdd
            fourierMaxFrequencyOfInterest = self.resolSam / normalizedFreq
        else:
            fourierMaxFrequencyOfInterest = self.resolSam / self._fourierMaxFrequencyOfInterest[iterN]
        maxFreq = fourierMaxFrequencyOfInterest
        args += ' --max_sam %(maxFreq)s' % locals()
    
    self.runJob("xmipp_resolution_fsc", args, numberOfMpi=1, **kwargs)


def runStoreResolutionStep(self, resolIterMd, resolIterMaxMd, sampling):
    self._log.info("compute resolution 1")
    #compute resolution
    mdRsol = xmipp.MetaData(resolIterMd)
    mdResolOut = xmipp.MetaData()
    mdResolOut.importObjects(mdRsol, xmipp.MDValueLT(xmipp.MDL_RESOLUTION_FRC, 0.5))
    self._log.info("compute resolution 2")
    if mdResolOut.size()==0:
        mdResolOut.clear()
        mdResolOut.addObject()
        id=mdResolOut.firstObject()
        mdResolOut.setValue(xmipp.MDL_RESOLUTION_FREQREAL, sampling*2., id)
        mdResolOut.setValue(xmipp.MDL_RESOLUTION_FRC, 0.5, id)
    else:
        mdResolOut.sort()
    
    id = mdResolOut.firstObject()
    filterFrequence = mdResolOut.getValue(xmipp.MDL_RESOLUTION_FREQREAL, id)
    frc = mdResolOut.getValue(xmipp.MDL_RESOLUTION_FRC, id)
    
    md = xmipp.MetaData()
    id = md.addObject()
    md.setColumnFormat(False)
    
    md.setValue(xmipp.MDL_RESOLUTION_FREQREAL, filterFrequence, id)
    md.setValue(xmipp.MDL_RESOLUTION_FRC, frc, id)
    md.setValue(xmipp.MDL_SAMPLINGRATE, sampling, id)
    md.write(resolIterMaxMd, xmipp.MD_APPEND)


def insertFilterVolumeStep(self, iterN, refN, **kwargs):
    reconstructedVolume = self._getFileName('reconstructedFileNamesIters', iter=iterN, ref=refN)
    reconstructedFilteredVolume = self.reconstructedFilteredFileNamesIters[iterN][refN]
    
    if not self.doLowPassFilter:
        return self._insertCopyFileStep(reconstructedVolume, reconstructedFilteredVolume)
    else:
        return self._insertFunctionStep('filterVolumeStep', iterN, refN, self._constantToAddToFiltration[iterN])


def runFilterVolumeStep(self, iterN, refN, constantToAddToFiltration):
    reconstructedVolume = self._getFileName('reconstructedFileNamesIters', iter=iterN, ref=refN)
    reconstructedFilteredVolume = self.reconstructedFilteredFileNamesIters[iterN][refN]
    if self.useFscForFilter:
        if self._fourierMaxFrequencyOfInterest[iterN+1] == -1:
            fourierMaxFrequencyOfInterest = self.resolSam / self._getFourierMaxFrequencyOfInterest(iterN, refN)
            print "el valor de la resolucion es :", self._getFourierMaxFrequencyOfInterest(iterN, refN)
            filterInPxAt = fourierMaxFrequencyOfInterest + constantToAddToFiltration
        else:
            filterInPxAt = constantToAddToFiltration
    
    if (filterInPxAt > 0.5):
        copyFile(reconstructedVolume, reconstructedFilteredVolume)
    else:
        args = ' -i %(volume)s -o %(filteredVol)s --fourier low_pass %(filter)s'
        params = {'volume': reconstructedVolume,
                  'filteredVol': reconstructedFilteredVolume,
                  'filter' : filterInPxAt
                  }
        self.runJob("xmipp_transform_filter", args % params)


def runCreateOutputStep(self):
    ''' Create standard output results_images, result_classes'''
    #creating results files
    imgSet = self.inputParticles.get()
    lastIter = self.numberOfIterations.get()
    if self.numberOfReferences != 1:
        inDocfile = self._getFileName('docfileInputAnglesIters', iter=lastIter)
        ClassFnTemplate = '%(rootDir)s/reconstruction_Ref3D_%(ref)03d.vol'
        
        allExpImagesinDocfile = xmipp.FileName()
        all_exp_images="all_exp_images"
        allExpImagesinDocfile.compose(all_exp_images, inDocfile)
        
        dataClasses = self._getFileName('sqliteClasses')
        
        createClassesFromImages(imgSet, str(allExpImagesinDocfile), dataClasses, 
                                SetOfClasses3D, xmipp.MDL_REF3D, ClassFnTemplate, lastIter)
        
        classes = self._createSetOfClasses3D(imgSet)
        clsSet = SetOfClasses3D(dataClasses)
        classes.appendFromClasses(clsSet)
        
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(imgSet.getSamplingRate())
        
        for refN in self.allRefs():
            volFn = self._getFileName('reconstructedFileNamesIters', iter=lastIter, ref=refN)
            vol = Volume()
            vol.setFileName(volFn)
            volumes.append(vol)
    
        self._defineOutputs(outputVolumes=volumes)
        self._defineOutputs(outputClasses=classes)
        self._defineSourceRelation(self.inputParticles, volumes)
        self._defineSourceRelation(self.inputParticles, classes)
        self._defineSourceRelation(self.input3DReferences, volumes)
        self._defineSourceRelation(self.input3DReferences, classes)
    else:
        volFn = self._getFileName('reconstructedFileNamesIters', iter=lastIter, ref=1)
        vol = Volume()
        vol.setFileName(volFn)
        vol.setSamplingRate(imgSet.getSamplingRate())
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineSourceRelation(self.input3DReferences, vol)
        
        #create set of images
        imgSetOut = self._createSetOfParticles("_iter_%03d" %lastIter)
        self._fillParticlesFromIter(imgSetOut, lastIter)
        
        self._defineOutputs(outputParticles=imgSetOut)
        self._defineSourceRelation(self.inputParticles, imgSetOut)
        self._defineSourceRelation(self.input3DReferences, imgSetOut)
