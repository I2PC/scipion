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
import os, glob, math
from pyworkflow.em import *
import xmipp





# Functions outside th loop loop for xmipp_projection_matching
def insertExecuteCtfGroupsStep(self, **kwargs):
    #...
    self._insertRunJobStep('xmipp_ctf_group') #...

def insertInitAngularReferenceFileStep(self, **kwargs):
    #...
    self._insertRunJobStep('') #...

# Functions in loop for xmipp_projection_matching

def insertMaskReferenceStep(self, iterN, refN, **kwargs):
    #...
    maskedFileName = self.getFilename('MaskedFileNamesIters', iter=iterN, ref=refN)
    ReconstructedFilteredVolume = self.reconstructedFilteredFileNamesIters[iterN-1][refN]
    
    self._insertRunJobStep('xmipp_transform_mask') #...


def insertAngularProjectLibraryStep(self, iterN, refN, **kwargs):
    self._args = ' -i %(maskedFileNamesIter)s --experimental_images %(experimentalImages)s -o %(projectLibraryRootName)s --sampling_rate %(samplingRate)s --sym %(symmetry)s'
    self._args += 'h --compute_neighbors --method %(projectionMethod)s' 

    ###need one block per reference
    # Project all references
    print '* Create projection library'
    
    if isinstance(self.input3DReferences.get(), Volume):
        xDim, yDim, zDim, _ = self.input3DReferences.get().getDim()
    else:
        xDim, yDim, zDim, = self.input3DReferences.get().getDimensions()
    
    memoryUsed = (xDim * yDim * zDim * 8) / pow(2,20)
    projectionMethod = self.getEnumText('projectionMethod')
    expImages = self.blockWithAllExpImages + '@' + self.docFileInputAngles[iterN-1]
    
    params = {'maskedFileNamesIter' : self.maskedFileNamesIter,
              'experimentalImages' : expImages,
              'projectLibraryRootName' : self.projectLibraryRootName,
              'samplingRate' : self.angSamplingRateDeg[iterN],
              'symmetry' : self.symmetry[iterN],
              'projectionMethod' : projectionMethod,
              }
    
    if projectionMethod == 'fourier':
        memoryUsed = memoryUsed * 6
        
        if self.fourierMaxFrequencyOfInterest == -1:
            md = xmipp.MetaData(self._getFileName('resolutionXmdPrevIterMax'))
            id = md.firstObject()
            fourierMaxFrequencyOfInterest = Float(md.getValue(xmipp.MDL_RESOLUTION_FREQREAL, id))
            fourierMaxFrequencyOfInterest = self.resolSam / fourierMaxFrequencyOfInterest + self.constantToAddToFiltration
            
            if fourierMaxFrequencyOfInterest > 0.5:
                fourierMaxFrequencyOfInterest = 0.5
            elif fourierMaxFrequencyOfInterest < 0.:
                fourierMaxFrequencyOfInterest = 0.001
        
        params['paddingAngularProjection'] = self.PaddingAngularProjection.get()
        params['fourierMaxFrequencyOfInterest'] = fourierMaxFrequencyOfInterest
        params['kernelAngularProjection'] = self.getEnumText('kernelAngularProjection')
        self._args += ' %(paddingAngularProjection)s %(fourierMaxFrequencyOfInterest)s %(kernelAngularProjection)s'
        
    if self.maxChangeInAngles < 181:
        self._args += ' --near_exp_data --angular_distance %(maxChangeInAngles)s'
    else:
        self._args += ' --angular_distance -1'
    
    if self.perturbProjectionDirections[iterN]:
        self._args +=' --perturb %(perturb)s'
        self._params['perturb'] = math.sin(math.radians(self.angSamplingRateDeg[iterN])) / 4.

    if self.doRestricSearchbyTiltAngle:
        self._args += ' --min_tilt_angle %(tilt0)s --max_tilt_angle %(tiltF)s'
        self._params['perturb'] = math.sin(math.radians(self.angSamplingRateDeg[iterN])) / 4.
        parameters += \
              ' --min_tilt_angle ' + str(Tilt0) + \
              ' --max_tilt_angle ' + str(TiltF)
 
#     if (DoCtfCorrection):
#         parameters += \
#               ' --groups ' + CtfGroupSubsetFileName
#     processorsToUse=NumberOfMpi * NumberOfThreads
#     if processorsToUse>1:
#         memoryAvailable=getMemoryAvailable()
#         processorsToUse=min(processorsToUse,floor(memoryAvailable/memoryUsed))
#     if (DoParallel and processorsToUse>1):
#         parameters = parameters + ' --mpi_job_size ' + str(MpiJobSize)
#     if (len(SymmetryGroupNeighbourhood) > 1):
#         parameters += \
#           ' --sym_neigh ' + SymmetryGroupNeighbourhood + 'h'
#     if (OnlyWinner):
#         parameters += \
#               ' --only_winner '
# 
#     runJob(_log, 'xmipp_angular_project_library',
#                          parameters,
#                          processorsToUse)
#     if (not DoCtfCorrection):
#         src = ProjectLibraryRootName.replace(".stk", '_sampling.xmd')
#         dst = src.replace('sampling.xmd', 'group%06d_sampling.xmd' % 1)
#         copyFile(_log, src, dst)
        self._insertRunJobStep('xmipp_angular_project_library', kk, **kwargs)

def insertProjectionMatchingStep(self, iterN, refN, **kwargs):
    #...
    self._insertRunJobStep('xmipp_angular_projection_matching') #...


def insertAssignImagesToReferencesStep(self, iterN, **kwargs):
    #...
    self._insertRunJobStep('') #...


def insertAngularClassAverageStep(self, iterN, **kwargs):
    #...
    self._insertRunJobStep('xmipp_angular_class_average') #...


def insertReconstructionStep(self, iterN, refN, suffix='', **kwargs):
    #...
    self._insertRunJobStep('xmipp_reconstruct_fourier', **kwargs) #...


def insertComputeResolutionStep(self, iterN, refN, **kwargs):
    #...
    self._insertRunJobStep('xmipp_resolution_fsc') #...


def insertFilterVolumeStep(self, iterN, refN, **kwargs):
    #...
    self._insertRunJobStep('xmipp_transform_filter') #...

