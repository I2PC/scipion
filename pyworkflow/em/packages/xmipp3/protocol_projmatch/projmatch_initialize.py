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
This sub-module contains several initialization functions related to 
Projection Matching.
"""

from os.path import join

# taken from angular_project_library in xmipp

def createFilenameTemplates(self):
    """ Centralize how files are called for iterations and references. """
    # Setup some params to be used in command line formatting
    LibraryDir = "ReferenceLibrary"
    self._params = {
                    'ReferenceVolumeName': 'reference_volume.vol',
                    'LibraryDir': LibraryDir,
                    'ProjectLibraryRootName': join(LibraryDir, "gallery"),
                    'ProjMatchDir': "ProjMatchClasses",
                    'ProjMatchName': self.name,
                    'ClassAverageName': 'class_average',
                    #ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
                    'ForReconstructionSel': "reconstruction.sel",
                    'ForReconstructionDoc': "reconstruction.doc",
                    'MultiAlign2dSel': "multi_align2d.sel",
                    'BlockWithAllExpImages' : 'all_exp_images',
                    'DocFileWithOriginalAngles': 'original_angles.doc',
                    'Docfile_with_current_angles': 'current_angles',
                    'Docfile_with_final_results': 'results.xmd',
                    'FilteredReconstruction': "filtered_reconstruction",
                    'ReconstructedVolume': "reconstruction",
                    'MaskReferenceVolume': "masked_reference",
                    'OutputFsc': "resolution.fsc",
                    'CtfGroupDirectory': "CtfGroups",
                    'CtfGroupRootName': "ctf",
                    'CtfGroupSubsetFileName': "ctf_images.sel",
                    }
    # Also setup as protocol attributes
    for k, v in self._params.iteritems():
        setattr(self, k, v)
                          
    Iter = 'Iter_%(iter)03d'
    Ref3D = 'Ref3D_%(ref)03d'
    Ctf = 'CtfGroup_%(ctf)06d'
    IterDir = self._getExtraPath(Iter)
    
    #ProjMatchDirs = join(IterDir, '%(ProjMatchDir)s.doc')
    ProjMatchDirs = join(IterDir, '%(ProjMatchDir)s')
    _OutClassesXmd = join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.xmd')
    _OutClassesXmdS1 = join(ProjMatchDirs, '%(ProjMatchName)s_split_1_' + Ref3D + '.xmd')
    _OutClassesXmdS2 = join(ProjMatchDirs, '%(ProjMatchName)s_split_2_' + Ref3D + '.xmd')
    CtfGroupBase = join(self.workingDirPath(), self.ctfGroupDirectory, '%(CtfGroupRootName)s')
    ProjLibRootNames = join(IterDir, '%(ProjectLibraryRootName)s_' + Ref3D)
    myDict = {
            # Global filenames templates
            'IterDir': IterDir,
            'ProjMatchDirs': ProjMatchDirs,
            'DocfileInputAnglesIters': join(IterDir, '%(Docfile_with_current_angles)s.doc'),
            'LibraryDirs': join(IterDir, '%(LibraryDir)s'),
            'ProjectLibraryRootNames': ProjLibRootNames,
#                'ProjMatchRootNames': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.doc'),
            'ProjMatchRootNames': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.sqlite'),
            'ProjMatchRootNamesWithoutRef': join(ProjMatchDirs, '%(ProjMatchName)s.doc'),
            'OutClasses': join(ProjMatchDirs, '%(ProjMatchName)s'),
            'OutClassesXmd': _OutClassesXmd,
            'OutClassesStk': join(ProjMatchDirs, '%(ProjMatchName)s_' + Ref3D + '.stk'),
            'OutClassesDiscarded': join(ProjMatchDirs, '%(ProjMatchName)s_discarded.xmd'),
            'ReconstructionXmd': Ref3D + '@' +_OutClassesXmd,
            'ReconstructionXmdSplit1': Ref3D + '@' +_OutClassesXmdS1,
            'ReconstructionXmdSplit2': Ref3D + '@' +_OutClassesXmdS2,
            'MaskedFileNamesIters': join(IterDir, '%(MaskReferenceVolume)s_' + Ref3D + '.vol'),
            'ReconstructedFileNamesIters': join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '.vol'),
            'ReconstructedFileNamesItersSplit1': join(IterDir, '%(ReconstructedVolume)s_split_1_' + Ref3D + '.vol'),
            'ReconstructedFileNamesItersSplit2': join(IterDir, '%(ReconstructedVolume)s_split_2_' + Ref3D + '.vol'),
            'ReconstructedFilteredFileNamesIters': join(IterDir, '%(ReconstructedVolume)s_filtered_' + Ref3D + '.vol'),
            'ResolutionXmdFile': join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '_frc.xmd'),
            'ResolutionXmd': 'resolution@' + join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '_frc.xmd'),
            'ResolutionXmdMax': 'resolution_max@' + join(IterDir, '%(ReconstructedVolume)s_' + Ref3D + '_frc.xmd'),
            'MaskedFileNamesIters': join(IterDir, '%(MaskReferenceVolume)s_' + Ref3D + '.vol'),
            # Particular templates for executeCtfGroups  
            'ImageCTFpairs': CtfGroupBase + '_images.sel',
            'CTFGroupSummary': CtfGroupBase + 'Info.xmd',
            'StackCTFs': CtfGroupBase + '_ctf.stk',
            'StackWienerFilters': CtfGroupBase + '_wien.stk',
            'SplitAtDefocus': CtfGroupBase + '_split.doc',
            # Particular templates for angular_project_library 
            'ProjectLibraryStk': ProjLibRootNames + '.stk',
            'ProjectLibraryDoc': ProjLibRootNames + '.doc',
            'ProjectLibrarySampling': ProjLibRootNames + '_sampling.xmd',
            'ProjectLibraryGroupSampling': ProjLibRootNames + '_group%(group)06d_sampling.xmd',
            }

    self._updateFilenamesDict(myDict)
    
    
def initializeLists(self):
    """ Load lists with values for iterations...values are taken from string values.
    """
    # Construct special filename list with zero special case
    self.docFileInputAngles = [self.docFileWithOriginalAngles] + [self._getFileName('DocfileInputAnglesIters', iter=i) for i in self.allIters()]
    #print 'self.docFileInputAngles: ', self.docFileInputAngles
    self.reconstructedFileNamesIters = [[None] + self.referenceFileNames]
    for iterN in self.allIters():
        self.reconstructedFileNamesIters.append([None] + [self._getFileName('ReconstructedFileNamesIters', iter=iterN, ref=r) for r in self.allRefs()])
    
    self.reconstructedFilteredFileNamesIters = [[None] + self.referenceFileNames]
    for iterN in self.allIters():
        self.reconstructedFilteredFileNamesIters.append([None] + [self._getFileName('ReconstructedFilteredFileNamesIters', iter=iterN, ref=r) for r in self.allRefs()])
    
    maxFreq = self.fourierMaxFrequencyOfInterest
    n2 = self.numberOfIterations.get() + 2
    
    if self.doComputeResolution=='0' or self.doComputeResolution==0:
        self.fourierMaxFrequencyOfInterest = [maxFreq] * n2
    else:
        self.fourierMaxFrequencyOfInterest = [-1] * n2
        self.fourierMaxFrequencyOfInterest[1] = maxFreq
    
    #parameter for projection matching
    self.align2DIterNr = self.itersFloatValues('align2DIterNr')
    self.align2dMaxChangeOffset = self.itersFloatValues('align2dMaxChangeOffset')
    self.align2dMaxChangeRot = self.itersFloatValues('align2dMaxChangeRot')
    self.angSamplingRateDeg = self.itersFloatValues('angSamplingRateDeg')
    self.constantToAddToFiltration = self.itersFloatValues('constantToAddToFiltration')
    self.constantToAddToMaxReconstructionFrequency = self.itersFloatValues('constantToAddToMaxReconstructionFrequency')
    self.discardPercentage = self.itersFloatValues('discardPercentage')
    self.discardPercentagePerClass = self.itersFloatValues('discardPercentagePerClass')
    self.doAlign2D = self.itersBoolValues('doAlign2D')
    self.doComputeResolution = self.itersBoolValues('doComputeResolution')
    self.doSplitReferenceImages = self.itersBoolValues('doSplitReferenceImages')
    self.innerRadius = self.itersFloatValues('innerRadius', firstValue=False) # FIXME: merging bool and floats
    self.maxChangeInAngles = self.itersFloatValues('maxChangeInAngles')
    self.maxChangeOffset = self.itersFloatValues('maxChangeOffset')
    self.minimumCrossCorrelation = self.itersFloatValues('minimumCrossCorrelation')
    self.onlyWinner = self.itersBoolValues('onlyWinner')
    self.outerRadius = self.itersFloatValues('outerRadius', firstValue=False) # FIXME: merging bool and floats
    self.perturbProjectionDirections = self.itersBoolValues('perturbProjectionDirections')
    self.referenceIsCtfCorrected = self.itersFloatValues(str(self.ReferenceIsCtfCorrected.get()) + " True") # FIXME: using bool string and float list
    self.scaleNumberOfSteps = self.itersFloatValues('scaleNumberOfSteps')
    self.scaleStep = self.itersFloatValues('scaleStep')
    self.search5DShift = self.itersFloatValues('search5DShift')
    self.search5DStep = self.itersFloatValues('search5DStep')
    self.symmetryGroup = self.itersFloatValues('symmetryGroup')
