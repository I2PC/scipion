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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-module contains several initialization functions related to 
Projection Matching.
"""

from os.path import join, exists
from pyworkflow.object import String

# taken from angular_project_library in xmipp

def createFilenameTemplates(self):
    """ Centralize how files are called for iterations and references. """
    # Setup some params to be used in command line formatting
    LibraryDir = "ReferenceLibrary"
    projectLibraryRootName = join(LibraryDir, "gallery")
    
    self._params = {
                    'referenceVolumeName': 'reference_volume.vol',
                    'libraryDir': LibraryDir,
                    'projectLibraryRootName': projectLibraryRootName,
                    'projMatchDir': "ProjMatchClasses",
                    'projMatchName': self.getClassName(), 
                    'classAverageName': 'class_average',
                    #ProjMatchRootName = ProjMatchDir + "/" + ProjMatchName
                    'forReconstructionSel': "reconstruction.sel",
                    'forReconstructionDoc': "reconstruction.doc",
                    'multiAlign2dSel': "multi_align2d.sel",
                    'blockWithAllExpImages' : 'all_exp_images',
                    'docFileWithOriginalAngles': 'original_angles.doc',
                    'docfile_with_current_angles': 'current_angles',
                    'docfile_with_final_results': 'results.xmd',
                    'filteredReconstruction': "filtered_reconstruction",
                    'reconstructedVolume': "reconstruction",
                    'maskReferenceVolume': "masked_reference",
                    'outputFsc': "resolution.fsc",
                    'ctfGroupDirectory': self._getExtraPath("CtfGroups"),
                    'ctfGroupRootName': "ctf",
                    'ctfGroupSubsetFileName': "ctf_images.sel",
                    }
    # Also setup as protocol attributes
    for k, v in self._params.iteritems():
        setattr(self, k, v)
                          
    Iter = 'iter_%(iter)03d'
    Ref3D = 'Ref3D_%(ref)03d'
    Ctf = 'CtfGroup_%(ctf)06d'
    IterDir = self._getExtraPath(Iter)
    ProjMatchDirs = join(IterDir, self.projMatchDir)
    
    def iterFile(key, suffix=''):
        return join(IterDir, (key % self._params) + suffix)
    
    def projmatchFile(key, suffix=''):
        return join(ProjMatchDirs, (key % self._params) + suffix)    
    
    _OutClassesXmd = projmatchFile('%(projMatchName)s_', Ref3D + '.xmd')
    _OutClassesXmdS1 = projmatchFile('%(projMatchName)s_split_1_', Ref3D + '.xmd')
    _OutClassesXmdS2 = projmatchFile('%(projMatchName)s_split_2_', Ref3D + '.xmd')
    CtfGroupBase = join(self.ctfGroupDirectory, self.ctfGroupRootName)
    ProjLibRootNames = join(IterDir, projectLibraryRootName + '_' + Ref3D)
    
    
    myDict = {
            # Global filenames templates
            'iterDir': IterDir,
            'projMatchDirs': ProjMatchDirs,
            'docfileInputAnglesIters': iterFile('%(docfile_with_current_angles)s.doc'),
            'libraryDirs': iterFile('%(libraryDir)s'),
            'projectLibraryRootNames': ProjLibRootNames,
            'projMatchRootNames': projmatchFile('%(projMatchName)s_', Ref3D + '.sqlite'),
            'projMatchRootNamesWithoutRef': projmatchFile('%(projMatchName)s.doc'),
            'outClasses': projmatchFile('%(projMatchName)s'),
            'outClassesXmd': _OutClassesXmd,
            'outClassesStk': projmatchFile('%(projMatchName)s_', Ref3D + '.stk'),
            'outClassesDiscarded': projmatchFile('%(projMatchName)s_discarded.xmd'),
            'reconstructionXmd': Ref3D + '@' +_OutClassesXmd,
            'reconstructionXmdSplit1': Ref3D + '@' +_OutClassesXmdS1,
            'reconstructionXmdSplit2': Ref3D + '@' +_OutClassesXmdS2,
            'maskedFileNamesIters': iterFile('%(maskReferenceVolume)s_', Ref3D + '.vol'),
            'reconstructedFileNamesIters': iterFile('%(reconstructedVolume)s_', Ref3D + '.vol'),
            'reconstructedFileNamesItersSplit1': iterFile('%(reconstructedVolume)s_split_1_', Ref3D + '.vol'),
            'reconstructedFileNamesItersSplit2': iterFile('%(reconstructedVolume)s_split_2_', Ref3D + '.vol'),
            'reconstructedFilteredFileNamesIters': iterFile('%(reconstructedVolume)s_filtered_', Ref3D + '.vol'),
            'resolutionXmdFile': iterFile('%(reconstructedVolume)s_', Ref3D + '_frc.xmd'),
            'resolutionXmd': 'resolution@' + iterFile('%(reconstructedVolume)s_', Ref3D + '_frc.xmd'),
            'resolutionXmdMax': 'resolution_max@' + iterFile('%(reconstructedVolume)s_', Ref3D + '_frc.xmd'),
            'maskedFileNamesIters': iterFile('%(maskReferenceVolume)s_', Ref3D + '.vol'),
            # Particular templates for executeCtfGroups  
            'imageCTFpairs': CtfGroupBase + '_images.sel',
            'cTFGroupSummary': CtfGroupBase + 'Info.xmd',
            'stackCTFs': CtfGroupBase + '_ctf.stk',
            'stackWienerFilters': CtfGroupBase + '_wien.stk',
            'splitAtDefocus': CtfGroupBase + '_split.doc',
            'ctfGroupBase' : CtfGroupBase,
            # Particular templates for angular_project_library 
            'projectLibraryStk': ProjLibRootNames + '.stk',
            'projectLibraryDoc': ProjLibRootNames + '.doc',
            'projectLibrarySampling': ProjLibRootNames + '_sampling.xmd',
            'projectLibraryGroupSampling': ProjLibRootNames + '_group%(group)06d_sampling.xmd',
            # ADDED FOR SCIPION
            'inputParticlesXmd': self._getPath('input_particles.xmd'),
            'inputParticlesDoc': self._getExtraPath(self.docFileWithOriginalAngles),
            'sqliteClasses' : self._getPath('classes_scipion.sqlite'), 
            'particlesScipion' : self._getPath('particles_iter_%(iter)03d.sqlite')
            }
    myDict.update(self._params)
    self._updateFilenamesDict(myDict)
    self.selFileName = self._getFileName('inputParticlesXmd')
    
    if self.doCTFCorrection:
        if exists(self.setOfDefocus.get()):
            self.ctfDatName = self.setOfDefocus
        else:
            self.ctfDatName = self.selFileName
    
def initializeLists(self):
    """ Load lists with values for iterations...values are taken from string values.
    """
    # Construct special filename list with zero special case
    self.docFileInputAngles = [self._getFileName('inputParticlesDoc')] + [self._getFileName('docfileInputAnglesIters', iter=i) for i in self.allIters()]
    #print 'self.docFileInputAngles: ', self.docFileInputAngles
    self.reconstructedFileNamesIters = [[None] + self.referenceFileNames]
    for iterN in self.allIters():
        self.reconstructedFileNamesIters.append([None] + [self._getFileName('reconstructedFileNamesIters', iter=iterN, ref=r) for r in self.allRefs()])
    
    self.reconstructedFilteredFileNamesIters = [[None] + self.referenceFileNames]
    for iterN in self.allIters():
        self.reconstructedFilteredFileNamesIters.append([None] + [self._getFileName('reconstructedFilteredFileNamesIters', iter=iterN, ref=r) for r in self.allRefs()])
    
    maxFreq = self.fourierMaxFrequencyOfInterest.get()
    n2 = self.numberOfIterations.get() + 2
    
    if self.doComputeResolution:
        self._fourierMaxFrequencyOfInterest = [-1] * n2
        self._fourierMaxFrequencyOfInterest[1] = maxFreq
    else:
        self._fourierMaxFrequencyOfInterest = [maxFreq] * n2
    
    #parameter for projection matching
    self._align2DIterNr = self.itersFloatValues('align2DIterNr')
    self._align2dMaxChangeOffset = self.itersFloatValues('align2dMaxChangeOffset')
    self._align2dMaxChangeRot = self.itersFloatValues('align2dMaxChangeRot')
    self._angSamplingRateDeg = self.itersFloatValues('angSamplingRateDeg')
    self._artLambda = self.itersFloatValues('artLambda')
    self._constantToAddToFiltration = self.itersFloatValues('constantToAddToFiltration')
    self._constantToAddToMaxReconstructionFrequency = self.itersFloatValues('constantToAddToMaxReconstructionFrequency')
    self._discardPercentage = self.itersFloatValues('discardPercentage')
    self._discardPercentagePerClass = self.itersFloatValues('discardPercentagePerClass')
    self._doAlign2D = self.itersBoolValues('doAlign2D')
    self._doSplitReferenceImages = self.itersBoolValues('doSplitReferenceImages')
    self._innerRadius = self.itersFloatValues('innerRadius')
    self._maxChangeInAngles = self.itersFloatValues('maxChangeInAngles')
    self._maxChangeOffset = self.itersFloatValues('maxChangeOffset')
    self._minimumCrossCorrelation = self.itersFloatValues('minimumCrossCorrelation')
    self._onlyWinner = self.itersBoolValues('onlyWinner')
    self._outerRadius = self.itersFloatValues('outerRadius')
    self._perturbProjectionDirections = self.itersBoolValues('perturbProjectionDirections')
    #FIXME: this need to use later self.referenceIsCtfCorrected value or take it
    # from the input references objects
    self.referenceIsCtfCorrected = String('True')
    self._referenceIsCtfCorrected = self.itersBoolValues('referenceIsCtfCorrected') 
    
    self._scaleNumberOfSteps = self.itersFloatValues('scaleNumberOfSteps')
    self._scaleStep = self.itersFloatValues('scaleStep')
    self._search5DShift = self.itersFloatValues('search5DShift')
    self._search5DStep = self.itersFloatValues('search5DStep')
    self._symmetry = self.itersStringValues('symmetry')
