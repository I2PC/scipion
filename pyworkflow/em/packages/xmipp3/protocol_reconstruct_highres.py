# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
Protocol to perform high-resolution reconstructions
"""

from pyworkflow.object import Float
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam
from pyworkflow.utils.path import cleanPath, makePath, copyFile, moveFile, createLink
from pyworkflow.em.protocol import ProtRefine3D
from pyworkflow.em.data import SetOfVolumes
from pyworkflow.em.metadata.utils import getFirstRow
from convert import writeSetOfParticles
from os.path import join, exists
from pyworkflow.em.convert import ImageHandler

from xmipp import MetaData, MDL_RESOLUTION_FRC, MDL_RESOLUTION_FREQREAL, MDL_SAMPLINGRATE, MDL_WEIGHT_SSNR, MDL_WEIGHT_SIGNIFICANT, \
                  MDL_WEIGHT, MD_APPEND, MDL_XSIZE, MDL_WEIGHT_CONTINUOUS2, MDL_ANGLE_DIFF, MDL_IMAGE, MDL_IMAGE1, MDL_IMAGE_ORIGINAL, \
                  MDL_COUNT, MDL_SHIFT_X, MDL_CONTINUOUS_X, MDL_WEIGHT_JUMPER, MDL_CTF_DEFOCUSU, MDL_CTF_MODEL, MDL_PARTICLE_ID, \
                  MDL_ZSCORE_RESCOV, MDL_ZSCORE_RESVAR, MDL_ZSCORE_RESMEAN, Image
from xmipp3 import HelicalFinder

class XmippProtReconstructHighRes(ProtRefine3D, HelicalFinder):
    """Reconstruct a volume at high resolution"""
    _label = 'highres'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                      'run of type *%s* class and some of the input parameters'
                      'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam, label="Full-size Images", important=True, 
                      condition='not doContinue', pointerClass='SetOfParticles',
                      help='Select a set of images at full resolution')
        form.addParam('phaseFlipped', BooleanParam, label="Images have been phase flipped", default=True, 
                      condition='not doContinue', help='Choose this option if images have been phase flipped')
        form.addParam('inputVolumes', PointerParam, label="Initial volumes", important=True,
                      condition='not doContinue', pointerClass='Volume, SetOfVolumes',
                      help='Select a set of volumes with 2 volumes or a single volume')
        form.addParam('particleRadius', IntParam, default=-1, 
                     condition='not doContinue', label='Radius of particle (px)',
                     help='This is the radius (in pixels) of the spherical mask covering the particle in the input images')       

        form.addParam('continueRun', PointerParam, pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('numberOfIterations', IntParam, default=6, label='Number of iterations')
        form.addParam("saveSpace", BooleanParam, default=False, label="Remove intermediary files")
        
        form.addSection(label='Weights')
        form.addParam('weightSSNR', BooleanParam, label="Weight by SSNR?", default=True,
                      help='Weight input images by SSNR')
        form.addParam('weightSignificance', BooleanParam, label="Weight by Significance?", default=True,
                      help='Weight input images by angular assignment significance')
        form.addParam('weightContinuous', BooleanParam, label="Weight by Continuous cost?", default=True,
                      help='Weight input images by angular assignment cost')
        form.addParam('weightJumper', BooleanParam, label="Weight by angular stability?", default=True,
                      help='Weight input images by angular stability between iterations')
        form.addParam('weightAstigmatism', BooleanParam, label="Weight by astigmatism?", default=True,
                      help='Give lower weight to those images whose astigmatic CTF would not allow to reach high resolution.' \
                           'This weight is calculated by the continuous assignment.')
        form.addParam('weightAstigmatismSigma', FloatParam, label="Astigmatism sigma", default=120, expertLevel=LEVEL_ADVANCED, 
                        condition="weightAstigmatism", help="Sigma in degrees for the CTF phase");
        form.addParam('weightResiduals', BooleanParam, label="Weight by residuals?", default=True,
                      help="Analyze how different are the image residuals, it only works after running the continuous assignment");

        form.addSection(label='Next Reference')
        form.addParam('nextResolutionCriterion',FloatParam, label="FSC criterion", default=0.143, 
                      help='The resolution of the reconstruction is defined as the inverse of the frequency at which '\
                      'the FSC drops below this value. Typical values are 0.143 and 0.5')
        form.addParam('nextLowPass', BooleanParam, label="Low pass filter?", default=True,
                      help='Apply a low pass filter to the previous iteration whose maximum frequency is '\
                           'the current resolution(A) + resolutionOffset(A). If resolutionOffset>0, then fewer information' \
                           'is used (meant to avoid overfitting). If resolutionOffset<0, then more information is allowed '\
                           '(meant for a greedy convergence).')
        form.addParam('nextResolutionOffset', FloatParam, label="Resolution offset (A)", default=3, condition='nextLowPass')
        form.addParam('nextSpherical', BooleanParam, label="Spherical mask?", default=True,
                      help='Apply a spherical mask of the size of the particle')
        form.addParam('nextPositivity', BooleanParam, label="Positivity?", default=True,
                      help='Remove from the next reference all negative values')
        form.addParam('nextMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be between 0 (remove these pixels) and 1 (let them pass). Smooth masks are recommended.')
        form.addParam('nextReferenceScript', StringParam, label="Next reference command", default="", expertLevel=LEVEL_ADVANCED, 
                      help='A command template that is used to generate next reference. The following variables can be used ' 
                           '%(sampling)s %(dim)s %(volume)s %(iterDir)s. The command should read Spider volumes and modify the input volume.'
                           'the command should be accessible either from the PATH or provide the absolute path.\n'
                           'Examples: \n'
                           'xmipp_transform_filter -i %(volume)s --fourier low_pass 15 --sampling %(sampling)s\n' 
                           '/home/joe/myScript %(volume)s sampling=%(sampling)s dim=%(dim)s')
        form.addParam('nextRemove', BooleanParam, label="Remove reference to save space?", default=True, expertLevel=LEVEL_ADVANCED, 
                      help='Remove reference volumes once they are not needed any more.')

        form.addSection(label='Angular assignment')
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=10,
                      help='Maximum shift as a percentage of the image size')
        line=form.addLine('Tilt angle:', help='0 degrees represent top views, 90 degrees represent side views')
        line.addParam('angularMinTilt', FloatParam, label="Min.", default=0)
        line.addParam('angularMaxTilt', FloatParam, label="Max.", default=90)
        groupSignificant = form.addGroup('Global')
        groupSignificant.addParam('significantMaxResolution', FloatParam, label="Global assignment if resolution is worse than (A)", default=12,
                      help='Significant assignment is always performed on the first iteration. Starting from the second, you may '\
                      'decide whether to perform it or not. Note that the significant angular assignment is a robust angular assignment '\
                      'meant to avoid local minima, although it may take time to calculate.')
        groupSignificant.addParam('significantSignificance', FloatParam, label="Significance (%)", default=99.75)
        groupSignificant.addParam('significantGrayValues', BooleanParam, label="Optimize gray values?", default=True)
        groupContinuous = form.addGroup('Local')
        groupContinuous.addParam('continuousMinResolution', FloatParam, label="Continuous assignment if resolution is better than (A)", default=15,
                      help='Continuous assignment can produce very accurate assignments if the initial assignment is correct.')
        groupContinuous.addParam('contShift', BooleanParam, label="Optimize shifts?", default=True,
                      help='Optimize shifts within a limit')
        groupContinuous.addParam('contScale', BooleanParam, label="Optimize scale?", default=True,
                      help='Optimize scale within a limit')
        groupContinuous.addParam('contMaxScale', FloatParam, label="Max. scale variation", default=0.02, expertLevel=LEVEL_ADVANCED)
        groupContinuous.addParam('contAngles', BooleanParam, label="Optimize angles?", default=True,
                      help='Optimize angles within a limit')
        groupContinuous.addParam('contGrayValues', BooleanParam, label="Optimize gray values?", default=True,
                      help='Optimize gray values. Do not perform this unless the reconstructed volume is gray-compatible with the projections,'\
                      ' i.e., the volumes haven been produced from projections')
        groupContinuous.addParam('contDefocus', BooleanParam, label="Optimize defocus?", default=True)
        groupContinuous.addParam('contPadding', IntParam, label="Fourier padding factor", default=2, expertLevel=LEVEL_ADVANCED,
                      help='The volume is zero padded by this factor to produce projections')
        
        form.addSection(label='Post-processing')
        form.addParam('postBFactor', BooleanParam, label="Correct for B-factor?", default=True)
        form.addParam('postBmin', FloatParam, label="Min. Resolution", default=10, condition='postBFactor', help='In Angstroms')
        form.addParam('postBmax', FloatParam, label="Max. Resolution", default=5, condition='postBFactor', help='In Angstroms')
        form.addParam('postNonnegativity', BooleanParam, label="Remove negative values?", default=True)        
        form.addParam('postAdHocMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be between 0 (remove these pixels) and 1 (let them pass). Smooth masks are recommended.')
        groupMask = form.addGroup('Mask')
        groupMask.addParam('postMask', BooleanParam, label="Construct mask for the reconstructed volume?", default=True)
        groupMask.addParam('postMaskThreshold', FloatParam, label="Mask sigma threshold", default=1, expertLevel=LEVEL_ADVANCED, condition="postMask",
                           help="In standard deviation units")
        groupMask.addParam('postDoMaskRemoveSmall', BooleanParam, label="Remove small objects?", default=True, expertLevel=LEVEL_ADVANCED,
                           condition="postMask")
        groupMask.addParam('postMaskRemoveSmallThreshold', IntParam, label="Small size", default=50, expertLevel=LEVEL_ADVANCED,
                           condition="postMask and postDoMaskRemoveSmall", help="An object is small if it has fewer than this number of voxels")
        groupMask.addParam('postMaskKeepLargest', BooleanParam, label="Keep largest component", default=True, expertLevel=LEVEL_ADVANCED,
                           condition="postMask")
        groupMask.addParam('postDoMaskDilate', BooleanParam, label="Dilate mask", default=True, expertLevel=LEVEL_ADVANCED,
                           condition="postMask")
        groupMask.addParam('postMaskDilateSize', IntParam, label="Dilation size", default=2, expertLevel=LEVEL_ADVANCED,
                           condition="postMask and postDoMaskDilate", help="In voxels")
        groupMask.addParam('postDoMaskSmooth', BooleanParam, label="Smooth borders", default=True, expertLevel=LEVEL_ADVANCED,
                           condition="postMask")
        groupMask.addParam('postMaskSmoothSize', FloatParam, label="Smooth size", default=2, expertLevel=LEVEL_ADVANCED, 
                           condition="postMask and postDoMaskSmooth", help="In voxels")
        groupMask.addParam('postMaskSymmetry', StringParam, label="Mask symmetry", default="c1", expertLevel=LEVEL_ADVANCED,
                           condition="postMask",
                           help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                           'If no symmetry is present, give c1')
        groupSymmetry = form.addGroup('Symmetry')
        groupSymmetry.addParam('postSymmetryWithinMask', BooleanParam, label="Symmetrize volume within mask?", default=False)
        groupSymmetry.addParam('postSymmetryWithinMaskType', StringParam, label="Mask symmetry", default="i1", condition="postSymmetryWithinMask",
                           help='If See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                           'If no symmetry is present, give c1')
        groupSymmetry.addParam('postSymmetryWithinMaskMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                               help='The mask values must be between 0 (remove these pixels) and 1 (let them pass). Smooth masks are recommended.')
        groupSymmetry.addParam('postSymmetryHelical', BooleanParam, label="Apply helical symmetry?", default=False)
        groupSymmetry.addParam('postSymmetryHelicalRadius', IntParam, label="Radius", default=-1, condition='postSymmetryHelical',
                               help="In voxels")
        groupSymmetry.addParam('postSymmetryHelicalDihedral', BooleanParam, label="Dihedral symmetry", default=False,
                               condition='postSymmetryHelical')
        groupSymmetry.addParam('postSymmetryHelicalMinRot', FloatParam, label="Min. Rotation", default=0, condition='postSymmetryHelical',
                               help="In degrees")
        groupSymmetry.addParam('postSymmetryHelicalMaxRot', FloatParam, label="Max. Rotation", default=360, condition='postSymmetryHelical',
                               help="In degrees")
        groupSymmetry.addParam('postSymmetryHelicalMinZ', FloatParam, label="Min. Z shift", default=0, condition='postSymmetryHelical',
                               help="In angstroms")
        groupSymmetry.addParam('postSymmetryHelicalMaxZ', FloatParam, label="Max. Z shift", default=40, condition='postSymmetryHelical',
                               help="In angstroms")
        form.addParam('postDoPseudo', BooleanParam, label="Convert to pseudoatoms", default=False)
        form.addParam('postPseudoRadius', FloatParam, label="Pseudoatoms radius", default=1.2, condition="postDoPseudo",
                      expertLevel=LEVEL_ADVANCED, help="In voxels")
        form.addParam('postScript', StringParam, label="Post-processing command", default="", expertLevel=LEVEL_ADVANCED, 
                      help='A command template that is used to post-process the reconstruction. The following variables can be used ' 
                           '%(sampling)s %(dim)s %(volume)s %(iterDir)s. The command should read Spider volumes and modify the input volume.'
                           'the command should be accessible either from the PATH or provide the absolute path.\n'
                           'Examples: \n'
                           'xmipp_transform_filter -i %(volume)s --fourier low_pass 15 --sampling %(sampling)s\n' 
                           '/home/joe/myScript %(volume)s sampling=%(sampling)s dim=%(dim)s')

        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.imgsFn=self._getExtraPath('images.xmd')
        if self.doContinue:
            self.copyAttributes(self.continueRun.get(), 'inputParticles')
            self.copyAttributes(self.continueRun.get(), 'phaseFlipped')
            self.copyAttributes(self.continueRun.get(), 'particleRadius')
            self._insertFunctionStep('copyBasicInformation')
            firstIteration=self.getNumberOfPreviousIterations()+1
        else:
            self._insertFunctionStep('convertInputStep', self.inputParticles.getObjId())
            if self.weightSSNR:
                self._insertFunctionStep('doWeightSSNR')
            self._insertFunctionStep('doIteration000', self.inputVolumes.getObjId())
            firstIteration=1
        self.TsOrig=self.inputParticles.get().getSamplingRate()
        for self.iteration in range(firstIteration,firstIteration+self.numberOfIterations.get()):
            self.insertIteration(self.iteration)
    
    def insertIteration(self,iteration):
        self._insertFunctionStep('globalAssignment',iteration)
        self._insertFunctionStep('localAssignment',iteration)
        self._insertFunctionStep('weightParticles',iteration)
        self._insertFunctionStep('reconstruct',iteration)
        self._insertFunctionStep('postProcessing',iteration)
        self._insertFunctionStep('evaluateReconstructions',iteration)
        self._insertFunctionStep('cleanDirectory',iteration)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, inputParticlesId):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn)
        self.runJob('xmipp_metadata_utilities','-i %s --fill image1 constant noImage'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "image1=image"'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate rename_column "itemId particleId"'%self.imgsFn,numberOfMpi=1)
        imgsFnId=self._getExtraPath('imagesId.xmd')
        self.runJob('xmipp_metadata_utilities','-i %s --operate keep_column particleId -o %s'%(self.imgsFn,imgsFnId),numberOfMpi=1)

    def getNumberOfPreviousIterations(self):
        from glob import glob
        fnDirs=sorted(glob(self.continueRun.get()._getExtraPath("Iter???")))
        lastDir=fnDirs[-1]
        return int(lastDir[-3:])

    def copyBasicInformation(self):
        previousRun=self.continueRun.get()
        copyFile(previousRun._getExtraPath('images.xmd'),self._getExtraPath('images.xmd'))
        copyFile(previousRun._getExtraPath('imagesId.xmd'),self._getExtraPath('imagesId.xmd'))
        if previousRun.weightSSNR:
            copyFile(previousRun._getExtraPath('ssnrWeights.xmd'),self._getExtraPath('ssnrWeights.xmd'))
        
        lastIter=self.getNumberOfPreviousIterations()
        for i in range(0,lastIter+1):
            createLink(previousRun._getExtraPath("Iter%03d"%i),join(self._getExtraPath("Iter%03d"%i)))

    def doWeightSSNR(self):
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        self.runJob("xmipp_image_ssnr", "-i %s -R %d --sampling %f --normalizessnr"%\
                    (self.imgsFn,R,self.inputParticles.get().getSamplingRate()))
        self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate keep_column "particleId weightSSNR" '%\
                    (self.imgsFn,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
        
    def doIteration000(self, inputVolumesId):
        fnDirCurrent=self._getExtraPath('Iter000')
        makePath(fnDirCurrent)
        
        # Split data
        self.runJob("xmipp_metadata_split","-i %s --oroot %s/images -n 2"%(self.imgsFn,fnDirCurrent),numberOfMpi=1)
        for i in range(1,3):
            moveFile("%s/images%06d.xmd"%(fnDirCurrent,i),"%s/images%02d.xmd"%(fnDirCurrent,i))
        
        # Get volume sampling rate
        TsCurrent=self.inputVolumes.get().getSamplingRate()
        self.writeInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE,TsCurrent)

        # Copy reference volumes and window if necessary
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        self.writeInfoField(fnDirCurrent,"size",MDL_XSIZE,newXdim)
        
        img = ImageHandler()
        if isinstance(self.inputVolumes.get(),SetOfVolumes):
            i=1
            for vol in self.inputVolumes.get():
                fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
                img.convert(vol, fnVol)
                if newXdim!=vol.getDim()[0]:
                    self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnVol,newXdim),numberOfMpi=1)
                i+=1
        else:
            fnVol1=join(fnDirCurrent,"volume%02d.vol"%1)
            fnVol2=join(fnDirCurrent,"volume%02d.vol"%2)
            vol=self.inputVolumes.get()
            img.convert(vol, fnVol1)
            if newXdim!=vol.getDim()[0]:
                self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnVol1,newXdim),numberOfMpi=1)
            self.runJob('xmipp_transform_randomize_phases',"-i %s -o %s --freq discrete 0.25"%(fnVol1,fnVol2),numberOfMpi=1)
        
        # Compare both reconstructions
        self.evaluateReconstructions(0)

    def evaluateReconstructions(self,iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        fnVol1=join(fnDirCurrent,"volume%02d.vol"%1)
        fnVol2=join(fnDirCurrent,"volume%02d.vol"%2)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE)
        
        # Align volumes
        letter=self.symmetryGroup.get()[0]
        fnVolAvg=join(fnDirCurrent,"volumeAvg.mrc")
        self.runJob('xmipp_image_operate','-i %s --plus %s -o %s'%(fnVol1,fnVol2,fnVolAvg),numberOfMpi=1)
        self.runJob('xmipp_image_operate','-i %s --mult 0.5'%fnVolAvg,numberOfMpi=1)
        if letter!='i' and letter!='d':
            self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnVolAvg,fnVol1),numberOfMpi=1)
            self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnVolAvg,fnVol2),numberOfMpi=1)
     
        # Estimate resolution
        fnFsc=join(fnDirCurrent,"fsc.xmd")
        self.runJob('xmipp_resolution_fsc','--ref %s -i %s -o %s --sampling_rate %f'%(fnVol1,fnVol2,fnFsc,TsCurrent),numberOfMpi=1)
        fnBeforeVol1=join(fnDirCurrent,"volumeBeforePostProcessing%02d.vol"%1)
        fnBeforeVol2=join(fnDirCurrent,"volumeBeforePostProcessing%02d.vol"%2)
        if exists(fnBeforeVol1) and exists(fnBeforeVol2):
            fnBeforeFsc=join(fnDirCurrent,"fscBeforePostProcessing.xmd")
            self.runJob('xmipp_resolution_fsc','--ref %s -i %s -o %s --sampling_rate %f'%(fnBeforeVol1,fnBeforeVol2,fnBeforeFsc,TsCurrent),
                        numberOfMpi=1)
        md = MetaData(fnFsc)
        resolution=2*TsCurrent
        for objId in md:
            fsc = md.getValue(MDL_RESOLUTION_FRC,objId)
            if fsc<self.nextResolutionCriterion.get():
                resolution=md.getValue(MDL_RESOLUTION_FREQREAL,objId)
                break
        self.writeInfoField(fnDirCurrent,"resolution",MDL_RESOLUTION_FREQREAL,resolution)
        
        # Filter the average to that resolution
        self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%(fnVolAvg,resolution,TsCurrent),numberOfMpi=1)
        self.runJob('xmipp_image_header','-i %s --sampling_rate %f'%(fnVolAvg,TsCurrent),numberOfMpi=1)
        
        # A little bit of statistics (accepted and rejected particles, number of directions, ...)
        if iteration>0:
            from xmipp import AGGR_MAX
            for i in range(1,3):
                fnAnglesi = join(fnDirCurrent,"angles%02d.xmd"%i)
                mdAngles = MetaData(fnAnglesi)
                mdUnique    = MetaData()
                mdUnique.aggregateMdGroupBy(mdAngles, AGGR_MAX, [MDL_PARTICLE_ID], MDL_WEIGHT, MDL_WEIGHT) 
                mdUnique.sort(MDL_PARTICLE_ID)
                fnAnglesUnique = join(fnDirCurrent,"imagesUsed%02d.xmd"%i)
                mdUnique.write(fnAnglesUnique)
    
            fnUsed=join(fnDirCurrent,"imagesUsed.xmd")
            fnUsed1=join(fnDirCurrent,"imagesUsed01.xmd")
            fnUsed2=join(fnDirCurrent,"imagesUsed02.xmd")
            self.runJob('xmipp_metadata_utilities',"-i %s --set union_all %s -o %s"%(fnUsed1,fnUsed2,fnUsed),numberOfMpi=1)
            cleanPath(fnUsed1)
            cleanPath(fnUsed2)
            fnAngles=join(fnDirCurrent,"angles.xmd")
            fnUsedId=join(fnDirCurrent,"imagesUsedId.xmd")
            self.runJob('xmipp_metadata_utilities',"-i %s --operate keep_column particleId -o %s"%(fnUsed,fnUsedId),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities',"-i %s --set natural_join %s"%(fnUsed,fnAngles),numberOfMpi=1)
    
            fnImages=self._getExtraPath("images.xmd")
            fnImagesId=self._getExtraPath('imagesId.xmd')
            fnImagesRejected=join(fnDirCurrent,"imagesRejected.xmd")
            self.runJob('xmipp_metadata_utilities',"-i %s --set subtraction %s particleId -o %s"%(fnImagesId,fnUsedId,fnImagesRejected),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities',"-i %s --set natural_join %s"%(fnImagesRejected,fnImages),numberOfMpi=1)
            cleanPath(fnUsedId)
    
            from pyworkflow.em.metadata.utils import getSize
            Nimages=getSize(fnImages)
            Nrepeated=getSize(join(fnDirCurrent,"angles.xmd"))
            Nunique=getSize(fnUsed)
            Nrejected=getSize(fnImagesRejected)
            
            fh=open(join(fnDirCurrent,"statistics.txt"),'w')
            fh.write("Number of input    images: %d\n"%Nimages)
            fh.write("Number of used     images: %d\n"%Nunique)
            fh.write("Number of rejected images: %d\n"%Nrejected)
            fh.write("Average number of directions per used image: %f\n"%(float(Nrepeated)/Nunique))
            fh.close()
    
    def readInfoField(self,fnDir,block,label):
        md = MetaData("%s@%s"%(block,join(fnDir,"iterInfo.xmd")))
        return md.getValue(label,md.firstObject())

    def writeInfoField(self,fnDir,block,label, value):
        md = MetaData()
        objId=md.addObject()
        md.setValue(label,value,objId)
        md.write("%s@%s"%(block,join(fnDir,"iterInfo.xmd")),MD_APPEND)
    
    def prepareImages(self,fnDirPrevious,fnDir,TsCurrent):
        print "Preparing images to sampling rate=",TsCurrent
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        if newXdim<40:
            newXdim=long(40)
            TsCurrent=Xdim*(self.TsOrig/newXdim)
        self.writeInfoField(fnDir,"sampling",MDL_SAMPLINGRATE,TsCurrent)
        self.writeInfoField(fnDir,"size",MDL_XSIZE,newXdim)
        
        # Prepare particles
        fnDir0=self._getExtraPath("Iter000")
        fnNewParticles=join(fnDir,"images.stk")
        if newXdim!=Xdim:
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(self.imgsFn,fnNewParticles,newXdim))
        else:
            self.runJob("xmipp_image_convert","-i %s -o %s --save_metadata_stack %s"%(self.imgsFn,fnNewParticles,join(fnDir,"images.xmd")),
                        numberOfMpi=1)
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        R=min(round(R*self.TsOrig/TsCurrent*(1+self.angularMaxShift.get()*0.01)),newXdim/2)
        self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnNewParticles,R))
        fnSource=join(fnDir,"images.xmd")
        for i in range(1,3):
            fnImagesi=join(fnDir,"images%02d.xmd"%i)
            self.runJob('xmipp_metadata_utilities','-i %s --set intersection %s/images%02d.xmd particleId particleId -o %s'%\
                        (fnSource,fnDir0,i,fnImagesi),numberOfMpi=1)
        cleanPath(fnSource)
        
    def prepareReferences(self,fnDirPrevious,fnDir,TsCurrent,targetResolution):
        print "Preparing references to sampling rate=",TsCurrent
        fnMask=''
        newXdim=self.readInfoField(fnDir,"size",MDL_XSIZE)
        if self.nextMask.hasValue():
            fnMask=join(fnDir,"mask.vol")
            self.prepareMask(self.nextMask.get(), fnMask, TsCurrent, newXdim)
        TsPrevious=self.readInfoField(fnDirPrevious,"sampling",MDL_SAMPLINGRATE)
        for i in range(1,3):
            fnPreviousVol=join(fnDirPrevious,"volume%02d.vol"%i)
            fnReferenceVol=join(fnDir,"volumeRef%02d.vol"%i)
            if TsPrevious!=TsCurrent:
                self.runJob("xmipp_image_resize","-i %s -o %s --dim %d"%(fnPreviousVol,fnReferenceVol,newXdim),numberOfMpi=1)
            else:
                if self.nextLowPass or self.nextSpherical or self.nextPositivity or fnMask!='':
                    copyFile(fnPreviousVol, fnReferenceVol)
                else:
                    createLink(fnPreviousVol, fnReferenceVol)
            self.runJob('xmipp_transform_filter','-i %s --fourier fsc %s --sampling %f'%(fnReferenceVol,join(fnDirPrevious,"fsc.xmd"),TsCurrent))
            if self.nextLowPass:
                self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%\
                            (fnReferenceVol,targetResolution+self.nextResolutionOffset.get(),TsCurrent),numberOfMpi=1)
            if self.nextSpherical:
                R=self.particleRadius.get()
                if R<=0:
                    R=self.inputParticles.get().getDimensions()[0]/2
                self.runJob('xmipp_transform_mask','-i %s --mask circular -%d'%\
                            (fnReferenceVol,round(R*self.TsOrig/TsCurrent)),numberOfMpi=1)
            if self.nextPositivity:
                self.runJob('xmipp_transform_threshold','-i %s --select below 0 --substitute value 0'%fnReferenceVol,numberOfMpi=1)
            if fnMask!='':
                self.runJob('xmipp_image_operate','-i %s --mult %s'%(fnReferenceVol,fnMask),numberOfMpi=1)
            if self.nextReferenceScript!="":
                scriptArgs = {'volume': fnReferenceVol,
                              'sampling': TsCurrent,
                              'dim': newXdim,
                              'iterDir': fnDir}
                cmd = self.nextReferenceScript % scriptArgs
                self.runJob(cmd, '', numberOfMpi=1)
            
        if fnMask!='':
            cleanPath(fnMask)

    def prepareMask(self,maskObject,fnMask,TsMaskOut,XdimOut):
        img=ImageHandler()
        img.convert(maskObject, fnMask)
        self.runJob('xmipp_image_resize',"-i %s --factor %f"%(fnMask,maskObject.getSamplingRate()/TsMaskOut),numberOfMpi=1)
        maskXdim, _, _, _ =img.getDimensions((1,fnMask))
        if XdimOut!=maskXdim:
            self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnMask,XdimOut),numberOfMpi=1)

    def calculateAngStep(self,newXdim,TsCurrent,ResolutionAlignment):
        k=newXdim*TsCurrent/ResolutionAlignment # Freq. index
        from math import atan2,pi
        return 2*atan2(1,k)*180.0/pi # Corresponding angular step

    def globalAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        makePath(fnDirCurrent)
        previousResolution=self.readInfoField(fnDirPrevious,"resolution",MDL_RESOLUTION_FREQREAL)

        if iteration==1 or previousResolution>=self.significantMaxResolution.get():
            fnGlobal=join(fnDirCurrent,"globalAssignment")
            makePath(fnGlobal)
    
            targetResolution=max(previousResolution,self.significantMaxResolution.get())
            TsCurrent=max(self.TsOrig,targetResolution/3)
            self.prepareImages(fnDirPrevious,fnGlobal,TsCurrent)
            self.prepareReferences(fnDirPrevious,fnGlobal,TsCurrent,targetResolution)

            # Calculate angular step at this resolution
            ResolutionAlignment=previousResolution
            if self.nextLowPass:
                ResolutionAlignment+=self.nextResolutionOffset.get()
            newXdim=self.readInfoField(fnGlobal,"size",MDL_XSIZE)
            angleStep=self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)
            self.writeInfoField(fnGlobal,"angleStep",MDL_ANGLE_DIFF,float(angleStep))
            
            # Significant alignment
            alpha=1-0.01*self.significantSignificance.get()
            for i in range(1,3):
                fnDirSignificant=join(fnGlobal,"significant%02d"%i)
                fnImgs=join(fnGlobal,"images%02d.xmd"%i)
                makePath(fnDirSignificant)

                # Create defocus groups
                self.runJob("xmipp_ctf_group","--ctfdat %s -o %s/ctf:stk --pad 2.0 --sampling_rate %f --phase_flipped  --error 0.1 --resol %f"%\
                            (fnImgs,fnDirSignificant,TsCurrent,targetResolution),numberOfMpi=1)
                moveFile("%s/ctf_images.sel"%fnDirSignificant,"%s/ctf_groups.xmd"%fnDirSignificant)
                cleanPath("%s/ctf_split.doc"%fnDirSignificant)
                md = MetaData("numberGroups@%s"%join(fnDirSignificant,"ctfInfo.xmd"))
                fnCTFs="%s/ctf_ctf.stk"%fnDirSignificant
                numberGroups=md.getValue(MDL_COUNT,md.firstObject())

                # Generate projections
                fnReferenceVol=join(fnGlobal,"volumeRef%02d.vol"%i)
                fnGallery=join(fnDirSignificant,"gallery%02d.stk"%i)
                fnGalleryMd=join(fnDirSignificant,"gallery%02d.xmd"%i)
                args="-i %s -o %s --sampling_rate %f --sym %s --min_tilt_angle %f --max_tilt_angle %f"%\
                     (fnReferenceVol,fnGallery,angleStep,self.symmetryGroup,self.angularMinTilt.get(),self.angularMaxTilt.get())
                self.runJob("xmipp_angular_project_library",args)
                cleanPath(join(fnDirSignificant,"gallery_angles%02d.doc"%i))
                moveFile(join(fnDirSignificant,"gallery%02d.doc"%i), fnGalleryMd)

                fnAngles=join(fnGlobal,"anglesDisc%02d.xmd"%i)
                for j in range(1,numberGroups+1):
                    fnAnglesGroup=join(fnDirSignificant,"angles_group%02d.xmd"%j)
                    if not exists(fnAnglesGroup):
                        fnGroup="ctfGroup%06d@%s/ctf_groups.xmd"%(j,fnDirSignificant)
                        fnGalleryGroup=join(fnDirSignificant,"gallery_group%06d.stk"%j)
                        fnGalleryGroupMd=join(fnDirSignificant,"gallery_group%06d.xmd"%j)
                        self.runJob("xmipp_transform_filter",
                                    "-i %s -o %s --fourier binary_file %d@%s --save_metadata_stack %s --keep_input_columns"%\
                                    (fnGalleryMd,fnGalleryGroup,j,fnCTFs,fnGalleryGroupMd))
                        args='-i %s --initgallery %s --odir %s --sym %s --iter 1 --alpha0 %f --alphaF %f --angularSampling %f --maxShift %d '\
                             '--minTilt %f --maxTilt %f --useImed --angDistance %f --dontReconstruct'%\
                             (fnGroup,fnGalleryGroupMd,fnDirSignificant,self.symmetryGroup,alpha,alpha,angleStep,\
                              round(self.angularMaxShift.get()*newXdim/100),self.angularMinTilt.get(),self.angularMaxTilt.get(),2*angleStep)
                        self.runJob('xmipp_reconstruct_significant',args)
                        moveFile(join(fnDirSignificant,"angles_iter001_00.xmd"),join(fnDirSignificant,"angles_group%02d.xmd"%j))
                        self.runJob("rm -f",fnDirSignificant+"/images_*iter00?_*.xmd",numberOfMpi=1)
                        if j==1:
                            copyFile(fnAnglesGroup, fnAngles)
                        else:
                            self.runJob("xmipp_metadata_utilities","-i %s --set union %s"%(fnAngles,fnAnglesGroup),numberOfMpi=1)
                if self.saveSpace:
                    self.runJob("rm -f",fnDirSignificant+"/gallery*",numberOfMpi=1)
                
                if self.significantGrayValues and previousResolution>self.continuousMinResolution.get():
                    fnRefinedStack=join(fnDirSignificant,"imagesRefined%02d.stk"%i)
                    self.runJob("xmipp_angular_continuous_assign2","-i %s -o %s --ref %s --optimizeGray"%\
                                (fnAngles,fnRefinedStack,fnReferenceVol))
                    fnRefinedXmd=join(fnDirSignificant,"imagesRefined%02d.xmd"%i)
                    moveFile(fnRefinedXmd,fnAngles)

    def adaptShifts(self, fnSource, TsSource, fnDest, TsDest):
        K=TsSource/TsDest
        copyFile(fnSource,fnDest)
        row=getFirstRow(fnDest)
        if row.containsLabel(MDL_SHIFT_X):
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "shiftX=%f*shiftX"'%(fnDest,K),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "shiftY=%f*shiftY"'%(fnDest,K),numberOfMpi=1)
        if row.containsLabel(MDL_CONTINUOUS_X):
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "continuousX=%f*continuousX"'%(fnDest,K),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "continuousY=%f*continuousY"'%(fnDest,K),numberOfMpi=1)

    def localAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        previousResolution=self.readInfoField(fnDirPrevious,"resolution",MDL_RESOLUTION_FREQREAL)
        if previousResolution<=self.continuousMinResolution.get():
            fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
            fnDirLocal=join(fnDirCurrent,"localAssignment")
            makePath(fnDirLocal)

            targetResolution=previousResolution
            TsCurrent=max(self.TsOrig,targetResolution/3)
            self.writeInfoField(fnDirLocal,"sampling",MDL_SAMPLINGRATE,TsCurrent)
            TsCurrent=self.readInfoField(fnDirLocal,"sampling",MDL_SAMPLINGRATE) # Write and read to guarantee consistency with previous directories 
            
            # Prepare images and references
            produceNewReferences=True
            fnDirGlobal=join(fnDirCurrent,"globalAssignment")
            if exists(fnDirGlobal):
                TsGlobal=self.readInfoField(fnDirGlobal,"sampling",MDL_SAMPLINGRATE)
                if TsGlobal==TsCurrent:
                    produceNewReferences=False
            if produceNewReferences:
                self.prepareImages(fnDirPrevious,fnDirLocal,TsCurrent)
                self.prepareReferences(fnDirPrevious,fnDirLocal,TsCurrent,targetResolution)
            else:
                newXdim=self.readInfoField(fnDirGlobal,"size",MDL_XSIZE)
                self.writeInfoField(fnDirLocal,"size",MDL_XSIZE,newXdim)
                for i in range(1,3):
                    createLink(join(fnDirGlobal,"images%02d.xmd"%i),join(fnDirLocal,"images%02d.xmd"%i))
                    createLink(join(fnDirGlobal,"volumeRef%02d.vol"%i),join(fnDirLocal,"volumeRef%02d.vol"%i))

            # Compute maximum angular deviation
            ResolutionAlignment=previousResolution
            if self.nextLowPass:
                ResolutionAlignment+=self.nextResolutionOffset.get()
            newXdim=self.readInfoField(fnDirLocal,"size",MDL_XSIZE)
            maxAngle=3*self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)

            for i in range(1,3):
                fnLocalXmd=join(fnDirLocal,"anglesCont%02d.xmd"%i)
                if not exists(fnLocalXmd):
                    fnLocalImages=join(fnDirLocal,"images%02d.xmd"%i)
    
                    # Starting angles
                    fnLocalAssignment=join(fnDirLocal,"anglesDisc%02d.xmd"%i)
                    if exists(fnDirGlobal):
                        fnGlobalAssignment=join(fnDirGlobal,"anglesDisc%02d.xmd"%i)
                        TsGlobal=self.readInfoField(fnDirGlobal,"sampling",MDL_SAMPLINGRATE)
                        if TsGlobal==TsCurrent:
                            copyFile(fnGlobalAssignment,fnLocalAssignment)
                        else:
                            self.adaptShifts(fnGlobalAssignment,TsGlobal,fnLocalAssignment,TsCurrent)
                    else:
                        TsPrevious=self.readInfoField(fnDirPrevious,"sampling",MDL_SAMPLINGRATE)
                        self.adaptShifts(join(fnDirPrevious,"angles%02d.xmd"%i),TsPrevious,fnLocalAssignment,TsCurrent)
                    self.runJob("xmipp_metadata_utilities","-i %s --operate drop_column image"%fnLocalAssignment,numberOfMpi=1)
                    self.runJob("xmipp_metadata_utilities","-i %s --set join %s particleId"%(fnLocalAssignment,fnLocalImages),numberOfMpi=1)
    
                    fnVol=join(fnDirLocal,"volumeRef%02d.vol"%i)
                    fnLocalStk=join(fnDirLocal,"anglesCont%02d.stk"%i)
                    
                    R=self.particleRadius.get()
                    if R<=0:
                        R=self.inputParticles.get().getDimensions()[0]/2
                    R=round(R*self.TsOrig/TsCurrent)
                    args="-i %s -o %s --sampling %f --Rmax %d --padding %d --ref %s --max_resolution %f --applyTo image1"%\
                       (fnLocalAssignment,fnLocalStk,TsCurrent,R,self.contPadding.get(),fnVol,previousResolution)
                    if self.contShift:
                        args+=" --optimizeShift --max_shift %d"%round(self.angularMaxShift.get()*newXdim*0.01)
                    if self.contScale:
                        args+=" --optimizeScale --max_scale %f"%self.contMaxScale.get() 
                    if self.contAngles:
                        args+=" --optimizeAngles --max_angular_change %f"%maxAngle
                    if self.contGrayValues:
                        args+=" --optimizeGray"
                    if self.contDefocus:
                        args+=" --optimizeDefocus"
                        if self.phaseFlipped:
                            args+=" --phaseFlipped"
                    if self.weightResiduals:
                        args+=" --oresiduals %s"%join(fnDirLocal,"residuals%02i.stk"%i)
                    self.runJob("xmipp_angular_continuous_assign2",args)
                    self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnLocalStk,R))

    def weightParticles(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        from math import exp
        for i in range(1,3):
            # Grab file
            fnDirGlobal=join(fnDirCurrent,"globalAssignment")
            fnDirLocal=join(fnDirCurrent,"localAssignment")
            fnAnglesCont=join(fnDirLocal,"anglesCont%02d.xmd"%i)
            fnAnglesDisc=join(fnDirGlobal,"anglesDisc%02d.xmd"%i)
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            if exists(fnAnglesCont):
                copyFile(fnAnglesCont, fnAngles)
                TsCurrent=self.readInfoField(fnDirLocal,"sampling",MDL_SAMPLINGRATE)
                Xdim=self.readInfoField(fnDirLocal,"size",MDL_XSIZE)
            else:
                if exists(fnAnglesDisc):
                    copyFile(fnAnglesDisc, fnAngles)
                    TsCurrent=self.readInfoField(fnDirGlobal,"sampling",MDL_SAMPLINGRATE)
                    Xdim=self.readInfoField(fnDirGlobal,"size",MDL_XSIZE)
                else:
                    raise Exception("Angles for iteration "+str(iteration)+" not found")
            self.writeInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE,TsCurrent)
            self.writeInfoField(fnDirCurrent,"size",MDL_XSIZE,Xdim)
                
            if self.weightSSNR:
                self.runJob("xmipp_metadata_utilities","-i %s --set join %s particleId"%\
                            (fnAngles,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
            if self.weightJumper and iteration>1:
                fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
                fnPreviousAngles=join(fnDirPrevious,"angles%02d.xmd"%i)
                self.runJob("xmipp_angular_distance","--ang1 %s --ang2 %s --compute_weights --oroot %s"%\
                            (fnPreviousAngles,fnAngles,fnDirCurrent+"/jumper"),numberOfMpi=1)
                moveFile(fnDirCurrent+"/jumper_weights.xmd", fnAngles)

            if self.weightResiduals and exists(fnAnglesCont):
                fnCovariance=join(fnDirLocal,"covariance%02d.stk"%i)
                self.runJob("xmipp_image_residuals","-i %s -o %s --normalizeDivergence"%(fnAngles,fnCovariance),numberOfMpi=1)
                moveFile(join(fnDirLocal,"covariance%02d.xmd"%i),fnAngles)
            
            md=MetaData(fnAngles)
            doWeightSignificance=self.weightSignificance and md.containsLabel(MDL_WEIGHT_SIGNIFICANT)
            for objId in md:
                weight=1
                if self.weightSSNR:
                    aux=md.getValue(MDL_WEIGHT_SSNR,objId)
                    weight*=aux
                if doWeightSignificance:
                    aux=md.getValue(MDL_WEIGHT_SIGNIFICANT,objId)
                    weight*=aux
                if self.weightContinuous and exists(fnAnglesCont):
                    aux=md.getValue(MDL_WEIGHT_CONTINUOUS2,objId)
                    weight*=aux
                if self.weightResiduals and exists(fnAnglesCont):
                    aux=md.getValue(MDL_ZSCORE_RESCOV,objId)
                    aux/=3
                    weight*=exp(-0.5*aux*aux)
                    aux=md.getValue(MDL_ZSCORE_RESMEAN,objId)
                    aux/=3
                    weight*=exp(-0.5*aux*aux)
                    aux=md.getValue(MDL_ZSCORE_RESVAR,objId)
                    aux/=3
                    weight*=exp(-0.5*aux*aux)
                if self.weightJumper and iteration>1:
                    aux=md.getValue(MDL_WEIGHT_JUMPER,objId)
                    weight*=aux
                md.setValue(MDL_WEIGHT,weight,objId)
            md.write(fnAngles)
            
            # Weight by Astigmatism
            if self.weightAstigmatism:
                self.runJob("xmipp_transform_filter","-i %s --fourier astigmatism %f --sampling %f"%\
                            (fnAngles,self.weightAstigmatismSigma.get(),TsCurrent))
        fnAngles=join(fnDirCurrent,"angles.xmd")
        fnAngles1=join(fnDirCurrent,"angles01.xmd")
        fnAngles2=join(fnDirCurrent,"angles02.xmd")
        self.runJob('xmipp_metadata_utilities',"-i %s --set union_all %s -o %s"%(fnAngles1,fnAngles2,fnAngles),numberOfMpi=1)

    def reconstruct(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE)
        for i in range(1,3):
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
            # Reconstruct Fourier
            args="-i %s -o %s --sym %s --weight"%(fnAngles,fnVol,self.symmetryGroup)
            row=getFirstRow(fnAngles)
            if row.containsLabel(MDL_CTF_DEFOCUSU) or row.containsLabel(MDL_CTF_MODEL):
                args+=" --useCTF --sampling %f"%TsCurrent
                if self.phaseFlipped:
                    args+=" --phaseFlipped"
            self.runJob("xmipp_reconstruct_fourier",args)
            # Reconstruct ADMM
#            args="-i %s -o %s --sym %s"%(fnAngles,fnVol,self.symmetryGroup)
#            row=getFirstRow(fnAngles)
#            if row.containsLabel(MDL_CTF_DEFOCUSU) or row.containsLabel(MDL_CTF_MODEL):
#                if not self.phaseFlipped:
#                    args+=" --dontUseCTF"
#            self.runJob("xmipp_reconstruct_admm",args)
    
    def postProcessing(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE)
        if not self.postBFactor and not self.postNonnegativity and not self.postMask and not self.postSymmetryWithinMask \
           and not self.postSymmetryHelical and not self.postDoPseudo and self.postScript=="":
            return
        for i in range(1,3):
            fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
            fnBeforeVol=join(fnDirCurrent,"volumeBeforePostProcessing%02d.vol"%i)
            copyFile(fnVol,fnBeforeVol)
            volXdim = self.readInfoField(fnDirCurrent, "size", MDL_XSIZE)
            
            if self.postBFactor:
                fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
                previousResolution=self.readInfoField(fnDirPrevious, "resolution", MDL_RESOLUTION_FREQREAL)
                if previousResolution<self.postBmin.get():
                    fnAuxVol=join(fnDirCurrent,"volumeAux%02d.vol"%i)
                    args="-i %s --sampling %f --auto --fit_minres %f --fit_maxres %f --maxres %f -o %s"%\
                       (fnVol,TsCurrent,self.postBmin.get(),self.postBmax.get(),2*TsCurrent,fnAuxVol)
                    self.runJob("xmipp_volume_correct_bfactor",args,numberOfMpi=1)
                    moveFile(fnAuxVol,fnVol)
                    moveFile(fnAuxVol+".guinier",fnVol+".guinier")
            
            if self.postNonnegativity:
                self.runJob("xmipp_transform_threshold","-i %s --select below 0 --substitute value 0"%fnVol,numberOfMpi=1)
            
            if self.postAdHocMask.hasValue():
                fnMask=join(fnDirCurrent,"mask.vol")
                self.prepareMask(self.postAdHocMask.get(), fnMask, TsCurrent, volXdim)
                self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnVol,fnMask),numberOfMpi=1)
                cleanPath(fnMask)

            if self.postMask:
                fnMask=join(fnDirCurrent,"mask%02d.vol"%i)
                V=Image(fnVol)
                [_,std,_,_]=V.computeStats()
                threshold=self.postMaskThreshold.get()*std
                
                self.runJob("xmipp_transform_threshold","-i %s -o %s --select below %f --substitute binarize"%(fnVol,fnMask,threshold),
                            numberOfMpi=1)
                if self.postDoMaskRemoveSmall:
                    self.runJob("xmipp_transform_morphology","-i %s --binaryOperation removeSmall %d"%\
                                (fnMask,self.postMaskRemoveSmallThreshold.get()),numberOfMpi=1)
                if self.postMaskKeepLargest:
                    self.runJob("xmipp_transform_morphology","-i %s --binaryOperation keepBiggest"%fnMask,numberOfMpi=1)
                if self.postDoMaskDilate:
                    self.runJob("xmipp_transform_morphology","-i %s --binaryOperation dilation --size %d"%\
                                (fnMask,self.postMaskDilateSize.get()),numberOfMpi=1)
                if self.postDoMaskSmooth:
                    self.runJob("xmipp_transform_filter","-i %s --fourier real_gaussian %f"%\
                                (fnMask,self.postMaskSmoothSize.get()),numberOfMpi=1)
                if self.postMaskSymmetry!="c1":
                    self.runJob("xmipp_transform_symmetrize","-i %s --sym %s"%\
                                (fnMask,self.postMaskSymmetry.get()),numberOfMpi=1)
                self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnVol,fnMask),numberOfMpi=1)
                cleanPath(fnMask)

        if self.postSymmetryWithinMask:
            if self.postMaskSymmetry!="c1":
                fnMask=join(fnDirCurrent,"mask%02d.vol"%i)
                self.prepareMask(self.postSymmetryWithinMaskMask.get(),fnMask,TsCurrent,volXdim)
                self.runJob("xmipp_transform_symmetrize","-i %s --sym %s --mask_in %s"%\
                            (fnVol,self.postSymmetryWithinMaskType.get(),fnMask),numberOfMpi=1)
                cleanPath(fnMask)
        
        if self.postSymmetryHelical:
            z0=float(self.postSymmetryHelicalMinZ.get())
            zF=float(self.postSymmetryHelicalMaxZ.get())
            zStep=(zF-z0)/10
            rot0=float(self.postSymmetryHelicalMinRot.get())
            rotF=float(self.postSymmetryHelicalMaxRot.get())
            rotStep=(rotF-rot0)/10
            fnCoarse=join(fnDirCurrent,"coarseHelical%02d.xmd"%i)
            fnFine=join(fnDirCurrent,"fineHelical%02d.xmd"%i)
            radius=int(self.postSymmetryHelicalRadius.get())
            height=int(volXdim)
            self.runCoarseSearch(fnVol, z0, zF, zStep, rot0, rotF, rotStep, 1, fnCoarse, radius, height)
            self.runFineSearch(fnVol, fnCoarse, fnFine, z0, zF, rot0, rotF, radius, height)
            cleanPath(fnCoarse)
            self.runSymmetrize(fnVol, fnFine, fnVol, radius, height)
            if self.postSymmetryHelicalDihedral:
                self.runApplyDihedral(fnVol, fnFine, join(fnDirCurrent,"rotatedHelix.vol"), radius, height)

        if self.postDoPseudo:
            fnPseudo=join(fnDirCurrent,"pseudo")
            self.runJob("xmipp_volume_to_pseudoatoms","-i %s -o %s --sigma %f --sampling_rate %f --verbose 2 --targetError 1"%\
                        (fnVol,fnPseudo,self.postPseudoRadius.get()*TsCurrent,TsCurrent),numberOfMpi=1)
            moveFile(fnPseudo+"_approximation.vol", fnVol)
            self.runJob("rm","-f "+fnPseudo+"*")
        
        if self.postScript!="":
            img = ImageHandler()
            volXdim, _, _, _ =img.getDimensions((1,fnVol))
            scriptArgs = {'volume': fnVol,
                          'sampling': TsCurrent,
                          'dim': volXdim,
                          'iterDir': fnDirCurrent}
            cmd = self.postScript % scriptArgs
            self.runJob(cmd, '', numberOfMpi=1)

    def cleanDirectory(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        if self.saveSpace:
            fnGlobal=join(fnDirCurrent,"globalAssignment")
            fnLocal=join(fnDirCurrent,"localAssignment")
            if exists(fnGlobal):
                cleanPath(join(fnGlobal,"images.stk"))
            for i in range(1,3):
                if exists(fnGlobal):
                    cleanPath(join(fnGlobal,"images%02d.xmd"%i))
                    cleanPath(join(fnGlobal,"significant%02d"%i))
                    cleanPath(join(fnGlobal,"volumeRef%02d.vol"%i))
                if exists(fnLocal):
                    cleanPath(join(fnLocal,"images%02d.xmd"%i))
                    cleanPath(join(fnLocal,"anglesCont%02d.stk"%i))
                    cleanPath(join(fnLocal,"anglesDisc%02d.xmd"%i))
                    cleanPath(join(fnLocal,"volumeRef%02d.vol"%i))
                    if self.weightResiduals:
                        cleanPath(join(fnDirLocal,"covariance%02d.stk"%i))
                        cleanPath(join(fnDirLocal,"residuals%02i.stk"%i))

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if isinstance(self.inputVolumes.get(),SetOfVolumes) and self.inputVolumes.get().getSize()!=2:
            errors.append("The set of input volumes should have exactly 2 volumes")
        if self.postSymmetryWithinMask and not self.postSymmetryWithinMaskMask.hasValue():
            errors.append("Symmetrize within mask requires a mask")
        if self.significantMaxResolution.get()>=self.continuousMinResolution.get():
            errors.append("There is a gap in resolution for which no angular assignment is performed")
        return errors    
    
    def _summary(self):
        summary = []
        summary.append("Images: %s" % self.inputParticles.getNameId())
        summary.append("Volume(s): %s" % self.inputVolumes.getNameId())
        summary.append("symmetry: %s" % self.symmetryGroup.get())
        return summary
    
    def _methods(self):
        methods = []
        return methods
    
