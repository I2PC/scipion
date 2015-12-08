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

from glob import glob

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam, EnumParam
from pyworkflow.utils.path import cleanPath, makePath, copyFile, moveFile, createLink
from pyworkflow.em.protocol import ProtRefine3D
from pyworkflow.em.data import SetOfVolumes, Volume
from pyworkflow.em.metadata.utils import getFirstRow, getSize
from convert import writeSetOfParticles
from os.path import join, exists, split
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles, setXmippAttributes
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md
import pyworkflow.em as em

import xmipp

from xmipp3 import HelicalFinder



class XmippProtReconstructHighRes(ProtRefine3D, HelicalFinder):
    """This is a 3D refinement protocol whose main input is a volume and a set of particles.
       The set of particles has to be at full size (the finer sampling rate available), but
       the rest of inputs (reference volume and masks) can be at any downsampling factor.
       The protocol scales the input images and volumes to a reasonable size depending on
       the resolution of the previous iteration.
       
       The protocol works with any input volume, whichever its resolution, as long as it
       is a reasonable initial volume for the set of particles. The protocol does not
       resolve the heterogeneous problem (it assumes an homogeneous population),
       although it is somewhat tolerant through the use of particle weights in the
       reconstruction process.
       
       It is recommended to perform several global alignment iterations before entering
       into the local iterations. The switch from global to local should be performed when
       a substantial percentage of the particles do not move from one iteration to the next.
       
       The algorithm reports the cross correlation (global alignment) or cost (local) function
       per defocus group, so that we can see which was the percentile of each particle in its
       defocus group. You may want to perform iterations one by one, and remove from one
       iteration to the next, those particles that worse fit the model."""
    _label = 'highres'
    
    SPLIT_STOCHASTIC = 0
    SPLIT_FIXED = 1
    
    GLOBAL_ALIGNMENT = 0
    LOCAL_ALIGNMENT = 1
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                      'run of type *%s* class and some of the input parameters'
                      'will be taken from it.' % self.getClassName())
        form.addParam('inputParticles', PointerParam, label="Full-size Images", important=True, 
                      pointerClass='SetOfParticles', allowsNull=True,
                      help='Select a set of images at full resolution')
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
        
        form.addSection(label='Next Reference')
        form.addParam('nextResolutionCriterion',FloatParam, label="FSC criterion", default=0.143, 
                      help='The resolution of the reconstruction is defined as the inverse of the frequency at which '\
                      'the FSC drops below this value. Typical values are 0.143 and 0.5')
        form.addParam('nextLowPass', BooleanParam, label="Low pass filter?", default=True,
                      help='Apply a low pass filter to the previous iteration whose maximum frequency is '\
                           'the current resolution(A) + resolutionOffset(A). If resolutionOffset>0, then fewer information' \
                           'is used (meant to avoid overfitting). If resolutionOffset<0, then more information is allowed '\
                           '(meant for a greedy convergence).')
        form.addParam('nextResolutionOffset', FloatParam, label="Resolution offset (A)", default=2, condition='nextLowPass')
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
        form.addParam('splitMethod', EnumParam, label='Image split method', choices=['Stochastic','Fixed'], default=self.SPLIT_FIXED, expertLevel=LEVEL_ADVANCED)
        form.addParam('multiresolution', BooleanParam, label='Multiresolution approach', default=True, expertLevel=LEVEL_ADVANCED,
                      help="In the multiresolution approach the sampling rate of the images is adapted to the current resolution")
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=10,
                      help='Maximum shift as a percentage of the image size')
        line=form.addLine('Tilt angle:', help='0 degrees represent top views, 90 degrees represent side views', expertLevel=LEVEL_ADVANCED)
        line.addParam('angularMinTilt', FloatParam, label="Min.", default=0, expertLevel=LEVEL_ADVANCED)
        line.addParam('angularMaxTilt', FloatParam, label="Max.", default=90, expertLevel=LEVEL_ADVANCED)
        form.addParam('alignmentMethod', EnumParam, label='Image alignment', choices=['Global','Local'], default=self.GLOBAL_ALIGNMENT)

        form.addParam('shiftSearch5d', FloatParam, label="Shift search", default=7.0, condition='alignmentMethod==0',
                  expertLevel=LEVEL_ADVANCED, help="In pixels. The next shift is searched from the previous shift plus/minus this amount.")
        form.addParam('shiftStep5d', FloatParam, label="Shift step", default=2.0, condition='alignmentMethod==0', 
	              expertLevel=LEVEL_ADVANCED, help="In pixels")

        form.addParam('contShift', BooleanParam, label="Optimize shifts?", default=True, condition='alignmentMethod==1',
                      help='Optimize shifts within a limit')
        form.addParam('contMaxShiftVariation', FloatParam, label="Max. shift variation", default=2, condition='alignmentMethod==1', expertLevel=LEVEL_ADVANCED,
                                 help="Percentage of the image size")
        form.addParam('contScale', BooleanParam, label="Optimize scale?", default=True, condition='alignmentMethod==1',
                      help='Optimize scale within a limit')
        form.addParam('contMaxScale', FloatParam, label="Max. scale variation", default=0.02, condition='alignmentMethod==1', expertLevel=LEVEL_ADVANCED)
        form.addParam('contAngles', BooleanParam, label="Optimize angles?", default=True, condition='alignmentMethod==1',
                      help='Optimize angles within a limit')
        form.addParam('contGrayValues', BooleanParam, label="Optimize gray values?", default=False, condition='alignmentMethod==1',
                      help='Optimize gray values. Do not perform this unless the reconstructed volume is gray-compatible with the projections,'\
                      ' i.e., the volumes haven been produced from projections')
        form.addParam('contMaxGrayScale', FloatParam, label="Max. gray scale variation", default=5, condition='alignmentMethod==1', expertLevel=LEVEL_ADVANCED)
        form.addParam('contMaxGrayShift', FloatParam, label="Max. gray shift variation", default=5, condition='alignmentMethod==1', expertLevel=LEVEL_ADVANCED,
                                 help='As a factor of the image standard deviation')
        form.addParam('contDefocus', BooleanParam, label="Optimize defocus?", condition='alignmentMethod==1', default=True)
        form.addParam('contMaxDefocus', FloatParam, label="Max. defocus variation", default=200, condition='alignmentMethod==1', expertLevel=LEVEL_ADVANCED,
                                 help="In Angstroms")
        form.addParam('contPadding', IntParam, label="Fourier padding factor", default=2, condition='alignmentMethod==1', expertLevel=LEVEL_ADVANCED,
                      help='The volume is zero padded by this factor to produce projections')
        
        form.addSection(label='Weights')
        form.addParam('weightSSNR', BooleanParam, label="Weight by SSNR?", default=True,
                      help='Weight input images by SSNR')
        form.addParam('weightContinuous', BooleanParam, label="Weight by Continuous cost?", default=True, condition='alignmentMethod==1',
                      help='Weight input images by angular assignment cost')
        form.addParam('weightJumper', BooleanParam, label="Weight by angular stability?", default=True,
                      help='Weight input images by angular stability between iterations')
        form.addParam('weightCC', BooleanParam, label="Weight by CC percentile?", default=True,
                      help='Weight input images by their fitness (cross correlation) percentile in their defocus group')
        form.addParam('weightCCmin', FloatParam, label="Minimum CC weight", default=0.1, expertLevel=LEVEL_ADVANCED,
                      help='Weights are between this value and 1')
        form.addParam('minCTF', FloatParam, label="Minimum CTF value", default=0.1, expertLevel=LEVEL_ADVANCED,
                      help='A Fourier coefficient is not considered if its CTF is below this value. Note that setting a too low value for this parameter amplifies noise.')
        
        form.addSection(label='Post-processing')
        form.addParam('postAdHocMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be between 0 (remove these pixels) and 1 (let them pass). Smooth masks are recommended.')
        groupSymmetry = form.addGroup('Symmetry')
        groupSymmetry.addParam('postSymmetryWithinMask', BooleanParam, label="Symmetrize volume within mask?", default=False)
        groupSymmetry.addParam('postSymmetryWithinMaskType', StringParam, label="Mask symmetry", default="i1", condition="postSymmetryWithinMask",
                           help='If See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                           'If no symmetry is present, give c1')
        groupSymmetry.addParam('postSymmetryWithinMaskMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True, condition="postSymmetryWithinMask",
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
        form.addParam('postScript', StringParam, label="Post-processing command", default="", expertLevel=LEVEL_ADVANCED, 
                      help='A command template that is used to post-process the reconstruction. The following variables can be used ' 
                           '%(sampling)s %(dim)s %(volume)s %(iterDir)s. The command should read Spider volumes and modify the input volume.'
                           'the command should be accessible either from the PATH or provide the absolute path.\n'
                           'Examples: \n'
                           'xmipp_transform_filter -i %(volume)s --fourier low_pass 15 --sampling %(sampling)s\n' 
                           '/home/joe/myScript %(volume)s sampling=%(sampling)s dim=%(dim)s')

        form.addParallelSection(threads=1, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.imgsFn=self._getExtraPath('images.xmd')
        if self.doContinue:
            self.copyAttributes(self.continueRun.get(), 'particleRadius')
            self.copyAttributes(self.continueRun.get(), 'inputVolumes')
            if not self.inputParticles.hasValue():
                self.copyAttributes(self.continueRun.get(), 'inputParticles')
            else:
                self._insertFunctionStep('convertInputStep', self.inputParticles.getObjId())
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
        self._insertFunctionStep("createOutput")
    
    def insertIteration(self,iteration):
        if self.alignmentMethod==self.GLOBAL_ALIGNMENT:
            self._insertFunctionStep('globalAssignment',iteration)
        else:
            self._insertFunctionStep('localAssignment',iteration)
        self._insertFunctionStep('weightParticles',iteration)
        self._insertFunctionStep('qualifyParticles',iteration)
        self._insertFunctionStep('reconstruct',iteration)
        self._insertFunctionStep('postProcessing',iteration)
        self._insertFunctionStep('evaluateReconstructions',iteration)
        self._insertFunctionStep('cleanDirectory',iteration)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, inputParticlesId):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn)
        self.runJob('xmipp_metadata_utilities','-i %s --fill image1 constant noImage'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "image1=image"'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --fill particleId constant 1'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "particleId=itemId"'%self.imgsFn,numberOfMpi=1)
        imgsFnId=self._getExtraPath('imagesId.xmd')
        self.runJob('xmipp_metadata_utilities','-i %s --operate keep_column particleId -o %s'%(self.imgsFn,imgsFnId),numberOfMpi=1)

    def createOutput(self):
        # get last iteration
        fnIterDir=glob(self._getExtraPath("Iter*"))
        lastIter=len(fnIterDir)-1
        fnLastDir=self._getExtraPath("Iter%03d"%lastIter)
        fnLastVol=join(fnLastDir,"volumeAvg.mrc")
        Ts=self.readInfoField(fnLastDir,"sampling",xmipp.MDL_SAMPLINGRATE)
        if exists(fnLastVol):
            volume=Volume()
            volume.setFileName(fnLastVol)
            volume.setSamplingRate(Ts)
            self._defineOutputs(outputVolume=volume)
            self._defineSourceRelation(self.inputParticles.get(),volume)
            #self._defineSourceRelation(self.inputVolumes.get(),volume)

        fnLastAngles=join(fnLastDir,"angles.xmd")
        if exists(fnLastAngles):
            fnAngles=self._getPath("angles.xmd")
            self.iterMd = md.iterRows(fnAngles, sortByLabel=md.MDL_ITEM_ID)
            self.lastRow = next(self.iterMd)
            self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate modify_values "image=image1"'%(fnLastAngles,fnAngles),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate sort particleId'%fnAngles,numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate drop_column image1'%fnAngles,numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "itemId=particleId"'%fnAngles,numberOfMpi=1)
            imgSet = self.inputParticles.get()
            imgSetOut = self._createSetOfParticles()
            imgSetOut.copyInfo(imgSet)
            imgSetOut.setAlignmentProj()
            imgSetOut.copyItems(imgSet,
                                updateItemCallback=self._createItemMatrix)#,
                                #itemDataIterator=md.iterRows(fnAngles, sortByLabel=md.MDL_ITEM_ID))
            self._defineOutputs(outputParticles=imgSetOut)
            self._defineSourceRelation(self.inputParticles, imgSetOut)

    

    def _createItemMatrix(self, particle, row):
        # We are using here an special type of iteration over the metadata
        # in this case the metadata may contains less elements than the 
        # input set of particles, so row=None here
        row = self.lastRow
        if row is None or particle.getObjId() != row.getValue(xmipp.MDL_ITEM_ID):
            particle._appendItem = False            
        else:
            from pyworkflow.em.packages.xmipp3.convert import createItemMatrix
            
            createItemMatrix(particle, row, align=em.ALIGN_PROJ)
            setXmippAttributes(particle, row, xmipp.MDL_SHIFT_X, xmipp.MDL_SHIFT_Y, xmipp.MDL_ANGLE_TILT, xmipp.MDL_SCALE, xmipp.MDL_MAXCC, xmipp.MDL_MAXCC_PERCENTILE, xmipp.MDL_WEIGHT)
        if row.containsLabel(xmipp.MDL_ANGLE_DIFF0):
            setXmippAttributes(particle, row, xmipp.MDL_ANGLE_DIFF0, xmipp.MDL_WEIGHT_JUMPER0)
        if row.containsLabel(xmipp.MDL_CONTINUOUS_X):
            setXmippAttributes(particle, row, xmipp.MDL_CONTINUOUS_X, xmipp.MDL_CONTINUOUS_Y, xmipp.MDL_COST, xmipp.MDL_WEIGHT_CONTINUOUS2, 
                               xmipp.MDL_CONTINUOUS_SCALE_X, xmipp.MDL_CONTINUOUS_SCALE_Y, xmipp.MDL_COST_PERCENTILE,
                               xmipp.MDL_CONTINUOUS_GRAY_A, xmipp.MDL_CONTINUOUS_GRAY_B)
            if row.containsLabel(xmipp.MDL_ANGLE_DIFF):
                setXmippAttributes(particle, row, xmipp.MDL_ANGLE_DIFF, xmipp.MDL_WEIGHT_JUMPER)
            if row.containsLabel(xmipp.MDL_ANGLE_DIFF2):
                setXmippAttributes(particle, row, xmipp.MDL_ANGLE_DIFF2, xmipp.MDL_WEIGHT_JUMPER2)
            if row.containsLabel(xmipp.MDL_WEIGHT_SSNR):
                setXmippAttributes(particle, row, xmipp.MDL_WEIGHT_SSNR)
                
            try:
                self.lastRow = next(self.iterMd)
            except StopIteration:
                self.lastRow = None

    def getLastFinishedIter(self):
        fnFscs=sorted(glob(self._getExtraPath("Iter???/fsc.xmd")))
        lastDir=split(fnFscs[-1])[0]
        return int(lastDir[-3:])
    
    def getNumberOfPreviousIterations(self):
        fnDirs=sorted(glob(self.continueRun.get()._getExtraPath("Iter???")))
        lastDir=fnDirs[-1]
        return int(lastDir[-3:])

    def copyBasicInformation(self):
        previousRun=self.continueRun.get()
        if not self.inputParticles.hasValue():
            copyFile(previousRun._getExtraPath('images.xmd'),self._getExtraPath('images.xmd'))
            copyFile(previousRun._getExtraPath('imagesId.xmd'),self._getExtraPath('imagesId.xmd'))
        if previousRun.weightSSNR:
            copyFile(previousRun._getExtraPath('ssnrWeights.xmd'),self._getExtraPath('ssnrWeights.xmd'))
        elif self.weightSSNR:
            self.doWeightSSNR()
        
        lastIter=self.getNumberOfPreviousIterations()
        for i in range(0,lastIter+1):
            createLink(previousRun._getExtraPath("Iter%03d"%i),join(self._getExtraPath("Iter%03d"%i)))

    def doWeightSSNR(self):
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        self.runJob("xmipp_image_ssnr", "-i %s -R %d --sampling %f --normalizessnr"%\
                    (self.imgsFn,R,self.inputParticles.get().getSamplingRate()),numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate keep_column "particleId weightSSNR" '%\
                    (self.imgsFn,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
        
    def doIteration000(self, inputVolumesId):
        fnDirCurrent=self._getExtraPath('Iter000')
        makePath(fnDirCurrent)
        
        # Split data
        if self.splitMethod == self.SPLIT_FIXED:
            self.runJob("xmipp_metadata_split","-i %s --oroot %s/images -n 2"%(self.imgsFn,fnDirCurrent),numberOfMpi=1)
            for i in range(1,3):
                moveFile("%s/images%06d.xmd"%(fnDirCurrent,i),"%s/images%02d.xmd"%(fnDirCurrent,i))
        
        # Get volume sampling rate
        TsCurrent=self.inputVolumes.get().getSamplingRate()
        self.writeInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)

        # Copy reference volumes and window if necessary
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        self.writeInfoField(fnDirCurrent,"size",xmipp.MDL_XSIZE,newXdim)
        
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
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
        
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
        mdFSC = xmipp.MetaData(fnFsc)
        resolution=2*TsCurrent
        for objId in mdFSC:
            fsc = mdFSC.getValue(xmipp.MDL_RESOLUTION_FRC,objId)
            if fsc<self.nextResolutionCriterion.get():
                resolution=mdFSC.getValue(xmipp.MDL_RESOLUTION_FREQREAL,objId)
                break
        self.writeInfoField(fnDirCurrent,"resolution",xmipp.MDL_RESOLUTION_FREQREAL,resolution)
        
        # Filter the average to that resolution
        self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%(fnVolAvg,resolution,TsCurrent),numberOfMpi=1)
        self.runJob('xmipp_image_header','-i %s --sampling_rate %f'%(fnVolAvg,TsCurrent),numberOfMpi=1)
        
        # A little bit of statistics (accepted and rejected particles, number of directions, ...)
        if iteration>0:
            from xmipp import AGGR_MAX
            for i in range(1,3):
                fnAnglesi = join(fnDirCurrent,"angles%02d.xmd"%i)
                mdAngles = xmipp.MetaData(fnAnglesi)
                mdUnique    = xmipp.MetaData()
                mdUnique.aggregateMdGroupBy(mdAngles, AGGR_MAX, [xmipp.MDL_PARTICLE_ID], xmipp.MDL_WEIGHT, xmipp.MDL_WEIGHT) 
                mdUnique.sort(xmipp.MDL_PARTICLE_ID)
                fnAnglesUnique = join(fnDirCurrent,"imagesUsed%02d.xmd"%i)
                mdUnique.write(fnAnglesUnique)
    
            fnUsed=join(fnDirCurrent,"imagesUsed.xmd")
            fnUsed1=join(fnDirCurrent,"imagesUsed01.xmd")
            fnUsed2=join(fnDirCurrent,"imagesUsed02.xmd")
            self.runJob('xmipp_metadata_utilities',"-i %s --set union %s -o %s"%(fnUsed1,fnUsed2,fnUsed),numberOfMpi=1)
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
        mdInfo = xmipp.MetaData("%s@%s"%(block,join(fnDir,"iterInfo.xmd")))
        return mdInfo.getValue(label,mdInfo.firstObject())

    def writeInfoField(self,fnDir,block,label, value):
        mdInfo = xmipp.MetaData()
        objId=mdInfo.addObject()
        mdInfo.setValue(label,value,objId)
        mdInfo.write("%s@%s"%(block,join(fnDir,"iterInfo.xmd")),xmipp.MD_APPEND)
    
    def prepareImages(self,fnDirPrevious,fnDir,TsCurrent,getShiftsFrom=''):
        print "Preparing images to sampling rate=",TsCurrent
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        if newXdim<40:
            newXdim=long(40)
            TsCurrent=Xdim*(self.TsOrig/newXdim)
        self.writeInfoField(fnDir,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
        self.writeInfoField(fnDir,"size",xmipp.MDL_XSIZE,newXdim)
        
        # Prepare particles
        fnDir0=self._getExtraPath("Iter000")
        fnNewParticles=join(fnDir,"images.stk")
        if newXdim!=Xdim:
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(self.imgsFn,fnNewParticles,newXdim),numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        else:
            self.runJob("xmipp_image_convert","-i %s -o %s --save_metadata_stack %s"%(self.imgsFn,fnNewParticles,join(fnDir,"images.xmd")),
                        numberOfMpi=1)
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        R=min(round(R*self.TsOrig/TsCurrent*(1+self.angularMaxShift.get()*0.01)),newXdim/2)
        self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnNewParticles,R),numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        fnSource=join(fnDir,"images.xmd")
        if self.splitMethod==self.SPLIT_STOCHASTIC:
            self.runJob('xmipp_metadata_utilities','-i %s --set intersection %s particleId particleId -o %s/all_images.xmd'%\
                        (fnSource,self._getExtraPath('images.xmd'),fnDir),numberOfMpi=1)
            self.runJob("xmipp_metadata_split","-i %s/all_images.xmd --oroot %s/images -n 2"%(fnDir,fnDir),numberOfMpi=1)
            cleanPath("%s/all_images.xmd"%fnDir)
            for i in range(1,3):
                moveFile("%s/images%06d.xmd"%(fnDir,i),"%s/images%02d.xmd"%(fnDir,i))
        else:
            for i in range(1,3):
                fnImagesi=join(fnDir,"images%02d.xmd"%i)
                self.runJob('xmipp_metadata_utilities','-i %s --set intersection %s/images%02d.xmd particleId particleId -o %s'%\
                            (fnSource,fnDir0,i,fnImagesi),numberOfMpi=1)
        cleanPath(fnSource)
        
        if getShiftsFrom!="":
            fnPreviousAngles=join(getShiftsFrom,"angles.xmd")
            TsPrevious=self.readInfoField(getShiftsFrom,"sampling",xmipp.MDL_SAMPLINGRATE)
            fnAux=join(fnDir,"aux.xmd")
            for i in range(1,3):
                fnImagesi=join(fnDir,"images%02d.xmd"%i)
                self.runJob('xmipp_metadata_utilities','-i %s --set join %s particleId particleId -o %s'%\
                            (fnImagesi,fnPreviousAngles,fnAux),numberOfMpi=1)
                self.adaptShifts(fnAux, TsPrevious, fnImagesi, TsCurrent)
            cleanPath(fnAux)
        
    def prepareReferences(self,fnDirPrevious,fnDir,TsCurrent,targetResolution):
        print "Preparing references to sampling rate=",TsCurrent
        fnMask=''
        newXdim=self.readInfoField(fnDir,"size",xmipp.MDL_XSIZE)
        if self.nextMask.hasValue():
            fnMask=join(fnDir,"mask.vol")
            self.prepareMask(self.nextMask.get(), fnMask, TsCurrent, newXdim)
        TsPrevious=self.readInfoField(fnDirPrevious,"sampling",xmipp.MDL_SAMPLINGRATE)
        for i in range(1,3):
            fnPreviousVol=join(fnDirPrevious,"volume%02d.vol"%i)
            fnReferenceVol=join(fnDir,"volumeRef%02d.vol"%i)
            if TsPrevious!=TsCurrent:
                self.runJob("xmipp_image_resize","-i %s -o %s --dim %d"%(fnPreviousVol,fnReferenceVol,newXdim),numberOfMpi=1)
            else:
                copyFile(fnPreviousVol, fnReferenceVol)
            self.runJob('xmipp_transform_filter','-i %s --fourier fsc %s --sampling %f'%(fnReferenceVol,join(fnDirPrevious,"fsc.xmd"),TsCurrent),numberOfMpi=1)
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
        return atan2(1,k)*180.0/pi # Corresponding angular step

    def globalAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        makePath(fnDirCurrent)
        previousResolution=self.readInfoField(fnDirPrevious,"resolution",xmipp.MDL_RESOLUTION_FREQREAL)

        if self.alignmentMethod==self.GLOBAL_ALIGNMENT:
            fnGlobal=join(fnDirCurrent,"globalAssignment")
            makePath(fnGlobal)
    
            targetResolution=previousResolution*0.8
            if self.multiresolution:
                TsCurrent=max(self.TsOrig,targetResolution/3)
            else:
                TsCurrent=self.TsOrig
            getShiftsFrom=''
            if iteration>1:
                getShiftsFrom=fnDirPrevious
            self.prepareImages(fnDirPrevious,fnGlobal,TsCurrent,getShiftsFrom)
            self.prepareReferences(fnDirPrevious,fnGlobal,TsCurrent,targetResolution)

            # Calculate angular step at this resolution
            ResolutionAlignment=previousResolution
            if self.nextLowPass:
                ResolutionAlignment+=self.nextResolutionOffset.get()
            newXdim=self.readInfoField(fnGlobal,"size",xmipp.MDL_XSIZE)
            angleStep=self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)
            angleStep=max(angleStep,3.0)
            self.writeInfoField(fnGlobal,"angleStep",xmipp.MDL_ANGLE_DIFF,float(angleStep))
            
            # Global alignment
            for i in range(1,3):
                fnDirSignificant=join(fnGlobal,"significant%02d"%i)
                fnImgs=join(fnGlobal,"images%02d.xmd"%i)
                makePath(fnDirSignificant)

                # Create defocus groups
                row=getFirstRow(fnImgs)
                if row.containsLabel(xmipp.MDL_CTF_MODEL) or row.containsLabel(xmipp.MDL_CTF_DEFOCUSU):
                    self.runJob("xmipp_ctf_group","--ctfdat %s -o %s/ctf:stk --pad 2.0 --sampling_rate %f --phase_flipped  --error 0.1 --resol %f"%\
                                (fnImgs,fnDirSignificant,TsCurrent,targetResolution),numberOfMpi=1)
                    moveFile("%s/ctf_images.sel"%fnDirSignificant,"%s/ctf_groups.xmd"%fnDirSignificant)
                    cleanPath("%s/ctf_split.doc"%fnDirSignificant)
                    mdInfo = xmipp.MetaData("numberGroups@%s"%join(fnDirSignificant,"ctfInfo.xmd"))
                    fnCTFs="%s/ctf_ctf.stk"%fnDirSignificant
                    numberGroups=mdInfo.getValue(xmipp.MDL_COUNT,mdInfo.firstObject())
                    ctfPresent=True
                else:
                    numberGroups=1
                    ctfPresent=False
                    fnCTFs=""

                # Generate projections
                fnReferenceVol=join(fnGlobal,"volumeRef%02d.vol"%i)
                fnGallery=join(fnDirSignificant,"gallery%02d.stk"%i)
                fnGalleryMd=join(fnDirSignificant,"gallery%02d.xmd"%i)
                args="-i %s -o %s --sampling_rate %f --sym %s --min_tilt_angle %f --max_tilt_angle %f"%\
                     (fnReferenceVol,fnGallery,angleStep,self.symmetryGroup,self.angularMinTilt.get(),self.angularMaxTilt.get())
                args+=" --compute_neighbors --angular_distance -1 --experimental_images %s"%self._getExtraPath("images.xmd")
                self.runJob("xmipp_angular_project_library",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
                cleanPath(join(fnDirSignificant,"gallery_angles%02d.doc"%i))
                moveFile(join(fnDirSignificant,"gallery%02d.doc"%i), fnGalleryMd)
                fnAngles=join(fnGlobal,"anglesDisc%02d.xmd"%i)
                for j in range(1,numberGroups+1):
                    fnAnglesGroup=join(fnDirSignificant,"angles_group%03d.xmd"%j)
                    if not exists(fnAnglesGroup):
                        if ctfPresent:
                            fnGroup="ctfGroup%06d@%s/ctf_groups.xmd"%(j,fnDirSignificant)                            
                            fnGalleryGroup=fnGallery
                            fnGalleryGroupMd=fnGalleryMd
                        else:
                            fnGroup=fnImgs
                            fnGalleryGroupMd=fnGalleryMd
                        if getSize(fnGroup)==0: # If the group is empty
                            continue
                        maxShift=round(self.angularMaxShift.get()*newXdim/100)
                        R=self.particleRadius.get()
                        if R<=0:
                            R=self.inputParticles.get().getDimensions()[0]/2
                        R=R*self.TsOrig/TsCurrent
                        args='-i %s -o %s --ref %s --Ri 0 --Ro %d --max_shift %d --search5d_shift %d --search5d_step %f --mem 2 --append --pad 2.0'%\
                             (fnGroup,join(fnDirSignificant,"angles_group%03d.xmd"%j),fnGalleryGroup,R,maxShift,self.shiftSearch5d.get(),self.shiftStep5d.get())
                        if ctfPresent:
                            args+=" --ctf %d@%s"%(j,fnCTFs)
                        if self.numberOfMpi>1:
                            args+=" --mpi_job_size 2"
                        self.runJob('xmipp_angular_projection_matching',args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
                        if j==1:
                            copyFile(fnAnglesGroup, fnAngles)
                        else:
                            self.runJob("xmipp_metadata_utilities","-i %s --set union %s"%(fnAngles,fnAnglesGroup),numberOfMpi=1)
                self.runJob("xmipp_metadata_utilities","-i %s --set join %s image"%(fnAngles,fnImgs),numberOfMpi=1)
                if self.saveSpace and ctfPresent:
                    self.runJob("rm -f",fnDirSignificant+"/gallery*",numberOfMpi=1)
                
    def adaptShifts(self, fnSource, TsSource, fnDest, TsDest):
        K=TsSource/TsDest
        copyFile(fnSource,fnDest)
        row=getFirstRow(fnDest)
        if row.containsLabel(xmipp.MDL_SHIFT_X):
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "shiftX=%f*shiftX"'%(fnDest,K),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "shiftY=%f*shiftY"'%(fnDest,K),numberOfMpi=1)
        if row.containsLabel(xmipp.MDL_CONTINUOUS_X):
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "continuousX=%f*continuousX"'%(fnDest,K),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "continuousY=%f*continuousY"'%(fnDest,K),numberOfMpi=1)

    def localAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        if self.alignmentMethod==self.LOCAL_ALIGNMENT:
            fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
            fnDirLocal=join(fnDirCurrent,"localAssignment")
            makePath(fnDirLocal)

            previousResolution=self.readInfoField(fnDirPrevious,"resolution",xmipp.MDL_RESOLUTION_FREQREAL)
            targetResolution=previousResolution*0.8
            if self.multiresolution:
                TsCurrent=max(self.TsOrig,targetResolution/3)
            else:
                TsCurrent=self.TsOrig
            self.writeInfoField(fnDirLocal,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
            TsCurrent=self.readInfoField(fnDirLocal,"sampling",xmipp.MDL_SAMPLINGRATE) # Write and read to guarantee consistency with previous directories 
            
            # Prepare images and references
            produceNewReferences=True
            fnDirGlobal=join(fnDirCurrent,"globalAssignment")
            if exists(fnDirGlobal):
                TsGlobal=self.readInfoField(fnDirGlobal,"sampling",xmipp.MDL_SAMPLINGRATE)
                if TsGlobal==TsCurrent:
                    produceNewReferences=False
            if produceNewReferences:
                self.prepareImages(fnDirPrevious,fnDirLocal,TsCurrent,fnDirPrevious)
                self.prepareReferences(fnDirPrevious,fnDirLocal,TsCurrent,targetResolution)
            else:
                newXdim=self.readInfoField(fnDirGlobal,"size",xmipp.MDL_XSIZE)
                self.writeInfoField(fnDirLocal,"size",xmipp.MDL_XSIZE,newXdim)
                for i in range(1,3):
                    createLink(join(fnDirGlobal,"images%02d.xmd"%i),join(fnDirLocal,"images%02d.xmd"%i))
                    createLink(join(fnDirGlobal,"volumeRef%02d.vol"%i),join(fnDirLocal,"volumeRef%02d.vol"%i))

            # Compute maximum angular deviation
            ResolutionAlignment=previousResolution
            if self.nextLowPass:
                ResolutionAlignment+=self.nextResolutionOffset.get()
            newXdim=self.readInfoField(fnDirLocal,"size",xmipp.MDL_XSIZE)
            maxAngle=3*self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)

            for i in range(1,3):
                fnLocalXmd=join(fnDirLocal,"anglesCont%02d.xmd"%i)
                if not exists(fnLocalXmd):
                    fnLocalImages=join(fnDirLocal,"images%02d.xmd"%i)
    
                    # Starting angles
                    fnLocalAssignment=join(fnDirLocal,"anglesDisc%02d.xmd"%i)
                    if exists(fnDirGlobal):
                        fnGlobalAssignment=join(fnDirGlobal,"anglesDisc%02d.xmd"%i)
                        TsGlobal=self.readInfoField(fnDirGlobal,"sampling",xmipp.MDL_SAMPLINGRATE)
                        if TsGlobal==TsCurrent:
                            copyFile(fnGlobalAssignment,fnLocalAssignment)
                        else:
                            self.adaptShifts(fnGlobalAssignment,TsGlobal,fnLocalAssignment,TsCurrent)
                    else:
                        TsPrevious=self.readInfoField(fnDirPrevious,"sampling",xmipp.MDL_SAMPLINGRATE)
                        fnAux=join(fnDirLocal,"aux.xmd")
                        self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s particleId particleId -o %s"%\
                                    (join(fnDirPrevious,"angles.xmd"),fnLocalImages,fnAux),numberOfMpi=1)
                        self.adaptShifts(fnAux,TsPrevious,fnLocalAssignment,TsCurrent)
                        cleanPath(fnAux)
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
                        args+=" --optimizeShift --max_shift %f"%(self.contMaxShiftVariation.get()*newXdim*0.01)
                    if self.contScale:
                        args+=" --optimizeScale --max_scale %f"%self.contMaxScale.get() 
                    if self.contAngles:
                        args+=" --optimizeAngles --max_angular_change %f"%maxAngle
                    if self.contGrayValues:
                        args+=" --optimizeGray --max_gray_scale %f --max_gray_shift %f"%(self.contMaxGrayScale.get(),self.contMaxGrayShift.get())
                    if self.contDefocus:
                        args+=" --optimizeDefocus --max_defocus_change %f"%self.contMaxDefocus.get()
                    if self.inputParticles.get().isPhaseFlipped():
                        args+=" --phaseFlipped"
                    #if self.weightResiduals:
                    #    args+=" --oresiduals %s"%join(fnDirLocal,"residuals%02i.stk"%i)
                    self.runJob("xmipp_angular_continuous_assign2",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
                    self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnLocalStk,R),numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())

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
                TsCurrent=self.readInfoField(fnDirLocal,"sampling",xmipp.MDL_SAMPLINGRATE)
                Xdim=self.readInfoField(fnDirLocal,"size",xmipp.MDL_XSIZE)
            else:
                if exists(fnAnglesDisc):
                    copyFile(fnAnglesDisc, fnAngles)
                    TsCurrent=self.readInfoField(fnDirGlobal,"sampling",xmipp.MDL_SAMPLINGRATE)
                    Xdim=self.readInfoField(fnDirGlobal,"size",xmipp.MDL_XSIZE)
                else:
                    raise Exception("Angles for iteration "+str(iteration)+" not found")
            self.writeInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
            self.writeInfoField(fnDirCurrent,"size",xmipp.MDL_XSIZE,Xdim)
                
            if self.weightSSNR:
                row=getFirstRow(fnAngles)
                if row.containsLabel(xmipp.MDL_WEIGHT_SSNR):
                    self.runJob("xmipp_metadata_utilities","-i %s --operate drop_column weightSSNR"%fnAngles,numberOfMpi=1)
                self.runJob("xmipp_metadata_utilities","-i %s --set join %s particleId"%\
                            (fnAngles,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
            if self.weightJumper and iteration>1:
                fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
                if self.splitMethod == self.SPLIT_FIXED:
                    fnPreviousAngles=join(fnDirPrevious,"angles%02d.xmd"%i)
                else:
                    fnPreviousAngles=join(fnDirCurrent,"aux.xmd")
                    self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s particleId particleId -o %s"%\
                                (join(fnDirPrevious,"angles.xmd"),fnAngles,fnPreviousAngles),numberOfMpi=1)
                self.runJob("xmipp_angular_distance","--ang1 %s --ang2 %s --compute_weights --oroot %s --sym %s --check_mirrors"%\
                            (fnPreviousAngles,fnAngles,fnDirCurrent+"/jumper",self.symmetryGroup),numberOfMpi=1)
                moveFile(fnDirCurrent+"/jumper_weights.xmd", fnAngles)
                if self.splitMethod == self.SPLIT_STOCHASTIC:
                    cleanPath(fnPreviousAngles)
                if iteration>2:
                    fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-2))
                    if self.splitMethod == self.SPLIT_FIXED:
                        fnPreviousAngles=join(fnDirPrevious,"angles%02d.xmd"%i)
                    else:
                        fnPreviousAngles=join(fnDirCurrent,"aux.xmd")
                        self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s particleId particleId -o %s"%\
                                    (join(fnDirPrevious,"angles.xmd"),fnAngles,fnPreviousAngles),numberOfMpi=1)
                    self.runJob("xmipp_angular_distance","--ang1 %s --ang2 %s --compute_weights --oroot %s --set 2 --sym %s --check_mirrors"%\
                                (fnPreviousAngles,fnAngles,fnDirCurrent+"/jumper",self.symmetryGroup),numberOfMpi=1)
                    moveFile(fnDirCurrent+"/jumper_weights.xmd", fnAngles)
                    if self.splitMethod == self.SPLIT_STOCHASTIC:
                        cleanPath(fnPreviousAngles)

            #if self.weightResiduals and exists(fnAnglesCont):
            #    fnCovariance=join(fnDirLocal,"covariance%02d.stk"%i)
            #    self.runJob("xmipp_image_residuals","-i %s -o %s --normalizeDivergence"%(fnAngles,fnCovariance),numberOfMpi=1)
            #    moveFile(join(fnDirLocal,"covariance%02d.xmd"%i),fnAngles)
            
            mdAngles=xmipp.MetaData(fnAngles)
            for objId in mdAngles:
                weight=1.0
                if self.weightSSNR:
                    aux=mdAngles.getValue(xmipp.MDL_WEIGHT_SSNR,objId)
                    weight*=aux
                if self.weightContinuous and exists(fnAnglesCont) and self.alignmentMethod==self.LOCAL_ALIGNMENT:
                    aux=mdAngles.getValue(xmipp.MDL_WEIGHT_CONTINUOUS2,objId)
                    weight*=aux
                #if self.weightResiduals and exists(fnAnglesCont):
                #    aux=mdAngles.getValue(xmipp.MDL_ZSCORE_RESCOV,objId)
                #    aux/=3
                #    weight*=exp(-0.5*aux*aux)
                #    aux=mdAngles.getValue(xmipp.MDL_ZSCORE_RESMEAN,objId)
                #    aux/=3
                #    weight*=exp(-0.5*aux*aux)
                #    aux=mdAngles.getValue(xmipp.MDL_ZSCORE_RESVAR,objId)
                #    aux/=3
                #    weight*=exp(-0.5*aux*aux)
                if self.weightJumper and iteration>1:
                    w1=mdAngles.getValue(xmipp.MDL_WEIGHT_JUMPER,objId)
                    w2=1.0
                    if iteration>2:
                        w2=mdAngles.getValue(xmipp.MDL_WEIGHT_JUMPER2,objId)
                    weight*=w1*w2
                
                mdAngles.setValue(xmipp.MDL_WEIGHT,weight,objId)
            mdAngles.write(fnAngles)
            
        fnAngles=join(fnDirCurrent,"angles.xmd")
        fnAngles1=join(fnDirCurrent,"angles01.xmd")
        fnAngles2=join(fnDirCurrent,"angles02.xmd")
        self.runJob('xmipp_metadata_utilities',"-i %s --set union %s -o %s"%(fnAngles1,fnAngles2,fnAngles),numberOfMpi=1)

    def qualifyParticles(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnAngles=join(fnDirCurrent,"angles.xmd")
        fnAnglesQualified=join(fnDirCurrent,"angles_qualified.xmd")
        
        # Qualify according to CC and COST by defocus groups
        row=getFirstRow(fnAngles)
        if row.containsLabel(xmipp.MDL_CTF_MODEL) or row.containsLabel(xmipp.MDL_CTF_DEFOCUSU):
            previousResolution=self.readInfoField(fnDirPrevious,"resolution",xmipp.MDL_RESOLUTION_FREQREAL)
            TsCurrent=self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
            numberGroups=50
            self.runJob("xmipp_ctf_group","--ctfdat %s -o %s/ctf:stk --simple %d"%\
                        (fnAngles,fnDirCurrent,numberGroups),numberOfMpi=1)
            moveFile("%s/ctf_images.sel"%fnDirCurrent,"%s/ctf_groups.xmd"%fnDirCurrent)
            ctfPresent=True
        else:
            numberGroups=1
            ctfPresent=False

        for j in range(1,numberGroups+2):
            fnAnglesGroup=join(fnDirCurrent,"angles_group%03d.xmd"%j)
            if ctfPresent:
                fnGroup="ctfGroup%06d@%s/ctf_groups.xmd"%(j,fnDirCurrent)
            else:
                fnGroup=fnAngles
            if getSize(fnGroup)>0:
                if row.containsLabel(xmipp.MDL_MAXCC):
                    self.runJob("xmipp_metadata_utilities","-i %s --operate percentile maxCC maxCCPerc -o %s"%(fnGroup,fnAnglesGroup),numberOfMpi=1)
                    fnGroup=fnAnglesGroup    
                if row.containsLabel(xmipp.MDL_COST):
                    self.runJob("xmipp_metadata_utilities","-i %s --operate percentile cost costPerc -o %s"%(fnGroup,fnAnglesGroup),numberOfMpi=1)          
                if not exists(fnAnglesQualified):
                    copyFile(fnAnglesGroup, fnAnglesQualified)
                else:
                    self.runJob("xmipp_metadata_utilities","-i %s --set union %s"%(fnAnglesQualified,fnAnglesGroup),numberOfMpi=1)
                cleanPath(fnAnglesGroup)
        if ctfPresent:
            cleanPath("%s/ctf_groups.xmd"%fnDirCurrent)
        moveFile(fnAnglesQualified, fnAngles)

        if self.weightCC:
            mdAngles=xmipp.MetaData(fnAngles)
            weightCCmin=float(self.weightCCmin.get())
            for objId in mdAngles:
                if self.alignmentMethod==self.LOCAL_ALIGNMENT:
                    w=mdAngles.getValue(xmipp.MDL_COST_PERCENTILE,objId)
                else:
                    w=mdAngles.getValue(xmipp.MDL_MAXCC_PERCENTILE,objId)
                weight=mdAngles.getValue(xmipp.MDL_WEIGHT,objId)
                weight*=weightCCmin+w*(1-weightCCmin)
                mdAngles.setValue(xmipp.MDL_WEIGHT,weight,objId)
            mdAngles.write(fnAngles)

    def reconstruct(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
        for i in range(1,3):
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
            # Reconstruct Fourier
            args="-i %s -o %s --sym %s --weight --thr %d"%(fnAngles,fnVol,self.symmetryGroup,self.numberOfThreads.get())
            row=getFirstRow(fnAngles)
            if row.containsLabel(xmipp.MDL_CTF_DEFOCUSU) or row.containsLabel(xmipp.MDL_CTF_MODEL):
                args+=" --useCTF --sampling %f --minCTF %f"%(TsCurrent,self.minCTF.get())
                if self.inputParticles.get().isPhaseFlipped():
                    args+=" --phaseFlipped"
            self.runJob("xmipp_reconstruct_fourier",args,numberOfMpi=self.numberOfMpi.get()+1)
            # Reconstruct ADMM
#            args="-i %s -o %s --sym %s"%(fnAngles,fnVol,self.symmetryGroup)
#            row=getFirstRow(fnAngles)
#            if row.containsLabel(xmipp.MDL_CTF_DEFOCUSU) or row.containsLabel(xmipp.MDL_CTF_MODEL):
#                if not self.inputParticles.get().isPhaseFlipped():
#                    args+=" --dontUseCTF"
#            self.runJob("xmipp_reconstruct_admm",args)
    
    def postProcessing(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
        if not self.postSymmetryWithinMask \
           and not self.postSymmetryHelical and self.postScript=="" and not self.postAdHocMask.hasValue():
            return
        for i in range(1,3):
            fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
            fnBeforeVol=join(fnDirCurrent,"volumeBeforePostProcessing%02d.vol"%i)
            copyFile(fnVol,fnBeforeVol)
            volXdim = self.readInfoField(fnDirCurrent, "size", xmipp.MDL_XSIZE)
            
            if self.postAdHocMask.hasValue():
                fnMask=join(fnDirCurrent,"mask.vol")
                self.prepareMask(self.postAdHocMask.get(), fnMask, TsCurrent, volXdim)
                self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnVol,fnMask),numberOfMpi=1)
                cleanPath(fnMask)

            if self.postSymmetryWithinMask:
                if self.postMaskSymmetry!="c1":
                    fnMask=join(fnDirCurrent,"mask.vol")
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
                    cleanPath(join(fnGlobal,"volumeRef%02d.vol"%i))
                if exists(fnLocal):
                    cleanPath(join(fnLocal,"images%02d.xmd"%i))
                    cleanPath(join(fnLocal,"anglesCont%02d.stk"%i))
                    cleanPath(join(fnLocal,"anglesDisc%02d.xmd"%i))
                    cleanPath(join(fnLocal,"volumeRef%02d.vol"%i))
                    #if self.weightResiduals:
                    #    cleanPath(join(fnLocal,"covariance%02d.stk"%i))
                    #    cleanPath(join(fnLocal,"residuals%02i.stk"%i))

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if isinstance(self.inputVolumes.get(),SetOfVolumes) and self.inputVolumes.get().getSize()!=2:
            errors.append("The set of input volumes should have exactly 2 volumes")
        if self.postSymmetryWithinMask and not self.postSymmetryWithinMaskMask.hasValue():
            errors.append("Symmetrize within mask requires a mask")
        if not self.doContinue and not self.inputParticles.hasValue():
            errors.append("You must provide input particles")
        return errors    
    
    def _summary(self):
        summary = []
        summary.append("Symmetry: %s" % self.symmetryGroup.get())
        summary.append("Number of iterations: "+str(self.numberOfIterations))
        if self.alignmentMethod==self.GLOBAL_ALIGNMENT:
            summary.append("Global alignment, shift search: %f in steps of %f"%(self.shiftSearch5d.get(), self.shiftStep5d.get()))
        else:
            auxStr="Local alignment, refining: "
            if self.contShift:
                auxStr+="shifts "
            if self.contScale:
                auxStr+="scale "
            if self.contAngles:
                auxStr+="angles "
            if self.contGrayValues:
                auxStr+="gray "
            if self.contDefocus:
                auxStr+="defocus"
            summary.append(auxStr)
        auxStr="Weights: "
        if self.weightSSNR:
            auxStr+="SSNR "
        if self.weightContinuous and self.alignmentMethod==self.LOCAL_ALIGNMENT:
            auxStr+="Continuous "
        if self.weightJumper:
            auxStr+="Jumper"
        summary.append(auxStr)
        if self.postSymmetryWithinMask:
            summary.append("Symmetrizing within mask: "+self.postMaskSymmetry)
        if self.postSymmetryHelical:
            summary.append("Looking for helical symmetry")
        return summary
    
    def _methods(self):
        strline = ''
        if hasattr(self, 'outputVolume') or True:
            strline += 'We processed %d particles from %s ' % (self.inputParticles.get().getSize(), 
                                                                self.getObjectTag('inputParticles'))
            strline += 'using %s as reference and Xmipp highres procedure. ' % (self.getObjectTag('inputVolumes'))
            if self.symmetryGroup!="c1":
                strline+="We imposed %s symmetry. "%self.symmetryGroup
            strline += "We performed %d iterations of "%self.numberOfIterations.get()
            if self.alignmentMethod==self.GLOBAL_ALIGNMENT:
                strline+=" global alignment (shift search: %f in steps of %f pixels)"%(self.shiftSearch5d.get(), self.shiftStep5d.get())
            else:
                strline+=" local alignment, refining "
                if self.contShift:
                    strline+="shifts "
                if self.contScale:
                    strline+="scale "
                if self.contAngles:
                    strline+="angles "
                if self.contGrayValues:
                    strline+="gray "
                if self.contDefocus:
                    strline+="defocus"
            strline+=". "
            if self.weightSSNR or (self.weightContinuous and self.alignmentMethod==self.LOCAL_ALIGNMENT) or self.weightJumper:
                strline+="For reconstruction, we weighted the images according to "
                if self.weightSSNR:
                    strline+="their SSNR "
                if self.weightContinuous and self.alignmentMethod==self.LOCAL_ALIGNMENT:
                    strline+=", their correlation in the continuous alignment "
                if self.weightJumper:
                    strline+=", and their angular stability"
                strline+=". "
            if self.postAdHocMask.hasValue():
                strline+="We masked the reconstruction with %s. "%self.getObjectTag('postAdHocMask')
                if self.postSymmetryWithinMask:
                    strline+="We imposed %s symmetry within the mask %s. "%(self.postSymmetryWithinMaskType.get(),self.getObjectTag('postSymmetryWithinMaskMask'))
            if self.postSymmetryHelical:
                strline+="Finally, we imposed helical symmetry. "
        return [strline]
