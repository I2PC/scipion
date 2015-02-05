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
This sub-package contains wrapper around Projection Outliers Xmipp program
"""

from pyworkflow.object import Float
from pyworkflow.protocol.constants import LEVEL_EXPERT
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam
from pyworkflow.utils.path import cleanPath, makePath, copyFile
from pyworkflow.em.protocol import ProtRefine3D
from convert import writeSetOfParticles
from os.path import join

from xmipp import MetaData, MDL_RESOLUTION_FRC, MDL_RESOLUTION_FREQREAL, MDL_SAMPLINGRATE, MD_APPEND
        
class XmippProtReconstructHighRes(ProtRefine3D):
    """Reconstruct a volume at high resolution"""
    _label = 'reconstruct highres'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, label="Images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select a set of images at full resolution')
        form.addParam('inputVolumes', PointerParam, label="Initial volumes", important=True,
                      pointerClass='Volume',
                      help='Select a set of volumes with 2 volumes')

        form.addParam('particleRadius', IntParam, default=-1, 
                     label='Radius of particle (px)',
                     help='This is the radius (in pixels) of the spherical mask covering the particle')       
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        
        form.addSection(label='Weights')
        form.addParam('weightSSNR', BooleanParam, label="Weight by SSNR?", default=True,
                      help='Weight input images by SSNR')
        form.addParam('weightSignificance', BooleanParam, label="Weight by Significance?", default=True,
                      help='Weight input images by angular assignment significance')

        form.addSection(label='Next Reference')
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
        # Falta una mascara mas ajustada
        form.addParam('nextRemove', BooleanParam, label="Remove reference to save space?", default=True,
                      help='Remove reference volumes once they are not needed any more.')

        form.addSection(label='Angular assignment')
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=5,
                      help='Maximum shift as a percentage of the image size')
        form.addParam('angularMinTilt', FloatParam, label="Min. Tilt", default=0,
                      help='Minimum tilt. 0 degrees represent top views, 90 degrees represent side views')
        form.addParam('angularMaxTilt', FloatParam, label="Max. Tilt", default=90,
                      help='Maximum tilt. 0 degrees represent top views, 90 degrees represent side views')
        groupSignificant = form.addGroup('Significant')
        groupSignificant.addParam('significantPerform', BooleanParam, label="Perform significant from 2nd iteration?", default=True,
                      help='Significant assignment is always performed on the first iteration. Starting from the second, you may '\
                      'decide whether to perform it or not. Note that the significant angular assignment is a robust angular assignment '\
                      'meant to avoid local minima, although it may take time to calculate.')
        groupSignificant.addParam('significantSignificance', FloatParam, label="Significance (%)", default=99.5)
        groupSignificant.addParam('significantStrict', BooleanParam, label="Strict directions", default=True,
                      help='If directions are strict, some images may not participate at all in the 3D reconstruction because they are '\
                      'not significant enough for any of the volume directions')

        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.imgsFn=self._getExtraPath('images.xmd')
        self.TsOrig=self.inputParticles.get().getSamplingRate()
        self._insertFunctionStep('convertInputStep', self.inputParticles)
        if self.weightSSNR:
            self._insertFunctionStep('doWeightSSNR')
        self._insertFunctionStep('doIteration000', self.inputVolumes)
        self.iteration=1
        self.insertIteration(self.iteration)
    
    def insertIteration(self,iteration):
        self._insertFunctionStep('determineSamplingRate',iteration)
        self._insertFunctionStep('calculateReference',iteration)
        self._insertFunctionStep('globalAssignment',iteration)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, inputParticles):
        writeSetOfParticles(inputParticles.get(),self.imgsFn)

    def doWeightSSNR(self):
        R=self.particleRadius.get()
        if R<=0:
            R=self.particleRadius.get().getDimensions()[0]/2
        self._insertRunJobStep("xmipp_image_ssnr", "-i %s -R %d --sampling %f --normalizessnr"%\
                               (self.imgsFn,R,self.inputParticles.get().getSamplingRate()))
        
    def doIteration000(self, inputVolumes):
        fnDir0=self._getExtraPath('Iter000')
        makePath(fnDir0)
        
        # Split data
        self.runJob("xmipp_metadata_split","-i %s --oroot %s/images -n 2"%(self.imgsFn,fnDir0),numberOfMpi=1)
        
        # Get volume sampling rate
        TsCurrent=inputVolumes.get().getSamplingRate()
        self.writeInfoField(fnDir0,"sampling",MDL_SAMPLINGRATE,TsCurrent)

        # Copy reference volumes and window if necessary
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=round(Xdim*self.TsOrig/TsCurrent)
        
        from pyworkflow.em.convert import ImageHandler
        img = ImageHandler()
        i=1
        for vol in inputVolumes.get():
            fnVol=join(fnDir0,"volume%06d.vol"%i)
            img.convert(vol, fnVol)
            if newXdim!=vol.getDim()[0]:
                self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnVol,newXdim),numberOfMpi=1)
            i+=1
        
        self.alignVolumes(fnDir0)
        self.estimateResolution(fnDir0,TsCurrent)
        
    def alignVolumes(self,fnDir):
        letter=self.symmetryGroup.get()[0]
        if letter=='i' or letter=='d':
            return
        fnVol1=join(fnDir,"volume%06d.vol"%1)
        fnVol2=join(fnDir,"volume%06d.vol"%2)
        fnVolAvg=join(fnDir,"volumeAvg.vol")
        self.runJob('xmipp_image_operate','-i %s --plus %s -o %s'%(fnVol1,fnVol2,fnVolAvg),numberOfMpi=1)
        self.runJob('xmipp_image_operate','-i %s --mult 0.5'%fnVolAvg,numberOfMpi=1)
        self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnVolAvg,fnVol1),numberOfMpi=1)
        self.runJob('xmipp_volume_align','--i1 %s --i2 %s --local --apply'%(fnVolAvg,fnVol2),numberOfMpi=1)
        cleanPath(fnVolAvg)

    def estimateResolution(self,fnDir,Ts):
        fnVol1=join(fnDir,"volume%06d.vol"%1)
        fnVol2=join(fnDir,"volume%06d.vol"%2)
        fnFsc=join(fnDir,"fsc.xmd")
        self.runJob('xmipp_resolution_fsc','--ref %s -i %s -o %s --sampling_rate %f'%(fnVol1,fnVol2,fnFsc,Ts),numberOfMpi=1)
        md = MetaData(fnFsc)
        resolution=2*Ts
        for objId in md:
            fsc = md.getValue(MDL_RESOLUTION_FRC,objId)
            if fsc<0.5:
                resolution=md.getValue(MDL_RESOLUTION_FREQREAL,objId)
                break
        self.writeInfoField(fnDir,"resolution",MDL_RESOLUTION_FREQREAL,resolution)
        return resolution
    
    def readInfoField(self,fnDir,block,label):
        md = MetaData("%s@%s"%(block,join(fnDir,"iterInfo.xmd")))
        return md.getValue(label,md.firstObject())

    def writeInfoField(self,fnDir,block,label, value):
        md = MetaData()
        objId=md.addObject()
        md.setValue(label,value,objId)
        md.write("%s@%s"%(block,join(fnDir,"iterInfo.xmd")),MD_APPEND)
    
    def determineSamplingRate(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        makePath(fnDirCurrent)
        
        TsPrevious=self.readInfoField(fnDirPrevious,"sampling",MDL_SAMPLINGRATE)
        ResolutionPrevious=self.readInfoField(fnDirPrevious,"resolution",MDL_RESOLUTION_FREQREAL)
        TsCurrent=max(self.TsOrig,ResolutionPrevious/3)
        self.writeInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE,TsCurrent)
        
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=round(Xdim*self.TsOrig/TsCurrent)
        
        # Prepare particles
        fnDir0=self._getExtraPath("Iter000")
        if TsCurrent!=self.TsOrig:
            fnNewParticles=join(fnDirCurrent,"images.stk")
            self.runJob("xmipp_image_resize","-i %s -o %s --dim %d"%\
                        (self.imgsFn,fnNewParticles,newXdim))
            self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%\
                        (fnNewParticles,ResolutionPrevious+self.nextResolutionOffset.get(),TsCurrent))
            fnSource=join(fnDirCurrent,"images.xmd")
            for i in range(1,3):
                self.runJob('xmipp_metadata_utilities','-i %s --set intersection %s/images%06d.xmd itemId itemId -o %s'%\
                            (fnSource,fnDir0,i,join(fnDirCurrent,"images%06d.xmd"%i)),numberOfMpi=1)
            cleanPath(fnSource)
        else:
            for i in range(1,3):
                copyFile(join(fnDir0,"images%06d.xmd"%i),join(fnDirCurrent,"images%06d.xmd"%i))
        
        # Prepare volumes
        if TsCurrent!=TsPrevious:
            self.runJob("xmipp_image_resize","-i %s -o %s --dim %d"%\
                        (join(fnDirPrevious,"volume%06d.vol"%1),join(fnDirCurrent,"volumeRef%06d.vol"%1),newXdim),numberOfMpi=1)
            self.runJob("xmipp_image_resize","-i %s -o %s --dim %d"%\
                        (join(fnDirPrevious,"volume%06d.vol"%1),join(fnDirCurrent,"volumeRef%06d.vol"%2),newXdim),numberOfMpi=1)
    
    def calculateReference(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE)
        ResolutionPrevious=self.readInfoField(fnDirPrevious,"resolution",MDL_RESOLUTION_FREQREAL)
        for i in range(1,3):
            fnReferenceVol=join(fnDirCurrent,"volumeRef%06d.vol"%i)
            if self.nextLowPass:
                self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%\
                            (fnReferenceVol,ResolutionPrevious+self.nextResolutionOffset.get(),TsCurrent),numberOfMpi=1)
            if self.nextSpherical:
                self.runJob('xmipp_transform_mask','-i %s --mask circular -%d'%\
                            (fnReferenceVol,round(self.particleRadius.get()*self.TsOrig/TsCurrent)),numberOfMpi=1)
            if self.nextPositivity:
                self.runJob('xmipp_transform_threshold','-i %s --select below 0 --substitute value 0'%fnReferenceVol,numberOfMpi=1)

    def globalAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE)
        if iteration==1 or self.significantPerform:
            ResolutionAlignment=self.readInfoField(fnDirPrevious,"resolution",MDL_RESOLUTION_FREQREAL)
            if self.nextLowPass:
                ResolutionAlignment+=self.nextResolutionOffset.get()
            
            # Calculate frequency index of this resolution 
            Xdim=self.inputParticles.get().getDimensions()[0]
            newXdim=round(Xdim*self.TsOrig/TsCurrent)
            k=newXdim*TsCurrent/ResolutionAlignment # Freq. index
            from math import atan2,pi
            angleStep=max(2.5,2*atan2(1,k)*180.0/pi) # Corresponding angular step
            
            alpha=1-0.01*self.significantSignificance.get()
            for i in range(1,3):
                fnReferenceVol=join(fnDirCurrent,"volumeRef%06d.vol"%i)
                fnParticles=join(fnDirCurrent,"images%06d.xmd"%i)
                fnDirSignificant=join(fnDirCurrent,"significant%d"%i)
                makePath(fnDirSignificant)
                args='-i %s --initvolumes %s --odir %s --sym %s --iter 1 --alpha0 %f --alphaF %f --angularSampling %f --maxShift %d '\
                     '--minTilt %f --maxTilt %f --useImed --angDistance %f'%(fnParticles,fnReferenceVol,fnDirSignificant,self.symmetryGroup,\
                                                                             alpha,alpha,angleStep,\
                                                                             round(self.angularMaxShift.get()*newXdim/100),\
                                                                             self.angularMinTilt.get(),self.angularMaxTilt.get(),2*angleStep)
                if self.significantStrict:
                    args+=" --strictDirection"
                self.runJob('xmipp_reconstruct_significant',args)
        else:
            TsPrevious=self.readInfoField(fnDirCurrent,"sampling",MDL_SAMPLINGRATE)
            K=TsPrevious/TsCurrent
            for i in range(1,3):
                fnPrevious=join(fnDirPrevious,"images%06d.xmd"%i)
                fnCurrent=join(fnDirCurrent,"images%06d.xmd"%i)
                self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate modify_values "shiftX=%f*shiftX"'%\
                            (fnPrevious,fnCurrent,K),numberOfMpi=1)
                self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate modify_values "shiftY=%f*shiftY"'%\
                            (fnPrevious,fnCurrent,K),numberOfMpi=1)
                

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors    
    
    def _summary(self):
        summary = []
        summary.append("Images: %i" % self.inputParticles.getNameId())
        summary.append("Volume: %s" % self.inputVolume.getNameId())
        summary.append("symmetry: %s" % self.symmetryGroup.get())
        return summary
    
    def _methods(self):
        methods = []
        return methods
    
