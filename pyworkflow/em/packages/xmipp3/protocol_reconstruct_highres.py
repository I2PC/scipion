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
                  MDL_COUNT, MDL_SHIFT_X, MDL_CONTINUOUS_X, Image

#Continuar un procesamiento anterior
#Criterio de convergencia
#Phase flipping interno para poder corregir por el desenfoque
#Hacer local a partir de una resolucion
#Global hasta una resolucion
#FSC as a filter to avoid overfitting
#Jumper weight
        
class XmippProtReconstructHighRes(ProtRefine3D):
    """Reconstruct a volume at high resolution"""
    _label = 'highres'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, label="Images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select a set of images at full resolution')
        form.addParam('inputVolumes', PointerParam, label="Initial volumes", important=True,
                      pointerClass='Volume, SetOfVolumes',
                      help='Select a set of volumes with 2 volumes or a single volume')

        form.addParam('particleRadius', IntParam, default=-1, 
                     label='Radius of particle (px)',
                     help='This is the radius (in pixels) of the spherical mask covering the particle in the input images')       
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('numberOfIterations', IntParam, default=6, label='Max. number of iterations')
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
        # COSS: Falta un script de nextReference
        form.addParam('nextRemove', BooleanParam, label="Remove reference to save space?", default=True, expertLevel=LEVEL_ADVANCED, 
                      help='Remove reference volumes once they are not needed any more.')

        form.addSection(label='Angular assignment')
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=10,
                      help='Maximum shift as a percentage of the image size')
        line=form.addLine('Tilt angle:', help='0 degrees represent top views, 90 degrees represent side views')
        line.addParam('angularMinTilt', FloatParam, label="Min.", default=0)
        line.addParam('angularMaxTilt', FloatParam, label="Max.", default=90)
        groupSignificant = form.addGroup('Significant')
        groupSignificant.addParam('significantMaxResolution', FloatParam, label="Global assignment if resolution is worse than (A)", default=12,
                      help='Significant assignment is always performed on the first iteration. Starting from the second, you may '\
                      'decide whether to perform it or not. Note that the significant angular assignment is a robust angular assignment '\
                      'meant to avoid local minima, although it may take time to calculate.')
        groupSignificant.addParam('significantSignificance', FloatParam, label="Significance (%)", default=99.5)
        groupContinuous = form.addGroup('Continuous')
        groupContinuous.addParam('continuousMinResolution', FloatParam, label="Continuous assignment if resolution is better than (A)", default=15,
                      help='Continuous assignment can produce very accurate assignments if the initial assignment is correct.')
        groupContinuous.addParam('contShift', BooleanParam, label="Optimize shifts?", default=True,
                      help='Optimize shifts within a limit')
        groupContinuous.addParam('contScale', BooleanParam, label="Optimize scale?", default=True,
                      help='Optimize scale within a limit')
        groupContinuous.addParam('contMaxScale', FloatParam, label="Max. scale variation", default=0.02, expertLevel=LEVEL_ADVANCED)
        groupContinuous.addParam('contAngles', BooleanParam, label="Optimize angles?", default=True,
                      help='Optimize angles within a limit')
        groupContinuous.addParam('contGrayValues', BooleanParam, label="Optimize gray values?", default=False,
                      help='Optimize gray values. Do not perform this unless the reconstructed volume is gray-compatible with the projections,'\
                      ' i.e., the volumes haven been produced from projections')
        groupContinuous.addParam('contPadding', IntParam, label="Fourier padding factor", default=2, expertLevel=LEVEL_ADVANCED,
                      help='The volume is zero padded by this factor to produce projections')
        
        # COSS: Falta un script de postprocessing

        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.imgsFn=self._getExtraPath('images.xmd')
        self.TsOrig=self.inputParticles.get().getSamplingRate()
        self._insertFunctionStep('convertInputStep', self.inputParticles.getObjId())
        if self.weightSSNR:
            self._insertFunctionStep('doWeightSSNR')
        self._insertFunctionStep('doIteration000', self.inputVolumes.getObjId())
        self.iteration=1
        self.insertIteration(self.iteration)
    
    def insertIteration(self,iteration):
        self._insertFunctionStep('globalAssignment',iteration)
        self._insertFunctionStep('localAssignment',iteration)
#        self._insertFunctionStep('weightParticles',iteration)
#        self._insertFunctionStep('reconstruct',iteration)
        # COSS: Falta postprocessing
#        self._insertFunctionStep('evaluateReconstructions',iteration)
#        self._insertFunctionStep('cleanDirectory',iteration)
#        self._insertFunctionStep('decideNextIteration',iteration)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, inputParticlesId):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn)
        self.runJob('xmipp_metadata_utilities','-i %s --fill image1 constant noImage'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "image1=image"'%self.imgsFn,numberOfMpi=1)

    def doWeightSSNR(self):
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        self.runJob("xmipp_image_ssnr", "-i %s -R %d --sampling %f --normalizessnr"%\
                    (self.imgsFn,R,self.inputParticles.get().getSamplingRate()))
        self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate keep_column "itemId weightSSNR" '%\
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
            copyFile(fnVol1, fnVol2)
        
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
        md = MetaData(fnFsc)
        resolution=2*TsCurrent
        for objId in md:
            fsc = md.getValue(MDL_RESOLUTION_FRC,objId)
            if fsc<self.nextResolutionCriterion.get():
                resolution=md.getValue(MDL_RESOLUTION_FREQREAL,objId)
                break
        self.writeInfoField(fnDirCurrent,"resolution",MDL_RESOLUTION_FREQREAL,resolution)
        
        # Filter the average to that resolution
        self.runJob('xmipp_transform_filter','-i %s --fourier low_pass %f --sampling %f'%(fnVolAvg,resolution,TsCurrent))
        self.runJob('xmipp_image_header','-i %s --sampling_rate %f'%(fnVolAvg,TsCurrent),numberOfMpi=1)
        
        # COSS: Falta algo de estadistica: particulas aceptadas, numero de repeticiones, particulas rechazadas
    
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
        
        # Prepare particlesfnMask
        fnDir0=self._getExtraPath("Iter000")
        fnNewParticles=join(fnDir,"images.stk")
        self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(self.imgsFn,fnNewParticles,newXdim))
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        R=min(round(R*self.TsOrig/TsCurrent*(1+self.angularMaxShift.get()*0.01)),newXdim/2)
        self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnNewParticles,R))
        self.runJob('xmipp_transform_filter','-i %s --fourier fsc %s --sampling %f'%(fnNewParticles,join(fnDirPrevious,"fsc.xmd"),TsCurrent))
        fnSource=join(fnDir,"images.xmd")
        for i in range(1,3):
            fnImagesi=join(fnDir,"images%02d.xmd"%i)
            self.runJob('xmipp_metadata_utilities','-i %s --set intersection %s/images%02d.xmd itemId itemId -o %s'%\
                        (fnSource,fnDir0,i,fnImagesi),numberOfMpi=1)
            self.runJob('xmipp_metadata_utilities','-i %s --operate rename_column "itemId particleId"'%fnImagesi,numberOfMpi=1)
        cleanPath(fnSource)
        
    def prepareReferences(self,fnDirPrevious,fnDir,TsCurrent,targetResolution):
        print "Preparing references to sampling rate=",TsCurrent
        fnMask=''
        newXdim=self.readInfoField(fnDir,"size",MDL_XSIZE)
        if self.nextMask.hasValue():
            img = ImageHandler()
            fnMask=join(fnDir,"mask.vol")
            img.convert(self.nextMask.get(), fnMask)
            TsMask=self.nextMask.get().getSamplingRate()
            self.runJob('xmipp_image_resize',"-i %s --factor %f"%(fnMask,TsMask/TsCurrent),numberOfMpi=1)
            maskXdim, _, _, _ =img.getDimensions((1,fnMask))
            if newXdim!=maskXdim:
                self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnMask,newXdim),numberOfMpi=1)
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
        if fnMask!='' and self.saveSpace:
            cleanPath(fnMask)

    def calculateAngStep(self,newXdim,TsCurrent,ResolutionAlignment):
        k=newXdim*TsCurrent/ResolutionAlignment # Freq. index
        from math import atan2,pi
        return max(3,2*atan2(1,k)*180.0/pi) # Corresponding angular step

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

                fnAngles=join(fnGlobal,"angles%02d.xmd"%i)
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
            self.writeInfoField(fnDirLocal,"resolution",MDL_RESOLUTION_FREQREAL,TsCurrent)
            
            # Prepare images and references
            self.prepareImages(fnDirPrevious,fnDirLocal,TsCurrent)
            produceNewReferences=True
            fnDirGlobal=join(fnDirCurrent,"globalAssignment")
            if exists(fnDirGlobal):
                TsGlobal=self.readInfoField(fnDirGlobal,"sampling",MDL_SAMPLINGRATE)
                if TsGlobal==TsCurrent:
                    produceNewReferences=False
            if produceNewReferences:
                self.prepareReferences(fnDirPrevious,fnDirLocal,TsCurrent,targetResolution)
            else:
                createLink(join(fnDirGlobal,"volumeRef%02d.vol"%1),join(fnDirLocal,"volumeRef%02d.vol"%1))
                createLink(join(fnDirGlobal,"volumeRef%02d.vol"%2),join(fnDirLocal,"volumeRef%02d.vol"%2))

            # Compute maximum angular deviation
            ResolutionAlignment=previousResolution
            if self.nextLowPass:
                ResolutionAlignment+=self.nextResolutionOffset.get()
            newXdim=self.readInfoField(fnDirLocal,"size",MDL_XSIZE)
            maxAngle=3*self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)

            for i in range(1,3):
                fnLocalImages=join(fnDirLocal,"images%02d.xmd"%i)

                # Starting angles
                fnLocalAssignment=join(fnDirLocal,"angles%02d.xmd"%i)
                if exists(fnDirGlobal):
                    fnGlobalAssignment=join(fnDirGlobal,"angles%02d.xmd"%i)
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
                self.runJob("xmipp_angular_continuous_assign2",args)
                self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnLocalStk,R))
                # COSS: Falta continuous con CTF 

    def weightParticles(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        # COSS: Falta outliers
        # COSS: Falta clusterability
        for i in range(1,3):
            fnAngles=join(fnDirCurrent,"anglesCont%02d.xmd"%i)
            if self.weightSSNR:
                self.runJob("xmipp_metadata_utilities","-i %s --set join %s itemId"%\
                            (fnAngles,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
            if self.weightJumper and iteration>1:
                # COSS: Falta hacerlo bien
                self.runJob("xmipp_angular_distance","-i %s --compute_weights --oroot %s"%\
                            (fnAngles,fnDirCurrent),numberOfMpi=1)
            md=MetaData(fnAngles)
            for objId in md:
                weight=1
                if self.weightSSNR:
                    aux=md.getValue(MDL_WEIGHT_SSNR,objId)
                    weight*=aux
                if self.weightSignificance and (iteration==1 or iteration<=self.significantIterations.get()):
                    aux=md.getValue(MDL_WEIGHT_SIGNIFICANT,objId)
                    weight*=aux
                if self.weightContinuous:
                    aux=md.getValue(MDL_WEIGHT_CONTINUOUS2,objId)
                    weight*=aux
                md.setValue(MDL_WEIGHT,weight,objId)
            md.write(fnAngles)

    def reconstruct(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        for i in range(1,3):
            fnAngles=join(fnDirCurrent,"anglesCont%02d.xmd"%i)
            fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
            # COSS: Falta Fourier con CTF
            self.runJob("xmipp_reconstruct_fourier","-i %s -o %s --sym %s --weight"%(fnAngles,fnVol,self.symmetryGroup))
    
    def cleanDirectory(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        # COSS: Falta clean
    
    def decideNextIteration(self, iteration):
        # COSS: Falta un criterio para decidir si otra iteracion
        if iteration<self.numberOfIterations.get():
            self.iteration+=1
            self.insertIteration(self.iteration)
        
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if isinstance(self.inputVolumes.get(),SetOfVolumes) and self.inputVolumes.get().getSize()!=2:
            errors.append("The set of input volumes should have exactly 2 volumes")
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
    
