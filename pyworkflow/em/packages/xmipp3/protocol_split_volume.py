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
Protocol to split a volume in two volumes based on a set of images
"""

from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam, NumericListParam
from pyworkflow.utils.path import cleanPath, makePath, copyFile, moveFile
from pyworkflow.utils.process import runJob
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.em.data import Volume
from pyworkflow.em.metadata.utils import getFirstRow
from convert import writeSetOfParticles
from os.path import join, exists
from pyworkflow.em.convert import ImageHandler

import xmipp


class XmippProtSplitvolume(ProtClassify3D):
    """Split volume in two"""
    _label = 'split volume'
    
    def __init__(self, **args):
        ProtClassify3D.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, label="Full-size Images", important=True, 
                      pointerClass='SetOfParticles',
                      help='Select a set of images at full resolution')
        form.addParam('phaseFlipped', BooleanParam, label="Images have been phase flipped", default=True, 
                      help='Choose this option if images have been phase flipped')
        form.addParam('inputVolume', PointerParam, label="Initial volume", important=True,
                      pointerClass='Volume')
        form.addParam('particleRadius', IntParam, default=-1, 
                     label='Radius of particle (px)',
                     help='This is the radius (in pixels) of the spherical mask covering the particle in the input images')       
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('nextMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be binary: 0 (remove these voxels) and 1 (let them pass).')
        form.addParam('splitFraction', FloatParam, label="Split fraction", default=0.02, expertLevel=LEVEL_ADVANCED, 
                        help="The input data set will be split in groups of this size");
        form.addParam('splitNumber', IntParam, label="Number of splits", default=500, expertLevel=LEVEL_ADVANCED, 
                        help="Number of random splits to analyze");
        form.addParam('splitPercentiles', NumericListParam, label="Splitted percentiles", default="5 95", expertLevel=LEVEL_ADVANCED, 
                        help="Split percentiles");
        
        form.addSection(label='Weights')
        form.addParam('weightSSNR', BooleanParam, label="Weight by SSNR?", default=True,
                      help='Weight input images by SSNR')
        form.addParam('weightSignificance', BooleanParam, label="Weight by Significance?", default=True,
                      help='Weight input images by angular assignment significance')
        form.addParam('weightAstigmatism', BooleanParam, label="Weight by astigmatism?", default=True,
                      help='Give lower weight to those images whose astigmatic CTF would not allow to reach high resolution.' \
                           'This weight is calculated by the continuous assignment.')
        form.addParam('weightAstigmatismSigma', FloatParam, label="Astigmatism sigma", default=120, expertLevel=LEVEL_ADVANCED, 
                        condition="weightAstigmatism", help="Sigma in degrees for the CTF phase");

        form.addSection(label='Angular assignment')
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=10,
                      help='Maximum shift as a percentage of the image size')
        line=form.addLine('Tilt angle:', help='0 degrees represent top views, 90 degrees represent side views')
        line.addParam('angularMinTilt', FloatParam, label="Min.", default=0)
        line.addParam('angularMaxTilt', FloatParam, label="Max.", default=90)
        groupSignificant = form.addGroup('Global')
        groupSignificant.addParam('significantMaxResolution', FloatParam, label="Perform the split at this resolution (A)", default=12,
                      help='Images will be rescaled to adapt to this resolution')
        groupSignificant.addParam('significantSignificance', FloatParam, label="Significance (%)", default=99.75)

        form.addParallelSection(threads=4, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.imgsFn=self._getExtraPath('images.xmd')
        self._insertFunctionStep('convertInputStep', self.inputParticles.getObjId())
        if self.weightSSNR:
            self._insertFunctionStep('doWeightSSNR')
        self.TsOrig=self.inputParticles.get().getSamplingRate()
        self._insertFunctionStep('globalAssignment')
        self._insertFunctionStep('weightParticles')
        stepId=self._insertFunctionStep('reconstruct')
        listOfIds=[]
        for i in range(self.splitNumber.get()):
            listOfIds.append(self._insertFunctionStep('split',i,prerequisites=[stepId]))
        self._insertFunctionStep('generateSplittedVolumes',prerequisites=listOfIds)
        self._insertFunctionStep('cleanDirectory')
        self._insertFunctionStep('createOutput')

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, inputParticlesId):
        writeSetOfParticles(self.inputParticles.get(),self.imgsFn)
        self.runJob('xmipp_metadata_utilities','-i %s --fill image1 constant noImage'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate modify_values "image1=image"'%self.imgsFn,numberOfMpi=1)
        self.runJob('xmipp_metadata_utilities','-i %s --operate rename_column "itemId particleId"'%self.imgsFn,numberOfMpi=1)

    def createOutput(self):
        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate(self.inputParticles.get().getSamplingRate())
        Nvols = len(self.splitPercentiles.get().split())
        fnStack = self._getPath("splittedVolumes.stk")
        for i in range(Nvols):
            vol = Volume()
            vol.setLocation(i+1, fnStack)
            volumesSet.append(vol)
        
        self._defineOutputs(outputVolumes=volumesSet)
        self._defineSourceRelation(self.inputParticles.get(), volumesSet)
        
    def doWeightSSNR(self):
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        self.runJob("xmipp_image_ssnr", "-i %s -R %d --sampling %f --normalizessnr"%\
                    (self.imgsFn,R,self.inputParticles.get().getSamplingRate()))
        self.runJob('xmipp_metadata_utilities','-i %s -o %s --operate keep_column "particleId weightSSNR" '%\
                    (self.imgsFn,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
        
    def readInfoField(self,fnDir,block,label):
        md = xmipp.MetaData("%s@%s"%(block,join(fnDir,"iterInfo.xmd")))
        return md.getValue(label,md.firstObject())

    def writeInfoField(self,fnDir,block,label, value):
        md = xmipp.MetaData()
        objId=md.addObject()
        md.setValue(label,value,objId)
        md.write("%s@%s"%(block,join(fnDir,"iterInfo.xmd")),xmipp.MD_APPEND)
    
    def prepareImages(self,fnDir,TsCurrent):
        print "Preparing images to sampling rate=",TsCurrent
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        if newXdim<40:
            newXdim=long(40)
            TsCurrent=Xdim*(self.TsOrig/newXdim)
        self.writeInfoField(fnDir,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
        self.writeInfoField(fnDir,"size",xmipp.MDL_XSIZE,newXdim)
        
        # Prepare particles
        fnNewParticles=join(fnDir,"images.stk")
        if newXdim!=Xdim:
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(self.imgsFn,fnNewParticles,newXdim),
                        numberOfMpi=self.numberOfThreads.get())
        else:
            self.runJob("xmipp_image_convert","-i %s -o %s --save_metadata_stack %s"%(self.imgsFn,fnNewParticles,join(fnDir,"images.xmd")),
                        numberOfMpi=1)
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        R=min(round(R*self.TsOrig/TsCurrent*(1+self.angularMaxShift.get()*0.01)),newXdim/2)
        self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnNewParticles,R))
        
    def prepareReferences(self,fnDir,TsCurrent,targetResolution):
        print "Preparing references to sampling rate=",TsCurrent
        fnMask=''
        newXdim=self.readInfoField(fnDir,"size",xmipp.MDL_XSIZE)
        if self.nextMask.hasValue():
            fnMask=join(fnDir,"mask.vol")
            self.prepareMask(self.nextMask.get(), fnMask, TsCurrent, newXdim)

        fnReferenceVol=join(fnDir,"volumeRef.vol")
        img=ImageHandler()
        img.convert(self.inputVolume.get(), fnReferenceVol)
        Xdim=self.inputVolume.get().getDim()[0]
        if Xdim!=newXdim:
            self.runJob("xmipp_image_resize","-i %s --fourier %d"%(fnReferenceVol,newXdim),numberOfMpi=1)
        
        if fnMask!='':
            self.runJob('xmipp_image_operate','-i %s --mult %s'%(fnReferenceVol,fnMask),numberOfMpi=1)

    def prepareMask(self,maskObject,fnMask,TsMaskOut,XdimOut):
        img=ImageHandler()
        img.convert(maskObject, fnMask)
        self.runJob('xmipp_image_resize',"-i %s --factor %f"%(fnMask,maskObject.getSamplingRate()/TsMaskOut),numberOfMpi=1)
        maskXdim, _, _, _ =img.getDimensions((1,fnMask))
        if XdimOut!=maskXdim:
            self.runJob('xmipp_transform_window',"-i %s --size %d"%(fnMask,XdimOut),numberOfMpi=1)
        self.runJob('xmipp_transform_threshold',"-i %s --select below 0.5 --substitute binarize"%fnMask,numberOfMpi=1)

    def calculateAngStep(self,newXdim,TsCurrent,ResolutionAlignment):
        k=newXdim*TsCurrent/ResolutionAlignment # Freq. index
        from math import atan2,pi
        return 2*atan2(1,k)*180.0/pi # Corresponding angular step

    def globalAssignment(self):
        iteration=1
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        makePath(fnDirCurrent)

        fnGlobal=join(fnDirCurrent,"globalAssignment")
        makePath(fnGlobal)

        targetResolution=self.significantMaxResolution.get()
        TsCurrent=max(self.TsOrig,targetResolution/3)
        self.prepareImages(fnGlobal,TsCurrent)
        self.prepareReferences(fnGlobal,TsCurrent,targetResolution)

        # Calculate angular step at this resolution
        ResolutionAlignment=targetResolution
        newXdim=self.readInfoField(fnGlobal,"size",xmipp.MDL_XSIZE)
        angleStep=self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)
        self.writeInfoField(fnGlobal,"angleStep",xmipp.MDL_ANGLE_DIFF,float(angleStep))
        
        # Significant alignment
        alpha=1-0.01*self.significantSignificance.get()
        fnDirSignificant=join(fnGlobal,"significant")
        fnImgs=join(fnGlobal,"images.xmd")
        makePath(fnDirSignificant)

        # Create defocus groups
        row=getFirstRow(fnImgs)
        if row.containsLabel(xmipp.MDL_CTF_MODEL) or row.containsLabel(xmipp.MDL_CTF_DEFOCUSU):
            self.runJob("xmipp_ctf_group","--ctfdat %s -o %s/ctf:stk --pad 2.0 --sampling_rate %f --phase_flipped  --error 0.1 --resol %f"%\
                        (fnImgs,fnDirSignificant,TsCurrent,targetResolution),numberOfMpi=1)
            moveFile("%s/ctf_images.sel"%fnDirSignificant,"%s/ctf_groups.xmd"%fnDirSignificant)
            cleanPath("%s/ctf_split.doc"%fnDirSignificant)
            md = xmipp.MetaData("numberGroups@%s"%join(fnDirSignificant,"ctfInfo.xmd"))
            fnCTFs="%s/ctf_ctf.stk"%fnDirSignificant
            numberGroups=md.getValue(xmipp.MDL_COUNT,md.firstObject())
            ctfPresent=True
        else:
            numberGroups=1
            ctfPresent=False

        # Generate projections
        fnReferenceVol=join(fnGlobal,"volumeRef.vol")
        fnGallery=join(fnDirSignificant,"gallery.stk")
        fnGalleryMd=join(fnDirSignificant,"gallery.xmd")
        args="-i %s -o %s --sampling_rate %f --sym %s --min_tilt_angle %f --max_tilt_angle %f"%\
             (fnReferenceVol,fnGallery,angleStep,self.symmetryGroup,self.angularMinTilt.get(),self.angularMaxTilt.get())
        self.runJob("xmipp_angular_project_library",args)
        cleanPath(join(fnDirSignificant,"gallery_angles.doc"))
        moveFile(join(fnDirSignificant,"gallery.doc"), fnGalleryMd)

        fnAngles=join(fnGlobal,"anglesDisc.xmd")
        for j in range(1,numberGroups+1):
            fnAnglesGroup=join(fnDirSignificant,"angles_group%02d.xmd"%j)
            if not exists(fnAnglesGroup):
                if ctfPresent:
                    fnGroup="ctfGroup%06d@%s/ctf_groups.xmd"%(j,fnDirSignificant)
                    fnGalleryGroup=join(fnDirSignificant,"gallery_group%06d.stk"%j)
                    fnGalleryGroupMd=join(fnDirSignificant,"gallery_group%06d.xmd"%j)
                    self.runJob("xmipp_transform_filter",
                                "-i %s -o %s --fourier binary_file %d@%s --save_metadata_stack %s --keep_input_columns"%\
                                (fnGalleryMd,fnGalleryGroup,j,fnCTFs,fnGalleryGroupMd))
                else:
                    fnGroup=fnImgs
                    fnGalleryGroupMd=fnGalleryMd
                args='-i %s --initgallery %s --odir %s --sym %s --iter 1 --alpha0 %f --alphaF %f --angularSampling %f --maxShift %d '\
                     '--minTilt %f --maxTilt %f --useImed --angDistance %f --dontReconstruct'%\
                     (fnGroup,fnGalleryGroupMd,fnDirSignificant,self.symmetryGroup,alpha,alpha,angleStep,\
                      round(self.angularMaxShift.get()*newXdim/100),self.angularMinTilt.get(),self.angularMaxTilt.get(),2*angleStep)
                self.runJob('xmipp_reconstruct_significant',args,numberOfMpi=self.numberOfThreads.get())
                moveFile(join(fnDirSignificant,"angles_iter001_00.xmd"),join(fnDirSignificant,"angles_group%02d.xmd"%j))
                self.runJob("rm -f",fnDirSignificant+"/images_*iter00?_*.xmd",numberOfMpi=1)
                if j==1:
                    copyFile(fnAnglesGroup, fnAngles)
                else:
                    self.runJob("xmipp_metadata_utilities","-i %s --set union %s"%(fnAngles,fnAnglesGroup),numberOfMpi=1)
        if ctfPresent:
            self.runJob("rm -f",fnDirSignificant+"/gallery*",numberOfMpi=1)
        
    def weightParticles(self):
        iteration=1
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        from math import exp
        # Grab file
        fnDirGlobal=join(fnDirCurrent,"globalAssignment")
        fnAnglesDisc=join(fnDirGlobal,"anglesDisc.xmd")
        fnAngles=join(fnDirCurrent,"angles.xmd")

        copyFile(fnAnglesDisc, fnAngles)
        TsCurrent=self.readInfoField(fnDirGlobal,"sampling",xmipp.MDL_SAMPLINGRATE)
        Xdim=self.readInfoField(fnDirGlobal,"size",xmipp.MDL_XSIZE)
        self.writeInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
        self.writeInfoField(fnDirCurrent,"size",xmipp.MDL_XSIZE,Xdim)
            
        if self.weightSSNR:
            self.runJob("xmipp_metadata_utilities","-i %s --set join %s particleId"%\
                        (fnAngles,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)
        
        md=xmipp.MetaData(fnAngles)
        doWeightSignificance=self.weightSignificance and md.containsLabel(xmipp.MDL_WEIGHT_SIGNIFICANT)
        for objId in md:
            weight=1
            if self.weightSSNR:
                aux=md.getValue(xmipp.MDL_WEIGHT_SSNR,objId)
                weight*=aux
            if doWeightSignificance:
                aux=md.getValue(xmipp.MDL_WEIGHT_SIGNIFICANT,objId)
                weight*=aux
            md.setValue(xmipp.MDL_WEIGHT,weight,objId)
        md.write(fnAngles)
        
        # Weight by Astigmatism
        if self.weightAstigmatism:
            self.runJob("xmipp_transform_filter","-i %s --fourier astigmatism %f --sampling %f"%\
                        (fnAngles,self.weightAstigmatismSigma.get(),TsCurrent))

    def reconstruct(self):
        iteration=1
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        TsCurrent=self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
        fnAngles=join(fnDirCurrent,"angles.xmd")
        fnVol=join(fnDirCurrent,"volume.vol")

        args="-i %s -o %s --sym %s --weight"%(fnAngles,fnVol,self.symmetryGroup)
        row=getFirstRow(fnAngles)
        if row.containsLabel(xmipp.MDL_CTF_DEFOCUSU) or row.containsLabel(xmipp.MDL_CTF_MODEL):
            args+=" --useCTF --sampling %f"%TsCurrent
            if self.phaseFlipped:
                args+=" --phaseFlipped"
        self.runJob("xmipp_reconstruct_fourier",args)
    
    def split(self,i):
        from pyworkflow.em.metadata.utils import getSize
        from math import ceil
        from random import randint
        from time import sleep
        sleep(randint(1,10))
        fnDirCurrent=self._getExtraPath("Iter001")
        fnAngles=self._getExtraPath("Iter001/angles.xmd")
        fnSplit=self._getExtraPath('split%05d.xmd'%i)
        Nimages=getSize(fnAngles)
        Nsplit=ceil(Nimages*self.splitFraction.get())
        self.runJob("xmipp_metadata_utilities","-i %s -o %s --operate random_subset %d"%(fnAngles,fnSplit,Nsplit),numberOfMpi=1)
    
        fnVol=self._getExtraPath('split%05d.vol'%i)
        args="-i %s -o %s --sym %s --weight"%(fnSplit,fnVol,self.symmetryGroup)
        row=getFirstRow(fnAngles)
        if row.containsLabel(xmipp.MDL_CTF_DEFOCUSU) or row.containsLabel(xmipp.MDL_CTF_MODEL):
            TsCurrent=self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
            args+=" --useCTF --sampling %f"%TsCurrent
            if self.phaseFlipped:
                args+=" --phaseFlipped"
        self.runJob("xmipp_reconstruct_fourier",args,numberOfMpi=1)
        self.runJob("xmipp_image_operate","-i %s --minus %s"%(fnVol,join(fnDirCurrent,"volume.vol")))
    
    def generateSplittedVolumes(self):
        fnSplitMD=self._getExtraPath("split.xmd")
        fnBasis=self._getPath("basis.stk")
        fnAvg=self._getExtraPath("Iter001/volume.vol")
        fnOut=self._getPath("splittedVolumes.stk")

        fnMask=self._getExtraPath("Iter001/globalAssignment/mask.vol")
        self.runJob("xmipp_metadata_selfile_create",'-p "%s" -o %s'%(self._getExtraPath("split*.vol"),fnSplitMD))
        args="-i %s --saveBasis %s --generatePCAVolumes %s --avgVolume %s --opca %s"%\
                    (fnSplitMD,fnBasis,self.splitPercentiles.get(),fnAvg,fnOut)
        if exists(fnMask):
            args+=" --mask binary_file %s"%fnMask
        self.runJob("xmipp_volume_pca",args)

    def cleanDirectory(self):
        iteration=1
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        fnGlobal=join(fnDirCurrent,"globalAssignment")
        if exists(fnGlobal):
            cleanPath(join(fnGlobal,"images.stk"))
        if exists(fnGlobal):
            cleanPath(join(fnGlobal,"images.xmd"))
            cleanPath(join(fnGlobal,"significant"))
            cleanPath(join(fnGlobal,"volumeRef.vol"))

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.numberOfMpi.get()!=1:
            errors.append("The number of MPI processors have to be 1")
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
    
