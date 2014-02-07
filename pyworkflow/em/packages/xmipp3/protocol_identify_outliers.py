# **************************************************************************
# *
# * Authors:     Javier Vargas and Adrian Quintana (jvargas@cnb.csic.es aquintana@cnb.csic.es)
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
#from pyworkflow.em.packages.xmipp3.convert import locationToXmipp,\
#    writeSetOfVolumes
from pyworkflow.em.packages.xmipp3.convert import createXmippInputImages
"""
This sub-package contains the protocol projection outliers
"""
from pyworkflow.em import *  
from xmipp3 import ProjMatcher

class XmippProtIdentifyOutliers(ProtClassify3D, ProjMatcher):
    """ Protocol for screening a number of classes comparing them to a volume. """
    _label = 'Identify Outliers'
    
    def __init__(self, **args):
        ProtClassify3D.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputImages', PointerParam, label="Input particles", important=True, 
              pointerClass='SetOfParticles',
              help='Select the input iamges from the project.'
                   'It should be a SetOfParticles class') 
        form.addParam('inputVolume', PointerParam, label="Volume to compare images to", important=True, 
#              pointerClass='Volume',
              pointerClass='SetOfVolumes',
              help='Provide a volume against which images will be compared'
                   'It should be a Volume class')
        form.addParam('volumeIsCTFCorrected', BooleanParam, label="Volume has been corrected by CTF", default=False,
                      help='Volume has been corrected by CTF')
        form.addParam('symmetryGroup', StringParam, default="c1",
              label='Symmetry group', 
              help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=3, expertLevel=LEVEL_EXPERT,
                      label="Angular sampling",
                      help='In degrees. This sampling defines how fine the projection gallery from the volume is explored.')
         
        form.addParallelSection(mpi=8)
        
        
    def _insertAllSteps(self):
        
        
        
         # Projection matching
        self.fnAngles = self._getTmpPath('angles.xmd')
#        self.images = locationToXmipp(*self.inputImages.get().getLocation())
        self.images = createXmippInputImages(self, self.inputImages.get())
        print "taka", self.images

        self._insertFunctionStep("projMatchStep", self.Volume)        

        # Prepare output
        fnOutputImages=self._getPath('images.xmd')
        self.insertRunJobStep("xmipp_metadata_utilities","-i %s --set join %s -o %s" % (self.fnAngles, self.images, fnOutputImages), numberOfMpi=1)
        
        # Produce difference images
        fnDiff=self._getExtraPath("diff.stk")
        fnAligned=self._getExtraPath("images_aligned.xmd")
        self._insertFunctionStep("produceAlignedImages", fnIn=fnOutputImages, fnOut=fnAligned, 
                        volumeIsCTFCorrected=self.VolumeIsCTFCorrected)
        
        # Evaluate each image
        fnAutoCorrelations=self.extraPath("autocorrelations.xmd")
        self.insertRunJobStep("xmipp_image_residuals", "-i %s -o %s --save_metadata_stack %s"%
                              (fnDiff,self.extraPath("autocorrelations.stk"),fnAutoCorrelations),
                              [fnAutoCorrelations],NumberOfMpi=1)
        self.insertRunJobStep("xmipp_metadata_utilities", "-i %s --set merge %s"%(fnAligned,fnAutoCorrelations),NumberOfMpi=1)
        self.insertStep("deleteFile",filename=fnAutoCorrelations)

        # Prepare output
        self.insertRunJobStep("xmipp_metadata_utilities","-i %s --set join %s"%(fnOutputImages,fnAligned),NumberOfMpi=1)
        self.insertRunJobStep("xmipp_metadata_utilities","-i %s --operate sort zScoreResCov desc"%fnOutputImages,NumberOfMpi=1)  

            
        self._insertFunctionStep('createOutputStep', prerequisites=alignSteps)
        
    def _getMaskArgs(self):
        maskArgs = ''
        if self.applyMask:
            if self.maskType.get() == ALIGN_MASK_CIRCULAR:
                maskArgs+=" --mask circular -%d" % self.maskRadius.get()
            else:
                maskArgs+=" --mask binary_file %s" % self.volMask
        return maskArgs
    
    def _getAlignArgs(self):
        alignArgs = ''
        if self.alignmentAlgorithm.get() == ALIGN_ALGORITHM_FAST_FOURIER:
            alignArgs += " --frm"
        elif self.alignmentAlgorithm.get() == ALIGN_ALGORITHM_LOCAL:
            alignArgs += " --local --rot %f %f 1 --tilt %f %f 1 --psi %f %f 1 -x %f %f 1 -y %f %f 1 -z %f %f 1 --scale %f %f 0.005" %\
               (self.initialRotAngle.get(), self.initialRotAngle.get(),\
                self.initialTiltAngle.get(), self.initialTiltAngle.get(),\
                self.initialInplaneAngle.get(), self.initialInplaneAngle.get(),\
                self.initialShiftX.get(), self.initialShiftX.get(),\
                self.initialShiftY.get(), self.initialShiftY.get(),\
                self.initialShiftZ.get(),self.initialShiftZ.get(),\
                self.initialScale.get(), self.initialScale.get())
        else:
            alignArgs += " --rot %f %f %f --tilt %f %f %f --psi %f %f %f -x %f %f %f -y %f %f %f -z %f %f %f --scale %f %f %f" %\
               (self.minRotationalAngle.get(), self.maxRotationalAngle.get(), self.stepRotationalAngle.get(),\
                self.minTiltAngle.get(), self.maxTiltAngle.get(), self.stepTiltAngle.get(),\
                self.minInplaneAngle.get(), self.maxInplaneAngle.get(), self.stepInplaneAngle.get(),\
                self.minimumShiftX.get(), self.maximumShiftX.get(), self.stepShiftX.get(),\
                self.minimumShiftY.get(), self.maximumShiftY.get(), self.stepShiftY.get(),\
                self.minimumShiftZ.get(), self.maximumShiftZ.get(), self.stepShiftZ.get(),\
                self.minimumScale.get(), self.maximumScale.get(), self.stepScale.get())
               
        return alignArgs
        
    def alignVolumeStep(self, refVolFn, volFn, maskArgs, alignArgs):
        args = "--i1 %s --i2 %s --apply" % (refVolFn, volFn)
        args += maskArgs
        args += alignArgs
        
        self.runJob(None, "xmipp_volume_align", args)
        if self.alignmentAlgorithm.get() == ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL:
            args = "--i1 %s --i2 %s --apply --local" % (refVolFn, volFn)
            self.runJob(None, "xmipp_volume_align", args)
      
    def createOutputStep(self):
#        print "summary"
        volumesSet = self._createSetOfVolumes()
        readSetOfVolumes(self.alignedMd, volumesSet)
        volumesSet.copyInfo(self.inputVols)
        volumesSet.write()
        
        self._defineOutputs(outputVolumes=volumesSet)
        self._defineTransformRelation(self.inputVols, volumesSet)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Reference volume: [%s] " % self.inputReferenceVolume.get().getFirstItem().getNameId())
            summary.append("Input volume: [%s] " % self.inputVolumes.get().getNameId())
            summary.append("Alignment method: %s" % self.alignmentAlgorithm.get())
                
            return summary
        
    def validate(self):
        errors = []
        (Xdim,Ydim,Zdim,Ndim)=getImageSize(self.Volume)
        if Xdim!=self.Xdim:
            errors.append("Make sure that the volume and the images have the same size")
        return errors  
        
    def _citations(self):
        if self.alignmentAlgorithm.get() == ALIGN_ALGORITHM_FAST_FOURIER:
            return ['Chen2013']
            