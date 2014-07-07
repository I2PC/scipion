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
"""
This sub-package contains the XmippRCT protocol
"""
from pyworkflow.em import ProtInitialVolume, String, PointerParam, IntParam, StringParam, BooleanParam, FloatParam, STEPS_PARALLEL
from pyworkflow.protocol.constants import  LEVEL_EXPERT
from convert import readSetOfVolumes, getImageLocation
from pyworkflow.utils import exists
from pyworkflow.em.data import Volume
from xmipp3 import XmippMdRow

from itertools import izip
import xmipp


class XmippProtRCT(ProtInitialVolume):
    """ Computes from a set of projections/classes using RCT algorithm """
    _label = 'random conical tilt'
    
    def __init__(self, **args):
        ProtInitialVolume.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.summaryInfo = String()

    #--------------------------- DEFINE param functions --------------------------------------------        
    def _defineParams(self, form):
        form.addSection(label='Input')
                
        #TODO: Input can be either a SetOfParticles or a SetOfClasses2D
        form.addParam('inputParticlesTiltPair', PointerParam, label="Input particles tilt pair", important=True, 
                      pointerClass='ParticlesTiltPair',
                      help='Select the input particles tilt pair from the project.')   
        
        form.addParam('inputImages', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfImages',
                      help='Select the input images or classes from the project.')    

        form.addSection(label='Alignment parameters')

        form.addParam('thinObject', BooleanParam, default=False, 
                      label='Thin Object', 
                      help=' If the object is thin, then the tilted projections can be stretched to match the untilted projections')
                       
        form.addParam('maxShift', IntParam, default="10", expertLevel=LEVEL_EXPERT,
                      label="Maximum allowed shift for tilted particles (pixels)", 
                      help='Particles that shift more will be discarded. A value larger than the '
                      'image size will not discard any particle.')
        
        form.addParam('skipTranslation', BooleanParam, default=False, expertLevel=LEVEL_EXPERT,
                      label='Skip tilted translation alignment', 
                      help=' If the tilted image quality is very low, then this alignment might result in poor estimates.')

        form.addSection(label='Reconstruction')

        form.addParam('additionalParams', StringParam, default="-n 5 -l 0.01",
                      label='Additional reconstruction parameters', 
                      help='See: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Reconstruct_art_v31')
        
        form.addParam('doFilter', BooleanParam, default=True, 
                      label='Filter reconstructed volumes?', 
                      help='Filtering may be useful to remove noise, especially when few particles '
                      'contribute to the reconstruction.')
        
        form.addParam('resoLowPassFilter', FloatParam, default=0.2,
                      label='Resolution of the low-pass fitler (dig.freq)',
                      help='Resolution of the low-pass fitler (dig.freq)')        

        form.addParallelSection(mpi=2)
            
         
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):

        rctImagesFn = self._getExtraPath('rct_images.xmd')
        
        #We are now using setofimages but in the future it should work with setofclasses
        classNo = 1
        blockMd = "class%06d_images@%s" % (classNo, rctImagesFn)
        self.appendRctImages(blockMd)
                                        
        classNameIn = blockMd
        classNameOut = self._getExtraPath("rct_images_%06d.xmd" % classNo)
        classVolumeOut = self._getPath("rct_%06d.vol" % classNo)
        
        if self.inputImages.get().hasRepresentative():
            classImage = self.inputImages.get().getRepresentative()
        else:
            classImage = ""
        
        self._insertFunctionStep('reconstructClass', classNameIn, classNameOut, classImage, classVolumeOut)
        
        self._insertFunctionStep('createOutputStep', classVolumeOut)
       
    def appendRctImages(self, blockMd):
        
        rctImagesMd = xmipp.MetaData()
        
        uImages = self.inputImages.get()
        tImages = self.inputParticlesTiltPair.get().getTilted()
        sangles = self.inputParticlesTiltPair.get().getCoordsPair().get().getAngles()
        uMics = uImages.getCoordinates().getMicrographs()
        tMics = tImages.getCoordinates().getMicrographs()
        print "ANGLES"
        sangles.printAll()
                    
        for uImg in uImages:
            imgId = uImg.getObjId()
            tImg = tImages[imgId]
            objId = rctImagesMd.addObject()
            pairRow = XmippMdRow()
            pairRow.setValue(xmipp.MDL_IMAGE, getImageLocation(uImg))
            uCoord = uImg.getCoordinate()
            micId = uCoord.getMicId()
            print "MICID=%s" % micId
            uMic = uMics[micId]
            angles = sangles[micId]
            pairRow.setValue(xmipp.MDL_MICROGRAPH, uMic.getFileName())
            pairRow.setValue(xmipp.MDL_XCOOR, uCoord.getX())
            pairRow.setValue(xmipp.MDL_YCOOR, uCoord.getY())
            pairRow.setValue(xmipp.MDL_ENABLED, 1)
            pairRow.setValue(xmipp.MDL_ITEM_ID, long(imgId))
            #TODO: WHERE DO WE GET ALL THISINFO FROM??? (ask J.M.)
            pairRow.setValue(xmipp.MDL_REF, 1)
            pairRow.setValue(xmipp.MDL_FLIP, bool(1))
            pairRow.setValue(xmipp.MDL_SHIFT_X, float(1))
            pairRow.setValue(xmipp.MDL_SHIFT_Y, float(1))
            pairRow.setValue(xmipp.MDL_ANGLE_PSI, float(1))
            
            pairRow.setValue(xmipp.MDL_IMAGE_TILTED, getImageLocation(tImg))
            tMic = tMics[micId]
            pairRow.setValue(xmipp.MDL_MICROGRAPH_TILTED, tMic.getFileName())
            (angleY, angleY2, angleTilt) = angles.getAngles()
            pairRow.setValue(xmipp.MDL_ANGLE_Y, float(angleY))
            pairRow.setValue(xmipp.MDL_ANGLE_Y2, float(angleY2))
            pairRow.setValue(xmipp.MDL_ANGLE_TILT, float(angleTilt))
            
            #pairRow.appendToMd(rctImagesMd, objId)
            pairRow.writeToMd(rctImagesMd, objId)
            rctImagesMd.write(blockMd)
            
    def reconstructClass(self, classIn, classOut, classImage, classVolumeOut):

    # If class image doesn't exists, generate it by averaging
        if len(classImage) == 0:
            classRootOut = classOut.replace(".xmd", "") + "_"
            statsFn = self._getExtraPath('stats.xmd')
            args = "-i %(classIn)s --save_image_stats %(classRootOut)s -o %(statsFn)s"
            self.runJob("xmipp_image_statistics", args % locals())

            classImage = classRootOut + "average.xmp"
            
        centerMaxShift = self.maxShift.get()
        args = "-i %(classIn)s -o %(classOut)s --ref %(classImage)s --max_shift %(centerMaxShift)d" % locals()
        
        if self.thinObject.get():
            args += " --do_stretch"
        
        if self.skipTranslation.get():
            args += " --do_not_align_tilted"
        
        self.runJob("xmipp_image_align_tilt_pairs", args)
        
        reconstructAdditionalParams = self.additionalParams.get()

        args = "-i %(classOut)s -o %(classVolumeOut)s %(reconstructAdditionalParams)s" % locals()
        self.runJob("xmipp_reconstruct_art", args)
        
        if exists(classVolumeOut):
            mdFn = self._getPath('volumes.xmd')
            md = xmipp.MetaData()
            
            if exists(mdFn):
                md.read(mdFn)
            objId = md.addObject()
            md.setValue(xmipp.MDL_IMAGE, classVolumeOut, objId)
                        
            if self.doFilter.get():
                filteredVolume = classVolumeOut.replace('.vol','_filtered.vol')
                lowPassFilter = self.resoLowPassFilter.get()
                args = "-i %(classVolumeOut)s -o %(filteredVolume)s --fourier low_pass %(lowPassFilter)f" % locals()
                self.runJob("xmipp_transform_filter", args)
                objId = md.addObject()
                md.setValue(xmipp.MDL_IMAGE, filteredVolume, objId)
            md.write(mdFn)
                    
         
    def createOutputStep(self, classVolumeOut):
        if self.doFilter.get():
            classVolumeOut = classVolumeOut.replace('.vol','_filtered.vol')
                    
        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate(self.inputParticlesTiltPair.get().getUntilted().getSamplingRate())
        
        #TODO: Change when we have more than one class
        #for k in range(1, self.numberOfClasses):
            
        vol = Volume()
        vol.setFileName(classVolumeOut)
        volumesSet.append(vol)

        self._defineOutputs(outputVolumes=volumesSet)
        self._defineSourceRelation(self.inputParticlesTiltPair.get(), volumesSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Summary not implemented yet.")
        return summary
        
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
            
    #--------------------------- UTILS functions --------------------------------------------           

