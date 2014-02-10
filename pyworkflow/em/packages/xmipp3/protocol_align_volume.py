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
from pyworkflow.em.packages.xmipp3.convert import locationToXmipp, writeSetOfVolumes
"""
This sub-package contains the  protocol
"""
from pyworkflow.em import *  
from convert import readSetOfClasses2D, createXmippInputClasses2D, createXmippInputVolumes, readSetOfVolumes
from math import floor
from xmipp import MetaData, MD_APPEND, MDL_MAXCC, MDL_WEIGHT, MDL_IMAGE, \
    MDL_VOLUME_SCORE_SUM, MDL_VOLUME_SCORE_SUM_TH, MDL_VOLUME_SCORE_MEAN, MDL_VOLUME_SCORE_MIN
#    removeFilenamePrefix
from pyworkflow.utils.path import moveFile, cleanPath, copyFile, removeExt
from protlib_xmipp import getMdSize


ALIGN_MASK_CIRCULAR = 0
ALIGN_MASK_BINARY_FILE = 1

ALIGN_ALGORITHM_EXHAUSTIVE = 0
ALIGN_ALGORITHM_LOCAL = 1
ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL = 2
ALIGN_ALGORITHM_FAST_FOURIER = 3

class XmippProtAlignVolume(ProtAlignVolume):
    """ Protocol to align a set of volumes. """
    _label = 'align volume'
    
    def __init__(self, **args):
        ProtAlignVolume.__init__(self, **args)
        
    def _defineParams(self, form):
        form.addSection(label='Volume parameters')
        form.addParam('inputReferenceVolume', PointerParam, label="Reference volume", important=True, 
#                      pointerClass='Volume', 
                    pointerClass='SetOfVolumes',
                      help='Reference volume to be used for the alignment')    
        form.addParam('inputVolumes', PointerParam, label="Input volume/s", important=True, 
                      pointerClass='SetOfVolumes,Volume',  
                      help='The input volume will be aligned to the reference volume')
        
        form.addSection(label='Mask')
        form.addParam('applyMask', BooleanParam, default=False, 
                      label='Apply mask?',
                      help='Apply a 3D Binary mask to the volumes')
        form.addParam('maskType', EnumParam, choices=['circular','binary file'], default=ALIGN_MASK_CIRCULAR, 
                      label='Mask type', display=EnumParam.DISPLAY_COMBO, condition='applyMask',
                      help='Select the type of mask you want to apply')
        form.addParam('maskRadius', IntParam, default=-1, condition='applyMask and maskType==%d' % ALIGN_MASK_CIRCULAR,
                      label='Mask radius', 
                      help='Insert the radius for the mask')
        form.addParam('maskFile', PointerParam, condition='applyMask and maskType==%d' % ALIGN_MASK_BINARY_FILE,
                      pointerClass='VolumeMask', label='Mask file', 
                      help='Select the volume mask object')
        
        form.addSection(label='Search method and space')
        form.addParam('alignmentAlgorithm', EnumParam, choices=['exhaustive','local', 'exhaustive + local', 'fast fourier'], default=ALIGN_ALGORITHM_EXHAUSTIVE, 
                      label='Alignment algorithm', display=EnumParam.DISPLAY_COMBO,
                      help='Exhaustive searches all possible combinations within a search space.'
                            'Local searches around a given position.'
                            'Be aware that the Fast Fourier algorithm requires a special compilation'
                            'of Xmipp (DO_FRM=True in install.sh). It performs the same job as the exhaustive method but much faster.')
        
        form.addParam('minRotationalAngle', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Minimum rotational angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER), 
                      help='')
        form.addParam('maxRotationalAngle', FloatParam, default=360, expertLevel=LEVEL_EXPERT,
                      label='Maximum rotational angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('stepRotationalAngle', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Step rotational angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        
        form.addParam('minTiltAngle', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Minimum tilt angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')        
        form.addParam('maxTiltAngle', FloatParam, default=180, expertLevel=LEVEL_EXPERT,
                      label='Maximum tilt angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('stepTiltAngle', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Step tilt angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
                
        form.addParam('minInplaneAngle', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Minimum in-plane angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('maxInplaneAngle', FloatParam, default=360, expertLevel=LEVEL_EXPERT,
                      label='Maximum in-plane angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('stepInplaneAngle', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Step in-plane angle', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        
        form.addParam('minimumShiftX', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Minimum shiftX', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('maximumShiftX', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Maximum shiftX', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('stepShiftX', FloatParam, default=1, expertLevel=LEVEL_EXPERT,
                      label='Step shiftX', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('minimumShiftY', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Minimum shiftY', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('maximumShiftY', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Maximum shiftY', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('stepShiftY', FloatParam, default=1, expertLevel=LEVEL_EXPERT,
                      label='Step shiftY', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('minimumShiftZ', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Minimum shiftZ', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('maximumShiftZ', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Maximum shiftZ', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        form.addParam('stepShiftZ', FloatParam, default=1, expertLevel=LEVEL_EXPERT,
                      label='Step shiftZ', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL, ALIGN_ALGORITHM_FAST_FOURIER),
                      help='')
        
        form.addParam('minimumScale', FloatParam, default=1.0, expertLevel=LEVEL_EXPERT,
                      label='Minimum scale', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL),
                      help='')
        form.addParam('maximumScale', FloatParam, default=1.0, expertLevel=LEVEL_EXPERT,
                      label='Maximum scale', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL),
                      help='')
        form.addParam('stepScale', FloatParam, default=0.005, expertLevel=LEVEL_EXPERT,
                      label='Step scale', condition='alignmentAlgorithm==%d or alignmentAlgorithm==%d' % (ALIGN_ALGORITHM_EXHAUSTIVE, ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL),
                      help='')
        
        form.addParam('initialRotAngle', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Initial rotational angle', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')
        form.addParam('initialTiltAngle', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Initial tilt angle', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')        
        form.addParam('initialInplaneAngle', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Initial in-plane angle', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')
        form.addParam('initialShiftX', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Initial shiftX', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')
        form.addParam('initialShiftY', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Initial shiftY', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')  
        form.addParam('initialShiftZ', FloatParam, default=0, expertLevel=LEVEL_EXPERT,
                      label='Initial shiftZ', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')
        form.addParam('initialScale', FloatParam, default=1, expertLevel=LEVEL_EXPERT,
                      label='Initial scale', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL,
                      help='')  
        
        form.addParallelSection()
        
        
    def _insertAllSteps(self):
        
        # Check volsMd is a volume or a stack
        inputRefVol = self.inputReferenceVolume.get().getFirstItem()
        inputVols = self.inputVolumes.get()
        
        if isinstance(inputVols, Volume):
            inputVols = [inputVols]
        refVolFn = locationToXmipp(*inputRefVol.getLocation())
        
        maskArgs = self._getMaskArgs()
        alignArgs = self._getAlignArgs()
             
        volsFn = []

        self.alignedMd = self._getPath("volume.xmd")
        writeSetOfVolumes(inputVols, alignedMd)
        md = xmipp.MetaData(alignedMd)
        for objId in md:
            volFileName = md.getValue(xmipp.MDL_IMAGE, objId)
            volFn = self._getPath(basename(volFileName))
            print "taka", volFileName, volFn
            
            volsFn.append(volFn)
            
            ImageHandler().convert(volFileName, volFn)
            md.setValue(xmipp.MDL_IMAGE, volFn, objId)
         
            
#        md.write(alignedMd)    
        
        
        alignSteps = []
        for item in volsFn:
#           volFn = self._getPath(basename(vol.getFileName()))
#           copyFile(vol.getFileName(), volFn)
            
            
#            volFn = locationToXmipp(*vol.getLocation())
            
            stepId = self._insertFunctionStep('alignVolumeStep', refVolFn, item, 
                                              maskArgs, alignArgs, prerequisites=[])
            
            alignSteps.append(stepId)
            
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
        
        self.runJob("xmipp_volume_align", args)
        if self.alignmentAlgorithm.get() == ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL:
            args = "--i1 %s --i2 %s --apply --local" % (refVolFn, volFn)
            self.runJob("xmipp_volume_align", args)
      
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
        
    def _citations(self):
        if self.alignmentAlgorithm.get() == ALIGN_ALGORITHM_FAST_FOURIER:
            return ['Chen2013']
            