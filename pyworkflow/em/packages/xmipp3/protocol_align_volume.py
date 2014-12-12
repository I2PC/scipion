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
This sub-package contains the  protocol
"""

import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.em.packages.xmipp3.convert import getImageLocation


ALIGN_MASK_CIRCULAR = 0
ALIGN_MASK_BINARY_FILE = 1

ALIGN_ALGORITHM_EXHAUSTIVE = 0
ALIGN_ALGORITHM_LOCAL = 1
ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL = 2
ALIGN_ALGORITHM_FAST_FOURIER = 3


class XmippProtAlignVolume(em.ProtAlignVolume):
    """ Aligns a set of volumes using cross correlation """
    _label = 'align volume'
    
    def __init__(self, **args):
        em.ProtAlignVolume.__init__(self, **args)
        self.stepsExecutionMode = em.STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Volume parameters')
        form.addParam('inputReference', params.PointerParam, pointerClass='Volume', 
                      label="Reference volume", important=True, 
                      help='Reference volume to be used for the alignment.')    
        form.addParam('inputVolumes', params.MultiPointerParam, pointerClass='SetOfVolumes,Volume',  
                      label="Input volume(s)", important=True, 
                      help='Select one or more volumes (Volume or SetOfVolumes)\n'
                           'to be aligned againt the reference volume.')
        
        group1 = form.addGroup('Mask')
        group1.addParam('applyMask', params.BooleanParam, default=False, 
                      label='Apply mask?',
                      help='Apply a 3D Binary mask to the volumes')
        group1.addParam('maskType', params.EnumParam, choices=['circular','binary file'], default=ALIGN_MASK_CIRCULAR, 
                      label='Mask type', display=params.EnumParam.DISPLAY_COMBO, condition='applyMask',
                      help='Select the type of mask you want to apply')
        group1.addParam('maskRadius', params.IntParam, default=-1, condition='applyMask and maskType==%d' % ALIGN_MASK_CIRCULAR,
                      label='Mask radius', 
                      help='Insert the radius for the mask')
        group1.addParam('maskFile', params.PointerParam, condition='applyMask and maskType==%d' % ALIGN_MASK_BINARY_FILE,
                      pointerClass='VolumeMask', label='Mask file', 
                      help='Select the volume mask object')
        
        form.addSection(label='Search strategy')
        form.addParam('alignmentAlgorithm', params.EnumParam, choices=['exhaustive','local', 'exhaustive + local', 'fast fourier'], default=ALIGN_ALGORITHM_EXHAUSTIVE, 
                      label='Alignment algorithm', display=params.EnumParam.DISPLAY_COMBO,
                      help='Exhaustive searches all possible combinations within a search space.'
                            'Local searches around a given position.'
                            'Be aware that the Fast Fourier algorithm requires a special compilation'
                            'of Xmipp (DO_FRM=True in install.sh). It performs the same job as the exhaustive method but much faster.')
        
        anglesCond = 'alignmentAlgorithm!=%d' % ALIGN_ALGORITHM_LOCAL
        
        group = form.addGroup('Angles range', condition=anglesCond, expertLevel=params.LEVEL_ADVANCED)
        
        line = group.addLine('Rotational angle (deg)')
        line.addParam('minRotationalAngle', params.FloatParam, default=0, label='Min')
        line.addParam('maxRotationalAngle', params.FloatParam, default=360, label='Max')
        line.addParam('stepRotationalAngle', params.FloatParam, default=5, label='Step')
        
        line = group.addLine('Tilt angle (deg)', expertLevel=params.LEVEL_ADVANCED)        
        line.addParam('minTiltAngle', params.FloatParam, default=0, label='Min')        
        line.addParam('maxTiltAngle', params.FloatParam, default=180, label='Max')
        line.addParam('stepTiltAngle', params.FloatParam, default=5, label='Step')
        
        line = group.addLine('Inplane angle (deg)', expertLevel=params.LEVEL_ADVANCED)        
        line.addParam('minInplaneAngle', params.FloatParam, default=0, label='Min')        
        line.addParam('maxInplaneAngle', params.FloatParam, default=360, label='Max')
        line.addParam('stepInplaneAngle', params.FloatParam, default=5, label='Step')       
        
        group = form.addGroup('Shifts range', condition=anglesCond, expertLevel=params.LEVEL_ADVANCED)
        line = group.addLine('Shift X (px)')        
        line.addParam('minimumShiftX', params.FloatParam, default=0, label='Min')        
        line.addParam('maximumShiftX', params.FloatParam, default=0, label='Max')
        line.addParam('stepShiftX', params.FloatParam, default=1, label='Step') 
        
        line = group.addLine('Shift Y (px)', expertLevel=params.LEVEL_ADVANCED)        
        line.addParam('minimumShiftY', params.FloatParam, default=0, label='Min')        
        line.addParam('maximumShiftY', params.FloatParam, default=0, label='Max')
        line.addParam('stepShiftY', params.FloatParam, default=1, label='Step')        
        
        line = group.addLine('Shift Z (px)', expertLevel=params.LEVEL_ADVANCED)        
        line.addParam('minimumShiftZ', params.FloatParam, default=0, label='Min')        
        line.addParam('maximumShiftZ', params.FloatParam, default=0, label='Max')
        line.addParam('stepShiftZ', params.FloatParam, default=1, label='Step')         
        
        line = form.addLine('Scale ', expertLevel=params.LEVEL_ADVANCED, condition=anglesCond)        
        line.addParam('minimumScale', params.FloatParam, default=1, label='Min')        
        line.addParam('maximumScale', params.FloatParam, default=1, label='Max')
        line.addParam('stepScale', params.FloatParam, default=0.005, label='Step')          
                        
        group = form.addGroup('Initial values', condition='alignmentAlgorithm==%d' % ALIGN_ALGORITHM_LOCAL)
        line = group.addLine('Initial angles', expertLevel=params.LEVEL_ADVANCED,
                            )        
        line.addParam('initialRotAngle', params.FloatParam, default=0, label='Rot')        
        line.addParam('initialTiltAngle', params.FloatParam, default=0, label='Tilt')
        line.addParam('initialInplaneAngle', params.FloatParam, default=0, label='Psi') 

        line = group.addLine('Initial shifts ', expertLevel=params.LEVEL_ADVANCED)        
        line.addParam('initialShiftX', params.FloatParam, default=0, label='X')        
        line.addParam('initialShiftY', params.FloatParam, default=0, label='Y')
        line.addParam('initialShiftZ', params.FloatParam, default=0, label='Z')    
             
        group.addParam('initialScale', params.FloatParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                      label='Initial scale')  
        
        form.addParallelSection()
        
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _iterInputVolumes(self):
        for pointer in self.inputVolumes:
            item = pointer.get()
            itemId = item.getObjId()
            if isinstance(item, em.Volume):
                item.outputName = self._getExtraPath('output_vol%06d.vol' % itemId)
                yield item
            elif isinstance(item, em.SetOfVolumes):
                for vol in item:
                    vol.outputName = self._getExtraPath('output_vol%06d_%03d.vol' % (itemId, vol.getObjId()))
                    yield vol
        
    def _insertAllSteps(self):
        # Iterate throught all input volumes and align them 
        # againt the reference volume
        refFn = getImageLocation(self.inputReference.get())
        maskArgs = self._getMaskArgs()
        alignArgs = self._getAlignArgs()
        alignSteps = []
        
        for vol in self._iterInputVolumes():
            volFn = getImageLocation(vol)
            stepId = self._insertFunctionStep('alignVolumeStep', refFn, volFn, vol.outputName, 
                                              maskArgs, alignArgs, prerequisites=[])
            alignSteps.append(stepId)
            
        self._insertFunctionStep('createOutputStep', prerequisites=alignSteps)
        
    #--------------------------- STEPS functions --------------------------------------------
    def alignVolumeStep(self, refFn, inVolFn, outVolFn, maskArgs, alignArgs):
        
        args = "--i1 %s --i2 %s --apply %s" % (refFn, inVolFn, outVolFn)
        args += maskArgs
        args += alignArgs
        
        self.runJob("xmipp_volume_align", args)
        
        if self.alignmentAlgorithm == ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL:
            args = "--i1 %s --i2 %s --apply --local" % (refFn, outVolFn)
            self.runJob("xmipp_volume_align", args)
      
    def createOutputStep(self):
        volSet = self._createSetOfVolumes()
        volSet.setSampling(self.inputReference.get().getSamplingRate())
        
        for vol in self._iterInputVolumes():
            vol.setLocation(vol.outputName)
            volSet.append(vol) 
        
        self._defineOutputs(outputVolumes=volSet)
        #TODO: define the source relation from each input
    
    #--------------------------- INFO functions --------------------------------------------
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
        if self.alignmentAlgorithm == ALIGN_ALGORITHM_FAST_FOURIER:
            return ['Chen2013']
        
    #--------------------------- UTILS functions --------------------------------------------
    def _getMaskArgs(self):
        maskArgs = ''
        if self.applyMask:
            if self.maskType == ALIGN_MASK_CIRCULAR:
                maskArgs+=" --mask circular -%d" % self.maskRadius
            else:
                maskArgs+=" --mask binary_file %s" % self.volMask
        return maskArgs
    
    def _getAlignArgs(self):
        alignArgs = ''
        if self.alignmentAlgorithm == ALIGN_ALGORITHM_FAST_FOURIER:
            alignArgs += " --frm"
        elif self.alignmentAlgorithm == ALIGN_ALGORITHM_LOCAL:
            alignArgs += " --local --rot %f %f 1 --tilt %f %f 1 --psi %f %f 1 -x %f %f 1 -y %f %f 1 -z %f %f 1 --scale %f %f 0.005" %\
               (self.initialRotAngle, self.initialRotAngle,
                self.initialTiltAngle, self.initialTiltAngle,
                self.initialInplaneAngle, self.initialInplaneAngle,
                self.initialShiftX, self.initialShiftX,
                self.initialShiftY, self.initialShiftY,
                self.initialShiftZ,self.initialShiftZ,
                self.initialScale, self.initialScale)
        else:
            alignArgs += " --rot %f %f %f --tilt %f %f %f --psi %f %f %f -x %f %f %f -y %f %f %f -z %f %f %f --scale %f %f %f" %\
               (self.minRotationalAngle, self.maxRotationalAngle, self.stepRotationalAngle,
                self.minTiltAngle, self.maxTiltAngle, self.stepTiltAngle,
                self.minInplaneAngle, self.maxInplaneAngle, self.stepInplaneAngle,
                self.minimumShiftX, self.maximumShiftX, self.stepShiftX,
                self.minimumShiftY, self.maximumShiftY, self.stepShiftY,
                self.minimumShiftZ, self.maximumShiftZ, self.stepShiftZ,
                self.minimumScale, self.maximumScale, self.stepScale)
               
        return alignArgs
            