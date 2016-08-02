# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco     (jgomez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import os

from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume, VolumeMask
import pyworkflow.em.metadata as md

from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase



class ProtRelionSort(ProtRelionBase):
    """
    Relion post-processing protocol for automated masking,
    overfitting estimation, MTF-correction and B-factor sharpening.
    """
    _label = 'sort particles'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        
        form.addParam('inputParticles', PointerParam,
                      pointerClass="ProtRelionRefine3D, ProtRelionClassify2D, ProtRelionClassify3D",
                      label='Select a previous relion protocol',
                      help='Select a previous relion class2D, class3D or refine3D protocol.')
        #form.addParam('inputreferences', PointerParam,
        #              pointerClass="SetOfClasses2D, SetOfClasses3D, Volume",
        #              label='Select references',
        #              help='Select references: a set of classes or a 3D volume')
        form.addParam('initMaskThreshold', FloatParam, default=0.02,
                      condition='doAutoMask',
                      label='Initial binarisation threshold',
                      help='This threshold is used to make an initial binary mask from the'
                           "average of the two unfiltered half-reconstructions. If you don't"
                           "know what value to use, display one of the unfiltered half-maps in a 3D"
                           "surface rendering viewer, like Chimera, and find the lowest threshold"
                           " that gives no noise peaks outside the reconstruction.")
        form.addParam('maskDiameterA', IntParam, default=-1,
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a soft circular mask '
                           'with this <diameter>. '
                           'Make sure this diameter is not set too small because that may mask '
                           'away part of the signal! If set to a value larger than the image '
                           'size no masking will be performed.\n\n'
                           'The same diameter will also be used for a spherical mask of the '
                           'reference structures if no user-provided mask is specified.')
        form.addParam('doLowPass', IntParam, default=-1,
                      label='Low pass filter references to (A):',
                      help='Lowpass filter in Angstroms for the references (prevent Einstein-from-noise!)')

        form.addParam('doInvert', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Invert contrast of references?',
                      help='Density in particles is inverted w.r.t. density in template')
        form.addParam('doCTF', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be corrected inside the MAP refinement. '
                           'The resulting algorithm intrinsically implements the optimal linear, '
                           'or Wiener filter. Note that input particles should contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED, condition='doCTF',
                      label='Ignore CTFs until their first peak?',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. '
                           'Therefore, this option is not generally recommended.')
        form.addParam('minZ', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Min Z-value?',
                      help='Minimum Z-value to count in the sorting of outliers')

        form.addParallelSection(threads=0, mpi=1)
            
    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertSortingStep()
        self._insertFunctionStep('createOutputStep')
        
        
    def _initialize(self):
        # ToDo: implement a better way to get the pattern of unmasked maps.
        self.input = self.protRelionRefine.get()._getExtraPath('relion')
        print "INPUT: ", self.input
        
    def _insertSortingStep(self):
        
        output = self._getExtraPath('sorted')
        self.samplingRate = self.protRelionRefine.get().inputParticles.get().getSamplingRate()
        
        args = " --i %s --o %s --angpix %f" %(self.input, output, self.samplingRate)

        if self.doAutoMask:
            args += " --auto_mask --inimask_threshold %f" % self.initMaskThreshold.get()
            args += " --extend_inimask %d" % self.extendInitMask.get()
            args += " --width_mask_edge %d" % self.addMaskEdge.get()
        else:
            args += ' --mask %s' % self.mask.get().getFileName()
            
        mtfFile = self.mtf.get()
        if mtfFile:
            args += ' --mtf %s' % mtfFile
            
        if self.doAutoBfactor:
            args += ' --auto_bfac --autob_lowres %f' % self.bfactorLowRes.get()
            args += ' --autob_highres %f' % self.bfactorHighRes.get()
        else:
            args += ' --adhoc_bfac %f' % self.bfactor.get()
            
        if self.skipFscWeighting:
            args += ' --skip_fsc_weighting --low_pass %f' % self.lowRes.get()
            
        # Expert params
        args += ' --filter_edge_width %d' % self.filterEdgeWidth.get()
        args += ' --randomize_at_fsc %f' % self.randomizeAtFsc.get()
        
        self._insertFunctionStep('sortingStep', args)
        
    #--------------------------- STEPS functions -------------------------------
    def postProcessStep(self, params):
        self.runJob('relion_postprocess', params)
        
    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getExtraPath('postprocess.mrc'))
        volume.setSamplingRate(self.samplingRate)
        vol = self.protRelionRefine.get().outputVolume
        mask = VolumeMask()
        mask.setFileName(self._getExtraPath('postprocess_automask.mrc'))
        mask.setSamplingRate(self.samplingRate)

        self._defineOutputs(outputVolume=volume)
        self._defineOutputs(outputMask=mask)
        self._defineSourceRelation(vol, volume)
        self._defineSourceRelation(vol, mask)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        mtfFile = self.mtf.get()

        if mtfFile and not os.path.exists(mtfFile):
            errors.append("Missing MTF-file '%s'" % mtfFile)

        return errors
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        postStarFn = self._getExtraPath("postprocess.star")
        if os.path.exists(postStarFn):
            mdResol = md.RowMetaData(postStarFn)
            resol = mdResol.getValue(md.RLN_POSTPROCESS_FINAL_RESOLUTION)
            summary.append("Final resolution: *%0.2f*" % resol)
        
        return summary
    
    #--------------------------- UTILS functions -------------------------------
