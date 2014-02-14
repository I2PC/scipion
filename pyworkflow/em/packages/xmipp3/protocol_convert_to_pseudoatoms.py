# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Jan 2014
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

import math
from glob import glob

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.protocol.constants import LEVEL_EXPERT, LEVEL_ADVANCED
import xmipp
from pyworkflow.em.packages.xmipp3.convert import getImageLocation

#from xmipp3 import XmippProtocol
NMA_MASK_NONE = 0
NMA_MASK_THRE = 1
NMA_MASK_FILE = 2

# TODO: Move to 3D Tools
class XmippProtConvertToPseudoAtoms(ProtPreprocessVolumes):
    """ Protocol for converting an EM volume into pseudoatoms """
    _label = 'convert to pseudoatoms'
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, label="Input structure", important=True, 
                      pointerClass='Volume')  
        form.addParam('maskMode', EnumParam, choices=['none', 'threshold', 'file'], 
                      default=NMA_MASK_NONE, 
                      label='Mask mode', display=EnumParam.DISPLAY_COMBO,
                      help='')        
        form.addParam('maskThreshold', FloatParam, default=0.01, 
                      condition='maskMode==%d' % NMA_MASK_THRE,
                      label='Threshold value',
                      help='Gray values below this threshold are set to 0')
        form.addParam('volumeMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask volume', condition='maskMode==%d' % NMA_MASK_FILE,
                      )          
        form.addParam('pseudoAtomRadius', FloatParam, default=1, 
                      label='Pseudoatom radius (vox)',
                      help='Pseudoatoms are defined as Gaussians whose \n'
                           'standard deviation is this value in voxels') 
        form.addParam('pseudoAtomTarget', FloatParam, default=5,
                      expertLevel=LEVEL_ADVANCED, 
                      label='Volume approximation error(%)',
                      help='This value is a percentage (between 0.001 and 100) \n'
                           'specifying how fine you want to approximate the EM \n'
                           'volume by the pseudoatomic structure. Lower values \n'
                           'imply lower approximation error, and consequently, \n'
                           'more pseudoatoms.')        
              
        form.addParallelSection(threads=4, mpi=1)    
             
    def _insertAllSteps(self):
        inputStructure = self.inputStructure.get()
        fnMask = self._insertMaskStep()
        self.sampling = inputStructure.getSamplingRate()
        self.fnIn=getImageLocation(inputStructure)
        self._insertFunctionStep('convertToPseudoAtomsStep', self.fnIn, fnMask)
        self._insertFunctionStep('createChimeraScriptStep')
        self._insertFunctionStep('createOutputStep')
        
    def _insertMaskStep(self):
        """ Check the mask selected and insert the necessary steps.
        Return the mask filename if needed.
        """
        fnMask = ''
        if self.maskMode == NMA_MASK_THRE:
            fnMask = self._getExtraPath('mask.vol')
            maskParams = '-i %s -o %s --select below %f --substitute binarize' % (getImageLocation(self.inputStructure.get()),
                                                                                  fnMask, self.maskThreshold.get())
            self._insertRunJobStep('xmipp_transform_threshold', maskParams)
        elif self.maskMode == NMA_MASK_FILE:
            fnMask = getImageLocation(self.volumeMask.get())
        return fnMask
        
    def convertToPseudoAtomsStep(self, inputFn, fnMask):
        prefix = 'pseudoatoms'
        outputFn = self._getPath(prefix)
        sampling = self.sampling
        sigma = sampling * self.pseudoAtomRadius.get() 
        targetErr = self.pseudoAtomTarget.get()
        nthreads = self.numberOfThreads.get()
        params = "-i %(inputFn)s -o %(outputFn)s --sigma %(sigma)f --thr %(nthreads)d "
        params += "--targetError %(targetErr)f --sampling_rate %(sampling)f -v 2 --intensityColumn Bfactor"
        if fnMask:
            params += " --mask binary_file %(fnMask)s"
        print params%locals()
        self.runJob("xmipp_volume_to_pseudoatoms", params % locals())
        for suffix in ["_approximation.vol", "_distance.hist"]:
            moveFile(self._getPath(prefix+suffix), self._getExtraPath(prefix+suffix))
        cleanPattern(self._getPath(prefix+'_*'))
        
    def createChimeraScriptStep(self):
        radius = self.sampling * self.pseudoAtomRadius.get() 
        input = self.inputStructure.get()
        localInputFn = self._getBasePath(self.fnIn)
        createLink(self.fnIn, localInputFn)
        fhCmd = open(self._getPath("chimera.cmd"),'w')
        fhCmd.write("open pseudoatoms.pdb\n")
        fhCmd.write("rangecol bfactor,a 0 white 1 red\n")
        fhCmd.write("setattr a radius %f\n" % radius)
        fhCmd.write("represent sphere\n")
        fhCmd.write("open %s\n" % basename(localInputFn))
         
        threshold = 0.01
        if self.maskMode == NMA_MASK_THRE:
            self.maskThreshold.get()
        xdim, _, _, _ = input.getDim()
        origin = xdim / 2
        fhCmd.write("volume #1 level %f transparency 0.5 voxelSize %f originIndex %d\n" % (threshold, self.sampling, origin))
        fhCmd.close()
     
    def createOutputStep(self):
        pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), self.outputPdb)

    def _summary(self):
        summary = []
        summary.append('Pseudoatom radius (voxels): %f'%self.pseudoAtomRadius.get())
        summary.append('Approximation target error (%%): %f'%self.pseudoAtomTarget.get())
        return summary

    def _methods(self):
        summary = []
        summary.append('We converted the volume %s into a pseudoatomic representation with Gaussian atoms (sigma=%f A and a target error'\
                       ' of %f%%) [Nogales2013].'%(self.inputStructure.get().getNameId(),
                                     self.pseudoAtomRadius.get()*self.inputStructure.get().getSamplingRate(),
                                     self.pseudoAtomTarget.get()));
        if self.hasAttribute('outputPdb'):
            summary.append('We refer to the pseudoatomic model as %s.'%self.outputPdb.getNameId())
        return summary

    def _citations(self):
        return ['Nogales2013']
        