# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This module contains the protocol base class for Relion protocols
"""
from pyworkflow.em import *
from pyworkflow.utils.path import moveFile
from convert import writeSetOfParticles, readSetOfParticles
from pyworkflow.protocol.params import Positive

from protocol_base import ProtRelionBase

class ProtRelionPreprocessParticles(ProtProcessParticles, ProtRelionBase):
    """ Wrapper to Relion preprocess program.
    """
    _label = 'preprocess particles'
    
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
                   
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label="Input particles",  
                      help='Select the input images from the project.')   
        form.addParam('doRemoveDust', BooleanParam, default=False,
                      label='Remove dust from particles', 
                      help='Remove dust from particles.')
        form.addParam('whiteDust', IntParam, validators=[Positive],
                      condition='doRemoveDust',
                      label='White dust',
                      help='Sigma-values above which white dust will be removed.')
        form.addParam('blackDust', IntParam, validators=[Positive],
                      condition='doRemoveDust',
                      label='Black dust',
                      help='Sigma-values above which black dust will be removed.')   
        form.addParam('doNormalize', BooleanParam, default=False,
                      label='Normalize',
                      help='If set to True, particles will be normalized in the way RELION prefers it.')
        form.addParam('backRadius', IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.).')        
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        form.addParam('doScale', BooleanParam, default=False,
                      label='Scale particles', 
                      help='Re-scale the particles to this size (in pixels).')
        form.addParam('scaleSize', IntParam, validators=[Positive],
                      condition='doScale',
                      label='Scale size',
                      help='New particle size.')  
        form.addParam('doWindow', BooleanParam, default=False,
                      label='Window particles', 
                      help='Re-window the particles to this size (in pixels).')
        form.addParam('windowSize', IntParam, validators=[Positive],
                      condition='doWindow',
                      label='Window size',
                      help='New window size.')  
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        self.imgStar = self._getPath('input_particles.star')
        self.imgFn = self._getPath('input_particles.mrcs')
        self.outFn = self._getPath('particles.mrcs')
        self._insertFunctionStep("convertInputStep")
        print "IMAGES STAR=%s" % self.imgStar
        self._insertFunctionStep('processStep')
        
        self._insertFunctionStep('createOutputStep')
        
    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, self.imgStar, self.imgFn)
        
    def processStep(self):
        # Enter here to generate the star file or to preprocess the images
        
        size = self._getSize()
        
        imgFn = os.path.relpath(self.imgFn, self._getPath())
        
        params = ' --operate_on %(imgFn)s'
        if self.doNormalize:
            radius = self.backRadius.get()
            if radius <= 0:
                radius = size
            params = params + ' --norm --bg_radius %(radius)s'
            
        if self.doRemoveDust:
            wDust = self.whiteDust.get()
            bDust = self.blackDust.get()
            params = params + ' --white_dust %(wDust)s --black_dust %(bDust)s'
            
        if self.doInvert:
            params = params + ' --invert_contrast'            
        
        if self.doScale:
            eSize = self.scaleSize.get()
            params = params + ' --scale %(eSize)s'

        if self.doWindow:
            wSize = self.windowSize.get()
            params = params + ' --window %(wSize)s'

        self.runJob(self._getProgram('relion_preprocess'), params % locals(), cwd=self._getPath())
                             
        outputMrcs = glob(self._getPath('particles*.mrcs')) # In Relion 1.3 it is produces particles.mrcs.mrcs
        outFile = outputMrcs[0]
        if outFile != self.outFn:
            moveFile(outFile, self.outFn)
    
    def createOutputStep(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        outputStar = self._getPath('particles.star')
        readSetOfParticles(outputStar, imgSet, preprocessRow=self._preprocessRow)
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(self.inputParticles.get(), self.outputParticles)

#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        validateMsgs = []
        if self.doScale and self.scaleSize.get() % 2 != 0:
            validateMsgs += ["Only re-scaling to even-sized images is allowed in RELION."]
        if self.doWindow and self.windowSize.get() % 2 != 0:
            validateMsgs += ["Only re-windowing to even-sized images is allowed in RELION."]

        if self.doNormalize:
            size = self._getSize()
            if self.backRadius.get() > size:
                validateMsgs.append('Set a valid Background radius less than %d' % size)
        return validateMsgs
    
    def _summary(self):
        summary = []
        summary.append("Input particles: %s" % self.inputParticles.get().getName())
        
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Dust removal: %s" % self.doRemoveDust)
            summary.append("Normalize the background: %s" % self.doNormalize)
            summary.append("Invert contrast: %s" % self.doInvert)
            summary.append("Scaled particles to size: %s" % self.scaleSize.get())
            summary.append("Windowed particles to size: %s" % self.windowSize.get())
        return summary
                
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getSize(self):
        """ get the size of SetOfParticles object"""

        Xdim = self.inputParticles.get().getDimensions()[0]
        size = int(Xdim/2)
        return size
    
    def _preprocessRow(self, imgRow):
        from convert import setupCTF, prependToFileName
        prependToFileName(imgRow, self._getPath())
        setupCTF(imgRow, self.inputParticles.get().getSamplingRate())