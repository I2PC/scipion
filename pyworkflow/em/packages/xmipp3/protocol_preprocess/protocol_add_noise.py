# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *              Pablo Conesa (pconesa@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import basename

from pyworkflow import VERSION_1_1
from pyworkflow.utils import removeExt
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D
from pyworkflow.em.data import Volume
import pyworkflow.em as em
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles 



class XmippProtAddNoise(ProtRefine3D):
    """    
    Given a sets of volumes or particles the protocol adds noise to them 
    The types of noise are Uniform, Student and Gaussian.
    """
    GAUSSIAN_NOISE = 0
    STUDENT_NOISE = 1
    UNIFORM_NOISE = 2
    _lastUpdateVersion = VERSION_1_1
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        
        form.addParam('noiseType', EnumParam,
                      choices=['Gaussian', 'Student', 'Uniform'],
                      default = 0,
                      label="Noise Type")
        
        form.addParam('gaussianStd', FloatParam, default=0.08, 
                      condition='noiseType == %d' % self.GAUSSIAN_NOISE,
                      label="Standard Deviation", 
                      help='Please, introduce the standard deviation value.'
                      'Mean value can be changed in advanced mode.')
        
        form.addParam('gaussianMean', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      condition='noiseType == %d' % self.GAUSSIAN_NOISE,
                      label="Mean", 
                      help='Please, introduce the mean value (default = 0).')
        
        form.addParam('studentDf', FloatParam, default=1, 
                      condition='noiseType == %d' % self.STUDENT_NOISE,
                      label="Degree of Freedom", 
                      help='Please, introduce the Degree of Freedom.'
                      'Mean value can be changed in advanced mode.')
        
        form.addParam('studentStd', FloatParam, default=0.08, 
                      condition='noiseType == %d' % self.STUDENT_NOISE,
                      label="Standard Deviation", 
                      help='Please, introduce the standard deviation value.'
                      'Mean value can be changed in advanced mode.')
        
        form.addParam('studentMean', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      condition='noiseType == %d' % self.STUDENT_NOISE,
                      label="Mean", 
                      help='Please, introduce the mean value (default = 0).')
        
        form.addParam('uniformMin', FloatParam, default=0, 
                      condition='noiseType == %d' % self.UNIFORM_NOISE,
                      label="Minimum Value", 
                      help='Please, introduce the minimum value. (default = 0)')
        
        form.addParam('uniformMax', FloatParam, default=1, 
                      condition='noiseType == %d' % self.UNIFORM_NOISE,
                      label="Maximum Value", 
                      help='Please, introduce the maximum value (default = 1).')
        
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions ------------------------


    def _insertAllSteps(self):        
        self.micsFn = self._getPath()
        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('addNoiseStep')
        self._insertFunctionStep('createOutputStep')

    def _getTypeOfNoise(self):
        if self.noiseType == self.GAUSSIAN_NOISE:
            kindNoise = 'gaussian'
            noiseParams = '%f %f' % (self.gaussianStd, self.gaussianMean)
        if self.noiseType == self.STUDENT_NOISE:
            kindNoise = 'student'
            noiseParams = '%f %f %f' % (self.studentDf, self.studentStd
                                         , self.studentMean)
        if self.noiseType == self.UNIFORM_NOISE:
            kindNoise = 'uniform'
            noiseParams = '%f %f' % (self.uniformMin, self.uniformMax)
        return kindNoise, noiseParams

    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        
        validateMsgs = []
        if self.input and not self.input.hasValue():
            validateMsgs.append('Please provide input volume.')  
        return validateMsgs

    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputVolume')):
            summary.append("Output volume not ready yet.")
        else:
            summary.append("Volume with %s noise has been obtained" 
                           % (self.getEnumText("noiseType")) )
            #
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolume')):
            messages.append('Noisy volume has been obtained')
        return messages
    
    def _citations(self):
        return ['Do not apply']
    
    def getSummary(self):
        summary = []
        summary.append("Particles analyzed:")
        #summary.append("Particles picked: %d" %coordsSet.getSize())
        return "\n"#.join(summary)
    
    
class XmippProtAddNoiseVolumes(XmippProtAddNoise):
    """    
    Given a set of volumes, or a volume the protocol will add noise to them 
    The types of noise are Uniform, Student and Gaussian.
    """
    _label = 'add noise volume/s'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
       
        form.addParam('input', PointerParam,
                      pointerClass='SetOfVolumes, Volume',
                      label="Input Volume/s", 
                      help='Select a volume or Set of volumes.')
        
        XmippProtAddNoise._defineParams(self, form)
        
    def convertInputStep(self):
        pass
        
    def _getNoisyOutputPath(self, fnvol):
        fnNoisy = self._getExtraPath(removeExt(basename(fnvol)) + '_Noisy.vol')
        return fnNoisy

    def _addNoisetoVolumeStep(self, kindNoise, noiseParams, vol):
        fnvol = vol.getFileName()
        fnNoisy = self._getNoisyOutputPath(fnvol)
        params = " -i %s --type %s %s -o %s" % (fnvol, kindNoise, noiseParams, 
                                                    fnNoisy)
        self.runJob('xmipp_transform_add_noise', params, numberOfMpi=1)

    def addNoiseStep(self):
        kindNoise, noiseParams = self._getTypeOfNoise()
        
        inputSet = self.input.get()
        if isinstance(inputSet, em.Volume):
            self._addNoisetoVolumeStep(kindNoise, noiseParams, inputSet)
        else:
            for vol in self.input.get():
                self._addNoisetoVolumeStep(kindNoise, noiseParams, vol)

    def createOutputStep(self):
        #Output Volume/SetOfVolumes
        volInput = self.input.get()
        if self._isSingleVolume():
            # Create the output with the same class as
            # the input, that should be Volume or a subclass
            # of Volume like VolumeMask
            fnvol = volInput.getFileName()
            fnOutVol = self._getNoisyOutputPath(fnvol)
            
            volClass = volInput.getClass()
            vol = volClass() # Create an instance with the same class of input 
            vol.copyInfo(volInput)
            vol.setFileName(fnOutVol)
            self._defineOutputs(outputVolume=vol)
            self._defineSourceRelation(self.input.get(), vol)
        else:
            volumes = self._createSetOfVolumes()
            volumes.copyInfo(volInput)
            volumes.copyItems(volInput, updateItemCallback=self._updateNoisyPath)
            self._defineOutputs(outputVol=volumes)
            self._defineSourceRelation(self.input.get(), volumes)
            
#         self._defineTransformRelation(self.inputVolumes, self.outputVol)
    def _updateNoisyPath(self, vol, row):
        fnvol = vol.getFileName()
        fnOutVol = self._getNoisyOutputPath(fnvol)
        vol.setFileName(fnOutVol)
        
    def _isSingleVolume(self):
        return isinstance(self.input.get(), Volume)

class XmippProtAddNoiseParticles(XmippProtAddNoise):
    """    
    Given a set of particles, the protocol will add noise to them 
    The types of noise are Uniform, Student and Gaussian.
    """
    _label = 'add noise particles'
    
    #--------------------------- DEFINE param functions --------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
       
        form.addParam('input', PointerParam, pointerClass='SetOfParticles',
                      label="Input particles", 
                      help='Select a set of particles.')
        
        XmippProtAddNoise._defineParams(self, form)

    def convertInputStep(self):
        """ Read the input metadatata.
        """
        # Get the converted input micrographs in Xmipp format
        inputSet = self.input.get()
        inputPath = self._getExtraPath('inputSet')
        fnSet = inputPath+'.xmd'
        writeSetOfParticles(inputSet, fnSet)
        
    def addNoiseStep(self):
        kindNoise, noiseParams = self._getTypeOfNoise()
        params ='--save_metadata_stack'
        params += " -i %s --type %s %s -o %s" % (self._getExtraPath('inputSet.xmd'),
                                      kindNoise, noiseParams,
                                  self.getFileNameNoisyStk())
        
        self.runJob('xmipp_transform_add_noise', params, numberOfMpi=1)
       
    def createOutputStep(self):
        #Output Volume/SetOfVolumes
        particlesSet = self._createSetOfParticles()
        particlesSet.copyInfo(self.input.get())
        particlesSet.copyItems(self.input.get(),
                               updateItemCallback=self._updateParticle)
        self._defineOutputs(outputParticles=particlesSet)
        self._defineSourceRelation(self.input.get(), particlesSet)      

    def getFileNameNoisyStk(self):
        return self._getExtraPath('Noisy.stk')

    def _updateParticle(self, particle, row):
        #fn = particle.getFileName()
        fnOut = self.getFileNameNoisyStk()
        particle.setFileName(fnOut)

