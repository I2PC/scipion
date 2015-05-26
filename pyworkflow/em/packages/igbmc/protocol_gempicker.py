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

from os.path import join

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.utils.properties import Message 


MASK_CIRCULAR = 0
MASK_OBJECT = 1


class ProtGemPicker(em.ProtParticlePicking):
    """Protocol to pick particles in a set of micrographs using appion dogpicker"""
    _label = 'gempicker'
        
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        
        em.ProtParticlePicking._defineParams(self, form)
        
        form.addParam('inputReferences', params.PointerParam, 
                      pointerClass='SetOfAverages',
                      label='Input References', important=True,
                      help="Template averages to be used in picking.")
        
        form.addParam('refsHaveInvertedContrast', params.BooleanParam,
                      label='References have inverted contrast',
                      help='Set to Yes to indicate that the reference have inverted \n'
                           'contrast with respect to the particles in the micrographs.')
        
        form.addParam('rotAngle', params.IntParam, default=5,
                      label='Rotational angle search',
                      help='In-plane rotating angle in degrees (0 = no rotation)')
        
        line = form.addLine('Threshold (0, 1]', 
                            help='Threshold value for picking, select low and high values.')
        line.addParam('threshold', params.FloatParam, default=0.2,
                      label='Low')
        line.addParam('thresholdHigh', params.FloatParam, default=0.9,
                      label='High')
        
        form.addParam('maxPiks', params.IntParam, default=0,
                      label='Max particles per micrograph',
                      help='Maximum number of particles picked from each micrograph (0 = no limit)')
        form.addParam('sizePiks', params.IntParam, default=0,
                      label='Box size (pix)', 
                      help='Size of picked images (if not provided or = 0, use search image size)')
        
        form.addParam('maskType', params.EnumParam, 
                      choices=['circular', 'object'], default=0, 
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Mask type', 
                      help='Select which type of mask do you want to apply. '
                           'Only the pixels beneath this mask will be analyzed. '
                           'In the simplest case, a circular mask can be used. '
                           'Alternatively, a custom mask can be used '
                           'which follows the contour of the particle (but not too tightly).')
        form.addParam('maskRadius', params.IntParam, default=-1,
                      condition='maskType==%d' % MASK_CIRCULAR,
                      label='Mask radius (px)', 
                      help='If -1, the entire image (in pixels) will be considered.')
        form.addParam('inputMasks', params.MultiPointerParam, pointerClass='Mask', 
                      condition='maskType==%d' % MASK_OBJECT, 
                      label="Mask objects", 
                      help="Select a mask file")
               
        form.addParallelSection(threads=1, mpi=0)
        
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', 
                                 self.getInputMicrographs().strId(),
                                 self.inputReferences.get().strId())
        self._insertPickingSteps()
        self._insertFunctionStep('createOutputStep')
        
    def _insertPickingSteps(self):
        """ Insert the steps to launch the pikcer with different modes. 
        Prepare the command arguments that will be passed. 
        """
        args =  ' --dirTgt=%s' % self._getTmpPath('micrographs')
        args += ' --dirSch=%s' % self._getTmpPath('templates')
        #args += ' --dirMskRef=%s' % self._getTmpPath('maskRef')
        args += ' --dirMskSch=%s' % self._getTmpPath('maskSch')
        args += ' --dirRes=%s' % self._getExtraPath('') # put the output in the extra dir
        args += ' --angle2D=%d' % self.rotAngle
        args += ' --contrast=%d' % (0 if self.refsHaveInvertedContrast else 1) 
        args += ' --thresh=%0.3f' % self.threshold
        args += ' --threshHigh=%0.2f' % self.thresholdHigh
        args += ' --nPickMax=%d' % self.maxPiks

        for mode in [0, 1]:
            self._insertFunctionStep('runGemPickerStep', mode, args)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, micsId, refsId):
        """ This step will take of the convertions from the inputs.
        Micrographs: they will be linked if are in '.mrc' format, converted otherwise.
        References: will always be converted to '.mrc' format
        Mask: either converted ('.tif' format) or generated a circular one
        """
        ih = em.ImageHandler()
        micDir = self._getTmpPath('micrographs')
        pwutils.makePath(micDir)
        
        for mic in self.getInputMicrographs():
            # Create micrograph folder
            micName = mic.getFileName()
            # If micrographs are in .mrc format just link it
            # otherwise convert them
            outMic = join(micDir, pwutils.replaceBaseExt(micName, 'mrc'))
            
            if micName.endswith('.mrc'):
                pwutils.createLink(micName, outMic)
            else:
                ih.convert(mic, outMic)
                
        refDir = self._getTmpPath('templates')
        pwutils.makePath(refDir)
        # We will always convert the templates, since
        # they can be in an stack and link will not be possible sometimes
        inputRefs = self.inputReferences.get()
        
        for i, ref in enumerate(inputRefs):
            outRef = join(refDir, 'ref%02d.mrc' % (i+1))
            ih.convert(ref, outRef)
            
        self.createInputMasks(inputRefs)  
            
    def createInputMasks(self, inputRefs):
        """ Create the needed mask for picking.
        We should either generate a circular mask, 
        or convert the inputs (one or just one per reference 2d)
        """
        maskSchDir = self._getTmpPath('maskSch')
        pwutils.makePath(maskSchDir)
        ih = em.ImageHandler()
        
        if self.maskType == MASK_CIRCULAR:
            if self.maskRadius < 0: # usually -1
                radius = inputRefs.getDim()[0]/2 # use half of input dim
            else:
                radius = self.maskRadius.get()
            outMask = join(maskSchDir, 'ref01.tif')
            ih.createCircularMask(radius, inputRefs.getFirstItem(), outMask)
        else:
            for i, mask in enumerate(self.inputMasks.get()):
                outMask = join(maskSchDir, 'ref%02d.tif' % (i+1))
                ih.convert(mask.get(), outMask)
                
    def runGemPickerStep(self, mode, args):
        #TODO: really launch the gEMPicker
        args += ' --mode=%d' % mode
        args += ' --nGPU=%d' % 0 # Use 0 for now
        args += ' --nCPU=%d' % self.numberOfThreads
        
        self.runJob('gEMpicker', args)

    def createOutputStep(self):
        return #FIXME
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        #TODO: parse the .box files and add coordinates per micrograph
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(micSet, coordSet)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d" % self.getInputMicrographs().getSize())
        if(self.getOutputsSize() > 0):
            summary.append("Number of particles picked: %d" % self.getCoords().getSize())
            summary.append("Particle size: %d" % self.getCoords().getBoxSize())
            summary.append("Threshold: %0.2f" % self.threshold)
            if self.extraParams.hasValue():
                summary.append("And other parameters: %s" % self.extraParams)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s of size %d." % (self.getObjectTag(self.getInputMicrographs()), self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(), output.getBoxSize(), self.threshold.get()))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs
    
    def _citations(self):
        return ['Hoang2013']
    
    #--------------------------- UTILS functions --------------------------------------------------
    