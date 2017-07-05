# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (sharov@igbmc.fr) [2]
# *
# *
# * [1] Science for Life Laboratory, Stockholm University
# * [2] L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
# *
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

import os
from os.path import join, exists

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from convert import readSetOfCoordinates, runGempicker, getProgram


MASK_CIRCULAR = 0
MASK_OBJECT = 1


class ProtGemPicker(em.ProtParticlePickingAuto):
    """
    gEMpicker is a template-based cryo-EM particle picking program that use
    cross-correlation approach. The user may define a template particle in 
    one of several ways, and this is then used to pick other similar particles
    from the micrographs by using a fast Fourier transform (FFT) to calculate 
    the correlation score at each pixel between the template and the micrograph.
    Multiple micrographs may be processed in parallel, and the calculation may 
    be accelerated considerably by using one or more attached graphics 
    processors (GPUs).
    """
    _label = 'auto-picking'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        
        em.ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('inputReferences', params.PointerParam,
                      pointerClass='SetOfAverages',
                      label='Input References', important=True,
                      help="Template images (2D class averages or reprojections "
                           "from a reference volume) to be used in picking.")
        form.addParam('refsHaveInvertedContrast', params.BooleanParam,
                      default=False,
                      label='References have inverted contrast',
                      help='Set to Yes to indicate that the reference have '
                           'inverted contrast with respect to the particles '
                           'in the micrographs.')
        form.addParam('rotAngle', params.IntParam, default=5,
                      label='Rotational angle search',
                      help='In-plane rotating angle in degrees '
                           '(0 = no rotation)')
        
        line = form.addLine('Threshold in the range [0-1]',
                            help="Threshold value for picking, select low and "
                                 "high values for cut-off.")
        line.addParam('thresholdLow', params.FloatParam, default=0.1,
                      label='Low')
        line.addParam('thresholdHigh', params.FloatParam, default=0.5,
                      label='High')
        
        form.addParam('maxPeaks', params.IntParam, default=0,
                      label='Max particles per micrograph',
                      expertLevel=LEVEL_ADVANCED,
                      help="Maximum number of particles picked from each "
                           "micrograph (0 = no limit)")
        form.addParam('boxSize', params.IntParam, default=0,
                      label='Box size (pix)', expertLevel=LEVEL_ADVANCED,
                      help='Size of picked images (if 0, use search image size)')
        form.addParam('boxDist', params.IntParam, default=0,
                      label='Min distance between particles (pix)',
                      expertLevel=LEVEL_ADVANCED,
                      help="Minimal distance between the centers of picked "
                           "particles (if 0, use reference image half-size)")
        form.addParam('boxBorder', params.IntParam, default=0,
                      label='Min distance from micrograph border (pix)',
                      expertLevel=LEVEL_ADVANCED,
                      help="Minimal distance between box edge and "
                           "micrograph border")
        form.addParam('maskType', params.EnumParam,
                      choices=['circular', 'object'], default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Mask type',
                      help='Select which type of mask do you want to apply. '
                           'Only the pixels beneath this mask will be analyzed. '
                           'In the simplest case, a circular mask can be used. '
                           'Alternatively, a custom mask can be used '
                           'which follows the contour of the particle '
                           '(but not too tightly).')
        form.addParam('maskRadius', params.IntParam, default=-1,
                      condition='maskType==%d' % MASK_CIRCULAR,
                      label='Mask radius (px)',
                      help='If -1, the entire image (in pixels) will be '
                           'considered.')
        form.addParam('inputMasks', params.MultiPointerParam,
                      pointerClass='Mask',
                      condition='maskType==%d' % MASK_OBJECT,
                      label="Mask objects",
                      help="Select a mask file")
        form.addParam('useGPU', params.BooleanParam, default=True,
                      label='Use GPU',
                      help='Set to Yes to use GPU as well as CPU')
        form.addParam('numberOfGPUs', params.IntParam, default=1,
                      label='GPUs per process', condition='useGPU',
                      help='Select number of GPUs per process')
        form.addParallelSection(threads=1, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertInitialSteps(self):
        convId = self._insertFunctionStep('convertInputStep',
                                          self.getInputMicrographs().strId(),
                                          self.inputReferences.get().strId())
        return [convId]
    
    #--------------------------- STEPS functions -------------------------------
    
    def convertInputStep(self, micsId, refsId):
        """ This step will take of the conversions from the inputs.
        Micrographs: they will be linked if are in '.mrc' format, converted
        otherwise.
        References: will always be converted to '.mrc' format
        Mask: either converted (to '.tif' format) or generated a circular one
        """
        
        self.convertInputs(self._getExtraPath())
    
    def _pickMicrograph(self, mic, *args):
        micName = mic.getFileName()
        runGempicker(micName, self._getExtraPath(), self.useGPU.get(), args[0],
                     log=self._log)
    
    def createOutputStep(self):
        pass
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        useGPU = self.useGPU.get()
        if not exists(getProgram(useGPU)):
            errors.append("Binary '%s' does not exits. \n"
                          "Check configuration file: "
                          "~/.config/scipion/scipion.conf\n"
                          "and set GEMPICKER variablesproperly."
                          % getProgram(useGPU))
            print "os.environ['GEMPICKER_HOME']", os.environ['GEMPICKER_HOME']
            print "os.environ['GEMPICKER']", os.environ['GEMPICKER']
        # Check that the number of input masks (in case of non-circular mask)
        # should be the same of the number of references, if greater than one
        if self.maskType == MASK_OBJECT:
            n = len(self.inputMasks)
            if n > 1 and n != self.inputReferences.get().getSize():
                errors.append('If the number of input masks is greater than '
                              'one, it should be equal to the number of '
                              'references.')
        
        value1 = round(self.thresholdLow,1)
        value2 = round(self.thresholdHigh,1)
        
        if (self.thresholdLow < self.thresholdHigh and
                        0.0 <= value1 <= 1.0 and 0.0 <= value2 <= 1.0):
            pass
        else:
            errors.append('Wrong threshold values!')
        
        return errors
    
    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())
        if(self.getOutputsSize() > 0):
            summary.append("Number of particles picked: %d"
                           % self.getCoords().getSize())
            summary.append("Particle size: %d px"
                           % self.getCoords().getBoxSize())
            summary.append("Threshold range: %0.3f - " % self.thresholdLow +
                           "%0.3f" % self.thresholdHigh)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary
    
    def _methods(self):
        methodsMsgs = []
        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']
        methodsMsgs.append("Input micrographs %s."
                           % (self.getObjectTag(self.getInputMicrographs())))
        
        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append("%s: User picked %d particles with a particle "
                               "size of %d px and threshold range %0.3f - %0.3f."
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.thresholdLow.get(),
                                  self.thresholdHigh.get()))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)
        
        return methodsMsgs
    
    def _citations(self):
        return ['Hoang2013']
    
    #--------------------------- UTILS functions -------------------------------
    def _getPickArgs(self, threshold=True, workingDir=None):
        """ Return the Gempicker parameters for picking one micrograph.
         The command line will depends on the protocol selected parameters.
        """
        nGPUs = self.numberOfGPUs.get() if self.useGPU else 0
        nThreads = self.numberOfThreads.get()
        
        args = ' --nGPU=%d' % nGPUs
        args += ' --nCPU=%d' % nThreads
        # Put the output in the extra dir by default
        args += ' --dirRes=%s' % (workingDir or self._getExtraPath(''))
        args += ' --angle2D=%d' % self.rotAngle
        args += ' --contrast=%d' % (0 if self.refsHaveInvertedContrast else 1)
        if threshold:
            args += ' --thresh=%0.3f' % self.thresholdLow
            args += ' --threshHigh=%0.3f' % self.thresholdHigh
        args += ' --nPickMax=%d' % self.maxPeaks
        args += ' --boxSize=%d' % self.boxSize
        args += ' --boxDist=%d' % self.boxDist
        args += ' --boxBorder=%d' % self.boxBorder
        
        return [args]
    
    def convertInputs(self, workingDir):
        """ This step will take of the conversions from the inputs.
        Micrographs: they will be linked if are in '.mrc' format, converted otherwise.
        References: will always be converted to '.mrc' format
        Mask: either converted (to '.tif' format) or generated a circular one
        """
        # Create micrographs dir
        def makePath(d):
            p = join(workingDir, d)
            pwutils.cleanPath(p)
            pwutils.makePath(p)
            return p
        
        makePath('micrographs')
        
        refDir = makePath('templates')
        inputRefs = self.inputReferences.get()
        for i, ref in enumerate(inputRefs):
            outRef = join(refDir, 'ref%02d.mrc' % (i + 1))
            em.ImageHandler().convert(ref, outRef)
        
        maskSchDir = makePath('maskSch')
        self.createInputMasks(inputRefs, maskSchDir)
    
    def createInputMasks(self, inputRefs, maskSchDir):
        """ Create the needed mask for picking.
        We should either generate a circular mask,
        or convert the inputs (one or just one per reference 2d)
        """
        ih = em.ImageHandler()
        
        if self.maskType == MASK_CIRCULAR:
            if self.maskRadius < 0:  # usually -1
                radius = inputRefs.getDim()[0] / 2  # use half of input dim
            else:
                radius = self.maskRadius.get()
            outMask = join(maskSchDir, 'ref01.tif')
            ih.createCircularMask(radius, inputRefs.getFirstItem(), outMask)
        else:
            for i, mask in enumerate(self.inputMasks):
                outMask = join(maskSchDir, 'ref%02d.tif' % (i + 1))
                ih.convert(mask.get(), outMask)

    def readCoordsFromMics(self, workingDir, micDoneList,
                           outputCoords):
        if self.boxSize and self.boxSize > 0:
            outputCoords.setBoxSize(self.boxSize.get())
        else:
            outputCoords.setBoxSize(
                self.inputReferences.get().getDim()[0])

        readSetOfCoordinates(workingDir, micDoneList, outputCoords)

    def getCoordsDir(self):
        return self._getExtraPath()