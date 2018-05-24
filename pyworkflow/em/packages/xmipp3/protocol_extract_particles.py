# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *              Javier Vargas (jvargas@cnb.csic.es)
# *              Grigory Sharov (sharov@igbmc.fr)
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

from glob import glob
from os.path import exists, basename

import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER
from pyworkflow.protocol.constants import (STEPS_PARALLEL, LEVEL_ADVANCED,
                                           STATUS_FINISHED)
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtExtractParticles
from pyworkflow.em.data import Particle
from pyworkflow.em.constants import RELATION_CTF

from convert import (micrographToCTFParam, writeMicCoordinates, xmippToLocation,
                     setXmippAttributes)
from xmipp3 import XmippProtocol


class XmippProtExtractParticles(ProtExtractParticles, XmippProtocol):
    """Protocol to extract particles from a set of coordinates"""
    _label = 'extract particles'
    
    def __init__(self, **kwargs):
        ProtExtractParticles.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _definePreprocessParams(self, form):
        # downFactor should always be 1.0 or greater
        geOne = params.GE(1.0,
                          error='Value should be greater or equal than 1.0')

        form.addParam('downFactor', params.FloatParam, default=1.0,
                      validators=[geOne],
                      label='Downsampling factor',
                      help='Select a value greater than 1.0 to reduce the size '
                           'of micrographs before extracting the particles. '
                           'If 1.0 is used, no downsample is applied. '
                           'Non-integer downsample factors are possible. ')

        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='This is size of the boxed particles (in pixels). '
                           'Note that if you use downsample option, the '
                           'particles are boxed out after downsampling. '
                           'Use the wizard to check boxSize changes after '
                           'downsampling or using a different pixel size. ')

        form.addParam('doBorders', params.BooleanParam, default=True,
                      label='Fill pixels outside borders',
                      help='Xmipp by default skips particles whose boxes fall '
                           'outside of the micrograph borders. Set this '
                           'option to True if you want those pixels outside '
                           'the borders to be filled with the closest pixel '
                           'value available')
        
        form.addSection(label='Preprocess')

        form.addParam('doRemoveDust', params.BooleanParam, default=True,
                      label='Dust removal (Recommended)', important=True,
                      help='Sets pixels with unusually large values to random '
                           'values from a Gaussian with zero-mean and '
                           'unity-standard deviation.')

        form.addParam('thresholdDust', params.FloatParam, default=5,
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this '
                           'value times the standard deviation of the image '
                           'will be affected. For cryo, 3.5 is a good value. '
                           'For high-contrast negative stain, the signal '
                           'itself may be affected so that a higher value may '
                           'be preferable.')

        form.addParam('doInvert', params.BooleanParam, default=None,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')
        
        form.addParam('doFlip', params.BooleanParam, default=None,
                      label='Phase flipping',
                      help='Use the information from the CTF to compensate for '
                           'phase reversals.\n'
                           'Phase flip is recommended in Xmipp or Eman\n'
                           '(even Wiener filtering and bandpass filter are '
                           'recommended for obtaining better 2D classes)\n'
                           'Otherwise (Frealign, Relion, Spider, ...), '
                           'phase flip is not recommended.')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and '
                           'normalizes so that in the background there is 0 '
                           'mean and standard deviation 1.')
        form.addParam('normType', params.EnumParam,
                      choices=['OldXmipp','NewXmipp','Ramp'], default=2,
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1). \n'
                           'NewXmipp (mean(background)=0, '
                           'stddev(background)=1) \n  '
                           'Ramp (subtract background+NewXmipp).')
        form.addParam('backRadius', params.IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius (px)', expertLevel=LEVEL_ADVANCED,
                      help='Pixels outside this circle are assumed to be noise '
                           'and their stddev is set to 1. Radius for '
                           'background circle definition (in pix.). If this '
                           'value is 0, then half the box size is used.')
        
        form.addParam('patchSize', params.IntParam, default=-1, 
                      label='Patch size for the variance filter (px)', 
                      expertLevel=LEVEL_ADVANCED,
                      help='Windows size to make the variance filtter and '
                           'compute the Gini coeff. A twice of the particle '
                           'size is recommended. Set at -1 applies 1.5*BoxSize.')

        form.addParallelSection(threads=4, mpi=1)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertInitialSteps(self):
        # Just overwrite this function to load some info
        # before the actual processing
        self._setupBasicProperties()

        return []

    def _getExtractArgs(self):
        """ Should be implemented in sub-classes to define the argument
        list that should be passed to the picking step function.
        """
        return [self.doInvert.get(),
                self._getNormalizeArgs(),
                self.doBorders.get()]
    
    #--------------------------- STEPS functions -------------------------------
    def _extractMicrograph(self, mic, doInvert, normalizeArgs, doBorders):
        """ Extract particles from one micrograph """
        fnLast = mic.getFileName()
        baseMicName = pwutils.removeBaseExt(fnLast)
        outputRoot = str(self._getExtraPath(baseMicName))
        fnPosFile = self._getMicPos(mic)
        boxSize = self.boxSize.get()
        downFactor = self.downFactor.get()
        patchSize = self.patchSize.get() if self.patchSize.get() > 0 \
                    else int(boxSize*1.5*downFactor)

        particlesMd = 'particles@%s' % fnPosFile
        # If it has coordinates extract the particles
        if exists(fnPosFile):
            # Create a list with micrographs operations (programs in xmipp) and
            # the required command line parameters (except input/ouput files)
            micOps = []

            # Compute the variance and Gini coeff. of the part. and mic., resp.
            args  =  '--pos %s' % fnPosFile
            args += ' --mic %s' % fnLast
            args += ' --patchSize %d' % patchSize
            self.runJob('xmipp_coordinates_noisy_zones_filter', args)

            def getMicTmp(suffix):
                return self._getTmpPath(baseMicName + suffix)

            # Check if it is required to downsample our micrographs
            if self.notOne(downFactor):
                fnDownsampled = getMicTmp("_downsampled.xmp")
                args = "-i %s -o %s --step %f --method fourier"
                self.runJob('xmipp_transform_downsample',
                            args % (fnLast, fnDownsampled, downFactor))
                fnLast = fnDownsampled

            if self.doRemoveDust:
                fnNoDust = getMicTmp("_noDust.xmp")
                args = " -i %s -o %s --bad_pixels outliers %f"
                self.runJob('xmipp_transform_filter',
                            args % (fnLast, fnNoDust, self.thresholdDust))
                fnLast = fnNoDust

            if self._useCTF():
                # We need to write a Xmipp ctfparam file
                # to perform the phase flip on the micrograph
                fnCTF = self._getTmpPath("%s.ctfParam" % baseMicName)
                micrographToCTFParam(mic, fnCTF)
                # Insert step to flip micrograph
                if self.doFlip:
                    fnFlipped = getMicTmp('_flipped.xmp')
                    args = " -i %s -o %s --ctf %s --sampling %f"
                    self.runJob('xmipp_ctf_phase_flip',
                                args % (fnLast, fnFlipped, fnCTF,
                                        self._getNewSampling()))
                    fnLast = fnFlipped
            else:
                fnCTF = None

            args = " -i %s --pos %s" % (fnLast, particlesMd)
            args += " -o %s --Xdim %d" % (outputRoot, boxSize)

            if doInvert:
                args += " --invert"

            if fnCTF:
                args += " --ctfparam " + fnCTF
            
            if doBorders:
                args += " --fillBorders"

            self.runJob("xmipp_micrograph_scissor", args)

            # Normalize
            if normalizeArgs:
                self.runJob('xmipp_transform_normalize',
                            '-i %s.stk %s' % (outputRoot, normalizeArgs))
        else:
            self.warning("The micrograph %s hasn't coordinate file! "
                         % baseMicName)
            self.warning("Maybe you picked over a subset of micrographs")

        # Let's clean the temporary mrc micrographs
        if not pwutils.envVarOn("SCIPION_DEBUG_NOCLEAN"):
            pwutils.cleanPattern(self._getTmpPath(baseMicName) + '*')

    def _getNormalizeArgs(self):
        if not self.doNormalize:
            return ''

        normType = self.getEnumText("normType")
        args = "--method %s " % normType

        if normType != "OldXmipp":
            bgRadius = self.backRadius.get()
            if bgRadius <= 0:
                bgRadius = int(self.boxSize.get() / 2)
            args += " --background circle %d" % bgRadius

        return args

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        # doFlip can only be selected if CTF information
        # is available on picked micrographs
        if self.doFlip and not self._useCTF():
            errors.append('Phase flipping cannot be performed unless '
                          'CTF information is provided.')

        if self.doNormalize:
            if self.backRadius > int(self.boxSize.get() / 2):
                errors.append("Background radius for normalization should be "
                              "equal or less than half of the box size.")
        
        # We cannot check this if the protocol is in streaming.
        
        #self._setupCtfProperties() # setup self.micKey among others
        # if self._useCTF() and self.micKey is None:
        #     errors.append('Some problem occurs matching micrographs and CTF.\n'
        #                   'There were micrographs for which CTF was not found '
        #                   'either using micName or micId.\n')

        # Clear the CTFs if micrograph source is "same as picking" to avoid unconsistencies
        if not self._micsOther():
            self.inputMicrographs.set(None)

        return errors
    
    def _citations(self):
        return ['Vargas2013b']
    
    def _summary(self):
        summary = []
        summary.append("Micrographs source: %s"
                       % self.getEnumText("downsampleType"))
        summary.append("Particle box size: %d" % self.boxSize)
        
        if not hasattr(self, 'outputParticles'):

            summary.append("Output images not ready yet.") 
        else:
            summary.append("Particles extracted: %d" %
                           self.outputParticles.getSize())
        
        return summary
    
    def _methods(self):
        methodsMsgs = []

        if self.getStatus() == STATUS_FINISHED:
            msg = ("A total of %d particles of size %d were extracted"
                   % (self.getOutput().getSize(), self.boxSize))
            
            if self._micsOther():
                msg += (" from another set of micrographs: %s"
                        % self.getObjectTag('inputMicrographs'))
            
            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)

            if self.doRemoveDust:
                methodsMsgs.append("Removed dust over a threshold of %s."
                                   % self.thresholdDust)
            if self.doInvert:
                methodsMsgs.append("Inverted contrast on images.")
            if self._doDownsample():
                methodsMsgs.append("Particles downsampled by a factor of %0.2f."
                                   % self.downFactor)
            if self.doNormalize:
                methodsMsgs.append("Normalization: %s."
                                   % self.getEnumText('normType'))

        return methodsMsgs

    # --------------------------- UTILS functions ------------------------------
    def _convertCoordinates(self, mic, coordList):
        writeMicCoordinates(mic, coordList, self._getMicPos(mic),
                            getPosFunc=self._getPos)
    
    def _micsOther(self):
        """ Return True if other micrographs are used for extract. """
        return self.downsampleType == OTHER

    def _useCTF(self):
        return self.ctfRelations.hasValue()

    def _doDownsample(self):
        return self.downFactor > 1.0

    def notOne(self, value):
        return abs(value - 1) > 0.0001

    def _getNewSampling(self):
        newSampling = self.samplingMics

        if self._doDownsample():
            # Set new sampling, it should be the input sampling of the used
            # micrographs multiplied by the downFactor
            newSampling *= self.downFactor.get()

        return newSampling

    def _setupBasicProperties(self):
        # Set sampling rate (before and after doDownsample) and inputMics
        # according to micsSource type
        inputCoords = self.getCoords()
        mics = inputCoords.getMicrographs()
        self.samplingInput = inputCoords.getMicrographs().getSamplingRate()
        self.samplingMics = self.getInputMicrographs().getSamplingRate()
        self.samplingFactor = float(self.samplingMics / self.samplingInput)

        scale = self.getBoxScale()
        self.debug("Scale: %f" % scale)
        if self.notOne(scale):
            # If we need to scale the box, then we need to scale the coordinates
            getPos = lambda coord: (int(coord.getX() * scale),
                                    int(coord.getY() * scale))
        else:
            getPos = lambda coord: coord.getPosition()
        # Store the function to be used for scaling coordinates
        self._getPos = getPos

    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates or
        Other micrographs. """
        if not self._micsOther():
            return self.inputCoordinates.get().getMicrographs()
        else:
            return self.inputMicrographs.get()

    def _storeMethodsInfo(self, fnImages):
        """ Store some information when the protocol finishes. """
        mdImgs = md.MetaData(fnImages)
        total = mdImgs.size()
        mdImgs.removeDisabled()
        zScoreMax = mdImgs.getValue(md.MDL_ZSCORE, mdImgs.lastObject())
        numEnabled = mdImgs.size()
        numRejected = total - numEnabled
        msg = ""

        if self.doFlip:
            msg += "\nPhase flipping was performed."

        self.methodsVar.set(msg)

    def getCoords(self):
        return self.inputCoordinates.get()

    def getOutput(self):
        if (self.hasAttribute('outputParticles') and
            self.outputParticles.hasValue()):
            return self.outputParticles
        else:
            return None

    def getCoordSampling(self):
        return self.getCoords().getMicrographs().getSamplingRate()

    def getMicSampling(self):
        return self.getInputMicrographs().getSamplingRate()

    def getBoxScale(self):
        """ Computing the sampling factor between input and output.
        We should take into account the differences in sampling rate between
        micrographs used for picking and the ones used for extraction.
        The downsampling factor could also affect the resulting scale.
        """
        samplingPicking = self.getCoordSampling()
        samplingExtract = self.getMicSampling()
        f = samplingPicking / samplingExtract
        return f / self.downFactor.get() if self._doDownsample() else f

    def getBoxSize(self):
        # This function is needed by the wizard
        return int(self.getCoords().getBoxSize() * self.getBoxScale())

    def _getOutputImgMd(self):
        return self._getPath('images.xmd')

    def createParticles(self, item, row):
        from pyworkflow.em.packages.xmipp3.convert import rowToParticle
        
        particle = rowToParticle(row, readCtf=self._useCTF())
        coord = particle.getCoordinate()
        item.setY(coord.getY())
        item.setX(coord.getX())
        particle.setCoordinate(item)
        self.imgSet.append(particle)
        item._appendItem = False

    def readPartsFromMics(self, micList, outputParts):
        """ Read the particles extract for the given list of micrographs
        and update the outputParts set with new items.
        """
        p = Particle()
        for mic in micList:
            # We need to make this dict because there is no ID in the .xmd file
            coordDict = {}
            for coord in self.coordDict[mic.getObjId()]:
                pos = self._getPos(coord)
                if pos in coordDict:
                    print("WARNING: Ignoring duplicated coordinate: %s, id=%s" %
                          (coord.getObjId(), pos))
                coordDict[pos] = coord

            added = set() # Keep track of added coords to avoid duplicates
            for row in md.iterRows(self._getMicXmd(mic)):
                pos = (row.getValue(md.MDL_XCOOR), row.getValue(md.MDL_YCOOR))
                coord = coordDict.get(pos, None)
                if coord is not None and coord.getObjId() not in added:
                    # scale the coordinates according to particles dimension.
                    coord.scale(self.getBoxScale())
                    p.copyObjId(coord)
                    p.setLocation(xmippToLocation(row.getValue(md.MDL_IMAGE)))
                    p.setCoordinate(coord)
                    p.setMicId(mic.getObjId())
                    p.setCTF(mic.getCTF())
                    # adding the variance and Gini coeff. value of the mic zone
                    setXmippAttributes(p, row, md.MDL_SCORE_BY_VAR)
                    setXmippAttributes(p, row, md.MDL_SCORE_BY_GINI)

                    # disabled particles (in metadata) should not add to the
                    # final set
                    if row.getValue(md.MDL_ENABLED) > 0:
                        outputParts.append(p)
                        added.add(coord.getObjId())

            # Release the list of coordinates for this micrograph since it
            # will not be longer needed
            del self.coordDict[mic.getObjId()]

    def _getMicPos(self, mic):
        """ Return the corresponding .pos file for a given micrograph. """
        micBase = pwutils.removeBaseExt(mic.getFileName())
        return self._getExtraPath(micBase + ".pos")

    def _getMicXmd(self, mic):
        """ Return the corresponding .xmd with extracted particles
        for this micrograph. """
        micBase = pwutils.removeBaseExt(mic.getFileName())
        return self._getExtraPath(micBase + ".xmd")