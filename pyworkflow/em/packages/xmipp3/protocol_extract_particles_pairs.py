# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *              Javier Vargas (jvargas@cnb.csic.es)
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

from pyworkflow.em.packages.xmipp3.protocol_extract_particles import XmippProtExtractParticles
from pyworkflow.em.packages.xmipp3.constants import SAME_AS_PICKING, OTHER 
from pyworkflow.em import PointerParam, IntParam, EnumParam, BooleanParam, FloatParam
from pyworkflow.protocol.params import Positive
from convert import writeSetOfCoordinates, readSetOfParticles
from pyworkflow.utils.path import removeBaseExt, exists
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from itertools import izip
import xmipp
from pyworkflow.em.data_tiltpairs import ParticlesTiltPair, TiltPair
               
                            
class XmippProtExtractParticlesPairs(XmippProtExtractParticles):
    """Protocol to extract particles from a set of tilted pairs coordinates"""
    _label = 'extract particle pairs'
    
    def __init__(self, **args):
        XmippProtExtractParticles.__init__(self, **args)
        
    #--------------------------- DEFINE param functions --------------------------------------------   

    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputCoordinatesTiltedPairs', PointerParam, label="Coordinates tilted pairs", 
                      pointerClass='CoordinatesTiltPair',
                      help='Select the CoordinatesTiltPairs ')
        
        form.addParam('downsampleType', EnumParam, choices=['same as picking', 'other', 'original'], 
                      default=0, important=True, label='Downsampling type', display=EnumParam.DISPLAY_COMBO, 
                      help='Select the downsampling type.')
        
        form.addParam('downFactor', FloatParam, default=2, condition='downsampleType==1',
                      label='Downsampling factor',
                      help='This factor is always referred to the original sampling rate. '
                      'You may use independent downsampling factors for extracting the '
                      'particles, picking them and estimating the CTF. All downsampling '
                      'factors are always referred to the original sampling rate, and '
                      'the differences are correctly handled by Xmipp.')      
        
        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed particles, ' 
                      'actual particles may be smaller than this.')

        form.addSection(label='Preprocess')
        form.addParam('doRemoveDust', BooleanParam, default=True, important=True,
                      label='Dust removal (Recommended)', 
                      help='Sets pixels with unusually large values to random values from a Gaussian '
                      'with zero-mean and unity-standard deviation.')
        form.addParam('thresholdDust', FloatParam, default=3.5, 
                      condition='doRemoveDust', expertLevel=LEVEL_ADVANCED,
                      label='Threshold for dust removal',
                      help='Pixels with a signal higher or lower than this value times the standard '
                      'deviation of the image will be affected. For cryo, 3.5 is a good value.'
                      'For high-contrast negative stain, the signal itself may be affected so '
                      'that a higher value may be preferable.')
        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')
        
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize (Recommended)', 
                      help='It subtract a ramp in the gray values and normalizes so that in the '
                      'background there is 0 mean and standard deviation 1.')
        form.addParam('normType', EnumParam, choices=['OldXmipp','NewXmipp','Ramp'], default=2, 
                      condition='doNormalize', expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_COMBO,
                      label='Normalization type', 
                      help='OldXmipp (mean(Image)=0, stddev(Image)=1).  \n  '
                           'NewXmipp (mean(background)=0, stddev(background)=1)  \n  '
                           'Ramp (subtract background+NewXmipp).  \n  ')
        form.addParam('backRadius', IntParam, default=-1, condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.). '
                      'If this value is -1, then half the box size is used.', 
                      expertLevel=LEVEL_ADVANCED)
        
        form.addParallelSection(threads=4, mpi=1)


    def _insertAllSteps(self):
        """for each micrograph insert the steps to preprocess it
        """       
        self.uMics = self.inputCoordinatesTiltedPairs.get().getUntilted().getMicrographs()
        self.tMics = self.inputCoordinatesTiltedPairs.get().getTilted().getMicrographs()

        self.inputMics = self._createSetOfParticles('auxMics')
        self.inputMics.copyInfo(self.uMics)
        self.inputMics.setStore(False)
        
        for micU, micT in izip(self.uMics, self.tMics):
            micU.cleanObjId()
            micT.cleanObjId()
            self.inputMics.append(micU)
            self.inputMics.append(micT)

        self.samplingInput = self.uMics.getSamplingRate()
        

        if self.downsampleType.get() != OTHER:
            # If 'same as picking' or 'original' get sampling rate from input micrographs
            #TODO: Review this when downsampling before picking is possible
            self.samplingFinal = self.samplingInput
        else:
            # If 'other' multiply the input sampling rate by the factor provided
            self.samplingFinal = self.samplingInput*self.downFactor.get()
                        
        # Write pos files for each micrograph
        firstStepId = self._insertFunctionStep('writePosFilesStep')
           
        # For each micrograph insert the steps
        #run in parallel
        
        deps = []
        for mic in self.inputMics:
            localDeps = [firstStepId]
            micrographToExtract = mic.getFileName()
            micName = removeBaseExt(mic.getFileName())
            micId = mic.getObjId()

            # If downsample type is 'other' perform a downsample
            if self.downsampleType == OTHER:
                fnDownsampled = self._getTmpPath(micName+"_downsampled.xmp")
                downFactor = self.downFactor.get()
                args = "-i %(micrographToExtract)s -o %(fnDownsampled)s --step %(downFactor)f --method fourier"
                localDeps=[self._insertRunJobStep("xmipp_transform_downsample", args % locals(),prerequisites=localDeps)]
                micrographToExtract = fnDownsampled
                                                            
            # If remove dust 
            if self.doRemoveDust:
                fnNoDust = self._getTmpPath(micName+"_noDust.xmp")
                
                thresholdDust = self.thresholdDust.get() #TODO: remove this extra variable
                args=" -i %(micrographToExtract)s -o %(fnNoDust)s --bad_pixels outliers %(thresholdDust)f"
                localDeps=[self._insertRunJobStep("xmipp_transform_filter", args % locals(),prerequisites=localDeps)]
                micrographToExtract = fnNoDust
                                        
            #self._insertFunctionStep('getCTF', micId, micName, micrographToExtract)
            micName = removeBaseExt(mic.getFileName())
      
            # Actually extract
            deps.append(self._insertFunctionStep('extractParticlesStep', micId, micName, 
                                              None, micrographToExtract, prerequisites=localDeps))
        # TODO: Delete temporary files
                        
        # Insert step to create output objects      
        self._insertFunctionStep('createOutputStep', prerequisites=deps)

                
    #--------------------------- STEPS functions --------------------------------------------
    def writePosFilesStep(self):
        """ Write the pos file for each micrograph on metadata format (both untilted and tilted). """      
        
        writeSetOfCoordinates(self._getExtraPath(), self.inputCoordinatesTiltedPairs.get().getUntilted())
           
        writeSetOfCoordinates(self._getExtraPath(), self.inputCoordinatesTiltedPairs.get().getTilted())
                         
                
    def createOutputStep(self):
        # Create the SetOfImages objects on the database and the ImagesTiltPair

        mdUntilted = xmipp.MetaData()
        mdTilted = xmipp.MetaData()
        #for objId in mdPairs:
        for uMic, tMic in izip(self.uMics, self.tMics):
            umicName = removeBaseExt(uMic.getFileName())
            fnMicU = self._getExtraPath(umicName + ".xmd")
            fnPosU = self._getExtraPath(umicName + ".pos")
            # Check if there are picked particles in this micrographs
            if exists(fnMicU):
                mdMicU = xmipp.MetaData(fnMicU)
                mdPosU = xmipp.MetaData('particles@%s' % fnPosU)
                mdPosU.merge(mdMicU)
                mdUntilted.unionAll(mdPosU)
                tmicName = removeBaseExt(tMic.getFileName())
                fnMicT = self._getExtraPath(tmicName + ".xmd")
                fnPosT = self._getExtraPath(tmicName + ".pos")
                mdMicT = xmipp.MetaData(fnMicT)
                mdPosT = xmipp.MetaData('particles@%s' % fnPosT)
                mdPosT.merge(mdMicT)
                mdTilted.unionAll(mdPosT)

        # Write image metadatas (check if it is really necessary)
        fnTilted = self._getExtraPath("images_tilted.xmd")
        fnUntilted = self._getExtraPath("images_untilted.xmd")
        mdUntilted.write(fnUntilted)
        mdTilted.write(fnTilted)

        # Create outputs SetOfParticles both for tilted and untilted
        imgSetU = self._createSetOfParticles(suffix="Untilted")
        imgSetU.copyInfo(self.uMics)

        imgSetT = self._createSetOfParticles(suffix="Tilted")
        imgSetT.copyInfo(self.tMics)

        if self.downsampleType == OTHER:
            imgSetU.setSamplingRate(self.samplingFinal)
            imgSetT.setSamplingRate(self.samplingFinal)

        imgSetU.setCoordinates(self.inputCoordinatesTiltedPairs.get().getUntilted())
        imgSetT.setCoordinates(self.inputCoordinatesTiltedPairs.get().getTilted())

        #Read untilted and tilted particles on a temporary object (also disabled particles)
        imgSetAuxU = self._createSetOfParticles('auxU')
        imgSetAuxU.copyInfo(imgSetU)
        readSetOfParticles(fnUntilted, imgSetAuxU, removeDisabled=False)
        imgSetAuxU.write()

        imgSetAuxT = self._createSetOfParticles('auxT')
        imgSetAuxT.copyInfo(imgSetT)
        readSetOfParticles(fnTilted, imgSetAuxT, removeDisabled=False)
        imgSetAuxT.write()

        coordsT = self.inputCoordinatesTiltedPairs.get().getTilted()
        # For each untilted particle retrieve micId from SetOFCoordinates untilted
        for imgU, coordU in izip(imgSetAuxU, self.inputCoordinatesTiltedPairs.get().getUntilted()):
            #FIXME: REmove this check when sure that objIds are equal
            id = imgU.getObjId()
            if id != coordU.getObjId():
                raise Exception('ObjId in untilted is not equal!!!!')

            imgT = imgSetAuxT[id]
            coordT = coordsT[id]

            #If both particles are enabled append them
            if imgU.isEnabled() and imgT.isEnabled():
                imgU.setCoordinate(coordU)
                imgSetU.append(imgU)
                imgT.setCoordinate(coordT)
                imgSetT.append(imgT)

        # For each untilted particle retrieve micId from SetOFCoordinates tilted
        #for img in imgSetAuxU:
        # for img, coord in izip(imgSetAuxT, self.inputCoordinatesTiltedPairs.get().getTilted()):
        #     #FIXME: This can be slow to make a query to grab the coord, maybe use zip(imgSet, coordSet)???
        #     #FIXME: REmove this check when sure that objIds are equal
        #     if img.getObjId() != coord.getObjId():
        #         raise Exception('ObjId is not equal!!!!')
        #     #coord = self.inputCoordinatesTiltedPairs.get().getTilted()[img.getObjId()]
        #     img.setCoordinate(coord)
        #     #img.cleanObjId()
        #     imgSetT.append(img)

        imgSetU.write()
        imgSetT.write()

        self._storeMethodsInfo(fnUntilted)

        # Define output ParticlesTiltPair 
        outputset = ParticlesTiltPair(filename=self._getPath('particles_pairs.sqlite'))
        outputset.setTilted(imgSetT)
        outputset.setUntilted(imgSetU)
        for imgU, imgT in izip(imgSetU, imgSetT):
            outputset.append(TiltPair(imgU, imgT))

        outputset.setCoordsPair(self.inputCoordinatesTiltedPairs.get())
        self._defineOutputs(outputParticlesTiltPair=outputset)
        self._defineSourceRelation(self.inputCoordinatesTiltedPairs, outputset)
            
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
    def _citations(self):
        return ['Vargas2013b']
      
    #TODO: Refactor method below    
    def _summary(self):
        downsampleTypeText = {
                              SAME_AS_PICKING:'Same as picking',
                              OTHER: 'Other downsampling factor'}
        summary = []
        summary.append("Downsample type: %s" % downsampleTypeText.get(self.downsampleType.get()))
        if self.downsampleType == OTHER:
            summary.append("Downsampling factor: %d" % self.downFactor.get())
            
        summary.append("Particle box size: %d" % self.boxSize.get())
        
        if not hasattr(self, 'outputParticlesTiltPair'):
            summary.append("Output images not ready yet.") 
        else:
            summary.append("Particles extracted: %d" % (self.outputParticlesTiltPair.getTilted().getSize()))
            
        return summary

    def getCoords(self):
        if self.inputCoordinatesTiltedPairs.hasValue():
            return self.inputCoordinatesTiltedPairs.get().getUntilted()
        else:
            return None

    def _storeMethodsInfo(self, fnImages):
        """ Store some information when the protocol finishes. """
        self.methodsVar.set("")

    def getInputMicrographs(self):
        """ Return the micrographs associated to the SetOfCoordinates"""
        #return self.inputCoordinatesTiltedPairs.get().getMicsPair()
        if self.inputCoordinatesTiltedPairs.hasValue():
            return self.inputCoordinatesTiltedPairs.get().getUntilted().getMicrographs()
        else:
            return None
    
    def getOutput(self):
        if self.outputParticlesTiltPair.hasValue():
            return self.outputParticlesTiltPair
        else:
            return None
