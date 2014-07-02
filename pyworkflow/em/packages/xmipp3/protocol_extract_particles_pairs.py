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
"""
This sub-package contains the XmippProtExtractParticlesPairs protocol
"""


from pyworkflow.em.packages.xmipp3.protocol_extract_particles import XmippProtExtractParticles, ORIGINAL, SAME_AS_PICKING, OTHER 
from pyworkflow.em import CoordinatesTiltPair
from convert import writeSetOfCoordinates, readSetOfParticles
from pyworkflow.utils.path import removeBaseExt, exists
from pyworkflow.protocol.params import PointerParam
from itertools import izip
from glob import glob
import xmipp
from pyworkflow.em.data_tiltpairs import ParticlesTiltPair
               
class XmippProtExtractParticlesPairs(XmippProtExtractParticles):
    """Protocol to extract particles from a set of tilted pairs coordinates"""
    _label = 'extract particles pairs'
    
    def __init__(self, **args):
        XmippProtExtractParticles.__init__(self, **args)
        
    #--------------------------- DEFINE param functions --------------------------------------------   

    def _defineParamsInput1(self, form):
        form.addParam('inputCoordinatesTiltedPairs', PointerParam, label="Coordinates tilted pairs", 
                      pointerClass='CoordinatesTiltPair',
                      help='Select the CoordinatesTiltPairs ')
        
    def _defineParamsInput2(self, form):
        form.addParam('inputMicrographsTiltedPair', PointerParam, label="Micrographs tilted pair", 
                  condition='downsampleType != 1',
                  pointerClass='MicrographsTiltPair',
                  help='Select the original Micrographs tilted pair')

    def _setInputMicrographs(self):
        # Set sampling rate and inputMics according to downsample type
        
        self.samplingInput = self.inputCoordinatesTiltedPairs.get().getUntilted().getMicrographs().getSamplingRate()
        
        if self.downsampleType.get() == SAME_AS_PICKING:
            # If 'same as picking' get samplingRate from input micrographs  (both tilted and untilted)
            self.uMics = self.inputCoordinatesTiltedPairs.get().getUntilted().getMicrographs()
            self.tMics = self.inputCoordinatesTiltedPairs.get().getTilted().getMicrographs()
            
            self.samplingFinal = self.samplingInput
        else:
            self.uMics = self.inputMicrographsTiltedPair.get().getUntilted()
            self.tMics = self.inputMicrographsTiltedPair.get().getTilted()
                                            
            self.samplingOriginal = self.uMics.getSamplingRate()
            if self.downsampleType.get() == ORIGINAL:
                # If 'original' get sampling rate from original micrographs
                self.samplingFinal = self.samplingOriginal
            else:
                # IF 'other' multiply the original sampling rate by the factor provided
                self.samplingFinal = self.samplingOriginal*self.downFactor.get()

        self.inputMics = self._createSetOfParticles('auxMics')
        self.inputMics.copyInfo(self.uMics)
        self.inputMics.setStore(False)
        
#        self.inputMics.appendFromImages(self.tMics)
#        self.inputMics.appendFromImages(self.uMics)
        for micU, micT in izip(self.uMics, self.tMics):
            micU.cleanObjId()
            micT.cleanObjId()
            self.inputMics.append(micU)
            self.inputMics.append(micT)
        # Tilted pairs cannot be flipped
        self.doFlip.set(False)
        
        # Tilted pairs do not have ctfRelation
        self.ctfRelations.set(None)
                
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
            imgSetU.setSamplingRate(self.uMics.getSamplingRate()*self.downFactor.get())
            imgSetT.setSamplingRate(self.tMics.getSamplingRate()*self.downFactor.get())
        imgSetU.setCoordinates(self.inputCoordinatesTiltedPairs.get().getUntilted())
        imgSetT.setCoordinates(self.inputCoordinatesTiltedPairs.get().getTilted())
        
        imgSetAuxU = self._createSetOfParticles('auxU')
        imgSetAuxU.copyInfo(imgSetU)
        readSetOfParticles(fnUntilted, imgSetAuxU, False)
        imgSetAuxU.write()
        # For each untilted particle retrieve micId from SetOFCoordinates untilted
        for img, coord in izip(imgSetAuxU, self.inputCoordinatesTiltedPairs.get().getUntilted()):
            #FIXME: REmove this check when sure that objIds are equal
            if img.getObjId() != coord.getObjId(): 
                raise Exception('ObjId is not equal!!!!')
            img.setCoordinate(coord)
            #img.cleanObjId()
            imgSetU.append(img)

        imgSetAuxT = self._createSetOfParticles('auxT')
        imgSetAuxT.copyInfo(imgSetT)
        readSetOfParticles(fnTilted, imgSetAuxT, False)    
        imgSetAuxT.write()        
        # For each untilted particle retrieve micId from SetOFCoordinates tilted
        #for img in imgSetAuxU:
        for img, coord in izip(imgSetAuxT, self.inputCoordinatesTiltedPairs.get().getTilted()):
            #FIXME: This can be slow to make a query to grab the coord, maybe use zip(imgSet, coordSet)???
            #FIXME: REmove this check when sure that objIds are equal
            if img.getObjId() != coord.getObjId(): 
                raise Exception('ObjId is not equal!!!!')
            #coord = self.inputCoordinatesTiltedPairs.get().getTilted()[img.getObjId()]
            img.setCoordinate(coord)
            #img.cleanObjId()
            imgSetT.append(img)
            
        imgSetU.write()
        imgSetT.write()
        
        self._storeMethodsInfo(fnUntilted)
        
        # Define output ParticlesTiltPair 
        outputset = ParticlesTiltPair()
        outputset.setTilted(imgSetT)
        outputset.setUntilted(imgSetU)
        outputset.setCoordsPair(self.inputCoordinatesTiltedPairs)
        self._defineOutputs(outputParticlesTiltPair=outputset)
        self._defineSourceRelation(self.inputCoordinatesTiltedPairs, outputset)
            
    #--------------------------- INFO functions -------------------------------------------- 
    def _citations(self):
        return ['Vargas2013b']
      
    #TODO: Refactor method below    
    def _summary(self):
        downsampleTypeText = {
                              ORIGINAL:'Original micrographs',
                              SAME_AS_PICKING:'Same as picking',
                              OTHER: 'Other downsampling factor'}
        summary = []
        summary.append("_Downsample type_: %s" % downsampleTypeText.get(self.downsampleType.get()))
        if self.downsampleType == OTHER:
            summary.append("Downsampling factor: %d" % self.downFactor.get())
        summary.append("Particle box size: %d" % self.boxSize.get())
        
        if not hasattr(self, 'outputParticlesTiltPair'):
            summary.append("Output images not ready yet.") 
        else:
            summary.append("Particles extracted: %d" % (self.outputParticlesTiltPair.getTilted().getSize()))
            
        return summary
    
    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("Particle box size %d" % self.boxSize.get())

        if self.methodsInfo.hasValue():
            methodsMsgs.append(self.methodsInfo.get())
        
        methodsMsgs.append("Automatic Rejection method selected: %s" % (self.rejectionMethod))    

        methodsMsgs.append("Invert contrast performed?: %s" % (self.doInvert.get()))
        methodsMsgs.append("Normalize performed?: %s" % (self.doNormalize.get()))
        if self.doNormalize.get():
            methodsMsgs.append("Nomalization used: %s" % (self.getEnumText('normType')))
            methodsMsgs.append("Nomalization used: %s" % (self.backRadius.get()))
        methodsMsgs.append("Remove dust?: %s" % (self.doRemoveDust.get()))
        if self.doRemoveDust.get():
            methodsMsgs.append("Dust threshold: %s" % (self.thresholdDust.get()))            

        return methodsMsgs

    
