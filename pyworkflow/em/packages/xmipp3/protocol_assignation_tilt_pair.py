# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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

from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.utils.path import moveFile, replaceBaseExt

import xmipp
from convert import readSetOfCoordinates, izip


class XmippProtAssignationTiltPair(ProtParticlePicking):
    """    
    From two sets of points (tilted and untilted) the protocol determines the affine transformation between
    these sets.
    """
    _label = 'assignation tiltpair'
    
    def __init__(self, *args, **kwargs):
        ProtParticlePicking.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('untiltcoor', PointerParam, pointerClass='SetOfParticles',
                      label="Untilt coordinates",  
                      help='Select the metadata with untilt coordinates.')     

        form.addParam('tiltcoor', PointerParam, pointerClass='SetOfParticles', 
                      label="Tilt coordinates",  
                      help='Select the metadata with untilt coordinates.')
        
        form.addParam('tiltmic', PointerParam, pointerClass='SetOfMicrographs', 
                      label="Tilt micrography",  
                      help='Select a tilt micrography for getting its size.')  
            
        form.addParam('maxshift', IntParam, default='1000', expertLevel=LEVEL_ADVANCED,
                      label="Maximum shift (pixels)", 
                      help='Maximum allowed distance (in pixels) that the tilt micrograph can be shifted' 
                      'respect to the untilted micrograph') 
        
        form.addParam('particlesize', IntParam, default=100, 
                      label="Particle size (pixels)",  
                      help='It defines the size os the box which contains the particle')

        form.addParam('threshold', FloatParam, default=0.25, expertLevel=LEVEL_ADVANCED,
                      label="Threshold value",  
                      help='Parameter between 0 and 1 that allows to define if \n' 
                      'a tilt point can be matched with a certain untilt point. \n'
                      'The matching is performed only if the distance is lesser than \n'
                      'threshold * particlesize.')
                
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):        
        '''tiltMic = self.tiltcoor.get().getMicrographs()
        untiltMic = self.untiltcoor.get().getMicrographs()
        
        self._insertFunctionStep('convertInputStep')
        deps = [] # store volumes steps id to use as dependencies for last step
        for untiltMic, tiltMic in izip(untiltMic, tiltMic):
            self._insertFunctionStep('assignationStep', untiltMic.getFileName(), tiltMic.getFileName())
        
        self._insertFunctionStep('createOutputStep', prerequisites=deps)'''
        
        tiltMic = self.tiltmic.get()
         
        self._insertFunctionStep('convertInputStep')
        deps = [] # store volumes steps id to use as dependencies for last step
        for tiltMic in izip(tiltMic):
            self._insertFunctionStep('assignationStep', tiltMic.getFileName())
        
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
        
        
    def convertInputStep(self):
        """ Read the input metadatata.
        """
        inputCoordsUntilt = self.untiltcoor.get()
        inputCoordsTilt = self.tiltcoor.get()
        
        untiltMic = inputCoordsUntilt.getMicrographs()
        tiltMic = inputCoordsTilt.getMicrographs()

        
        readSetOfCoordinates(self._getExtraPath('tilt'), tiltMic, inputCoordsTilt)
        readSetOfCoordinates(self._getExtraPath('untilt'), untiltMic, inputCoordsUntilt)

    
    def assignationStep(self, fnuntilt, fntilt):

        tiltMic = self.tiltmic.get()
        params =  ' --untiltcoor %s' % join(self._getExtraPath('untilt'), replaceBaseExt(fnuntilt, 'pos'))        
        params += ' --tiltcoor %s' % join(self._getExtraPath('tilt'), replaceBaseExt(fntilt, 'pos'))
        params += ' --tiltmic %s' % tiltMic.getFileName()
        params += ' --maxshift %f' % self.maxshift.get()
        params += ' --particlesize %f' % self.particlesize.get()
        params += ' --threshold %f' % self.threshold.get()


        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        self.runJob('xmipp_assignation_tilt_pair', params, numberOfMpi=nproc,numberOfThreads=nT)
        
          
    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        
        for i, vol in enumerate(self._iterInputVols()):
            volume = vol.clone()
            volDir = self._getVolDir(i+1)
            volPrefix = 'vol%03d_' % (i+1)
            validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
            moveFile(join(volDir, 'validation.xmd'), 
                     validationMd)
            clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
            moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
            
            md = xmipp.MetaData(validationMd)
            weight = md.getValue(xmipp.MDL_WEIGHT, md.firstObject())
            volume.weight = Float(weight)
            volume.clusterMd = String(clusterMd)
            volume.cleanObjId() # clean objects id to assign new ones inside the set
            outputVols.append(volume)                
        
        outputVols.setSamplingRate(volume.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        self._defineTransformRelation(self.inputVolumes, outputVols)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        #if self.inputVolumes.get() and not self.inputVolumes.hasValue():
        #    validateMsgs.append('Please provide an input reference volume.')
        if self.untiltcoor.get() and not self.untiltcoor.hasValue():
            validateMsgs.append('Please provide input particles.')  
        if self.tiltcoor.get() and not self.tiltcoor.hasValue():
            validateMsgs.append('Please provide input particles.')         
        return validateMsgs
        
    
    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            size = 0
            for i, vol in enumerate(self._iterInputVols()):
                size +=1
            summary.append("Volumes to validate: *%d* " % size)
            summary.append("Angular sampling: %s" % self.angularSampling.get())
            summary.append("Significance value: %s" % self.alpha.get())

        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and significant value of %f' % (self.angularSampling.get(), self.alpha.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _defineMetadataRootName(self, mdrootname,volId):
        
        if mdrootname=='P':
            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'clusteringTendency.xmd')
        if mdrootname=='Volume':

            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'validation.xmd')
            
    def _definePName(self):
        fscFn = self._defineMetadataRootName('P')
        return fscFn
    
    def _defineVolumeName(self,volId):
        fscFn = self._defineMetadataRootName('Volume',volId)
        return fscFn

