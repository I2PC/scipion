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

from os.path import split, splitext

from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.utils.path import makePath
from pyworkflow.em.data_tiltpairs import TiltPair, CoordinatesTiltPair
from pyworkflow.em import ProtParticlePicking
from pyworkflow.em.packages.xmipp3 import XmippProtocol

import xmipp
from convert import readSetOfCoordinates, writeSetOfCoordinates, izip




class XmippProtAssignmentTiltPair(ProtParticlePicking, XmippProtocol):
    """    
    From two sets of points (tilted and untilted) the protocol determines the affine transformation between
    these sets.
    """
    _label = 'assignment tiltpair'
    
    def __init__(self, *args, **kwargs):
        ProtParticlePicking.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('tiltpair', PointerParam, pointerClass='MicrographsTiltPair',
                      label="Micrograph tilt pair",  
                      help='Select micrographs tilt pair.')
        
        form.addParam('untiltcoor', PointerParam, pointerClass='SetOfCoordinates', 
                      label="Untilt coordinates",  
                      help='Select the metadata with untilt coordinates.')
        
        form.addParam('tiltcoor', PointerParam, pointerClass='SetOfCoordinates', 
                      label="Tilt coordinates",  
                      help='Select the metadata with tilt coordinates.')     
        
        form.addParam('threshold', FloatParam, default=0.25, expertLevel=LEVEL_ADVANCED,
                      label="Threshold value",  
                      help='Parameter between 0 and 1 that allows to define if \n' 
                      'a tilt point can be matched with a certain untilt point. \n'
                      'The matching is performed only if the distance is lesser than \n'
                      'threshold * particlesize.')

        form.addParam('maxshift', FloatParam, default='1000', expertLevel=LEVEL_ADVANCED,
                      label="Maximum shift (pixels)", 
                      help='Maximum allowed distance (in pixels) that the tilt micrograph can be shifted' 
                      'respect to the untilted micrograph')         

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):        
        self.micsFn = self._getPath('input_micrographs.xmd')
        # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')
        deps = []
        
        for tiltPair in self.tiltpair.get():
            micUntilted = tiltPair.getUntilted()
            micTilted = tiltPair.getTilted()
            Unpath, Unname = split(micUntilted.getFileName())
            untiltpath, unna = split(self.untiltcoor.get().getFileName())
            tiltpath, tna = split(self.tiltcoor.get().getFileName())
            Unname, ext = splitext(Unname)
            Tpath, Tname = split(micTilted.getFileName())
            Tname, ext = splitext(Tname)
            fnUntilt = 'particles@'+self._getExtraPath("untilted/")+Unname+'.pos'
            print fnUntilt
            fnTilt = 'particles@'+self._getExtraPath("tilted/")+Tname+'.pos'
            fnmicsize = tiltPair.getTilted().getFileName()
            stepId=self._insertFunctionStep('assignmentStep',fnUntilt,fnTilt, fnmicsize, self._getExtraPath())
            deps.append(stepId)
            
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
        

    def convertInputStep(self):
        """ Read the input metadatata.
        """
        # Get the converted input micrographs in Xmipp format
        makePath(self._getExtraPath("untilted"))
        makePath(self._getExtraPath("tilted"))
        writeSetOfCoordinates(self._getExtraPath("untilted"),self.untiltcoor.get())
        writeSetOfCoordinates(self._getExtraPath("tilted"),self.tiltcoor.get())
    
    
    def assignmentStep(self,fnuntilt, fntilt, fnmicsize, Unpath):

        params =  ' --untiltcoor %s' % fnuntilt        
        params += ' --tiltcoor %s' % fntilt
        params += ' --tiltmicsize %s' % fnmicsize
        params += ' --maxshift %f' % self.maxshift.get()
        params += ' --particlesize %f' % self.untiltcoor.get().getBoxSize()
        params += ' --threshold %f' % self.threshold.get()
        params += ' --odir %s' % Unpath


        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        self.runJob('xmipp_image_assignment_tilt_pair', params, numberOfMpi=nproc,numberOfThreads=nT)
        
          
    def createOutputStep(self):
        
        extradir = self._getExtraPath()
        inputset = self.tiltpair.get()
        uSet = inputset.getUntilted()
        tSet = inputset.getTilted()
        
        
        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = self._createSetOfCoordinates(uSet, suffix='Untilted')
        uCoordSet.setBoxSize(self.untiltcoor.get().getBoxSize())
        tCoordSet = self._createSetOfCoordinates(tSet, suffix='Tilted')
        readSetOfCoordinates(extradir, uSet, uCoordSet)
        uCoordSet.write()
        tCoordSet.setBoxSize(self.tiltcoor.get().getBoxSize())
        readSetOfCoordinates(extradir, tSet, tCoordSet)
        tCoordSet.write()
        
        
        # Create CoordinatesTiltPair object
        outputset = CoordinatesTiltPair(filename=self._getPath('coordinates_pairs.sqlite'))
        outputset.setTilted(tCoordSet)
        outputset.setUntilted(uCoordSet)
        outputset.setMicsPair(inputset)



        for coordU, coordT in izip(uCoordSet, tCoordSet):
            outputset.append(TiltPair(coordU, coordT))

        self._defineOutputs(outputCoordinatesTiltPair=outputset)
        self._defineSourceRelation(self.tiltpair, outputset)
        outputset.setObjComment(self.getSummary(outputset))

        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        if self.untiltcoor.get() and not self.untiltcoor.hasValue():
            validateMsgs.append('Please provide input coordinates.')  
        if self.tiltcoor.get() and not self.tiltcoor.hasValue():
            validateMsgs.append('Please provide input coordinates.')         
        return validateMsgs
        
   
    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputCoordinatesTiltPair')):
            summary.append("Output tilpairs not ready yet.")
        else:
            #summary.append("Particles matched: " )
            summary.append("Particle box size: %d" %self.untiltcoor.get().getBoxSize())
        return summary
    
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputCoordinatesTiltPair')):
            messages.append('The assignment has been performed using and affinity transformation [Publication: Not yet]')
        return messages
    
    def _citations(self):
        return ['Not yet']
    
    def getSummary(self, coordsSet):
        summary = []
        summary.append("Particles picked:")
        #summary.append("Particles picked: %d" %coordsSet.getSize())
        return "\n"#.join(summary)
    
    
    