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
from pyworkflow.em.data import EMObject
from pyworkflow.protocol.params import (PointerParam, FloatParam, EnumParam, STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.utils.path import makePath
from pyworkflow.em.data_tiltpairs import TiltPair, CoordinatesTiltPair, Coordinate
from pyworkflow.em import ProtParticlePicking
from pyworkflow.em.packages.xmipp3 import XmippProtocol, readAnglesFromMicrographs

import xmipp
from convert import readSetOfCoordinates, writeSetOfCoordinates, writeSetOfMicrographsPairs, izip

TYPE_COORDINATES = 0
TYPE_PARTICLES = 1


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
        
        form.addParam('typeOfSet', EnumParam, choices=['Coordinates', 'Particles'], 
                      default=0, label='Input type', display=EnumParam.DISPLAY_COMBO, 
                      help='Select a Set of Coordinates or a Set or Particles.')
        
        form.addParam('untiltCoor', PointerParam, pointerClass='SetOfCoordinates', 
                      label="Untilted coordinates", condition='typeOfSet==%d' % TYPE_COORDINATES,  
                      help='Select the metadata with untilt coordinates.')
        
        form.addParam('tiltCoor', PointerParam, pointerClass='SetOfCoordinates', 
                      label="Tilted coordinates", condition='typeOfSet==%d' % TYPE_COORDINATES,   
                      help='Select the metadata with tilt coordinates.')
        
        form.addParam('untiltPar', PointerParam, pointerClass='SetOfParticles', 
                      label="Untilted particles", condition='typeOfSet==%d' % TYPE_PARTICLES,
                      help='Select the metadata with untilt particles.')
        
        form.addParam('tiltPar', PointerParam, pointerClass='SetOfParticles', 
                      label="Tilted particles", condition='typeOfSet==%d' % TYPE_PARTICLES,   
                      help='Select the metadata with tilt particles.')
        
#         form.addParam('boxsiz', FloatParam, default=100, 
#                       label="Box Size", condition='typeOfSet==%d' % TYPE_PARTICLES, 
#                       help='Sometimes, particles do not have information about the particle size\n,'
#                       'i.e. filtered by ZScore. In those cases, the particle size or box size is requiered \n'
#                       'By default: boxSize = 100') 
        
        form.addParam('tiltAngle', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label="Tilt angle",  
                      help='Tilt angle estimation, the method will look for the assignment in the\n'
                      ' interval of [tilt_angle-15, tilt_angle+15].\n'
                      'By default: tilt angle = -1, if there is not any information about the tilt angle')     
        
        form.addParam('threshold', FloatParam, default=0.25, expertLevel=LEVEL_ADVANCED,
                      label="Threshold value",  
                      help='Parameter between 0 and 1 that allows to define if \n' 
                      'a tilt point can be matched with a certain untilt point. \n'
                      'The matching is performed only if the distance is lesser than \n'
                      'threshold * particlesize.')

        form.addParam('maxShift', FloatParam, default=1500, expertLevel=LEVEL_ADVANCED,
                      label="Maximum shift (pixels)", 
                      help='Maximum allowed distance (in pixels) that the tilt micrograph can be shifted' 
                      'respect to the untilted micrograph')         

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):        
        self.micsFn = self._getPath()
        #self.micsFn = self._getPath('input_micrographs.xmd')
        # Convert input into xmipp Metadata format
        convertId=self._insertFunctionStep('convertInputStep')
        deps = []
        for tiltPair in self.tiltpair.get():
            micUntilted = tiltPair.getUntilted()
            micTilted = tiltPair.getTilted()
            Unpath, Unname = split(micUntilted.getFileName())
            Unname, ext = splitext(Unname)
            Tpath, Tname = split(micTilted.getFileName())
            Tname, ext = splitext(Tname)
            fnUntilt = 'particles@'+self._getExtraPath("untilted/")+Unname+'.pos'
            fnTilt = 'particles@'+self._getExtraPath("tilted/")+Tname+'.pos'
            fnmicsize = tiltPair.getTilted().getFileName()
            fnposUntilt = 'particles@'+self._getExtraPath("")+Unname+'.pos'
            fnposTilt = 'particles@'+self._getExtraPath("")+Tname+'.pos'
            stepId = self._insertFunctionStep('assignmentStep',fnUntilt, fnTilt, fnmicsize, self._getExtraPath(),
                                              fnposUntilt, fnposTilt, self._getExtraPath(Unname), prerequisites=[convertId])
             
#            self._insertFunctionStep('estimateTiltAxisStep', fnposUntilt, fnposTilt, self._getExtraPath(Unname), prerequisites=[convertId])
             
            deps.append(stepId)
            
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
        

    def convertInputStep(self):
        """ Read the input metadatata.
        """
        # Get the converted input micrographs in Xmipp format
        makePath(self._getExtraPath("untilted"))
        makePath(self._getExtraPath("tilted"))
        

            
        if self.typeOfSet.get() == TYPE_PARTICLES:
            U_set = self.untiltPar.get().getCoordinates()
            T_set = self.tiltPar.get().getCoordinates()

            if (U_set or T_set) is None:
                U_set = self._createSetOfCoordinates(self.tiltpair.get().getUntilted(), 
                                                     suffix='_untilted')
                T_set = self._createSetOfCoordinates(self.tiltpair.get().getTilted(), 
                                                     suffix='_tilted')
                
                untiltPar = self.untiltPar.get()
                for particle_u in untiltPar:  
                    newCoord = particle_u.getCoordinate().clone()
                    newCoord.copyObjId(particle_u)
                    U_set.append(newCoord)
#
                tiltPar = self.tiltPar.get()
                for particle_t in tiltPar:
                    newCoord = particle_t.getCoordinate().clone()
                    newCoord.copyObjId(particle_t)
                    T_set.append(newCoord)

                aux = self.untiltPar.get()
                aux2 = aux[1]
                bos, y_, z_ = aux2.getDim()
                U_set.setBoxSize(bos)
                T_set.setBoxSize(bos)
        else:
            U_set = self.untiltCoor.get()
            T_set = self.tiltCoor.get() 
            
        writeSetOfCoordinates(self._getExtraPath("untilted"), U_set)
        writeSetOfCoordinates(self._getExtraPath("tilted"), T_set)
#         writeSetOfMicrographsPairs(self.tiltpair.get().getUntilted(), 
#                                    self.tiltpair.get().getTilted(), 
#                                    self.micsFn)
    
    def assignmentStep(self,fnuntilt, fntilt, fnmicsize, Unpath, fnposUntilt, fnposTilt, opath):

        params =  ' --untiltcoor %s' % fnuntilt        
        params += ' --tiltcoor %s' % fntilt
        params += ' --tiltmicsize %s' % fnmicsize
        params += ' --maxshift %f' % self.maxShift.get()
        if self.typeOfSet.get() == TYPE_PARTICLES:
                aux = self.untiltPar.get()
                aux2 = aux[1]
                boxsize, y_, z_ = aux2.getDim()
                params += ' --particlesize %f' % boxsize         
        else:
            params += ' --particlesize %f' % self.untiltCoor.get().getBoxSize()
        params += ' --threshold %f' % self.threshold.get()
        params += ' --odir %s' % Unpath


        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        self.runJob('xmipp_image_assignment_tilt_pair', params, numberOfMpi=nproc,numberOfThreads=nT)
        

        
        self.estimateTiltAxis(fnposUntilt, fnposTilt)#, self.micsFn)

        
        
    def estimateTiltAxis(self, fnposUntilt, fnposTilt):#, opath):

        params =  ' --untilted %s' % fnposUntilt        
        params += ' --tilted %s' % fnposTilt
        params += ' -o %s' % self._getPath()+'/input_micrographs.xmd'

        self.runJob('xmipp_angular_estimate_tilt_axis', params)
        
          
    def createOutputStep(self):
        
        extradir = self._getExtraPath()
        inputset = self.tiltpair.get()
        uSet = inputset.getUntilted()
        tSet = inputset.getTilted()
        
        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = self._createSetOfCoordinates(uSet, suffix='Untilted')
        if self.typeOfSet.get() == TYPE_PARTICLES:
#             boxsize_u = self.untiltPar.get().getCoordinates().getBoxSize()
#             if boxsize_u is None:
                aux = self.untiltPar.get()
                aux2 = aux[1]
                boxsize_u, y_, z_ = aux2.getDim()
        else:
            boxsize_u = self.untiltCoor.get().getBoxSize()
        
        uCoordSet.setBoxSize(boxsize_u)    
        tCoordSet = self._createSetOfCoordinates(tSet, suffix='Tilted')
        readSetOfCoordinates(extradir, uSet, uCoordSet)
        uCoordSet.write()
        
        if self.typeOfSet.get() == TYPE_PARTICLES:
#             boxsize_t = self.tiltPar.get().getCoordinates().getBoxSize()
#             if boxsize_t is None:
                aux = self.tiltPar.get()
                aux2 = aux[1]
                boxsize_t, y_, z_ = aux2.getDim()
        else:
            boxsize_t = self.tiltCoor.get().getBoxSize()
        
        tCoordSet.setBoxSize(boxsize_t)
        readSetOfCoordinates(extradir, tSet, tCoordSet)
        tCoordSet.write()
        

        setAngles = self._createSetOfAngles()
        pathangles = self.micsFn + '/input_micrographs.xmd'
        readAnglesFromMicrographs(pathangles, setAngles)
        
        setAngles.write()
        # Create CoordinatesTiltPair object
        outputset = CoordinatesTiltPair(filename=self._getPath('coordinates_pairs.sqlite'))
        outputset.setTilted(tCoordSet)
        outputset.setUntilted(uCoordSet)
        outputset.setAngles(setAngles)
        outputset.setMicsPair(inputset)
        outputset.setObjComment(self.getSummary(outputset))
#         print '-----------------------------'
#         print setAngles
#         print '-----------------------------'
        for coordU, coordT in izip(uCoordSet, tCoordSet):
            outputset.append(TiltPair(coordU, coordT))

        self._defineOutputs(outputCoordinatesTiltPair=outputset)
        self._defineSourceRelation(self.tiltpair, outputset)
        

        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        
        if self.typeOfSet.get() == TYPE_PARTICLES:
            validateMsgs = []
            if self.untiltPar.get() and not self.untiltPar.hasValue():
                validateMsgs.append('Please provide input particles.')  
            if self.tiltPar.get() and not self.tiltPar.hasValue():
                validateMsgs.append('Please provide input particles.')         
            return validateMsgs
        else:
            validateMsgs = []
            if self.untiltCoor.get() and not self.untiltCoor.hasValue():
                validateMsgs.append('Please provide input coordinates.')  
            if self.tiltCoor.get() and not self.tiltCoor.hasValue():
                validateMsgs.append('Please provide input coordinates.')         
            return validateMsgs
         
        
   
    def _summary(self):
        summary = []
        if  (not hasattr(self,'outputCoordinatesTiltPair')):
            summary.append("Output tilpairs not ready yet.")
        else:
            if self.getOutputsSize() > 0:
                return self._summary_aux()
            else:
                return [self.getSummary(self.outputCoordinatesTiltPair.getUntilted())]   
        return summary
   
   
    def _summary_aux(self):
        summary = []
        if self.getInputMicrographs() is  not None:
            summary.append("Number of input micrographs: %d" % self.getInputMicrographs().getSize())

        if(self.getOutputsSize() > 1):
            for key, output in self.iterOutputAttributes(EMObject):
                label = output.getObjLabel() if output.getObjLabel() != "" else key
                summary.append("*%s:*%s"%(key, output.getObjComment()))
                summary.append("      Particles picked: %d" %output.getSize())
                summary.append("      Particle size: %d \n" %output.getBoxSize())
        elif(self.getOutputsSize() == 1):
            summary.append(self.getCoords().getObjComment())
        return summary
    
    def getSummary(self, coordsSet):
        summary = []
        summary.append("Particles picked: %d" %coordsSet.getSize())
        summary.append("Particle size: %d" %self.outputCoordinatesTiltPair.getUntilted().getBoxSize())
        return "\n".join(summary)
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputCoordinatesTiltPair')):
            messages.append('The assignment has been performed using and affinity transformation [Publication: Not yet]')
        return messages
    
    def _citations(self):
        return ['Not yet']
    
   
    def getCoords(self):
        count = self.getOutputsSize()
        suffix = str(count) if count > 1 else ''
        outputName = 'outputCoordinatesTiltPair' + suffix
        return getattr(self, outputName)
#     
    def registerCoords(self, coordsDir):
        #from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates, readAnglesFromMicrographs

        count = self.getOutputsSize()
        suffix = str(count + 1) if count > 0 else ''
        inputset = self.tiltpair.get()
        uSet = inputset.getUntilted()
        tSet = inputset.getTilted()
        outputName = 'outputCoordinatesTiltPair' + suffix
        uSuffix = 'Untilted' + suffix
        tSuffix = 'Tilted' + suffix
        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = self._createSetOfCoordinates(uSet, suffix=uSuffix)
        readSetOfCoordinates(coordsDir, uSet, uCoordSet)
        uCoordSet.write()
        tCoordSet = self._createSetOfCoordinates(tSet, suffix=tSuffix)
        readSetOfCoordinates(coordsDir, tSet, tCoordSet)
        tCoordSet.write()

        # Read Angles from input micrographs
        micsFn = self._getPath('input_micrographs.xmd')
        setAngles = self._createSetOfAngles(suffix=suffix)
        readAnglesFromMicrographs(micsFn, setAngles)
        setAngles.write()
        # Create CoordinatesTiltPair object
        outputset = CoordinatesTiltPair(filename=self._getPath('coordinates_pairs%s.sqlite' % suffix))
        outputset.setTilted(tCoordSet)
        outputset.setUntilted(uCoordSet)
        outputset.setAngles(setAngles)
        outputset.setMicsPair(inputset)
        for coordU, coordT in izip(uCoordSet, tCoordSet):
            outputset.append(TiltPair(coordU, coordT))

        summary = self.getSummary(outputset)
        outputset.setObjComment(summary)
        outputs = {outputName: outputset}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.tiltpair, outputset)
        self._store()
        
    def getInputMicrographs(self):
        return self.tiltpair.get().getTilted()


    