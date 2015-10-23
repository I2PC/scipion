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
from pyworkflow.protocol.params import (PointerParam, FloatParam, EnumParam, STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.utils.path import makePath
from pyworkflow.em.data_tiltpairs import TiltPair, CoordinatesTiltPair, Coordinate
from pyworkflow.em import ProtParticlePicking
from pyworkflow.em.packages.xmipp3 import XmippProtocol

import xmipp
from convert import readSetOfCoordinates, writeSetOfCoordinates, izip

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
        
        form.addParam('tiltAngle', FloatParam, default=-1, 
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

        form.addParam('maxShift', FloatParam, default=1000, expertLevel=LEVEL_ADVANCED,
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
            #untiltpath, unna = split(self.untiltCoor.get().getFileName())
            #tiltpath, tna = split(self.tiltCoor.get().getFileName())
            Unname, ext = splitext(Unname)
            Tpath, Tname = split(micTilted.getFileName())
            Tname, ext = splitext(Tname)
            
            fnUntilt = 'particles@'+self._getExtraPath("untilted/")+Unname+'.pos'
            #print fnUntilt
            fnTilt = 'particles@'+self._getExtraPath("tilted/")+Tname+'.pos'
            fnmicsize = tiltPair.getTilted().getFileName()
            stepId=self._insertFunctionStep('assignmentStep',fnUntilt, fnTilt, fnmicsize, self._getExtraPath())
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
                #U_set = self._createSetOfCoordinates(micUntilted,Unname)
                U_set = self._createSetOfCoordinates(self.tiltpair.get().getUntilted(), 
                                                     suffix='_untilted')
                #T_set = self._createSetOfCoordinates(micTilted,Tname)
                T_set = self._createSetOfCoordinates(self.tiltpair.get().getTilted(), 
                                                     suffix='_tilted')
                
                untiltPar = self.untiltPar.get()
                for particle_u in untiltPar:  
                    newCoord = particle_u.getCoordinate().clone()
                    newCoord.copyObjId(particle_u)
#                     print '------------------'
#                     print "particle_u.getLocation()", particle_u.getLocation()                  
#                     print "particle_u.getobjId():", particle_u.getObjId()
                    U_set.append(newCoord)
#
                tiltPar = self.tiltPar.get()
                for particle_t in tiltPar:
                    newCoord = particle_t.getCoordinate().clone()
                    newCoord.copyObjId(particle_t)
                    T_set.append(newCoord)

#                 for i, mic_ in enumerate(self.tiltpair.get()):
#                 print '---------------------------'
#                 print indx
#                 print '---------------------------'
#                 for particle_u in self.untiltPar.get():
#                     coordOrig_u = particle_u.getCoordinate()
#                     coordAux_u = []
#                     if (indx + 1 == coordOrig_u.getMicId()):
#                         coordAux_u = Coordinate()
#                         coordAux_u.setX(coordOrig_u.getX())
#                         coordAux_u.setY(coordOrig_u.getY())
#                         coordAux_u.setMicName(coordOrig_u.getMicName())
#                         coordAux_u.setMicId(coordOrig_u.getMicId())
#                         coordAux_u.setMicrograph(micUntilted)
#                         U_set.append(coordAux_u)
#                 print U_set
#                     
#                 for particle_t in self.tiltPar.get():
#                     coordOrig_t = particle_t.getCoordinate()
#                     coordAux_t =[]
#                     if (indx + 1 == coordOrig_t.getMicId()):
#                         coordAux_t = Coordinate()
#                         coordAux_t.setX(coordOrig_t.getX())
#                         coordAux_t.setY(coordOrig_t.getY())
#                         coordAux_t.setMicName(coordOrig_t.getMicName())
#                         coordAux_t.setMicId(coordOrig_t.getMicId())
#                         coordAux_t.setMicrograph(micTilted)
#                         T_set.append(coordAux_t)
                aux = self.untiltPar.get()
                aux2 = aux[1]
                bos, y_, z_ = aux2.getDim()
                U_set.setBoxSize(bos)
                T_set.setBoxSize(bos)
        else:
            U_set = self.untiltCoor.get()
            T_set = self.tiltCoor.get() 
            
        print "type(U_set.getMicrographs())", U_set.getMicrographs()
        print "type(T_set.getMicrographs())", T_set.getMicrographs()
        writeSetOfCoordinates(self._getExtraPath("untilted"), U_set)
        writeSetOfCoordinates(self._getExtraPath("tilted"), T_set)
    
    
    def assignmentStep(self,fnuntilt, fntilt, fnmicsize, Unpath):

        params =  ' --untiltcoor %s' % fnuntilt        
        params += ' --tiltcoor %s' % fntilt
        params += ' --tiltmicsize %s' % fnmicsize
        params += ' --maxshift %f' % self.maxShift.get()
        if self.typeOfSet.get() == TYPE_PARTICLES:
#             boxsize = self.untiltPar.get().getCoordinates().getBoxSize()
#             if boxsize is None:
                aux = self.untiltPar.get()
                aux2 = aux[1]
                boxsize, y_, z_ = aux2.getDim()
                params += ' --particlesize %f' % boxsize
#             else:
#                 params += ' --particlesize %f' % boxsize
                
        else:
            params += ' --particlesize %f' % self.untiltCoor.get().getBoxSize()
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
            #summary.append("Particles matched: " )
            summary.append("Particle box size: %d" %self.untiltCoor.get().getBoxSize())
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
    
    
    