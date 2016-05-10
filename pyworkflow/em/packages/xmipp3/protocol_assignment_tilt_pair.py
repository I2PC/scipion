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

from itertools import izip

import pyworkflow.protocol.params as params
from pyworkflow.utils.path import makePath, removeBaseExt
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.data_tiltpairs import TiltPair, CoordinatesTiltPair

from convert import readAnglesFromMicrographs
from protocol_particle_pick_pairs import XmippProtParticlePickingPairs

from convert import readSetOfCoordinates, writeSetOfCoordinates

TYPE_COORDINATES = 0
TYPE_PARTICLES = 1


class XmippProtAssignmentTiltPair(XmippProtParticlePickingPairs):
    """    
    From two sets of points (tilted and untilted) the protocol determines
    the affine transformation between these sets.
    """
    _label = 'assignment tiltpair'

    def __init__(self, *args, **kwargs):
        XmippProtParticlePickingPairs.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

        def onChangeInputType():
            """ Dynamically change the PointerClass depending on the
             selected input type (either Coordinates or Particles).
            """
            pointerClass = 'SetOf%s' % self.getEnumText('inputType')
            self.getParam('untiltedSet').setPointerClass(pointerClass)
            self.getParam('tiltedSet').setPointerClass(pointerClass)

        self.inputType.trace(onChangeInputType)

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMicrographsTiltedPair', params.PointerParam,
                      pointerClass='MicrographsTiltPair',
                      important=True,
                      label="Micrograph tilt pair",
                      help='Select micrographs tilt pair.')

        form.addParam('inputType', params.EnumParam,
                      choices=['Coordinates', 'Particles'], default=0,
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Input type',
                      help='Select a Set of Coordinates or a Set or Particles.')

        form.addParam('untiltedSet', params.PointerParam,
                      pointerClass='SetOfCoordinates,SetOfParticles',
                      label="Untilted input",
                      help='Select the untilted input set, it can be either '
                           'coordinates or particles (that contains coordinates.')

        form.addParam('tiltedSet', params.PointerParam,
                      pointerClass='SetOfCoordinates,SetOfParticles',
                      label="Tilted input",
                      help='Select the tilted input set, it can be either '
                           'coordinates or particles (that contains coordinates. '
                           'It should be of the same type of the input untilted.')

        form.addParam('tiltAngle', params.FloatParam, default=-1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Tilt angle",
                      help='Tilt angle estimation, the method will look for '
                           'the assignment in the interval of [tilt_angle-15, '
                           'tilt_angle+15].\n By default: tilt angle = -1, if '
                           'there is not any information about the tilt angle')

        form.addParam('threshold', params.FloatParam, default=0.25,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Threshold value",
                      help='Parameter between 0 and 1 that allows to define if \n'
                      'a tilt point can be matched with a certain untilt point.\n'
                      'The matching is performed only if the distance is lesser\n'
                      'than threshold * particlesize.')

        form.addParam('maxShift', params.FloatParam, default=1500,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Maximum shift (pixels)",
                      help='Maximum allowed distance (in pixels) that the tilt '
                           'micrograph can be shifted respect to the untilted '
                           'micrograph')

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _particlesPos(self, *parts):
        return 'particles@%s.pos' % self._getExtraPath(*parts)

    def _micBaseName(self, mic):
        return removeBaseExt(mic.getFileName())

    def _coordsFromParts(self, micSet, partSet, suffix):
        """ Create and fill a SetOfCoordinates from a given SetOfParticles. """
        coordSet = self._createSetOfCoordinates(micSet, suffix='_untilted')

        for particle in partSet:
            coord = particle.getCoordinate().clone()
            coord.copyObjId(particle)
            coordSet.append(coord)
        coordSet.setBoxSize(partSet.getXDim())

        return coordSet

    def _getBoxSize(self):
        untiltedSet = self.untiltedSet.get()
        if isinstance(untiltedSet, SetOfParticles):
            boxSize = untiltedSet.getXDim()
        else:
            boxSize = untiltedSet.getBoxSize() # coordinates
        return boxSize

    def _insertAllSteps(self):
        self.micsFn = self._getPath()
        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep')
        deps = []
        for tiltPair in self.inputMicrographsTiltedPair.get():
            uName = self._micBaseName(tiltPair.getUntilted())
            tName = self._micBaseName(tiltPair.getTilted())
            stepId = self._insertFunctionStep('assignmentStep',
                                              self._particlesPos("untilted", uName),
                                              self._particlesPos("tilted", tName),
                                              tiltPair.getTilted().getFileName(),
                                              self._particlesPos(uName),
                                              self._particlesPos(tName),
                                              prerequisites=[convertId])
            deps.append(stepId)

        self._insertFunctionStep('createOutputStep', prerequisites=deps)

    def convertInputStep(self):
        """ Read the input metadatata. """
        # Get the converted input micrographs in Xmipp format
        makePath(self._getExtraPath("untilted"))
        makePath(self._getExtraPath("tilted"))

        uSet = self.untiltedSet.get()
        tSet = self.tiltedSet.get()

        # Get the untilted and tilted coordinates, depending on the input type
        if isinstance(uSet, SetOfParticles):
            uCoords = uSet.getCoordinates()
            tCoords = tSet.getCoordinates()

            # If there are not Coordinates associated to particles
            # we need to create and fill the set of coordinates
            if uCoords is None or tCoords is None:
                micTiltedPairs = self.inputMicrographsTiltedPair.get()
                uCoords = self._coordsFromParts(micTiltedPairs.getUntilted(),
                                                uSet, '_untilted')
                tCoords = self._coordsFromParts(micTiltedPairs.getTilted(),
                                                tSet, '_tilted')
        else:
            uCoords = uSet
            tCoords = tSet

        writeSetOfCoordinates(self._getExtraPath("untilted"), uCoords)
        writeSetOfCoordinates(self._getExtraPath("tilted"), tCoords)

    def assignmentStep(self, fnuntilt, fntilt, fnmicsize, fnposUntilt, fnposTilt):
        params =  ' --untiltcoor %s' % fnuntilt
        params += ' --tiltcoor %s' % fntilt
        params += ' --tiltmicsize %s' % fnmicsize
        params += ' --maxshift %f' % self.maxShift
        params += ' --particlesize %d' % self._getBoxSize()
        params += ' --threshold %f' % self.threshold
        params += ' --odir %s' % self._getExtraPath()
        self.runJob('xmipp_image_assignment_tilt_pair', params)

        # Estimate the tilt axis
        params =  ' --untilted %s' % fnposUntilt
        params += ' --tilted %s' % fnposTilt
        params += ' -o %s' % self._getPath('input_micrographs.xmd')
        self.runJob('xmipp_angular_estimate_tilt_axis', params)

    def createOutputStep(self):
        self.registerCoords(self._getExtraPath(), store=False)

    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        uSet = self.untiltedSet.get()
        tSet = self.tiltedSet.get()

        if (uSet is not None and tSet is not None and
            uSet.getClassName() != tSet.getClassName()):
            errors.append('Both untilted and tilted inputs should be of the '
                          'same type. ')

        return errors

    def _methods(self):
        messages = []
        if hasattr(self,'outputCoordinatesTiltPair'):
            messages.append('The assignment has been performed using and '
                            'affinity transformation [Publication: Not yet]')
        return messages

    def _citations(self):
        return ['Not yet']

