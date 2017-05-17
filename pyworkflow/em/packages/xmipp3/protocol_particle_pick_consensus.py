# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Consensus picking protocol
"""
import os
from math import sqrt

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking
from pyworkflow.protocol.constants import *
from pyworkflow.em.data import Coordinate

import numpy as np


class XmippProtConsensusPicking(ProtParticlePicking):
    """
    Protocol to estimate the agreement between different particle picking algorithms. The protocol
    takes several Sets of Coordinates calculated by different programs and/or different parameter
    settings. Let's say     we consider N independent pickings. Then, a coordinate is considered
    to be a correct particle if M pickers have selected the same particle (within a radius in
    pixels specified in the form).
    
    If you want to be very strict, then set M=N; that is, a coordinate represents a particle if
    it has been selected by all particles (this is the default behaviour). Then you may relax
    this condition by setting M=N-1, N-2, ...
    
    If you want to be very flexible, set M=1, in this way it suffices that 1 picker has
    selected the coordinate to be considered as a particle. Note that in this way, the cleaning
    of the dataset has to be performed by other means (screen particles, 2D and 3D 
    classification, ...).
    """
    _label = 'consensus picking'
    
    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        self.stepsExecutionMode = STEPS_SERIAL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.MultiPointerParam, pointerClass='SetOfCoordinates',
                      label="Input coordinates",
                      help='Select the set of coordinates to compare')
        form.addParam('consensusRadius',params.IntParam, default=10, label="Radius",
                      help="All coordinates within this radius (in pixels) are presumed to correspond to the same particle")
        form.addParam('consensus',params.IntParam, default=-1, label="Consensus", expertLevel=LEVEL_ADVANCED,
                      help="How many times need a particle to be selected to be considered as a consensus particle. "\
                           "Set to -1 to indicate that it needs to be selected by all algorithms. Set to 1 to indicate that "\
                           "it suffices that only 1 algorithm selects the particle")

        form.addParallelSection(threads=4, mpi=0)
        
#--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        deps = []
        for micrograph in self.inputCoordinates[0].get().getMicrographs():
            stepId = self._insertFunctionStep("calculateConsensusStep", 
                                              micrograph.getObjId(), prerequisites=[])
            deps.append(stepId)
        self._insertFunctionStep("createOutputStep", prerequisites=deps)
    
    def getInputMicrographs(self):
        return self.inputCoordinates[0].get().getMicrographs()
    
    def _summary(self):
        message = []
        for i, coordinates in enumerate(self.inputCoordinates):
            protocol = self.getMapper().getParent(coordinates.get())
            message.append("Method %d %s" % (i+1, protocol.getClassLabel()))
        message.append("Radius = %d" % self.consensusRadius)
        message.append("Consensus = %d" % self.consensus)
        return message
    
    def _methods(self):
        return []    
    
    def calculateConsensusStep(self, micId):
        # Take the sampling rates
        Tm = []
        for coordinates in self.inputCoordinates:
            Tm.append(coordinates.get().getMicrographs().getSamplingRate())
        
        # Get all coordinates for this micrograph
        coords = []
        Ncoords = 0
        n=0
        for coordinates in self.inputCoordinates:
            coordArray = np.asarray([x.getPosition() 
                                     for x in coordinates.get().iterCoordinates(micId)])
            coordArray *= Tm[n]/Tm[0]
            coords.append(coordArray)
            Ncoords += coordArray.shape[0]
            n+=1
        
        allCoords = np.zeros([Ncoords,2])
        votes = np.zeros(Ncoords)
        
        # Add all coordinates in the first method
        N0 = coords[0].shape[0]
        inAllMicrographs = self.consensus <= 0 or self.consensus == len(self.inputCoordinates)
        if N0==0 and inAllMicrographs:
            return
        elif N0>0:
            allCoords[0:N0,:] = coords[0]
            votes[0:N0] = 1
        
        # Add the rest of coordinates
        Ncurrent = N0
        for n in range(1, len(self.inputCoordinates)):
            for coord in coords[n]:
                if Ncurrent>0:
                    dist = np.sum((coord - allCoords[0:Ncurrent])**2, axis=1)
                    imin = np.argmin(dist)
                    if sqrt(dist[imin]) < self.consensusRadius:
                        newCoord = (votes[imin]*allCoords[imin,]+coord)/(votes[imin]+1)
                        allCoords[imin,] = newCoord
                        votes[imin] += 1
                    else:
                        allCoords[Ncurrent,:] = coord
                        votes[Ncurrent] = 1
                        Ncurrent += 1
                else:
                    allCoords[Ncurrent, :] = coord
                    votes[Ncurrent] = 1
                    Ncurrent += 1

        # Select those in the consensus
        if self.consensus <= 0:
            consensus = len(self.inputCoordinates)
        else:
            consensus = self.consensus.get()
        consensusCoords = allCoords[votes>=consensus,:]
        jaccardIdx = float(len(consensusCoords))/(float(len(allCoords))/len(self.inputCoordinates))
        # COSS: Possible problem with concurrent writes
        with open(self._getExtraPath('jaccard.txt'), "a") as fhJaccard:
            fhJaccard.write("%d %f\n"%(micId,jaccardIdx))
        
        # Write the consensus file only if there
        # are some coordinates (size > 0)
        if consensusCoords.size:
            np.savetxt(self._getExtraPath('consensus_%06d.txt' % micId), consensusCoords)
    
    def createOutputStep(self):
        firstCoords = self.inputCoordinates[0].get()
        inputMics = firstCoords.getMicrographs()
        setOfCoordinates = self._createSetOfCoordinates(inputMics)
        setOfCoordinates.setBoxSize(firstCoords.getBoxSize())
        
        # Read all consensus particles
        for micrograph in inputMics:
            fnTmp = self._getExtraPath('consensus_%06d.txt' % micrograph.getObjId())
            if os.path.exists(fnTmp):
                coords = np.loadtxt(fnTmp)
                if coords.size == 2:  # special case with only one coordinate in consensus
                    coords = [coords]
                for coord in coords:
                    aux = Coordinate()
                    aux.setMicrograph(micrograph)
                    aux.setX(coord[0])
                    aux.setY(coord[1])
                    setOfCoordinates.append(aux)
                #cleanPath(fnTmp)

        # Set output
        self._defineOutputs(outputCoordinates=setOfCoordinates)
        
        for coordinates in self.inputCoordinates:
            self._defineSourceRelation(coordinates, self.outputCoordinates)
