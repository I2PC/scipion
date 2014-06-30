#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:    Airen Zaldivar         (airenzp@gmail.com) 
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

import sys
from pyworkflow.em.data import EMObject
from pyworkflow.em.protocol import getProtocolFromDb
from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates, readAnglesFromMicrographs
from pyworkflow.em.data_tiltpairs import CoordinatesTiltPair, SetOfAngles


if __name__ == '__main__':


    dbpath = sys.argv[1]
    protid = sys.argv[2]
    prot = getProtocolFromDb(dbpath, protid)
    extradir = prot._getExtraPath()
    count = 0
    
    for key, output in prot.iterOutputAttributes(EMObject):
        count += 1
    
    suffix = str(count + 1) if count > 0 else ''
    if prot.getClassName() == "XmippProtParticlePicking":
        inputset = prot.inputMicrographs.get()
        outputName = 'outputCoordinates' + suffix
        outputset = prot._createSetOfCoordinates(inputset, suffix=suffix)#micrographs are the input set if protocol is not finished
        readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    if prot.getClassName() == "XmippProtParticlePickingPairs":
        inputset = prot.inputMicrographsTiltedPair.get()
        uSet = inputset.getUntilted()
        tSet = inputset.getTilted()
        outputName = 'outputCoordinatesTiltPair' + suffix
        uSuffix = 'Untilted' + suffix
        tSuffix = 'Tilted' + suffix
        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = prot._createSetOfCoordinates(uSet, suffix=uSuffix)
        readSetOfCoordinates(extradir, uSet, uCoordSet)
        uCoordSet.write()
        tCoordSet = prot._createSetOfCoordinates(tSet, suffix=tSuffix)
        readSetOfCoordinates(extradir, tSet, tCoordSet)
        tCoordSet.write()
        
        # Read Angles from input micrographs
        micsFn = prot._getPath('input_micrographs.xmd')
        setAngles = prot._createSetOfAngles(suffix=suffix)
        readAnglesFromMicrographs(micsFn, setAngles)
        setAngles.write()
        # Create CoordinatesTiltPair object
        outputset = CoordinatesTiltPair()
        outputset.setTilted(tCoordSet)
        outputset.setUntilted(uCoordSet)
        outputset.setAngles(setAngles)

    outputs = {outputName: outputset}
    prot._defineOutputs(**outputs)
    prot._defineSourceRelation(inputset, outputset)
    prot._store()
    
        
    