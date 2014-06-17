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

import os, sys
from pyworkflow.manager import Manager

from pyworkflow.em import *
import pyworkflow.em.packages.xmipp3 as xmipp3
from xmipp import *
from pyworkflow.utils.path import moveTree
from pyworkflow.em.packages.xmipp3 import readSetOfCoordinates



if __name__ == '__main__':

    projectId = sys.argv[1]
    project = Manager().loadProject(projectId)
    protId = int(sys.argv[2])
    prot = project.getProtocol(protId)
    
    
    
    inputset = prot.inputMicrographs.get()
    extradir = prot._getExtraPath()
    count = 0;
    for key, output in prot.iterOutputAttributes(EMObject):
        count += 1
    
    suffix = str(count + 1) if count > 0 else ''
    outputName = 'outputCoordinates' + suffix
    outputset = prot._createSetOfCoordinates(inputset, suffix=suffix)#micrographs are the input set if protocol is not finished
    readSetOfCoordinates(extradir, outputset.getMicrographs(), outputset)
    outputs = {outputName: outputset}
    prot._defineOutputs(**outputs)
    prot._defineSourceRelation(inputset, outputset)
    
    project._updateProtocol(prot)
        
    