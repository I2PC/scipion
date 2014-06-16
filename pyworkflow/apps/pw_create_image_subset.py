#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:    Airen Zaldivar         (airenzp@gmail.com) 
#               J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.manager import Manager
from pyworkflow.em import ProtUserSubSet
from pyworkflow.utils import timeit


def runSubsetProtocol(projectId, inputId, sqliteFile, 
                      outputType, protocolLabel, other):
    """ Load the project and launch the protocol to
    create the subset.
    """
    # Retrieve project and input object from db
    project = Manager().loadProject(projectId)
    inputObject = project.mapper.selectById(int(inputId))

    # Create the new protocol instance and set the input values
    prot = project.newProtocol(ProtUserSubSet)
    prot.inputObject.set(inputObject)
    prot.setObjLabel(protocolLabel)
    prot.sqliteFile.set(sqliteFile)
    prot.outputClassName.set(outputType)
    
    # Launch the protocol
    project.launchProtocol(prot, wait=True)
    

if __name__ == '__main__':
    #TODO: REMOVE THIS AFTER DEBUGGING
    print "ARGS: ", sys.argv
    other = sys.argv[6:]
    runSubsetProtocol(projectId=sys.argv[1],
                      inputId=sys.argv[2],
                      sqliteFile=sys.argv[3],
                      outputType=sys.argv[4],
                      protocolLabel=sys.argv[5], 
                      other=other)
    sys.exit(0)

