#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:    Josue Gomez Blanco         (jgomez@gmail.com) 
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
from pyworkflow.em.packages.xmipp3 import XmippProtRecalculateCTF
from pyworkflow.em.packages.brandeis import ProtRecalculateCTFFind

@timeit
def runRecalcCftProtocol(projectId, inputObjId, sqliteFile, pathFile):
    """ Load the project and launch the protocol to
    create the subset.
    """
    # Retrieve project, input protocol and object from db
    project = Manager().loadProject(projectId)
    inputObj = project.mapper.selectById(int(inputObjId))
    parentProtId = inputObj.getObjParentId()
    parentProt = project.mapper.selectById(parentProtId)
    
    if parentProt.getClassName() == "XmippProtCTFMicrographs":
        # Create the new protocol
        prot = project.newProtocol(XmippProtRecalculateCTF)
    else:   
        # Create the new protocol
        prot = project.newProtocol(ProtRecalculateCTFFind)
    
    useQueue = parentProt.useQueue()
    Mpi = parentProt.numberOfMpi.get()
    Threads = parentProt.numberOfThreads.get()
    # Define the input params of the new protocol
    prot._useQUeue = useQueue
    prot.numberOfMpi.set(Mpi)
    prot.numberOfThreads.set(Threads)
    prot.sqliteFile.set(sqliteFile)
    prot.inputCtf.set(inputObj)
    prot.inputValues.set(pathFile)
    
    # Launch the protocol
    project.launchProtocol(prot, wait=True)

if __name__ == '__main__':
    for i, arg in enumerate(sys.argv):
        print "%02d: %s" % (i, arg)
    
    runRecalcCftProtocol(projectId=sys.argv[1],
                        inputObjId=sys.argv[2],
                        sqliteFile = sys.argv[3],
                        pathFile=sys.argv[4])

