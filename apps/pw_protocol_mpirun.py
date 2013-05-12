#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
This module is responsible for launching protocol executions.
"""
import sys
from pyworkflow.utils import runProtocol
from mpi4py import MPI


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    projName = sys.argv[1]
    protId = int(sys.argv[2])
    
    if rank == 0:
        print "="*100
        print "projName: ", projName
        print "protId: ", protId        
        runProtocol(projName, protId, comm)
        
    else:
        from pyworkflow.manager import Manager
        project = Manager().createProject(projName) # Create the project to change dir
        from pyworkflow.utils.mpi import runJobMPISlave
        runJobMPISlave(comm)
