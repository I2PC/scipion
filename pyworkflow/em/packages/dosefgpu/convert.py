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

import os
from os.path import join

from pyworkflow.utils import Environ
from pyworkflow.em.data import MovieAlignment



def getEnviron():
    """ Return the envirion settings to run dosefgpu programs. """
    """ Setup the environment variables needed to launch Relion. """
    environ = Environ(os.environ)
    environ.update({
            'PATH': join(os.environ['DOSEFGPU_HOME']),
            'LD_LIBRARY_PATH': join(os.environ['DOSEFGPU_CUDA_LIB'])
            }, position=Environ.BEGIN)
    return environ  


def parseMovieAlignment(logFile):
    """ Get the first and last frames together with the shifts
    between frames aligned.
    """
    f = open(logFile)
    first = None
    shifts = []
    
    for line in f:
        l = line.strip()
        if 'Add Frame #' in l:
            parts = l.split()
            if first is None: # read the first frame number
                first = int(parts[2][1:]) # take the id from #000 format
            # take the id from the last two colums of the line
            shifts.append(int(parts[-2]))
            shifts.append(int(parts[-1]))
            
    f.close()
    
    return MovieAlignment(first=first, last=first+len(shifts)-1, shifts=shifts)
    
    