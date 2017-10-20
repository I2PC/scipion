# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import join
from pyworkflow.utils import Environ


def getEnviron():
    """ Setup the environment variables needed to launch Eman. """
    environ = Environ(os.environ)
    EMAN2DIR = os.environ['EMAN2DIR']
    pathList = [os.path.join(EMAN2DIR, d)
                for d in ['lib', 'bin', 'extlib/site-packages']]

    # This environment variable is used to setup OpenGL (Mesa)
    # library in remote desktops
    if 'REMOTE_MESA_LIB' in os.environ:
        pathList.append(os.environ['REMOTE_MESA_LIB'])

    environ.update({
            'PATH': join(EMAN2DIR, 'bin'),
            'LD_LIBRARY_PATH': os.pathsep.join(pathList),
            'PYTHONPATH': os.pathsep.join(pathList),
            'EMAN_PYTHON': os.path.join(EMAN2DIR, 'Python/bin/python')
            }, position=Environ.REPLACE)
    return environ