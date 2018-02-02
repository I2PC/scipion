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
import sys
from os.path import join#, exist
from pyworkflow.utils import Environ

def getVersion():
    path = os.environ['LOCSCALE_HOME']
    for v in getSupportedVersions():
        if v in path:
            return v
    return ''


def getSupportedVersions():
    return ['0.1']


def getSupportedEmanVersions():
    return ['2.11, 2.12, 2.2']


def getEmanVersion():
    path = os.environ['EMAN2DIR']
    for v in getSupportedEmanVersions():
        if v in path:
            return v
    return ''

def validateEmanVersion(protocol, errors):
    """ Validate if eman version is set properly according
     to installed version and the one set in the config file.
     Params:
        protocol: the input protocol calling to validate
        errors: a list that will be used to add the error message.
    """
    if getEmanVersion == '':
        errors.append('Eman%s is needed to execute this protocol'
                      % getSupportedEmanVersions())

def getEmanEnviron():
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

def setEnviron():
    """ Setup the environment variables needed to import localrec classes. """
    os.environ.update(getEmanEnviron())
    sys.path.append(os.path.join(os.environ["LOCSCALE_HOME"], "source"))
    if not 'EMAN_PYTHON' in os.environ:
        os.environ['EMAN_PYTHON'] = os.path.join(os.environ['EMAN2DIR'],
                                                 'Python/bin/python')

def getProgram(program):
    # For localscale python scripts, join the path to source
    if not 'EMAN_PYTHON' in os.environ:
        setEnviron()
    program = os.path.join(os.environ['LOCSCALE_HOME'], 'source', program)
    #raise Exception('EMAN_PYTHON is not load in environment')
    python = os.environ['EMAN_PYTHON']
    return '%(python)s %(program)s ' % locals()
