# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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
from pyworkflow.utils import Environ, getEnvVariable
from pyworkflow.em.packages.eman2 import EMAN_DIR_VAR

SCRATCHDIR = getEnvVariable('EMAN2SCRATCHDIR')


def getEnviron():
    """ Setup the environment variables needed to launch Eman. """
    environ = Environ(os.environ)
    EMAN2DIR = os.environ[('%s' % EMAN_DIR_VAR)]
    pathList = [os.path.join(EMAN2DIR, d) for d in ['lib', 'bin']]

    if getVersion() in ['2.11', '2.12']:
        pathList.append(os.path.join(EMAN2DIR, 'extlib/site-packages'))

    # This environment variable is used to setup OpenGL (Mesa)
    # library in remote desktops
    if 'REMOTE_MESA_LIB' in os.environ:
        pathList.append(os.environ['REMOTE_MESA_LIB'])

    environ.update({'PATH': join(EMAN2DIR, 'bin')}, position=Environ.BEGIN)

    environ.update({
        'LD_LIBRARY_PATH': os.pathsep.join(pathList),
        'PYTHONPATH': os.pathsep.join(pathList),
        'SCIPION_MPI_FLAGS': os.environ.get('EMANMPIOPTS', '')
    }, position=Environ.REPLACE)

    if getVersion() in ['2.11', '2.12']:
        environ.update({'EMAN_PYTHON': os.path.join(EMAN2DIR, 'Python/bin/python')
                        }, position=Environ.END)
    return environ


def getVersion():
    path = os.environ['EMAN2DIR']
    for v in getSupportedVersions():
        if v in path:
            return v
    return ''


def getSupportedVersions():
    return ['2.11', '2.12', '2.21']


def isNewVersion():
    return not getVersion().startswith("2.1")


def validateVersion(protocol, errors):
    """ Validate if eman version is set properly according
     to installed version and the one set in the config file.
     Params:
        protocol: the input protocol calling to validate
        errors: a list that will be used to add the error message.
    """
    protocol.validatePackageVersion('EMAN2DIR', errors)


def getEmanProgram(program):
    if not 'EMAN_PYTHON' in os.environ:
        if getVersion() in ['2.11', '2.12']:
            pyPath = 'Python/bin/python'
        else:
            pyPath = 'bin/python'
        os.environ['EMAN_PYTHON'] = os.path.join(os.environ['EMAN2DIR'], pyPath)
    # For EMAN2 python scripts, join the path to bin
    program = os.path.join(os.environ['EMAN2DIR'], 'bin', program)
    python = os.environ['EMAN_PYTHON']
    return '%(python)s %(program)s ' % locals()


def getEmanCommand(program, args):
    return getEmanProgram(program) + args


def getBoxerCommand(emanVersion, boxerVersion='new'):
    # returns boxer program depending on Eman version
    if emanVersion in ['2.11', '2.12']:
        return 'e2boxer.py'
    else:
        if boxerVersion == 'new':
            return 'e2boxer.py'
        else:
            return 'e2boxer_old.py'
