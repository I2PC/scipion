# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
"""
This module contains converter functions that will serve to:
1. Write from base classes to 3DFSC specific files
2. Read from 3DFSC files to base classes
"""

import os
from pyworkflow.utils import Environ
from pyworkflow.utils.path import join

NYSBC_3DFSC_HOME_VAR = 'NYSBC_3DFSC_HOME'


def getEnviron():
    """ Setup the environment variables needed to launch 3DFSC. """
    environ = Environ(os.environ)
    NYSBC_3DFSC_HOME = os.environ[('%s' % NYSBC_3DFSC_HOME_VAR)]

    environ.update({
        'PATH': join(NYSBC_3DFSC_HOME, 'ThreeDFSC'),
    }, position=Environ.BEGIN)

    if 'PYTHONPATH' in environ:
        # this is required for python virtual env to work
        environ.set('PYTHONPATH', '', position=Environ.BEGIN)
    return environ


_environ = getEnviron()

SUPPORTED_VERSIONS = ['2.5']


def getVersion():
    path = os.environ[NYSBC_3DFSC_HOME_VAR]
    for v in SUPPORTED_VERSIONS:
        if v in path:
            return v
    return ''


def findSphericity(fn):
    f = open(fn, 'r')
    sph = 0.
    for line in f.readlines():
        if 'Sphericity is ' in line:
            sph = float(line.split()[2])
    f.close()

    return sph
