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
1. Write from base classes to cryoEF specific files
2. Read from cryoEF files to base classes
"""

import os
from numpy import rad2deg, deg2rad
from numpy.linalg import inv

from pyworkflow.utils import Environ
from pyworkflow.utils.path import join

CRYOEF_HOME_VAR = 'CRYOEF_HOME'


def getEnviron():
    """ Setup the environment variables needed to launch cryoEF. """
    environ = Environ(os.environ)
    CRYOEF_HOME = os.environ[('%s' % CRYOEF_HOME_VAR)]

    environ.update({
        'PATH': join(CRYOEF_HOME, 'bin'),
    }, position=Environ.BEGIN)
    return environ


_environ = getEnviron()

SUPPORTED_VERSIONS = ['1.1.0']


def getVersion():
    path = os.environ[CRYOEF_HOME_VAR]
    for v in SUPPORTED_VERSIONS:
        if v in path:
            return v
    return ''


def parseOutput(filename):
    """ Retrieve efficiency, mean PSF res, stdev, worst and best PSF res
    from the output log file of the cryoEF execution.
    """
    result = []
    if os.path.exists(filename):
        f = open(filename)
        for line in f:
            if 'Efficiency:' in line:
                result.append(float(line.split()[1]))
            if 'Mean PSF resolution:' in line:
                result.append(float(line.split()[3]))
            if 'Standard deviation:' in line:
                result.append(float(line.split()[2]))
            if 'Worst PSF resolution:' in line:
                result.append(float(line.split()[3]))
            if 'Best PSF resolution:' in line:
                result.append(float(line.split()[3]))
        f.close()
    return result


def iterAngles(fn):
    f = open(fn)
    for line in f:
        rot, tilt = map(float, line.split())
        yield rot, tilt
    f.close()


def writeAnglesFn(img, fn):
    # get alignment parameters for each particle
    shifts, angles = geometryFromMatrix(img.getTransform().getMatrix())
    rot, tilt, _ = angles
    fn.write("%0.6f %0.6f\n" % (rot, tilt))


def geometryFromMatrix(matrix, inverseTransform=True):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix

    if inverseTransform:
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles

