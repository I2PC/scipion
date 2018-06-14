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

from pyworkflow.utils import Environ, replaceBaseExt
from pyworkflow.em.convert import ImageHandler

# we declarate global constants to multiple usage
LOCSCALE_HOME_VAR = 'LOCSCALE_HOME'
EMAN2DIR_VAR = 'EMAN2DIR'
LOCSCALE_HOME = os.environ[LOCSCALE_HOME_VAR]
EMAN2DIR = os.environ[EMAN2DIR_VAR]

def getVersion():
    version = ''
    for v in getSupportedVersions():
        if v in LOCSCALE_HOME:
            version = v
    return version

def getSupportedVersions():
    return ['0.1']

def getSupportedEmanVersions():
    """ LocScale needs eman to work.
    """
    return ['2.11', '2.12', '2.2']

def getEmanVersion():
    """ Returns a valid eman version installed or an empty string.
    """
    emanVersion = ''
    if os.path.exists(EMAN2DIR):
        for supV in ['eman-'+v for v in getSupportedEmanVersions()]:
            if supV in EMAN2DIR:
                emanVersion = supV
    return emanVersion

def validateEmanVersion(errors):
    """ Validate if eman version is set properly according
     to installed version and the one set in the config file.
     Params:
        protocol: the input protocol calling to validate
        errors: a list that will be used to add the error message.
    """
    if getEmanVersion() == '':
        errors.append('Eman (v: %s) is needed to execute this protocol'
                      % ', '.join(getSupportedEmanVersions()))
    return errors

def setEmanEnviron():
    """ Setup the environment variables needed to launch Eman and 
        use its modules.
    """
    env = Environ(os.environ)
    pathList = [os.path.join(EMAN2DIR, d)
                for d in ['lib', 'bin', 'extlib/site-packages']]

    env.update({'PATH': os.path.join(EMAN2DIR, 'bin'),
                'LD_LIBRARY_PATH': os.pathsep.join(pathList),
                'PYTHONPATH': os.pathsep.join(pathList),
                'EMAN_PYTHON': os.path.join(EMAN2DIR, 'Python/bin/python')}, 
                position=Environ.BEGIN)

    os.environ.update(env)

def getEmanPythonProgram(program):
    if not 'EMAN_PYTHON' in os.environ:
        setEmanEnviron()

    # locscale scripts are in $LOCSCALE_HOME/source
    program = os.path.join(LOCSCALE_HOME, 'source', program)
    python = os.environ['EMAN_PYTHON']

    return python, program

def convertBinaryVol(vol, outputDir):
    """ Convert binary volume to a mrc format.
    Params:
        vol: input volume object to be converted.
        outputDir: where to put the converted file(s)
    Return:
        new file name of the volume (converted or not).
    """
    ih = ImageHandler()

    def convertToMrc(fn):
        """ Convert from a format that is not read by Relion
        to mrc format.
        """
        newFn = os.path.join(outputDir, replaceBaseExt(fn, 'mrc'))
        ih.convert(fn, newFn)
        return newFn

    volFn = vol.getFileName()
    if ':' in volFn:
        volFn = volFn.split(':')[0]
    
    if not volFn.endswith('.mrc'):
        volFn = convertToMrc(volFn)

    return volFn
    