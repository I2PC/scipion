# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
from os.path import join, dirname, abspath, isdir
import re

from pyworkflow.object import String
from pyworkflow.utils import runJob, Environ
from pyworkflow.utils.path import replaceBaseExt
from pyworkflow.em.data import EMObject

END_HEADER = 'END BATCH HEADER'
PATH = abspath(dirname(__file__))
SCRIPTS_DIR = 'scripts'

# Regular expression for parsing vars in script header
# Match strings of the type "# key=value # some comment"
REGEX_KEYVALUE = re.compile('(?P<var>\[?[a-zA-Z0-9_-]+\]?)=(?P<value>"\S+")(?P<suffix>\s+.*)')


def getEnviron():
    """ Load the environment variables needed for Imagic.
    IMAGIC_ROOT is set to IMAGIC_DIR
    MPI-related vars are set if MPI_BIN is not empty is scipion.config
    IMAGIC_BATCH is needed for batch files to work.
    """
    env = Environ(os.environ)
    imagicdir = env.get('IMAGIC_DIR', None)  # Scipion definition

    if imagicdir is None:
        print "ERROR: Missing IMAGIC_DIR variable in scipion.conf file"

    mpi_dir = imagicdir + '/openmpi'

    if isdir(mpi_dir):
        env.update({'MPIHOME': mpi_dir,
                    'MPIBIN': mpi_dir + '/bin',
                    'OPAL_PREFIX': mpi_dir
                    })
        env.set('PATH', mpi_dir + '/bin', env.BEGIN)

    else:
        print "ERROR: IMAGIC_ROOT directory (", imagicdir, ") does not contain openmpi folder"

    env.update({'IMAGIC_ROOT': imagicdir,
                'IMAGIC_BATCH': "1",
                'FFTW_IPATH': imagicdir + '/fftw/include',
                'FFTW_LPATH': imagicdir + '/fftw/lib',
                'FFTW_LIB': 'fftw3f'
                })

    env.set('LD_LIBRARY_PATH', imagicdir + '/fftw/lib', env.BEGIN)
    env.set('LD_LIBRARY_PATH', imagicdir + '/lib', env.BEGIN)

    return env


environment = getEnviron()


def _getFile(*paths):
    return join(PATH, *paths)


def getScript(*paths):
    return _getFile(SCRIPTS_DIR, *paths)


def __substituteVar(match, paramsDict, lineTemplate):
    if match and match.groupdict()['var'] in paramsDict:
        d = match.groupdict()
        d['value'] = paramsDict[d['var']]
        return lineTemplate % d
    return None


def writeScript(inputScript, outputScript, paramsDict):
    """ Create a new Imagic script by substituting
    params in the input 'paramsDict'.
    """
    fIn = open(getScript(inputScript), 'r')
    fOut = open(outputScript, 'w')
    inHeader = True  # After the end of header, no more value replacement

    for i, line in enumerate(fIn):
        if END_HEADER in line:
            inHeader = False
        if inHeader:
            try:
                newLine = __substituteVar(REGEX_KEYVALUE.match(line), paramsDict,
                                          "%(var)s=%(value)s%(suffix)s\n")
                if newLine:
                    line = newLine
            except Exception, ex:
                print ex, "on line (%d): %s" % (i + 1, line)
                raise ex
        fOut.write(line)
    fIn.close()
    fOut.close()


def runTemplate(inputScript, paramsDict, log=None, cwd=None):
    """ This function will create a valid Imagic script
    by copying the template and replacing the values in dictionary.
    After the new file is read, the Imagic interpreter is invoked.
    Usually the execution should be done where the results will
    be left.
    """
    outputScript = replaceBaseExt(inputScript, 'b')

    if cwd is not None:
        outputScript = join(cwd, outputScript)

    # First write the script from the template with the substitutions
    writeScript(inputScript, outputScript, paramsDict)
    # Then proceed to run the script
    runScript(outputScript, log, cwd)


def runScript(inputScript, log=None, cwd=None):
    args = " %s" % inputScript
    shellPath = '/bin/bash'
    runJob(log, shellPath, args, env=dict(environment), cwd=cwd)


class ImagicPltFile(object):
    """ Handler class to read/write imagic plt file. """
    pass

    def __init__(self, filename, mode='r'):
        self._file = open(filename, mode)
        self._count = 0

    def Values(self, classNum):
        values = []
        for line in self._file:
            line = line.strip()
            fields = line.split()
            if int(float(fields[1])) == classNum:
                values.append(int(float(fields[0])))
        return values

    def close(self):
        self._file.close()


class MsaFile(EMObject):
    """ This is a container of files produced by MSA Imagic protocol.
    It is possible to use *.lis and *.plt files.
    """

    def __init__(self, **args):
        EMObject.__init__(self, **args)

        self.filename = String()

    def getFileName(self):
        return self.filename.get()


class EigFile(EMObject):
    """ This is a container of files produced by MSA Imagic protocol.
    It is possible to use eigen_img.img files.
    """

    def __init__(self, **args):
        EMObject.__init__(self, **args)

        self.filename = String()

    def getFileName(self):
        return self.filename.get()
