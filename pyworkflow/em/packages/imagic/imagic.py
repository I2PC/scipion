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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import join, dirname, abspath, isdir, exists
import re

from pyworkflow.utils import runJob, Environ
from pyworkflow.utils.path import replaceBaseExt

END_HEADER = 'END BATCH HEADER'
PATH = abspath(dirname(__file__))
SCRIPTS_DIR = 'scripts'

# Regular expression for parsing vars in script header
# Match strings of the type "# key=value # some comment"
REGEX_KEYVALUE = re.compile('(?P<var>\[?[a-zA-Z0-9_-]+\]?)=(?P<value>"\S+")(?P<suffix>\s+.*)')


def getEnviron():
    """ Load the environment variables needed for Imagic.
    IMAGIC_ROOT is set to IMAGIC_DIR
    MPI-related vars are set if IMAGIC_DIR/openmpi path exists
    IMAGIC_BATCH is needed for batch files to work.
    """
    env = Environ(os.environ)
    imagicdir = env.get('IMAGIC_DIR', None)  # Scipion definition

    if imagicdir is None or not isdir(imagicdir):
        print "ERROR: Missing IMAGIC_DIR variable in scipion.conf file or path does not exist."

    else:
        env.update({'IMAGIC_ROOT': imagicdir,
                    'IMAGIC_BATCH': "1",
                    'FFTW_IPATH': imagicdir + '/fftw/include',
                    'FFTW_LPATH': imagicdir + '/fftw/lib',
                    'FFTW_LIB': 'fftw3f'
                    })

        env.set('LD_LIBRARY_PATH', imagicdir + '/fftw/lib', env.BEGIN)
        env.set('LD_LIBRARY_PATH', imagicdir + '/lib', env.BEGIN)

    mpidir = imagicdir + '/openmpi'

    if isdir(mpidir):
        env.update({'MPIHOME': mpidir,
                    'MPIBIN': mpidir + '/bin',
                    'OPAL_PREFIX': mpidir
                    })
        env.set('PATH', mpidir + '/bin', env.BEGIN)

    else:
        print "Warning: IMAGIC_ROOT directory (", imagicdir, ") does not contain openmpi folder.\n", \
              "No MPI support will be enabled."

    return env


def getVersion():
    imagicdir = os.environ['IMAGIC_DIR']
    for v in getSupportedVersions():
        versionFile = join(imagicdir, 'version_' + v)
        if exists(versionFile):
            return v
    return ''


def getSupportedVersions():
    return ['110308', '160418']


def _getFile(*paths):
    return join(PATH, *paths)


def getScript(*paths):
    return _getFile(SCRIPTS_DIR, getVersion(), *paths)


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
    runJob(log, shellPath, args, env=getEnviron(), cwd=cwd)


class ImagicPltFile(object):
    """ Handler class to read/write imagic plt file. """

    def __init__(self, filename):
        self._filename = filename

    def iterRows(self):
        f = open(self._filename)

        for line in f:
            line = line.strip()
            fields = map(float, line.split())
            # rows contains tuples (float) of
            # image_number, class_number
            yield int(fields[0]), fields[1]

        f.close()


class ImagicLisFile(object):
    """ Handler class to read imagic lis file. """

    def __init__(self, filename, clsNum):
        self._filename = filename
        self._clsNum = clsNum
        self.varianceDict = {}
        self.quality1Dict = {}
        self.quality2Dict = {}

    def parseFile(self, source):
        """ Collect info for all classes into a common list"""
        paramsList = []
        read = False
        for line in source:
            line = line.strip().lower()
            if line.startswith("classes sorted by intra-class variance"):
                read = True
                # read file starting from this line
            else:
                if read and '#' not in line and ':' not in line:
                    fields = line.split()
                    if len(fields) == 5:
                        # we have 5 columns in line
                        classId = int(fields[2])
                        value = float(fields[1])
                        paramsList.append([classId, value])
                continue

        return paramsList

    def getParams(self):
        f = open(self._filename)
        cls = self._clsNum
        paramsList = self.parseFile(f)
        for param in paramsList[0:cls]:
            self.varianceDict[param[0]] = param[1]
        for param in paramsList[cls:(2*cls)]:
            self.quality1Dict[param[0]] = param[1]
        for param in paramsList[(2*cls):]:
            self.quality2Dict[param[0]] = param[1]
        f.close()

        return self.varianceDict, self.quality1Dict, self.quality2Dict
