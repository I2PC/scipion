#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco     (jgomez@cnb.csic.es)
# *              I. Foche Perez         (ifoche@cnb.csic.es)
# *              J. Burguet Castell     (jburguet@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
Launch the build process for a fresh installation of Scipion.

It will first download and install SCons, and then let it execute the
SConscript, which will install all the libraries, python modules and
em packages that are needed.
"""

from __future__ import with_statement
# Needed in case we run with Python 2.5

import sys
import os
from os.path import join, exists, dirname, abspath

try:
    from urllib2 import urlopen
    from StringIO import StringIO
except ImportError:
    # This will make it work with Python 3, once SCons also supports it!
    from urllib.request import urlopen
    from io import StringIO
import subprocess
import tarfile


SCIPION_HOME = os.environ.get('SCIPION_HOME', dirname(dirname(abspath(__file__))))
SOFTWARE = join(SCIPION_HOME, 'software')
INSTALL = join(SOFTWARE, 'install')
LOGFILE = join(SOFTWARE, 'log', 'scons.log')

SCONS = 'scons-2.3.2'

if sys.version_info[:2] >= (3, 0):
    sys.exit('Running with Python 3, but Python 2 needed for the installation.')
# Once SCons works with Python 3, we can update and take this out.


def build(args):
    """ Download SCons if needed and run it (with args) """

    # Go to SCIPION_HOME if we are not there.
    os.chdir(SCIPION_HOME)

    # First make sure that we have SCons installed.
    if not exists(join(SOFTWARE, 'bin', 'scons')):
        # Download and extract SCons.
        url = 'http://scipionwiki.cnb.csic.es/files/scipion/software/python/%s.tgz' % SCONS
        sys.stdout.write(
           'Downloading, unpacking & installing scons from %s ...\n' % url)
        tarfile.open(fileobj=StringIO(urlopen(url).read())).extractall(INSTALL)

        # Install it.
        command = [sys.executable, 'setup.py', 'install', '--prefix=%s' % SOFTWARE]
        sys.stdout.write('Executing: %s\n' % ' '.join(command))
        with open(LOGFILE, 'w') as logFile:
            r = subprocess.call(command, cwd=join(INSTALL, SCONS),
                                stdout=logFile, stderr=logFile)
        if r != 0:
            return r

    # Call SCons and show the output on the screen and logfile.
    proc = subprocess.Popen([sys.executable, 'software/bin/scons'] + args,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with open(LOGFILE, 'a') as logFile:
        while True:
            ret = proc.poll()
            if ret is not None:
                return ret
            line = proc.stdout.readline()
            sys.stdout.write(line)
            logFile.write(line)
            logFile.flush()



if __name__ == '__main__':
    sys.exit(build(sys.argv[1:]))
