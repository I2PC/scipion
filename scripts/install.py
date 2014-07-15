#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco     (jgomez@cnb.csic.es)
# *              I. Foche Perez (ifoche@cnb.csic.es)
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
"""
This script will generate the pw.bashrc and pw.cshrc file to include
"""

import sys
import os
from os.path import join, exists, dirname, abspath
import platform

from urllib2 import urlopen
from StringIO import StringIO
import subprocess
import tarfile

# OS boolean vars
MACOSX = (platform.system() == 'Darwin')
WINDOWS = (platform.system() == 'Windows')
LINUX = (platform.system() == 'Linux')

SCIPION_HOME = os.environ.get('SCIPION_HOME', dirname(dirname(abspath(__file__))))
SCIPION_SOFTWARE_PATH = join(SCIPION_HOME, 'software')
SCIPION_INSTALL_PATH = join(SCIPION_SOFTWARE_PATH, 'install')
SCIPION_INSTALL_LOG = join(SCIPION_SOFTWARE_PATH, 'log', 'scons.log')

SCONS = 'scons-2.3.1'

SCIPION_DIRS = ['SCIPION_DATA', 'SCIPION_LOGS', 'SCIPION_TESTS', 'SCIPION_USER_DATA', 'SCIPION_TMP']



def build(args):

    # First make sure that we have SCons installed.
    if not exists(join(SCIPION_INSTALL_PATH, SCONS)):
        # Download and extract SCons.
        url = 'http://scipionwiki.cnb.csic.es/files/scipion/software/python/%s.tgz' % SCONS
        print "Downloading, unpacking & installing scons from %s..." % url
        tarfile.open(fileobj=StringIO(urlopen(url).read())).extractall(SCIPION_INSTALL_PATH)

        # Install it.
        command = ['python', 'setup.py', 'install',
                   '--prefix=%s' % SCIPION_SOFTWARE_PATH]
        print 'Executing: %s' % ' '.join(command)
        with open(SCIPION_INSTALL_LOG, 'w') as logFile:
            r = subprocess.call(command, cwd=join(SCIPION_INSTALL_PATH, SCONS),
                                stdout=logFile, stderr=logFile, env=os.environ)
        if r != 0:
            return r

    # Call SCons and show the output on the screen and logfile.
    proc = subprocess.Popen(['software/bin/scons'] + args,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with open(SCIPION_INSTALL_LOG, 'a') as logFile:
        while True:
            ret = proc.poll()
            if ret is not None:
                return ret
            line = proc.stdout.readline()
            sys.stdout.write(line)
            logFile.write(line)
            logFile.flush()



if __name__ == '__main__':
    returncode = build(sys.argv[1:])
    # This script will exit with the same exit code as scons did
    sys.exit(returncode)
