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
from os.path import join, exists
import platform

from urllib2 import urlopen
from subprocess import call
import tarfile

# OS boolean vars
MACOSX = (platform.system() == 'Darwin')
WINDOWS = (platform.system() == 'Windows')
LINUX = (platform.system() == 'Linux')

SCIPION_HOME = os.environ['SCIPION_HOME']
SCIPION_SOFTWARE_PATH = join(SCIPION_HOME, 
                                    'software')
SCIPION_INSTALL_PATH = join(SCIPION_SOFTWARE_PATH,
                                    'install')
SCIPION_INSTALL_LOG = join(SCIPION_SOFTWARE_PATH,
                                   'log',
                                   'scons.log')
SCONS_VERSION = 'scons-2.3.1'
SCIPION_DIRS = ['SCIPION_DATA', 'SCIPION_LOGS', 'SCIPION_TESTS', 'SCIPION_USER_DATA', 'SCIPION_TMP']


    
def downloadScons():
    """ Download the scons tgz file and extract it. """
    SCONS_URL = "http://scipionwiki.cnb.csic.es/files/scipion/software/python/"
    SCONS_VERSION = 'scons-2.3.1.tgz'
    INSTALL_PATH = join(os.environ['SCIPION_HOME'], 'software', 'install')
    print "Downloading scons from " + SCONS_URL + SCONS_VERSION
    try:
        data = urlopen("%s/%s" % (SCONS_URL, SCONS_VERSION)).read()
        open(join(INSTALL_PATH, SCONS_VERSION), 'w').write(data)
    except Exception as e:
        print "\tError downloading %s (%s)" % (SCONS_VERSION, e)
        print "destination: %s" % INSTALL_PATH
        print "URL: %s/%s" % (SCONS_URL, SCONS_VERSION)
        return -1
    print "Unpacking scons"
    sourceTar = tarfile.open(join(INSTALL_PATH, SCONS_VERSION), 'r')
    sourceTar.extractall(INSTALL_PATH)
    sourceTar.close()
    os.remove(join(INSTALL_PATH, SCONS_VERSION))


def build(args):
    print "Scipion Home in : ", SCIPION_HOME

    if not '--purge' in args:
        if not exists(join(SCIPION_INSTALL_PATH, SCONS_VERSION)):
            # Download and untar Scons
            downloadScons()
            #Compile scons
            setupArgs = [["install", "--prefix=%s" % SCIPION_SOFTWARE_PATH]]
            if exists(SCIPION_INSTALL_LOG):
                os.remove(SCIPION_INSTALL_LOG)
            for sargs in setupArgs:
                command = ['python', join(SCIPION_INSTALL_PATH, SCONS_VERSION, 'setup.py')] + sargs
                print 'Executing %s on scons' % " ".join(sargs)
                with open(SCIPION_INSTALL_LOG, 'a') as logFile:
                    r = call(command, cwd=join(SCIPION_INSTALL_PATH, SCONS_VERSION),
                             stdout=logFile, stderr=logFile, env=os.environ)
                if r != 0:
                    return r

    return call('software/bin/scons %s | tee -a %s' % (' '.join(args), SCIPION_INSTALL_LOG),
                shell=True, env=os.environ)



if __name__ == '__main__':
    returncode = build(sys.argv[1:])
    # This script will exit with the same exit code as scons did
    sys.exit(returncode)
