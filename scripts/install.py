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
import os, sys, platform

from pyworkflow import SETTINGS
from pyworkflow.utils.path import makePath, copyFile, join
from pyworkflow.apps.config import writeDefaults
from urllib2 import urlopen
import subprocess
import tarfile

# OS boolean vars
MACOSX = platform.system() == 'Darwin'
WINDOWS = platform.system() == 'Windows'
LINUX = platform.system() == 'Linux'


def updateProjectSettings():
    # Update the settings to all existing projects
    from pyworkflow.manager import Manager
    manager = Manager()
    projects = manager.listProjects()
    
    for p in projects:
        proj = manager.loadProject(p.getName())
        projSettings = proj.settingsPath
        print "Copying settings to: ", join(p.getName(), projSettings)
        copyFile(SETTINGS, projSettings)


def downloadScons():
    """ Download the scons tgz file and extract it. """
    SCONS_URL = "http://scipionwiki.cnb.csic.es/files/scipion/software/python/"
    SCONS_VERSION = 'scons-2.3.1.tgz'
    INSTALL_PATH = os.path.join(os.environ['SCIPION_HOME'],
                                'software',
                                'install')
    print "Downloading scons from " + SCONS_URL + SCONS_VERSION
    try:
        data = urlopen("%s/%s" % (SCONS_URL, SCONS_VERSION)).read()
        open(os.path.join(INSTALL_PATH, SCONS_VERSION), 'w').write(data)
    except Exception as e:
        print "\tError downloading %s (%s)" % (SCONS_VERSION, e)
        print "destination: %s" % INSTALL_PATH
        print "URL: %s/%s" % (SCONS_URL, SCONS_VERSION)
        return -1
    print "Unpacking scons"
    sourceTar = tarfile.open(os.path.join(INSTALL_PATH, SCONS_VERSION), 'r')
    sourceTar.extractall(INSTALL_PATH)
    sourceTar.close()
    os.remove(os.path.join(INSTALL_PATH, SCONS_VERSION))


if __name__ == '__main__':
#    parser = argsparse.ArgumentParser(description=__doc__)
#    add = parser.add_argument # shortcut
#    add('--update', action='store_true',
#        help='update local installation')
#    add(name='name', metavar='METAVAR', nargs='*', help='help')
#    args = parser.parse_args()
    
    SCIPION_HOME = os.environ['SCIPION_HOME']
    SCIPION_SOFTWARE_PATH = os.path.join(SCIPION_HOME, 
                                        'software')
    SCIPION_INSTALL_PATH = os.path.join(SCIPION_SOFTWARE_PATH,
                                        'install')
    SCIPION_INSTALL_LOG = os.path.join(SCIPION_SOFTWARE_PATH,
                                       'log',
                                       'scons.log')
    SCONS_VERSION = 'scons-2.3.1'
    SCIPION_DIRS = ['SCIPION_DATA', 'SCIPION_LOGS', 'SCIPION_TESTS', 'SCIPION_USER_DATA', 'SCIPION_TMP']

    #necessary to properly build scons
    myenv = {}
    #myenv = os.environ.copy()

    print "Scipion Home in : ", SCIPION_HOME

    if not '--purge' in sys.argv:
        # Create DATA folders
        for d in SCIPION_DIRS:
            path = os.environ[d]
            if not os.path.exists(path):
                print "  creating %s folder: %s" % (d, path)
                makePath(path)
        
         # Write default configurations
#         writeDefaults()
#         
#         updateProjectSettings()
        if not os.path.exists(os.path.join(SCIPION_INSTALL_PATH, SCONS_VERSION)):
             # Download and untar Scons
            downloadScons()
            #Compile scons
            setupArgs=[["clean"],
                       ["build"], 
                       ["install", "--prefix=%s" % SCIPION_SOFTWARE_PATH]]
            if os.path.exists(SCIPION_INSTALL_LOG):
                os.remove(SCIPION_INSTALL_LOG)
            logFile = open(SCIPION_INSTALL_LOG, 'w+').close()
            for arg in setupArgs:
                logFile = open(SCIPION_INSTALL_LOG, 'a')
                command = ['python', os.path.join(SCIPION_INSTALL_PATH, SCONS_VERSION, 'setup.py')] + arg
                
                print 'Executing %s on scons' % (arg)
                scons = subprocess.Popen(command, 
                                         cwd=os.path.join(SCIPION_INSTALL_PATH, SCONS_VERSION),
                                         stdout=logFile,
                                         stderr=logFile,
                                         env=myenv)
                scons.wait()
                output = scons.communicate()[0]
                logFile.close()
                if scons.returncode is not 0:
                    sys.exit(scons.returncode)

    #Setting environmental vars before calling scons
    myenv['PATH'] = '%s:%s' % (os.path.join(SCIPION_SOFTWARE_PATH, 'bin'), os.environ['PATH'])
    DYN = 'LD_LIBRARY_PATH'
    if MACOSX:
        DYN = 'DYLD_FALLBACK_LIBRARY_PATH'
    if DYN in myenv:
        myenv[DYN] = '%s:%s' % (os.path.join(SCIPION_SOFTWARE_PATH, 'lib'), myenv[DYN])
        myenv[DYN] = '%s:%s' % (os.path.join(SCIPION_SOFTWARE_PATH, 'lib64'), myenv[DYN])

    logFile = open(SCIPION_INSTALL_LOG, 'a')
    args = ' '.join(sys.argv[1:])
    install = subprocess.Popen(['scons']+sys.argv[1:],
                               cwd=SCIPION_HOME, 
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               env=myenv)
    for line in install.stdout:
        install.stdout.flush()
        sys.stdout.flush()
        sys.stdout.write(line)
        logFile.write(line)
    install.wait()
    output = install.communicate()[0]
    logFile.close()
    
    # This script will exit with the same exit code as scons did
    sys.exit(install.returncode)


