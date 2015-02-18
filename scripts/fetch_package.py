#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     I. Foche Perez         (ifoche@cnb.csic.es)
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
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************

"""
Downloads the given EM Package to allow its installation.

It will download the package from the Scipion server or an specific git
repository.
"""

import sys
import tarfile
from io import BytesIO
import os
try:
    from urllib2 import urlopen, urlretrieve
    from StringIO import StringIO
except ImportError:
    # This will make it work with Python 3, once SCons also supports it!
    from urllib import urlopen, urlretrieve
    from io import StringIO
    
# BASE PACKAGE URL
_BASE_URL = 'http://metamagical.org/scipion_soft'

# RETURN STATES
_RETURN_OK = 0
_RETURN_BAD_ARG_LENGTH = 1
_RETURN_BAD_DOWNLOAD = 2
_RETURN_BAD_UNTAR = 3

_RETURN_DICT = {
                 _RETURN_OK: 'OK', # 0 - OK
                 _RETURN_BAD_ARG_LENGTH: 'Bad argument number. This script accepts only 3 arguments (EMPackage, PackageHome, TarFile)',       # 1 - BAD
                 _RETURN_BAD_DOWNLOAD: 'Error during download process',
                 _RETURN_BAD_UNTAR: 'Error during untar process'
}


# AUXILIAR FUNCS

# Exit function, just to print the error message before exiting
def exitScript(returnCode):
    if returnCode != 0:
        print >> sys.stderr, 'Script finished with error %d (%s)' % (returnCode, _RETURN_DICT.get(returnCode))
    sys.exit(returnCode)

# SCRIPT ITSELF
if len(sys.argv) != 4:
    exitScript(_RETURN_BAD_ARG_LENGTH)

_PACKAGE_NAME = sys.argv[1]
_PACKAGE_HOME = sys.argv[2]
_PACKAGE_TAR = sys.argv[3]

print 'name ',_PACKAGE_NAME
print 'home ',_PACKAGE_HOME
print 'tar', _PACKAGE_TAR

url = '%s/em/%s' % (_BASE_URL, _PACKAGE_TAR)
dstFolder = 'software/em'

print 'url %s' % url
print 'dstFolder %s' % dstFolder

#try:
if not os.path.exists('software/tmp/%s' % _PACKAGE_TAR):
    #Download
    u = urlopen(url)
    file_name = 'software/tmp/%s' % _PACKAGE_TAR
    f = open(file_name, 'w+')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print "Downloading: %s Bytes: %s" % (file_name, file_size)
    
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8)*(len(status)+1)
        print status,
    f.close()
#except:
#    exitScript(_RETURN_BAD_DOWNLOAD)
#try:
print 'dstFolder %s' % dstFolder
print os.getcwd()
packageTar = tarfile.open('software/tmp/%s' % _PACKAGE_TAR, 'r')
packageTar.extractall(dstFolder)
#except:
#    exitScript(_RETURN_BAD_UNTAR)

exitScript(_RETURN_OK)
