#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This module serve to download the needed software packages.
"""
import os, sys


PACKAGES_URL = "http://scipionwiki.cnb.csic.es/files/scipion/software/em/"
PACKAGES_DICT = {
                 #'xmipp': "xmipp.tgz", 
                 'bsoft': 'bsoft1_8_8_Fedora_12.tgz',
                 'ctffind': 'ctffind_V3.5.tgz',
                 'eman': "eman2.pre2-1.linux64.tar.gz",
                 'frealign': 'frealign_v9.07.tgz',                 
                 'relion': "relion-1.2.tgz",
                 'spider': "spider-web-21.13.tgz"                 
                 }


def download(package):
    os.system("wget " + PACKAGES_URL + PACKAGES_DICT[package])
    
def downloadAll():
    for package in PACKAGES_DICT:
        download(package)


if __name__ == '__main__':
    # Enter to EM packages directory
    os.chdir(os.path.join(os.environ['SCIPION_HOME'], 'software', 'em'))
    
    n = len(sys.argv)
    if n > 1:
        package = sys.argv[1]
        download(package)
    else:
        downloadAll()
