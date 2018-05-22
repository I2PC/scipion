# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
"""
This module contains converter functions that will serve to:
1. define ccp4 environ
TODO:
2. Read/Write CCP4 specific files
"""

import collections
import os
from pyworkflow.em.convert import ImageHandler
import pyworkflow.utils as pwutils
import struct
import getpass
from math import isnan


def getEnviron(first=True):
    environ = pwutils.Environ(os.environ)
    pos = pwutils.Environ.BEGIN if first else pwutils.Environ.END
    _ccp4_home = os.environ['CCP4_HOME']
    _ccp4_master, _dir = os.path.split(_ccp4_home)
    _username = getpass.getuser()
    environ.update({
            'PATH': os.path.join(_ccp4_home, 'bin'),
            'LD_LIBRARY_PATH': os.path.join(_ccp4_home, 'lib'),
            # # CCP4_MASTER is the location of the top-level directory
            # # containing ccp4-N.N.N.
            # export CCP4_MASTER=/home/roberto
            # export CCP4=$CCP4_MASTER/ccp4-6.5
            # alias xtal='pushd $CCP4_MASTER>/dev/null'
            'CCP4_MASTER': _ccp4_master,
            # alias ccp4='pushd $CCP4>/dev/null'
            'CCP4': _ccp4_home,
            # # CCP4_SCR: a per-user directory for run-time-generated scratch
            # # files.
            # export CCP4_SCR=/tmp/`whoami`
            'CCP4_SCR': os.path.join("/tmp", _username),
            # # This variable is set to ensure that the logfile output from programs
            # # compiled with Gfortran is in the correct order.
            # export GFORTRAN_UNBUFFERED_PRECONNECTED=Y
            'GFORTRAN_UNBUFFERED_PRECONNECTED':"Y",
            # # CBIN: location of the executables -- must be on your path
            # # (see below)
            # export CBIN=$CCP4/bin
            # alias cbin='pushd $CBIN>/dev/null'
            'CBIN': os.path.join(_ccp4_home, 'bin'),
            # # CLIB: location of (binary) library files such as libccp4.a
            # # and libccp4.so
            # export CLIB=$CCP4/lib
            # alias clib='pushd $CLIB>/dev/null'
            'CLIB': os.path.join(_ccp4_home, 'lib'),
            # # CLIBD: platform-independent data files
            # export CLIBD=$CCP4/lib/data
            # alias clibd='pushd $CLIBD>/dev/null'
            'CLIBD': os.path.join(_ccp4_home, 'lib','data'),
            # # CETC: executable scripts (NOT configuration files)
            # export CETC=$CCP4/etc
            # alias cetc='pushd $CETC>/dev/null'
            'CETC': os.path.join(_ccp4_home, 'etc'),
            # # CINCL: headers and two *.def files for handling
            # # "logical names" in CCP4
            # export CINCL=$CCP4/include
            # alias cincl='pushd $CINCL>/dev/null'
            'CINCL': os.path.join(_ccp4_home, 'include'),
            # # CHTML: html documentation
            # export CHTML=$CCP4/html
            # alias chtml='pushd $CHTML>/dev/null'
            'CHTML': os.path.join(_ccp4_home, 'html'),
            # # CEXAM: examples and some tests
            # export CEXAM=$CCP4/examples
            # alias cexam='pushd $CEXAM>/dev/null'
            'CEXAM': os.path.join(_ccp4_home, 'examples'),
            # # CCP4I_TOP: the top directory of the interface
            # export CCP4I_TOP=$CCP4/share/ccp4i
            # # source code directories
            # #export CLIBS=$CCP4/lib/libccp4
            # #alias clibs='pushd $CLIBS>/dev/null'
            # #export CPROG=$CCP4/src
            # #alias cprog='pushd $CPROG>/dev/null'
            'CCP4I_TOP': os.path.join(_ccp4_home, 'share','ccp4i'),
            # # MMCIFDIC: platform-dependent (not in $CLIBD) data file for
            # # the ccif library
            # export MMCIFDIC=$CLIB/ccp4/cif_mmdic.lib
            'MMCIFDIC': os.path.join(_ccp4_home, 'lib', 'cif_mmdic.lib'),
            # # CLIBD_MON: dictionary files for REFMAC5 (keep trailing /)
            # export CLIBD_MON=$CCP4/lib/data/monomers/
            'CLIBD_MON': os.path.join(_ccp4_home, 'lib', 'data', 'monomers'),
            # # CRANK: location of Crank automation suite within ccp4i
            # export CRANK=$CCP4I_TOP/crank
            'CRANK': os.path.join(_ccp4_home, 'crank'),
            # # CCP4_HELPDIR: location of the VMS-style help file used
            # # by (ip)mosflm
            # export CCP4_HELPDIR=$CCP4/help/            # NB trailing /
            'CCP4_HELPDIR': os.path.join(_ccp4_home, 'help'),
            }, position=pos)
    return environ

def runPhenixProgram(program, args="", extraEnvDict=None):
    """ Internal shortcut function to launch a Phenix program. """
    env = getEnviron()
    if extraEnvDict is not None:
        env.update(extraEnvDict)
    pwutils.runJob(None, program, args, env=env)


def getProgram(progName):
    """ Return the program binary that will be used. """
    if 'PHENIX_HOME' not in os.environ:
        return None
    return os.path.join(os.environ['PHENIX_HOME'], 'bin',
                        os.path.basename(progName))
