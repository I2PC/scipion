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

import os
import pyworkflow.em as em
import pyworkflow.utils as pwutils


#
chimeraPdbTemplateFileName= "scipionOut%04d.pdb"
chimeraScriptFileName = "cootScript.py"


def getEnviron(ccp4First=True):
    environ = pwutils.Environ(os.environ)
    pos = pwutils.Environ.BEGIN if ccp4First else pwutils.Environ.END
    environ.update({
            'PATH': os.path.join(os.environ['CCP4_HOME'], 'bin'),
            'LD_LIBRARY_PATH': os.path.join(os.environ['CCP4_HOME'], 'lib'),
            }, position=pos)
    return environ

def runCCP4Program(program, args=""):
    """ Internal shortcut function to launch a CCP4 program. """
    env=getEnviron()
    #env.update(_envDict)
    pwutils.runJob(None, program, args, env=env)

def adapBinFileToCCP4(inFileName,outFileName):
    if inFileName.endswith('.mrc'):
        pwutils.createLink(inFileName, outFileName)
    else:
        em.ImageHandler().convert(inFileName, outFileName)

def getProgram(progName):
    """ Return the program binary that will be used. """
    if (not 'CHIMERA_HOME' in os.environ):
        return None
    return os.path.join(os.environ['CHIMERA_HOME'], 'bin',
                           os.path.basename(progName))

