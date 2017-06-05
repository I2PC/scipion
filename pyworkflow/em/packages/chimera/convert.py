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
from pyworkflow.em.viewer import getChimeraEnviron
from pyworkflow.em.constants import SYM_I222, SYM_I222r, SYM_In25, SYM_In25r

#
chimeraPdbTemplateFileName= "scipionOut%04d.pdb"
chimeraMapTemplateFileName= "scipionOut%04d.mrc"
chimeraScriptFileName = "chimeraScript.py"

symMapperScipionchimera = {}
symMapperScipionchimera[SYM_I222]="222"
symMapperScipionchimera[SYM_I222r]="222r"
symMapperScipionchimera[SYM_In25]="n25"
symMapperScipionchimera[SYM_In25r]="n25r"

def getEnviron(ccp4First=True):
    return getChimeraEnviron()

def runChimeraProgram(program, args=""):
    """ Internal shortcut function to launch chimera program. """
    env=getEnviron()
    pwutils.runJob(None, program, args, env=env)

def getProgram(progName="chimera"):
    """ Return the program binary that will be used. """
    if (not 'CHIMERA_HOME' in os.environ):
        return None
    return os.path.join(os.environ['CHIMERA_HOME'], 'bin',
                           os.path.basename(progName))

