# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
import pyworkflow.utils as pwutils
from pyworkflow.em.constants import SYM_CYCLIC, SYM_DIHEDRAL, \
    SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222, SYM_I222r, SYM_In25, SYM_In25r
from pyworkflow.utils import  Environ

#
chimeraPdbTemplateFileName = "chimeraOut%04d.pdb"
chimeraMapTemplateFileName = "chimeraOut%04d.mrc"
chimeraScriptFileName = "chimeraScript.py"
sessionFile = "SESSION.py"

symMapperScipionchimera = {}
symMapperScipionchimera[SYM_CYCLIC] = "Cn"
symMapperScipionchimera[SYM_DIHEDRAL] = "Dn"
symMapperScipionchimera[SYM_TETRAHEDRAL] = "T"
symMapperScipionchimera[SYM_OCTAHEDRAL] = "O"
symMapperScipionchimera[SYM_I222] = "222"
symMapperScipionchimera[SYM_I222r] = "222r"
symMapperScipionchimera[SYM_In25] = "n25"
symMapperScipionchimera[SYM_In25r] = "n25r"


def getEnviron():
    return getChimeraEnviron()


def runChimeraProgram(program, args=""):
    """ Internal shortcut function to launch chimera program. """
    env = getEnviron()
    pwutils.runJob(None, program, args, env=env)


def getProgram(progName="chimera"):
    """ Return the program binary that will be used. """
    if 'CHIMERA_HOME' not in os.environ:
        return None
    return os.path.join(os.environ['CHIMERA_HOME'], 'bin',
                        os.path.basename(progName))


def createCoordinateAxisFile(dim, bildFileName="/tmp/axis.bild",
                             sampling=1, r1=0.1):
    ff = open(bildFileName, "w")
    arrowDict = {}
    arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = \
        sampling * dim * 3. / 4.
    arrowDict["r1"] = r1 * dim / 50.
    arrowDict["r2"] = 4 * r1
    arrowDict["rho"] = 0.75  # axis thickness

    ff.write(".color 1 0 0\n"
             ".arrow 0 0 0 %(x)d 0 0 %(r1)f %(r2)f %(rho)f\n"
             ".color 1 1 0\n"
             ".arrow 0 0 0 0 %(y)d 0 %(r1)f %(r2)f %(rho)f\n"
             ".color 0 0 1\n"
             ".arrow 0 0 0 0 0 %(z)d %(r1)f %(r2)f %(rho)f\n" %
             arrowDict)
    ff.close()

#def adaptOriginFromCCP4ToChimera(origin):
#    print "Do not use this"
#    exit(0)
#    return tuple(-1.0*x for x in origin)


def getChimeraEnviron():
    """ Return the proper environ to launch chimera.
    CHIMERA_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = Environ(os.environ)
    environ.set('PATH', os.path.join(os.environ['CHIMERA_HOME'], 'bin'),
                position=Environ.BEGIN)

    if "REMOTE_MESA_LIB" in os.environ:
        environ.set('LD_LIBRARY_PATH', os.environ['REMOTE_MESA_LIB'],
                    position=Environ.BEGIN)
    return environ


