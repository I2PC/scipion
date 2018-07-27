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
1. define phenix environ
TODO:
2. Read/Write PHENIX specific files
"""

import os
import pyworkflow.utils as pwutils

PHENIX_PYTHON = 'phenix.python '  # keep the ending space
PHENIX_SCRIPT_PATH = 'modules/cctbx_project/mmtbx/command_line'


def getEnviron(first=True):
    environ = pwutils.Environ(os.environ)
    pos = pwutils.Environ.BEGIN if first else pwutils.Environ.END
    _phenix_home = os.environ['PHENIX_HOME']

    # add to variable
    environ.update({
            'PATH': os.path.join(_phenix_home, 'build','bin'),
    }, position=pos)

    # replace variable value
    environ.update({
        'LIBTBX_BUILD': os.path.join(_phenix_home, 'build'),
        'LIBTBX_OPATH': os.environ['PATH'],
    }, position=pwutils.Environ.REPLACE)  # replace

    return environ


def runPhenixProgram(program, args=None, extraEnvDict=None, cwd=None):
    """ Internal shortcut function to launch a Phenix program. """
    env = getEnviron()
    if extraEnvDict is not None:
        env.update(extraEnvDict)
    program = PHENIX_PYTHON + program
    pwutils.runJob(None, program, args, env=env, cwd=cwd)


def getProgram(progName):
    """ Return the program binary that will be used. """
    if 'PHENIX_HOME' not in os.environ:
        return None

    return os.path.join(os.environ['PHENIX_HOME'],
                        PHENIX_SCRIPT_PATH,
                        os.path.basename(progName))
