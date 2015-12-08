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
This sub-package will contains Eman3.0 specific protocols
"""

import os
from os.path import join
from pyworkflow.utils import Environ


def getEnviron():
    """ Setup the environment variables needed to launch Eman. """
    environ = Environ(os.environ)
    EMAN2DIR = os.environ['EMAN2DIR']
    pathList = [os.path.join(EMAN2DIR, d) for d in ['lib', 'bin', 'extlib/site-packages']]
    
    environ.update({
            'PATH': join(EMAN2DIR, 'bin'),
            'LD_LIBRARY_PATH': os.pathsep.join(pathList),
            'PYTHONPATH': os.pathsep.join(pathList),
            'EMAN_PYTHON': os.path.join(EMAN2DIR, 'Python/bin/python')
            }, position=Environ.REPLACE)
    return environ


def getVersion():
    path = os.environ['EMAN2DIR']
    if '2.11' in path:
        return '2.11'
    elif '2.12' in path:
        return '2.12'
    else:
        return ''


def getEmanProgram(program):
    if not 'EMAN_PYTHON' in os.environ:
        os.environ['EMAN_PYTHON'] = os.path.join(os.environ['EMAN2DIR'], 'Python/bin/python')
    # For EMAN2 python scripts, join the path to bin
    # Comment this out cause there are eman programs that do not start with e2
    #if program.startswith('e2'):
    program = os.path.join(os.environ['EMAN2DIR'], 'bin', program)
    #raise Exception('EMAN_PYTHON is not load in environment')
    python = os.environ['EMAN_PYTHON']
    return '%(python)s %(program)s ' % locals()
    
    
def getEmanCommand(program, args):    
    return getEmanProgram(program) + args
    
    
