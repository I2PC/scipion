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

import os, sys

from pyworkflow.em import ProtPreprocessMicrographs   



class EmanProtPreprocessMicrographs(ProtPreprocessMicrographs):
    pass


def loadEnvironment():
    """ Load the environment variables needed for use EMAN2 tools. """
    EMAN2DIR = os.environ['EMAN2DIR']
    os.environ['PATH'] = "%(EMAN2DIR)s/bin" % locals() + os.pathsep + os.environ['PATH']
    pathList = [os.path.join(EMAN2DIR, d) for d in ['lib', 'bin', 'extlib/site-packages']]
    os.environ['PYTHONPATH'] = os.pathsep.join(pathList) + os.pathsep + os.environ['PYTHONPATH']
    os.environ['EMAN_PYTHON'] = os.path.join(EMAN2DIR, 'Python/bin/python')
    # Undo the xmipp_python additions to LD_LIBRARY_PATH
    if 'OLD_LD_LIBRARY_PATH' in os.environ:
        print 'OLD_LD_LIBRARY_PATH', os.environ['OLD_LD_LIBRARY_PATH']
        os.environ['LD_LIBRARY_PATH'] = os.environ['OLD_LD_LIBRARY_PATH']
    #sys.path += pathList
    
def getEmanCommand(program, args):
    if not 'EMAN_PYTHON' in os.environ:
        os.environ['EMAN_PYTHON'] = os.path.join(os.environ['EMAN2DIR'], 'Python/bin/python')
        #raise Exception('EMAN_PYTHON is not load in environment')
    python = os.environ['EMAN_PYTHON']
    
    return '%(python)s %(program)s %(args)s' % locals()
    
    
