# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
This sub-package contains data and protocol classes
wrapping Pytom
"""
import os

from pyworkflow.utils import Environ



def getEnviron():
    """ Create the needed environment for Pytom programs. """
    environ = Environ(os.environ)
    if 'PYTOM_HOME' in os.environ:
        PYTOM_HOME = os.environ['PYTOM_HOME']
    else:
        print ('Warning: Cannot find environ variable PYTOM_HOME. Using '
               '$SCIPION_HOME/software/em/pytom instead.')
        PYTOM_HOME = '%s/software/em/pytom' % os.environ['SCIPION_HOME']
    environ.update({
                    'LD_LIBRARY_PATH': os.path.join(PYTOM_HOME,'pytomc/libs/libtomc/libs'),
                    'PATH': os.path.join(PYTOM_HOME, 'bin'),
                    'PYTHONPATH': os.pathsep.join([os.path.join(PYTOM_HOME, suffix) for suffix in [
                                                   'pytomc', 
                                                   'pytomc/swigModules', 
                                                   'pytomc/sh_alignment/frm/swig', 
                                                   'external/lib',
                                                   'external/lib/python/site-packages']] + [os.path.dirname(PYTOM_HOME)])
                    }, 
                   position=Environ.BEGIN)
    return environ
