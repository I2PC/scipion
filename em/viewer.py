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
This module is mainly for the Viewer class, which 
serve as base for implementing visualization tools(Viewer sub-classes).
"""

from os.path import join

class Viewer():
    """ All visualization wrappers should user the Viewer class
    as base and provide the implementation to launch specific 
    command line tools in order to visualize objects.
    
    The _target class property should contains a list of string
    with the class names that this viewer is able to visualize.
    For example: _target = ['Image', 'SetOfImages']
    """
    _target = [] 
    
    def __init__(self, tmpPath='.', **args):
        self._tmpPath = tmpPath
        
    def _getTmpPath(self, *paths):
        return join(self._tmpPath, *paths)
    
    def visualize(self, obj):
        """ This method should make the necessary convertions
        and call the command line utilities to visualize this
        particular object.
        """
        pass



