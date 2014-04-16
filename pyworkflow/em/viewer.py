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
This module implement viewers for some type of common objects.
"""

import os
from pyworkflow.viewer import View, MessageView, CommandView, Viewer, DESKTOP_TKINTER
from pyworkflow.em.data import *
from pyworkflow.em.protocol import *
from showj import *


#------------------------ Some common Views ------------------

class DataView(View):
    """ Wrapper the arguments to showj (either web or desktop). """
    def __init__(self, path, **kwargs):
        View.__init__(self)
        self._loadPath(path)
        self._env = kwargs.get('env', None)
        
    def _loadPath(self, path):
        if '@' in path:
            self._tableName, self._path = path.split('@')
        else:
            self._tableName, self._path = None, path
            
    def show(self):
        runShowJ(self._path, env=self._env)
        
    def getPath(self):
        return self._path
    
    def getTableName(self):
        return self._tableName
        
        
class ObjectView(DataView):
    """ Wrapper to View but for displaying Scipion objects. """
    def __init__(self, path, *args, **kwargs):
        DataView.__init__(self, path, **kwargs)
        self._args = args
        
    def show(self):
        runScipionShowJ(self._path, *self._args, env=self._env)   
        
        
#------------------------ Some viewers ------------------------

class ChimeraViewer(Viewer):
    """ Wrapper to visualize PDB object with Chimera. """
    _environments = [DESKTOP_TKINTER]
    _targets = [PdbFile]
    
    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):        
        cls = type(obj)
        
        if issubclass(cls, PdbFile):
            fn = obj.getFileName()
            from protlib_gui_ext import chimera
            if obj.getPseudoAtoms():
                # Write an script for this case to use colored spheres
                pdb = basename(fn)
                fn += '_tmp_chimera.cmd'
                f = open(fn, 'w')
                f.write("open %s\n" % pdb)
                f.write("rangecol bfactor,a 0 white 1 red \n")
                f.write("setattr a radius 1 \n")
                f.write("represent sphere \n")
                f.close()
            chimera(fn)
        else:
            raise Exception('XmippViewer.visualize: can not visualize class: %s' % obj.getClassName())
 
class VmdViewer(Viewer):
    """ Wrapper to visualize PDB objects with VMD viewer. """
    _environments = [DESKTOP_TKINTER]
    _targets = [PdbFile]
    
    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):        
        cls = type(obj)
        
        if issubclass(cls, PdbFile):
            fn = obj.getFileName()
            os.system("vmd %s" % fn)
        else:
            raise Exception('XmippViewer.visualize: can not visualize class: %s' % obj.getClassName())     
        
