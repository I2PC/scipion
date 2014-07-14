# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
from pyworkflow.viewer import View, Viewer, DESKTOP_TKINTER
from pyworkflow.em.protocol import PdbFile
import pyworkflow as pw
from showj import (runJavaIJapp, ZOOM, ORDER, VISIBLE, 
                   MODE, PATH, TABLE_NAME, MODE_MD, RENDER)
# PATH is used by app/em_viewer.py



#------------------------ Some common Views ------------------

class DataView(View):
    """ Wrapper the arguments to showj (either web or desktop). """
    def __init__(self, path, viewParams={}, **kwargs):
        View.__init__(self)
        self._memory = '1g'
        self._loadPath(path)
        self._env = kwargs.get('env', os.environ.copy())
        self._viewParams = viewParams
        if self._path.endswith('.star'):
            from pyworkflow.em.packages.relion import addRelionLabelsToEnviron
            addRelionLabelsToEnviron(self._env)
            
    def _loadPath(self, path):
        # Check if there is a table name with @ in path
        # in that case split table name and path
        # table names can never starts with a number
        # this is considering an image inside an stack
        if '@' in path and path[0] not in '0123456789':
            self._tableName, self._path = path.split('@')
        else:
            self._tableName, self._path = None, path
            
    def show(self):        
        runJavaIJapp(self._memory, 'xmipp.viewer.Viewer', self.getShowJParams(), True, env=self._env)
    
    def getShowJParams(self):
        params = '-i ' + self._path
        for key, value in self._viewParams.items():
            params = "%s --%s %s"%(params, key, value)
        
        return params
    
    def getShowJWebParams(self):
    
    #=OLD SHOWJ WEB DOCUMENTATION===============================================
    # Extra parameters can be used to configure table layout and set render function for a column
    # Default layout configuration is set in ColumnLayoutProperties method in layout_configuration.py
    # 
    # Parameters are formed by: [label]___[property]: [value]. E.g.: id___visible:True or micrograph___renderFunc:"get_image_psd"
    # Properties to be configured are:
    #    visible: Defines if this column is displayed
    #    allowSetVisible: Defines if user can change visible property (show/hide this column).
    #    editable: Defines if this column is editable, ie user can change field value.
    #    allowSetEditable: Defines if user can change editable property (allow editing this column).
    #    renderable: Defines if this column is renderizable, ie it renders data column using renderFunc
    #    allowSetRenderable: Defines if user can change renderable property.
    #    renderFunc: Function to be used when this field is rendered. (it has to be inserted in render_column method)
    #    extraRenderFunc: Any extra parameters needed for rendering. Parameters are passed like in a url ie downsample=2&lowPass=3.5
    # 
    # Example:
    # extraParameters["id___visible"]=True
    # extraParameters["micrograph___renderFunc"]="get_image_psd"
    # extraParameters["micrograph___extraRenderFunc"]="downsample=2"
    #===========================================================================
    
        parameters = {
            'mode', # FOR MODE TABLE OR GALLERY
            'visible',
            'zoom',
            'order',
#            'allowSetVisible':'allowSetVisible',
#            'editable':'editable',
#            'allowSetEditable':'allowSetEditable',
#            'render': 'renderable',
#            'allowSetRenderable':'allowSetRenderable',
#            'renderFunc':'renderFunc',
#            'extraRenderFunc':'extraRenderFunc',
        }
        
        params = {}
        for key, value in self._viewParams.items():
            
            if key in parameters:
                if key in {'order','zoom', 'mode'}:
                    if key == 'mode' and value =='metadata':
                        value = 'table'
                    params[key] = value
                else:    
                    for val in value.split(' '):
                        params[val] = '%s___%s' % (val, key)
                     
        return params
        
    def getPath(self):
        return self._path
    
    def getTableName(self):
        return self._tableName
        
        
class ObjectView(DataView):
    """ Wrapper to DataView but for displaying Scipion objects. """
    def __init__(self, projectid, inputid, path, other='', viewParams={}, **kwargs):
        DataView.__init__(self, path, viewParams, **kwargs)
        self.python = pw.PYTHON
        self.scripts = pw.join('apps')
        self.type = type
        self.projectid = projectid
        self.inputid = inputid
        self.other = other
        
    def getShowJParams(self):
        # mandatory to provide scipion params
        params = DataView.getShowJParams(self) + ' --scipion %s %s \"%s\" %s %s'%(self.python, self.scripts,  self.projectid, self.inputid, self.other)
        return params
    
    def show(self):
        runJavaIJapp(self._memory, 'xmipp.viewer.scipion.ScipionViewer', self.getShowJParams(), True, env=self._env)
        
        
class ClassesView(ObjectView):
    """ Customized ObjectView for SetOfClasses. """
    def __init__(self, projectid, inputid, path, other='', viewParams={}, **kwargs):
        labels =  'enabled id _size _representative._filename'
        defaultViewParams = {ORDER:labels,
                             VISIBLE: labels, RENDER:'_representative._filename',
                             'sortby': '_size', 'label': '_size',
                             }
        defaultViewParams.update(viewParams)
        ObjectView.__init__(self, projectid, inputid, path, other, defaultViewParams, **kwargs)
        
        
class Classes3DView(ClassesView):
    """ Customized ObjectView for SetOfClasses. """
    def __init__(self, projectid, inputid, path, other='', viewParams={}, **kwargs):
        defaultViewParams = {ZOOM: '99', MODE: 'metadata'}
        defaultViewParams.update(viewParams)
        ClassesView.__init__(self, projectid, inputid, path, other, defaultViewParams, **kwargs)
            
            
class CoordinatesObjectView(DataView):
    """ Wrapper to View but for displaying Scipion objects. """
    def __init__(self, path, outputdir, viewParams={}, **kwargs):
        DataView.__init__(self, path, **kwargs)
        self.python = pw.PYTHON
        self.script = pw.join('apps', 'pw_create_coords_subset.py')
        self.outputdir = outputdir
        
    def getShowJParams(self):
        params = '--input %s --output %s --mode readonly'%(self._path, self.outputdir)
        return params
    
    def show(self):
        runJavaIJapp(self._memory, 'xmipp.viewer.particlepicker.training.SupervisedPickerRunner', self.getShowJParams(), True, env=self._env)
        
        
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
                pdb = os.path.basename(fn)
                fn += '_tmp_chimera.cmd'
                f = open(fn, 'w')
                f.write("open %s\n" % pdb)
                f.write("rangecol bfactor,a 0 white 1 red \n")
                f.write("setattr a radius 1 \n")
                f.write("represent sphere \n")
                f.close()
            chimera(fn)
        else:
            raise Exception('ChimeraViewer.visualize: can not visualize class: %s' % obj.getClassName())
 
 
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
            raise Exception('VmdViewer.visualize: can not visualize class: %s' % obj.getClassName())     
        
