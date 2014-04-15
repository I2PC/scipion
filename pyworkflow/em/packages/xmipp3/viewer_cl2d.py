# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_cl2d import XmippProtCL2D
from viewer import runScipionShowJ
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
import glob


CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2   
CLASS_CHOICES = ['Classes', 'Class Cores', 'Class Stable Cores']

LEVEL_LAST = 0
LEVEL_SEL = 1
LEVEL_CHOICES = ['last', 'selection']
        
        
class XmippCL2DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer cl2d'
    _targets = [XmippProtCL2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('classesToShow', EnumParam, choices=CLASS_CHOICES,
                      label="What to show", default=CLASS_CORES,
                      display=EnumParam.DISPLAY_LIST)
        form.addParam('doShowClassHierarchy', BooleanParam, default=True, 
                      label="Visualize class hierarchy.")      
        form.addParam('doShowLastLevel', EnumParam, default=LEVEL_LAST, 
                      choices=LEVEL_CHOICES,
                      display=EnumParam.DISPLAY_LIST,
                      label="Level to visualize")     
        form.addParam('showSeveralLevels', StringParam, default='',
              label='Levels selection', condition='doShowLastLevel==%d' % LEVEL_SEL,
              help='Specify a  list of levels like: 0,1,3 or 0-3 ')    
    
    def _getVisualizeDict(self):
        return {'doShowClassHierarchy': self._viewClassHierarchy,
                'doShowLastLevel': self._viewLevelFiles}        

    def _getSubset(self):
        fnSubset = ""
        classesToShow = self.classesToShow.get()
        if classesToShow == CLASSES:
            fnSubset = ""
        elif classesToShow == CLASS_CORES:
            fnSubset = "_core"
        elif classesToShow == CLASS_STABLE_CORES:
            fnSubset = "_stable_core" 
        return fnSubset
                        
    def _viewClassHierarchy(self, e=None):
        fnSubset = self._getSubset()
        fnHierarchy = self.protocol._getExtraPath("classes%s_hierarchy.txt" % fnSubset)
        if os.path.exists(fnHierarchy):
            return [self.textView([fnHierarchy])] 
            
    def _viewLevelFiles(self, e=None):
        fnSubset = self._getSubset()
        views = []
        errors = []
        obj = getattr(self.protocol, "outputClasses" + fnSubset)
        levelFiles = self.protocol._getLevelMdFiles(fnSubset)
        
        if levelFiles:
            levelFiles.sort()
            
            if self.doShowLastLevel == LEVEL_LAST:
                views.append(ObjectView("classes@"+levelFiles[-1], "Particles", 
                                        self._project.getName(), obj.strId(), 
                                        obj.getImages().strId()))
            else:
                if self.showSeveralLevels.empty():
                    errors.append('Please select the levels that you want to visualize.')
                else:
                    listOfLevels = []
                    try:
                        listOfLevels = self._getListFromRangeString(self.showSeveralLevels.get())
                    except Exception, ex:
                        errors.append('Invalid levels range.')
                        
                    files = []
                    for level in listOfLevels:
                        fn = self.protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                        if os.path.exists(fn):
                            files.append("classes_sorted@%s " % fn)
                        else:
                            errors.append('Level %s does not exist.' % level)
                    
                    inputImagesId = self.protocol.inputImages.get().strId()
                    for fn in files:                        
                        views.append(ObjectView(fn, "Particles", 
                                                self._project.getName(), 
                                                self.protocol.strId(), inputImagesId))
        
        if errors:
            views.append(self.errorMessage('\n'.join(errors), "Visualization errors"))
            
        return views        

