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
from viewer import runShowJ
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
import glob

CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2
   
        
class XmippCL2DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtCL2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer cl2d'
   

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('classesToShow', EnumParam, choices=['Classes', 'Class Cores', 'Class Stable Cores'],
                      label="What to show", default=CLASS_CORES,
                      display=EnumParam.DISPLAY_COMBO)
        form.addParam('doShowClassHierarchy', BooleanParam, label="Visualize class hierarchy.", default=True)      
        form.addParam('doShowLastLevel', BooleanParam, label="Visualize last level.", default=True)     
        form.addParam('showSeveralLevels', StringParam, default='',
              label='Visualize several levels', condition='not doShowLastLevel',
              help='Create a  list of levels like: 0,1,3 or 0-3 ')    
    
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
            showTextfileViewer(fnHierarchy,[fnHierarchy]) 
            
    def _viewLevelFiles(self, e=None):
        fnSubset = self._getSubset()
        levelFiles = self.protocol._getLevelMdFiles(fnSubset)
        if levelFiles:
            levelFiles.sort()
            lastLevelFile = levelFiles[-1]
            if self.doShowLastLevel:
                runShowJ("classes@"+lastLevelFile)
            else:
                if self.showSeveralLevels.empty():
                    self.formWindow.showError('Please select the levels that you want to visualize.')
                else:
                    listOfLevels = []
                    try:
                        listOfLevels = self._getListFromRangeString(self.showSeveralLevels.get())
                    except Exception, ex:
                        self.formWindow.showError('Invalid levels range.')
                        
                    #lastLevel = int(re.search('level_(\d\d)',lastLevelFile).group(1))
                    #if max(listOfLevels) <= lastLevel:
                    files = "";
                    for level in listOfLevels:
                        fn = self.protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                        if os.path.exists(fn):
                            files += "classes_sorted@"+fn+" "
                        else:
                            self.formWindow.showError('Level %s does not exist.' % level)
                    if files != "":
                        runShowJ(files)        
        
    def _viewAll(self, *args):

        self._viewLevelFiles()
            
        if self.doShowClassHierarchy:
            self._viewClassHierarchy()

    def getVisualizeDictWeb(self):
        return {'doShowClassHierarchy': "viewClassHierarchy",
                'doShowLastLevel': "viewLevelFiles"}
        
    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerCL2D"
        
