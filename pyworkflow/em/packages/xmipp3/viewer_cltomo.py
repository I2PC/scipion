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
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_cltomo import XmippProtCLTomo
from viewer import runShowJ
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
import glob

class XmippCLTomoViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer cltomo'
    _targets = [XmippProtCLTomo]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowLastLevel', BooleanParam, label="Visualize last level.", default=True)     
        form.addParam('showSeveralLevels', StringParam, default='',
              label='Visualize several levels', condition='not doShowLastLevel',
              help='Create a  list of levels like: 0,1,3 or 0-3 ')    
    
    def _getVisualizeDict(self):
        return {'doShowLastLevel': self._viewLevelFiles}        

    def _viewLevelFiles(self, e=None):
        import glob
        levelFiles=glob.glob(self.protocol._getExtraPath("results_classes_level*.xmd"))
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
                        fn = self.protocol._getExtraPath("results_classes_level_%02d.xmd"%level)
                        if os.path.exists(fn):
                            files += "classes@"+fn+" "
                        else:
                            self.formWindow.showError('Level %s does not exist.' % level)
                    if files != "":
                        runShowJ(files)        
        
    def _viewAll(self, *args):
        self._viewLevelFiles()

    def getVisualizeDictWeb(self):
        return {'doShowLastLevel': "viewLevelFiles"}
        
    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerCLTomo"
        
