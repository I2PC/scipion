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
import glob

CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2

class XmippDefCL2DViewer(Form):
    """Create the definition of parameters for
    the XmippProtCL2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Visualization')
        self.addParam('classesToShow', EnumParam, choices=['Classes', 'Class Cores', 'Class Stable Cores'],
                      label="What to show", default=CLASS_CORES,
                      display=EnumParam.DISPLAY_COMBO)
        self.addParam('doShowClassHierarchy', BooleanParam, label="Visualize class hierarchy.", default=True)      
        self.addParam('doShowLastLevel', BooleanParam, label="Visualize last level.", default=True)     
        self.addParam('showSeveralLevels', StringParam, default='',
              label='Visualize several levels', condition='not doShowLastLevel',
              help='Create a  list of levels like: 0,1,3 or 0-3 ')    

        
class XmippCL2DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtCL2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _definition = XmippDefCL2DViewer()
    _label = 'Xmipp Viewer CL2D'
    
        
    def _viewAll(self, *args):
        classesToShow = self.classesToShow.get()
        
        fnSubset = ""
        if classesToShow == CLASSES:
            fnSubset = ""
        elif classesToShow == CLASS_CORES:
            fnSubset = "_core"
        elif classesToShow == CLASS_STABLE_CORES:
            fnSubset = "_stable_core"
        
        levelFiles = self.protocol._getLevelMdFiles(fnSubset)
        if levelFiles:
            levelFiles.sort()
            lastLevelFile = levelFiles[-1]
            if self.doShowLastLevel.get():
                runShowJ("classes@"+lastLevelFile)
            else:
                listOfLevels = self.__getListFromRangeString(self.showSeveralLevels.get())
                lastLevel = int(re.search('level_(\d\d)',lastLevelFile).group(1))
                if max(listOfLevels) <= lastLevel:
                    files = "";
                    for level in listOfLevels:
                        print "en for con level=%s" % level
                        fn = self.protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                        if os.path.exists(fn):
                            files += "classes_sorted@"+fn+" "
                    if files != "":
                        runShowJ(files)
        if self.doShowClassHierarchy.get():
            fnHierarchy = self.protocol._getExtraPath("classes%s_hierarchy.txt" % fnSubset)
            #TODO: Adapt showTextFileViewer not to use Xmipp objects and to be scrollable
#            if os.path.exists(fnHierarchy):
#                showTextfileViewer(fnHierarchy,[fnHierarchy])

#TODO: This method should not be necessary, instead NumericListParam should return a list and not a String 
    def __getListFromRangeString(self, rangeStr):
        ''' Create a list of integer from a string with range definitions
        Examples:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]
        '''
        elements = rangeStr.split(',')
        values = []
        for e in elements:
            if '-' in e:
                limits = e.split('-')
                values += range(int(limits[0]), int(limits[1])+1)
            else:
                # If values are separated by comma also splitted 
                values += map(int, e.split())
        return values
            

        

