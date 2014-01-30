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
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_nma import XmippProtNMA
from protocol_nma_alignment import XmippProtAlignmentNMA
from viewer import runShowJ
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
# FIXME: REMOVE DEPENDENCIES FROM XMIPP
from protlib_gui_ext import chimera, showj
import glob

CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2
   
        
class XmippAlignmentNMAViewer(ProtocolViewer):
    """ Visualization of results from the NMA protocol
    """
    _label = 'viewer nma alignment'
    _targets = [XmippProtAlignmentNMA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
        
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayRawDeformation', StringParam, default='1',
                      label='Display raw deformation',
                      help='Type 1 to see the histogram of raw deformation number 1; \n'
                           'type 2 to see the histogram of raw deformation number 2, etc.\n'
                           'Type 1 2 to see the 2D plot of raw deformations number 1 vs 2.\n'
                           'Type 1 2 3 to see the 3D plot of raw deformations 1, 2 and 3; etc.'
                           )
        form.addParam('analyzeMatlab', BooleanParam, default=True, 
                      label="Analyze with Matlab?",
                      help='In MATLAB you may create image clusters and draw deformation trajectories.')
    
    def _getVisualizeDict(self):
        return {'displayRawDeformation': self._viewRawDeformation,
                'analyzeMatlab': self._viewWithMatlab
                } 
                        
    def _viewWithMatlab(self, paramName):
        os.system("matlab -r \"xmipp_nma_selection_tool(\'"+self._getPath()+"\')\"")
        
    def _viewRawDeformation(self, paramName):
        plotter = self._doViewRawDeformation(paramName)
        plotter.show()
        
    def _doViewRawDeformation(self, paramName):
        components = map(int, self.displayRawDeformation.get().split())
        dim = len(components)
        if dim > 0:
            modeList = []
            modeNameList = []
            # Get modes
            md = xmipp.MetaData(self.protocol.inputModes.get().getFileName())
            md.removeDisabled()
            
            for mode in components:
                if mode > md.size():
                    from protlib_gui_ext import showWarning
                    showWarning('Warning', "You don't have so many modes",parent=self.master)
                else:
                    mode -= 1
                    modeName = ""
                    for i, objId in enumerate(md):
                        modeName = md.getValue(xmipp.MDL_NMA_MODEFILE, objId)
                        if i == mode:
                            break
                    modeNameList.append(modeName)
                    modeList.append(mode)
            
            defFn = self.protocol._getExtraPath('deformations.txt')
            # Actually plot
            
            plotter = XmippNmaPlotter(defFn, dirname(modeNameList[0])) 
            baseList = [basename(n) for n in modeNameList]
            
            if dim == 1:
                plotter.plotArray1D("Histogram for mode: %s" % baseList[0], modeList[0], 
                                    "Deformation value", "Number of images")
            elif dim == 2:
                plotter.plotArray2D("Mode %s vs %s" % tuple(baseList), 
                                    modeList[0], modeList[1], *baseList)
            elif dim == 3:
                plotter.plotArray3D("Modes %s %s %s" % tuple(baseList), 
                                    modeList[0], modeList[1], modeList[2], *baseList)
            return plotter

    
    def getVisualizeDictWeb(self):
        return {'displayRawDeformation': 'doDisplayRawDeformation',
                'analyzeMatlab': 'doAnalyzeMatlab'
                } 
        
    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerNMAAlign"
        
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from plotter import XmippPlotter


class XmippNmaPlotter(XmippPlotter):
    """ Add some extra plot utilities to XmippPlotter class, mainly for
    NMA vectors plotting of the deformations.txt file.
    """
    def __init__(self, fnArray, vectorsPath):
        """ Load the data.
        Params:
            fnArray: path to deformations.txt data file.
            vectorsPath: the folder path in which the vectors are located.
        """
        self._data = np.loadtxt(fnArray)
        XmippPlotter.__init__(self, windowTitle=vectorsPath)
        
    def plotArray1D(self, title, col, xlabel, ylabel):
        data = self._data
        if len(data.shape) > 1:
            data = self._data[:,col]
        
        ax = self.createSubPlot(title, xlabel, ylabel)
        
        # histogram our data with numpy
        n, bins = np.histogram(data, 50)
        
        # get the corners of the rectangles for the histogram
        left = np.array(bins[:-1])
        right = np.array(bins[1:])
        bottom = np.zeros(len(left))
        top = bottom + n
        
        # we need a (numrects x numsides x 2) numpy array for the path helper
        # function to build a compound path
        XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
        
        # get the Path object
        barpath = path.Path.make_compound_path_from_polys(XY)
        
        # make a patch out of it
        patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
        ax.add_patch(patch)
        
        # update the view limits
        ax.set_xlim(left[0], right[-1])
        ax.set_ylim(bottom.min(), top.max())
    
    def plotArray2D(self, title, colX, colY, xlabel, ylabel):
        ax = self.createSubPlot(title, xlabel, ylabel)
        ax.plot(self._data[:,colX], self._data[:,colY], '.')
    
    def plotArray3D(self, title, colX, colY, colZ, xlabel, ylabel, zlabel):
        ax = p3.Axes3D(self.figure)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        a = self._data
        ax.scatter3D(a[:,colX], a[:,colY], a[:,colZ])
