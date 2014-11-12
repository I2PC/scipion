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

import numpy as np
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter


class XmippNmaPlotter(XmippPlotter):
    """ Add some extra plot utilities to XmippPlotter class, mainly for
    NMA vectors plotting of the deformations.txt file.
    """
    def __init__(self, fnArray, **kwargs):
        """ Load the data.
        Params:
            fnArray: path to deformations.txt data file.
            vectorsPath: the folder path in which the vectors are located.
        """
        self._data = np.loadtxt(fnArray)
        XmippPlotter.__init__(self, **kwargs)
        self.useLastPlot = False
        
    def createSubPlot(self, title, xlabel, ylabel):
        if self.useLastPlot and self.last_subplot:
            ax = self.last_subplot
            ax.cla()
            ax.set_title(title)
        else:
            ax = XmippPlotter.createSubPlot(self, title, xlabel, ylabel)
            
        return ax
            
    def plotArray1D(self, title, col, xlabel, ylabel):
        data = self._data
        if len(data.shape) > 1:
            data = self._data[:,col]
        
        ax = self.createSubPlot(title, xlabel, ylabel)
        ax.hist(data, 50)
        
        return
        
        # CHECK WITH COSS why not to use matplotlib histogram ???
#         # histogram our data with numpy
#         n, bins = np.histogram(data, 50)
#         
#         # get the corners of the rectangles for the histogram
#         left = np.array(bins[:-1])
#         right = np.array(bins[1:])
#         bottom = np.zeros(len(left))
#         top = bottom + n
#         
#         # we need a (numrects x numsides x 2) numpy array for the path helper
#         # function to build a compound path
#         XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
#         
#         import matplotlib.patches as patches
#         import matplotlib.path as path
#         # get the Path object
#         barpath = path.Path.make_compound_path_from_polys(XY)
#         
#         # make a patch out of it
#         patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
#         ax.add_patch(patch)
#         
#         # update the view limits
#         ax.set_xlim(left[0], right[-1])
#         ax.set_ylim(bottom.min(), top.max())
    
    def plotArray2D(self, title, colX, colY, xlabel, ylabel):
        ax = self.createSubPlot(title, xlabel, ylabel)
        #ax.plot(self._data[:,colX], self._data[:,colY], 'o')
        ax.scatter(self._data[:,colX], self._data[:,colY], c=self._data[:,colX])
    
    def plotArray3D(self, title, colX, colY, colZ, xlabel, ylabel, zlabel):
        import mpl_toolkits.mplot3d.axes3d as p3
        ax = p3.Axes3D(self.figure)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        a = self._data
        ax.scatter3D(a[:,colX], a[:,colY], a[:,colZ])
        # Disable tight_layout that is not available for 3D 
        self.tightLayoutOn = False
        

