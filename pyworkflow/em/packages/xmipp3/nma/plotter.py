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
    def __init__(self, **kwargs):
        """ Create the plotter, 'data' should be passed in **kwargs.
        """
        self._data = kwargs.get('data')
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
            
    def plotArray1D(self, title, xlabel, ylabel):
        xdata = self._data.getXData()
        ax = self.createSubPlot(title, xlabel, ylabel)
        ax.hist(xdata, 50)
    
    def plotArray2D(self, title, xlabel, ylabel):
        ax = self.createSubPlot(title, xlabel, ylabel)
        #ax.plot(self._data[:,colX], self._data[:,colY], 'o')
        xdata = self._data.getXData()
        ydata = self._data.getYData()
        weights = self._data.getWeights()
        cax = ax.scatter(xdata, ydata, c=weights)
        ax.figure.colorbar(cax)
    
    def plotArray3D(self, title, xlabel, ylabel, zlabel):
        import mpl_toolkits.mplot3d.axes3d as p3
        ax = p3.Axes3D(self.figure)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        xdata = self._data.getXData()
        ydata = self._data.getYData()
        zdata = self._data.getZData()
        weights = self._data.getWeights()
        cax = ax.scatter3D(xdata, ydata, zdata, c=weights)
        ax.figure.colorbar(cax)
        # Disable tight_layout that is not available for 3D 
        self.tightLayoutOn = False
        

