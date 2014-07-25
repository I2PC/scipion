# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import ObjectView, MODE, MODE_MD, ORDER, VISIBLE
from pyworkflow.em.plotter import EmPlotter

from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy
import collections
import numpy as np


class XmippCTFDiscrepancyViewer(Viewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer CTF Discrepancy'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtCTFDiscrepancy]
    
    def _visualize(self, obj):
        print "WARNING still under development"
        views = []
        fn = obj.outputCTF.getFileName()
    #            self._views.append(DataView(fn, viewParams={MODE: 'metadata'}))
        labels = 'id enabled _micObj._filename method1 method2 resolution _defocusU _defocusV _defocusAngle' 
        views.append(ObjectView(self._project.getName(), obj.strId(), fn, 
                                      viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        
        views.append(self._createPlots(obj))
        views.append(self._createMatrix(obj))
        
        return views 
    
    def _createPlots(self, obj):
        print "inside _createPlot"
        plotter = EmPlotter()
        plotter.createSubPlot("Resolution Discrepancies histogram", 
                      "Resolution", "# of Micrographs")
        resolution = [ctf.resolution.get() for ctf in obj.outputCTF]
        n = len(resolution)
        plotter.plotHist(resolution, nbins=min(n, 1+n/20))
        return plotter

    def _createMatrix(self,obj):
        inputCTFs=obj.inputCTFs
        _matrix = np.zeros(shape=(len(inputCTFs), len(inputCTFs)))
        _matrix[0][0]=1
        _matrix[1][0]=2
        _matrix[0][1]=3
        _matrix[1][1]=4
        
        ticksLablesMajor=[]

        self.methodNames = collections.OrderedDict()
        for i, ctf in enumerate(inputCTFs):
            protocol = obj.getMapper().getParent(ctf.get())
            name = "(%d) %s " % (i+1, protocol.getClassLabel())
            ticksLablesMajor.append(name)
        plotter = EmPlotter(fontsize=14)
        resolution=2
        plotter.createSubPlot("# micrographs with resolution\n lower than %d"%(resolution), 
                      "Resolution", "# of Micrographs")
#        plotter.plotMatrix(_matrix,cmap='seismic_r'
        plotter.plotMatrix(_matrix,cmap='Greens'
                        , xticksLablesMajor=ticksLablesMajor
                        , yticksLablesMajor=ticksLablesMajor)
        return plotter

