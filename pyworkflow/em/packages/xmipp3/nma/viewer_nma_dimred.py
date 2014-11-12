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

import os
from os.path import join, dirname, basename

from pyworkflow.viewer import (ProtocolViewer, CommandView,
                               DESKTOP_TKINTER, WEB_DJANGO)
from pyworkflow.protocol.params import StringParam, BooleanParam
from protocol_nma_dimred import XmippProtDimredNMA
import xmipp

from plotter import XmippNmaPlotter


        
class XmippDimredNMAViewer(ProtocolViewer):
    """ Visualization of results from the NMA protocol
    """
    _label = 'viewer nma dimred'
    _targets = [XmippProtDimredNMA]
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
        
    def _getVisualizeDict(self):
        return {'displayRawDeformation': self._viewRawDeformation,
                } 
                        
    def _viewRawDeformation(self, paramName):
        components = self.displayRawDeformation.get()
        return self._doViewRawDeformation(components)
        
    def _doViewRawDeformation(self, components):
#        components = map(int, self.displayRawDeformation.get().split())
        components = map(int, components.split())
        dim = len(components)
        views = []
        
        if dim > 0:
            modeList = [m - 1 for m in components]
            modeNameList = ['Mode %d' % m for m in components]
            missingList = []
                    
            if missingList:
                return [self.errorMessage("Invalid mode(s) *%s*\n." % (', '.join(missingList)), 
                              title="Invalid input")]
            
            defFn = self.protocol.getOutputMatrixFile()
            
            # Actually plot
            plotter = XmippNmaPlotter(defFn) 
            baseList = [basename(n) for n in modeNameList]
            
            if dim == 1:
                plotter.plotArray1D("Histogram for %s" % baseList[0], modeList[0], 
                                    "Deformation value", "Number of images")
            elif dim == 2:
                plotter.plotArray2D("%s vs %s" % tuple(baseList), 
                                    modeList[0], modeList[1], *baseList)
            elif dim == 3:
                plotter.plotArray3D("%s %s %s" % tuple(baseList), 
                                    modeList[0], modeList[1], modeList[2], *baseList)
            views.append(plotter)
            
        return views

    
