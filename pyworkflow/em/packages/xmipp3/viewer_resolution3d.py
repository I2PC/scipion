# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_resolution3d import *
from plotter import XmippPlotter
from xmipp import *
import numpy as np


FREQ_LABEL = 'frequency (1/A)'


class XmippResolution3DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer resolution3D'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtResolution3D]
    
    def _defineParams(self, form):
        form.addSection(label='Results')
        form.addParam('doShowFsc', BooleanParam, default=True, 
                      label="Display Fourier Shell Correlation?")
        form.addParam('doShowDpr', BooleanParam, default=True, 
                      label="Display Differential Phase Residual?")
        form.addParam('doShowStructureFactor', BooleanParam, default=True, 
                      label="Display Structure factor?")

    def _getVisualizeDict(self):
        return {'doShowFsc': self._viewFsc,
                'doShowDpr': self._viewDpr,
                'doShowStructureFactor': self._viewStructureFactor,
                }
    
    def _viewFsc(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = MetaData(fscFn)
        return [self._createPlot("Fourier Shell Correlation", FREQ_LABEL, 'FSC', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_FRC, color='r'),
                DataView(fscFn)]
        
    def _viewDpr(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = MetaData(fscFn)
        return [self._createPlot("Differential Phase Residual", FREQ_LABEL, 'DPR', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_DPR),
                DataView(fscFn)]
    
    def _createPlot(self, title, xTitle, yTitle, md, mdLabelX, mdLabelY, color='g'):        
        xplotter = XmippPlotter(1, 1, figsize=(4,4), windowTitle="Plot")
        xplotter.createSubPlot(title, xTitle, yTitle)
        xplotter.plotMdFile(md, mdLabelX, mdLabelY, color)
        return xplotter

    
    def _viewStructureFactor(self, e=None):
        strFactFn = self.protocol._defineStructFactorName()
        md = MetaData(strFactFn)
        return [self._createPlot("Structure Factor", 'frequency^2 (1/A^2)', 'log(Structure Factor)', 
                               md, xmipp.MDL_RESOLUTION_FREQ2, xmipp.MDL_RESOLUTION_LOG_STRUCTURE_FACTOR)]        
