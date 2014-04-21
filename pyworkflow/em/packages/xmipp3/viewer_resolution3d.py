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
from viewer import runShowJ
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
        form.addParam('doShowEstructureFactor', BooleanParam, default=True, 
                      label="Display Structure factor?")
#         form.addParam('doShowGuinier', BooleanParam, label="Display Guinier plot?", condition='self.protocol.doStructureFactor', default=True)
#         form.addParam('useMatlab', BooleanParam, label="Use Matlab for Guinier?", condition='self.protocol.doStructureFactor and doShowGuinier', default=True)
        form.addParam('doShowSsnr', BooleanParam, default=True, 
                      label="Display Spectral SNR?")
        form.addParam('doShowVssnr', BooleanParam, default=True, 
                      label="Display Volumetric Spectral SNR?")

    def _getVisualizeDict(self):
        return {'doShowFsc': self._viewFsc,
                'doShowDpr': self._viewDpr,
                'doShowEstructureFactor': self._viewEstructureFactor,
#                 'doShowGuinier': self._viewGuinier,
#                 'useMatlab': self._useMatlab,
                'doShowSsnr': self._viewSsnr,
                'doShowVssnr': self._viewVssnr
                }
    
    def _viewFsc(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = MetaData(fscFn)
        return [self._viewPlot("Fourier Shell Correlation", FREQ_LABEL, 'FSC', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_FRC, color='r'),
                self._showJ(fscFn)]
        
    def _viewDpr(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = MetaData(fscFn)
        return [self._viewPlot("Differential Phase Residual", FREQ_LABEL, 'DPR', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_DPR)]
    
    def _viewPlot(self, title, xTitle, yTitle, md, mdLabelX, mdLabelY, color='g'):        
        xplotter = XmippPlotter(1, 1, figsize=(4,4), windowTitle="Plot")
        xplotter.createSubPlot(title, xTitle, yTitle)
        xplotter.plotMdFile(md, mdLabelX, mdLabelY, color)
        return xplotter
        
    def _showJ(self, filename):
        return ProjectDataView(filename)
    
    def _viewEstructureFactor(self, e=None):
        strFactFn = self.protocol._defineStructFactorName()
        md = MetaData(strFactFn)
        return [self._viewPlot("Structure Factor", FREQ_LABEL, 'Structure Factor', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_STRUCTURE_FACTOR),
                self._viewPlot("Structure Factor", FREQ_LABEL, 'log(Structure Factor)', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_LOG_STRUCTURE_FACTOR),
                self._showJ(strFactFn)]        
    
    def _viewSsnr(self, e=None):
        pass
    
    def _viewVssnr(self, e=None):
        pass
    
