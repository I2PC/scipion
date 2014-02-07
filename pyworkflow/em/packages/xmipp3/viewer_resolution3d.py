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
"""
This module implement the wrappers aroung Xmipp ML3D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_resolution_3d import *
from viewer import runShowJ
from plotter import XmippPlotter
from xmipp import *

import numpy as np

class XmippResolution3DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtResolution3D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer resolution3D'

    
    def _defineParams(self, form):
        form.addSection(label='Results')
        form.addParam('doShowFsc', BooleanParam, label="Display Fourier Shell Correlation?", default=True)
        form.addParam('doShowDpr', BooleanParam, label="Display Differential Phase Residual?", default=True)
        form.addParam('doShowEstructureFactor', BooleanParam, label="Display Structure factor?", default=True)
#         form.addParam('doShowGuinier', BooleanParam, label="Display Guinier plot?", condition='self.protocol.doStructureFactor', default=True)
#         form.addParam('useMatlab', BooleanParam, label="Use Matlab for Guinier?", condition='self.protocol.doStructureFactor and doShowGuinier', default=True)
        form.addParam('doShowSsnr', BooleanParam, label="Display Spectral SNR?", default=True)
        form.addParam('doShowVssnr', BooleanParam, label="Display Volumetric Spectral SNR?", default=True)

    def _getVisualizeDict(self):
        return {'doShowFsc': self._viewFsc,
                'doShowDpr': self._viewDpr,
                'doShowEstructureFactor': self._viewEstructureFactor,
#                 'doShowGuinier': self._viewGuinier,
#                 'useMatlab': self._useMatlab,
                'doShowSsnr': self._viewSsnr,
                'doShowVssnr': self._viewVssnr
                }
        
    def _viewAll(self, *args):
        if self.doShowFsc and self.protocol.doFSC:
            self._viewFsc().show()
        if self.doShowDpr and self.protocol.doFSC:
            self._viewDpr().show()
        if self.doShowEstructureFactor and self.protocol.doStructureFactor:
            self._viewEstructureFactor().show()
#         if self.doShowGuinier and self.protocol.doStructureFactor:
#             self._viewGuinier()
        if self.doShowSsnr and self.protocol.doSSNR:
            self._viewSsnr()
        if self.doShowVssnr and self.protocol.doVSSNR:
            self._viewVssnr()
    
    def _viewFsc(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = md = MetaData(fscFn)
        self._viewPlot("Fourier Shell Correlation", 'frequency (1/A)', 'FSC', md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_FRC, color='r')
        self._showJ(fscFn)
        
    def _viewDpr(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = md = MetaData(fscFn)
        self._viewPlot("Differential Phase Residual", 'frequency (1/A)', 'DPR', md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_DPR)
    
    def _viewPlot(self, title, xTitle, yTitle, md, mdLabelX, mdLabelY, color='g'):
        
        xplotter = XmippPlotter(1, 1, figsize=(4,4), windowTitle="Plot")
        xplotter.createSubPlot(title, xTitle, yTitle)
        xplotter.plotMdFile(md, mdLabelX, mdLabelY, color)
        return xplotter
        
    def _showJ(self, filename):
        runShowJ(filename)
    
    def _viewEstructureFactor(self, e=None):
        
        strFactFn = self.protocol._defineStructFactorName()
        md = md = MetaData(strFactFn)
        self._viewPlot("Structure Factor", 'frequency (1/A)', 'Structure Factor', md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_STRUCTURE_FACTOR)
        self._viewPlot("Structure Factor", 'frequency (1/A)', 'log(Structure Factor)', md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_LOG_STRUCTURE_FACTOR)
        self._showJ(strFactFn)
        
    
    def _viewSsnr(self, e=None):
        pass
    
    def _viewVssnr(self, e=None):
        pass

    def getVisualizeDictWeb(self):
        return {'doShowFsc': '_viewFsc',
                'doShowDpr': '_viewDpr',
                'doShowEstructureFactor': '_viewEstructureFactor',
                'doShowSsnr': '_viewSsnr',
                'doShowVssnr': '_viewVssnr'
                }

    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerResolution3D"
    
