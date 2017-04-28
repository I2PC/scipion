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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from protocol_refmac import CCP4ProtRunRefmac
from pyworkflow.em import data, TextFileView
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em.viewer import ObjectView
from pyworkflow.em.showj import MODE, MODE_MD, ORDER, VISIBLE
from pyworkflow.protocol.params import FloatParam, IntParam, LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer, Viewer
from pyworkflow.utils.path import cleanPath
import numpy as np
import os


class CCP4ProtRunRefmacViewer(ProtocolViewer):
    """ This protocol computes the maximum resolution up to which two
     CTF estimations would be ``equivalent'', defining ``equivalent'' as having
      a wave aberration function shift smaller than 90 degrees
    """
    _label = 'viewer Refmac results'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CCP4ProtRunRefmac]
    _memory = False
    resolutionThresholdOLD = -1
    # temporary metadata file with ctf that has some resolution greathan than X
    tmpMetadataFile = 'viewersTmp.sqlite'

    def _defineParams(self, form):
        form.addSection(label='Visualization of Refmac results')
        # group = form.addGroup('Overall results')
        form.addParam('displayMask', LabelParam,
                      label="Masked volume",
                      help="Display FSC of masked map")
        form.addParam('showFinalResults', LabelParam,
                      label="Final Results Table",
                      help="Table of Final Results from refine.log file")
        form.addParam('showLogFile', LabelParam,
                      label = "Show log file",
                      help="refine.log file contains all results from Refmac")
        form.addParam('showLastIteration', LabelParam,
                      label = "Last Iteration Results Table",
                      help= "Table of Refmac results obtained from the last iteration")
        form.addParam('displayRFactorPlot', LabelParam,
                      label = "R-factor vs. iteration",
                      help= "Plot R-factor as a function of the iteration")
        form.addParam('displayFOMPlot', LabelParam,
                      label="FOM vs. iteration",
                      help= "Plot Figure Of Merit as a function of the iteration")
        form.addParam('displayLLPlot', LabelParam,
                      label="-LL vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayLLfreePlot', LabelParam,
                      label="-LLfree vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayGeometryPlot', LabelParam,
                      label="Geometry vs. iteration",
                      help="Plot Geometry as a function of the iteration:\n"  "Geometry includes rmsBOND (root mean square bond lengths),\n"
                           "zBOND (zscore of the deviation of bond lengths),\n rmsANGL (root mean square bond angles),\n" 
                           "zANGL (zscore of the deviation of bond angles),\n and rmsCHIRAL (root mean square of chiral index")

    def _getVisualizeDict(self):
        return {
            'displayMask': self._visualizeMask,
            'showFinalResults': self._visualizeFinalResults,
            'showLogFile': self._visualizeLogFile,
            'showLastIteration': self._visualizeLastIteration,
            'displayRFactorPlot': self._visualizeRFactorPlot,
            'displayFOMPlot': self._visualizeFOMPlot,
            'displayLLPlot': self._visualizeLLPlot,
            'displayLLfreePlot': self._visualizeLLfreePlot,
            'displayGeometryPlot': self._visualizeGeometryPlot
        }

    def _visualizeMask(self):
        pass

    def _visualizeFinalResults(self, e=None):
        views = []
        numberOfBins = self.visualizeHistogram.get()
        numberOfBins = min(numberOfBins, self.protocol.outputCTF.getSize())
        plotter = EmPlotter()
        plotter.createSubPlot("Resolution Discrepancies histogram",
                              "Resolution (A)", "# of Comparisons")
        resolution = [ctf._discrepancy_resolution.get() for ctf in self.protocol.outputCTF]
        # print "resolution", resolution
        # print "numberOfBins",numberOfBins,self.protocol.outputCTF.getSize()
        plotter.plotHist(resolution, nbins=numberOfBins)
        # ROB: why do I need this show?
        plotter.show()
        return views.append(plotter)

    def _visualizeLogFile(self, e=None):
        """ Check if the refine.log file is generated and if so, read the file."""
        refineLogFileName = self.protocol._getExtraPath(self.protocol.refineLogFileName)
        view = TextFileView(refineLogFileName)
        view.show()

    def _visualizeLastIteration(self, e=None):
        pass

    def _visualizeRFactorPlot(self, e=None):
        pass

    def _visualizeFOMPlot(self, e=None):
        pass

    def _visualizeLLPlot(self, e=None):
        pass

    def _visualizeLLfreePlot(self, e=None):
        pass

    def _visualizeGeometryPlot(self, e=None):
        pass

