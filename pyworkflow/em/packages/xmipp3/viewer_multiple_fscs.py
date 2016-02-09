# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
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

import os

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView, ChimeraView
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
from plotter import XmippPlotter

from protocol_multiple_fscs import XmippProtMultipleFSCs



class XmippMultipleFSCsViewer(XmippViewer):
    """ Visualize the output of protocol multiple fscs """
    _label = 'viewer multiple fscs'
    _targets = [XmippProtMultipleFSCs]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

    def _visualize(self, obj, **args):
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        xplotter = XmippPlotter(windowTitle="FSC")
        a = xplotter.createSubPlot("FSC", "Frequency (1/A)", "FSC")
        legends = []
        for i, vol in enumerate(self.protocol.inputVolumes):
            index = i + 1
            fnFSC = self.protocol._getExtraPath("volume_%02d.fsc" % index)
            if os.path.exists(fnFSC):
                label =  vol.get().getObjLabel() or 'Volume %d' % index
                self._plotFSC(a, fnFSC, label)
        a.plot([self.minInv, self.maxInv],[0.143, 0.143], color='black', linestyle='--')
        a.plot([self.minInv, self.maxInv],[0.5, 0.5], color='black', linestyle='--')
        a.grid(True)
        #xplotter.legend()
        a.legend()

        self._views.append(xplotter)

    def _plotFSC(self, a, fnFSC, label):
        from xmipp import MetaData, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_FRC
        md = MetaData(fnFSC)
        resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, f) for f in md]
        frc = [md.getValue(MDL_RESOLUTION_FRC, f) for f in md]
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, frc, label=label)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999
        if value:
            inv = 1/value
        return "1/%0.2f" % inv
