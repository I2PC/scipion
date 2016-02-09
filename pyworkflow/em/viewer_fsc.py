# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, Viewer
from plotter import EmPlotter
from data import FSC


class FscViewer(Viewer):
    """ Viewer to plot a FSC object. """
    _targets = [FSC]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)
        self.threshold = kwargs.get('threshold', 0.143)

    def _visualize(self, obj, **kwargs):
        x, y = obj.getData()
        plotter = EmPlotter(x=1, y=1, windowTitle='FSC')
        a = plotter.createSubPlot("FSC", "frequency (1/A)", "fsc")
        a.set_ylim([-0.1, 1.1])
        a.plot([0, x[-1]],
               [self.threshold, self.threshold],
               'b--', label=obj.getObjLabel())
        a.legend()
        plotter.plotData(x, y, '-')

        return [plotter]
