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
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.text import Text
from pyworkflow.em.protocol import ProtCreateFSC
from pyworkflow.utils.properties import Icon, Color

class FscViewer(Viewer):
    """ Viewer for plot a FSC object. """
    _targets = [FSC]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)
        self.threshold = kwargs.get('threshold', 0.143)
        self.legend=""

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
        return "1/%0.2f" % inv

    def _visualize(self, obj, figure='active', **kwargs):
        def createFSCObject(event):

            prot = self.project.newProtocol(ProtCreateFSC)
            prot.setObjLabel('FSC-%s'%self.label)#label)
            prot.setInputObj(self.protocol)#protocol may be not finished
            #prot.setParentObject(self.obj)#protocol may be not finished
            #prot.inputObj.set(self.protocol)#protocol may be not finished

            self.project.launchProtocol(prot)
            #prot.inputObject.set(self.obj)

        self.obj = obj
        x, y = obj.getData()
        plotter = EmPlotter(x=1, y=1, windowTitle='FSC', figure=figure)
        a = plotter.createSubPlot("FSC", "frequency (1/A)", "fsc")
        a.set_ylim([-0.1, 1.1])
        a.plot([0, x[-1]],
               [self.threshold, self.threshold],
               'k--')
        a.grid(True)
        _plotFormatter = FuncFormatter(self._formatFreq)
        a.xaxis.set_major_formatter(_plotFormatter)
        self.project = self.getProject()
        self.label = kwargs.get('label', self.protocol.getRunName())
        plotter.plotData(x, y, '-',label=self.label)
        a.legend()
        #####
        axcreateFSC = plt.axes([0.75, 0.02, 0.2, 0.050])
        #Button does not allow to define text color so
        #I write it directly
        axcreateFSC.text(0.5, 0.5, 'create FSC',
                             verticalalignment='center',
                             horizontalalignment='center',
                             transform=axcreateFSC.transAxes,color='white')
        bcreateFSC = Button(axcreateFSC, '',#leave label empty
                            color=Color.RED_COLOR,
                            hovercolor='maroon')
        bcreateFSC.on_clicked(createFSCObject)
        plt.show()
        ####
        return [plotter]
