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
This module implement the wrappers around Xmipp NMA protocol
visualization program.
"""

from pyworkflow.gui.project import ProjectWindow
from pyworkflow.protocol.params import LabelParam, IntParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import ObjectView, VmdView
from protocol_nma import XmippProtNMA
from plotter import XmippNmaPlotter
import xmipp

OBJCMD_NMA_PLOTDIST = "Plot distance profile"
OBJCMD_NMA_VMD = "Display VMD animation"


class XmippNMAViewer(ProtocolViewer):
    """ Visualization of results from the NMA protocol.    
        Normally, NMA modes with high collectivity and low NMA score are preferred.
    """
    _label = 'viewer nma'
    _targets = [XmippProtNMA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

#     def setProtocol(self, protocol):
#         ProtocolViewer.setProtocol(self, protocol)
#         inputPdb = protocol.inputStructure.get()
#         self.isEm.set(inputPdb.getPseudoAtoms())
        
    def _defineParams(self, form):
        form.addSection(label='Visualization')
  
        form.addParam('displayModes', LabelParam,
                      label="Display output Normal Modes?", important=True)
        form.addParam('displayMaxDistanceProfile', LabelParam,
                      label="Plot max distance profile?",
                      help="TODO: ADD HELP ABOUT MAX DISTANCE") 
        
        group = form.addGroup('Single mode')  
        group.addParam('modeNumber', IntParam, default=7,
              label='Mode number')
        group.addParam('displayVmd', LabelParam,
                       label='Display mode animation with VMD?') 
        group.addParam('displayDistanceProfile', LabelParam, default=False,
                      label="Plot mode distance profile?")
        
    def _getVisualizeDict(self):
        return {'displayModes': self._viewParam,
                'displayMaxDistanceProfile': self._viewParam,
                'displayVmd': self._viewSingleMode,
                'displayDistanceProfile': self._viewSingleMode,
                } 
                        
    def _viewParam(self, paramName):
        if paramName == 'displayModes':
            modes =  self.protocol.outputModes
            return [ObjectView(self._project, modes.strId(), modes.getFileName())]
        elif paramName == 'displayMaxDistanceProfile':
            fn = self.protocol._getExtraPath("maxAtomShifts.xmd")
            return [createShiftPlot(fn, "Maximum atom shifts", "maximum shift")]
    
    def _viewSingleMode(self, paramName):
        """ visualization for a selected mode. """
        modes =  self.protocol.outputModes
        modeNumber = self.modeNumber.get()
        mode = modes[modeNumber]
        
        if mode is None:
            return [self.errorMessage("Invalid mode number *%d*\n"
                                      "Display the output Normal Modes to see "
                                      "the availables ones." % modeNumber,
                                      title="Invalid input")]
        elif paramName == 'displayVmd':
            return [createVmdView(self.protocol, modeNumber)]
        elif paramName == 'displayDistanceProfile':
            return [createDistanceProfilePlot(self.protocol, modeNumber)]


def createShiftPlot(mdFn, title, ylabel):
    plotter = XmippNmaPlotter()
    plotter.createSubPlot(title, 'atom index', ylabel)
    plotter.plotMdFile(mdFn, None, xmipp.MDL_NMA_ATOMSHIFT)
    return plotter


def createDistanceProfilePlot(protocol, modeNumber):
    vectorMdFn = protocol._getExtraPath("distanceProfiles","vec%d.xmd"
                                        % modeNumber)
    plotter = createShiftPlot(vectorMdFn, "Atom shifts for mode %d"
                              % modeNumber, "shift")
    return plotter


def createVmdView(protocol, modeNumber):
    vmdFile = protocol._getExtraPath("animations", "animated_mode_%03d.vmd"
                                     % modeNumber)
    return VmdView('-e "%s"' % vmdFile)


def showDistanceProfilePlot(protocol, modeNumber):
    createDistanceProfilePlot(protocol, modeNumber).show()


def showVmdView(protocol, modeNumber):
    createVmdView(protocol, modeNumber).show()


ProjectWindow.registerObjectCommand(OBJCMD_NMA_PLOTDIST,
                                    showDistanceProfilePlot)
ProjectWindow.registerObjectCommand(OBJCMD_NMA_VMD,
                                    showVmdView)

