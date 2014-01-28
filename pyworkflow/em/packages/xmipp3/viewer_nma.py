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
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_nma import XmippProtNMA
from viewer import runShowJ
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
# FIXME: REMOVE DEPENDENCIES FROM XMIPP
from protlib_gui_ext import chimera, showj
import glob

CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2
   
        
class XmippNMAViewer(ProtocolViewer):
    """ Visualization of results from the NMA protocol
    """
    _label = 'viewer nma'
    _targets = [XmippProtNMA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def setProtocol(self, protocol):
        ProtocolViewer.setProtocol(self, protocol)
        isEm = not isinstance(protocol.inputStructure.get(), PdbFile)
        self.isEm.set(isEm)
        
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        #TODO: Just trying a trick to have hidden params
        form.addParam('isEm', BooleanParam, default=False, 
                      condition='False')
        form.addParam('displayPseudoAtom', BooleanParam, default=False, 
                      condition='isEm',
                      label="Display pseudoatoms representation?")
        form.addParam('displayPseudoAtomAproximation', BooleanParam, default=False,
                      condition='isEm',
                      label="Display pseudoatoms aproximation?")     
        form.addParam('displayModes', BooleanParam, default=False, 
                      label="Open list of modes?")
        form.addParam('displayMaxDistanceProfile', BooleanParam, default=False, 
                      label="Plot max distance profile?")     
        form.addParam('displayDistanceProfile', BooleanParam, default=False, 
                      label="Plot distance profile?")
                         
        form.addParam('singleMode', IntParam, default=7,
              label='Open specific mode', condition='True')    
    
    def _getVisualizeDict(self):
        return {'displayPseudoAtom': self._viewParam,
                'displayPseudoAtomAproximation': self._viewParam,
                'displayModes': self._viewParam,
                'displayMaxDistanceProfile': self._viewParam,
                'displayDistanceProfile': self._viewParam,
                'singleMode': self._viewParam,
                } 
                        
    def _viewParam(self, paramName):
        if paramName == 'displayPseudoAtom':
            chimera(self.protocol._getPath("chimera.cmd"))
        elif paramName == 'displayPseudoAtomAproximation':
            files = self.protocol.inputStructure.get().getFirstItem().getFileName()
            files += " " + self.protocol._getExtraPath('pseudoatoms_approximation.vol')
            showj(files)
        elif paramName == 'displayModes':
            showj(self.protocol._getPath('modes.xmd'))
        elif paramName == 'displayMaxDistanceProfile':
            fn = self.protocol._getExtraPath("maxAtomShifts.xmd")
            #TODO: Change this to real plot
            os.system('xmipp_metadata_plot -i %s -y nmaAtomShift --title "Maximum atom shifts" &' % fn)
        elif paramName == 'displayDistanceProfile':
            mode = self.singleMode.get()
            fn = self.protocol._getExtraPath("distanceProfiles","vec%d.xmd" % mode)
            #TODO: Change this to real plot
            os.system('xmipp_metadata_plot -i %s -y nmaAtomShift --title "Atom shifts for mode %d" &' % (fn, mode))
        elif paramName == 'singleMode':
            if self.singleMode.hasValue():
                vmdFile = self.protocol._getExtraPath("animations", "animated_mode_%03d.vmd" % self.singleMode.get())
                os.system("vmd -e %s" % vmdFile)
        
    def getVisualizeDictWeb(self):
        return {}
        
    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerNMA"
        
