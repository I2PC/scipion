# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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
This module implement the wrappers aroung Xmipp ML2D protocol
visualization program.
"""
from pyworkflow.viewer import Viewer
from pyworkflow.em import *
from protocol_ml2d import XmippProtML2D
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.gui.form import FormWindow


class XmippDefML2DViewer(Form):
    """Create the definition of parameters for
    the XmippProtML2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Visualization')
        self.addParam('lastIterRefs', BooleanParam, label="Visualize last iter references", default=True, 
                      help='Visualize last iteration references.')      
        self.addParam('itersLL', BooleanParam, label="Show Log-Likehood over iterations?", default=False, 
                      help='The Log-Likelihood value should increase.')      
        self.addParam('maxModelProb', BooleanParam, label="Show maximum model probability?", default=False, 
                      help='Show the maximum probability for a model, this should tend to be a deltha function.')      
        self.addParam('signalChange', BooleanParam, label="Show plot for signal change?", default=False, 
                      help='Should approach to zero when convergence.')      
        self.addParam('mirrorFracLastIter', BooleanParam, label="Show mirror fraction for last iteration?", default=True, 
                      help='he the mirror fraction of each reference in last iteration.')      
        
class XmippML2DViewer(Viewer, EMProtocol):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtML2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _definition = XmippDefML2DViewer()
    _label = 'Xmipp Viewer ML2D'
    
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        """Open the Protocol GUI Form given a Protocol instance"""
        #hosts = [host.getLabel() for host in self.settings.getHosts()]
        self.allowHeader.set(False)
        w = FormWindow("Protocol Run: " + self.getClassName(), self, 
                       self._viewAll, args['windows'], 
                       
                       visualizeDict={'lastIterRefs': self._viewParam,
                                      'itersLL': self._viewParam})
        w.visualizeMode = True
        w.show(center=True)
        
    def _viewAll(self, *args):
        print "HERE ALL SELECTED VISUALIZATORS WILL BE OPEN"
        
    def _viewParam(self, paramName):
        print "viewing param: ", paramName
        

