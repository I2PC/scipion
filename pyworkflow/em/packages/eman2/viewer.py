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
This module implement the first version of viewers using 
around xmipp_showj visualization program.
"""
import os
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER,
                               WEB_DJANGO, Viewer)
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
import pyworkflow.em as em
# from pyworkflow.em.plotter import EmPlotter
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (LabelParam, NumericRangeParam,
                                        EnumParam, BooleanParam, FloatParam)

from protocol_boxing import EmanProtBoxing
from protocol_initialmodel import EmanProtInitModel
from protocol_refineasy import EmanProtRefine

LAST_ITER = 0
ALL_ITERS = 1
SELECTED_ITERS = 2

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

FSC_UNMASK = 0
FSC_MASK = 1
FSC_MASKTIGHT = 2


class EmanViewer(XmippViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER]
    _targets = [EmanProtBoxing, EmanProtInitModel]
 
    def _visualize(self, obj, **args):
         
        if isinstance(obj, EmanProtBoxing):
            XmippViewer._visualize(self, obj.outputCoordinates)
             
        elif isinstance(obj, EmanProtInitModel):
            XmippViewer._visualize(self, obj.outputVolumes)


class RefineEasyViewer(ProtocolViewer):
    """ Visualization of Refine Easy."""
    
    _targets = [EmanProtRefine]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Refine Easy'
    
    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam, label="Which iteration do you want to visualize?", default=0, 
                      choices=['last','all','selection'], display=EnumParam.DISPLAY_LIST)
        form.addParam('iterSelection', NumericRangeParam, default='1',
              label='Selected iterations', condition='iterToShow==%d' % SELECTED_ITERS,
              help="""
*last*: only the last iteration will be visualized.
*all*: all iterations  will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                    
                   """)
        group = form.addGroup('Particles')
        
        group.addParam('showImagesAngularAssignment', LabelParam,
                       label='Particles angular assignment')
        
        group = form.addGroup('Volumes')
        
        group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'], 
              default=VOLUME_SLICES, display=EnumParam.DISPLAY_COMBO, 
              label='Display volume with',
              help='*slices*: display volumes as 2D slices along z axis.\n'
                   '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayAngDist', EnumParam, choices=['2D plot', 'chimera'], 
              default=ANGDIST_2DPLOT, display=EnumParam.DISPLAY_COMBO, 
              label='Display angular distribution',
              help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                   '*chimera*: display angular distribution using Chimera with red spheres.') 
        
        group = form.addGroup('Resolution plots')
        
        group.addParam('resolutionPlotsFSC', EnumParam, choices=['unmasked', 'masked', 'masked tight'], 
              default=FSC_UNMASK, display=EnumParam.DISPLAY_COMBO, 
              label='Display resolution plots (FSC)',
              help='*unmasked*: display FSC of unmasked maps.\n'
                   '*masked*: display FSC of masked maps.\n'
                   '*masked tight*: display FSC of masked tight maps.') 
        
        group.addParam('resolutionThresholdFSC', FloatParam, default=0.143, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold in resolution plots',
                      help='')
    
    def _getVisualizeDict(self):
        self._load()
        return {'showImagesAngularAssignment' : self._showImagesAngularAssignment,
#                 'displayVol': self._showVolumes,
#                 'displayAngDist': self._showAngularDistribution,
#                 'resolutionPlotsSSNR': self._showSSNR,
#                 'resolutionPlotsFSC': self._showFSC
                }
    
#===============================================================================
# showImagesAngularAssignment     
#===============================================================================

    def _showImagesAngularAssignment(self, paramName=None):
        views = []
        
        for it in self._iterations:
            print "iter: ", it, self._iterations
            fn = self.protocol._getIterData(it)
            v = self.createScipionPartView(fn)
            views.append(v)
        
        return views
    
    #--------------------------- UTILS functions --------------------------------------------
    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        
        labels =  'enabled id _size _filename _transform._matrix'
        viewParams = {em.ORDER:labels,
                      em.VISIBLE: labels, em.RENDER:'_filename',
                      'labels': 'id',
                      }
        return em.ObjectView(self._project, 
                          self.protocol.strId(), filename, other=inputParticlesId,
                          env=self._env, viewParams=viewParams)
    
    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self.protocol._createFilenameTemplates()
        self.protocol._createIterTemplates(self.protocol._getRun())
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()
        
        if self.iterToShow.get() == LAST_ITER:
            self._iterations = [self.lastIter]
        elif self.iterToShow.get() == ALL_ITERS:
            self._iterations = range(1, self.lastIter + 1)
        elif self.iterToShow.get() == SELECTED_ITERS:
            self._iterations = self._getListFromRangeString(self.iterSelection.get())
            
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq)
    
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
        return "1/%0.2f" % inv
