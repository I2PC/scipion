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
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_ml2d import XmippProtML2D
from viewer import runShowJ
from plotter import XmippPlotter
import numpy as np



class XmippDefML2DViewer(Form):
    """Create the definition of parameters for
    the XmippProtML2D protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Visualization')
        self.addParam('doShowClasses', BooleanParam, label="Visualize last iter references", default=True, 
                      help='Visualize last iteration references.')
        self.addParam('doShowPlots', BooleanParam, label="Show all plots per iteration?", default=True)
        
        self.addSection(label='Iteration plots')    
        self.addParam('doShowLL', BooleanParam, label="Show Log-Likehood over iterations?", default=False, 
                      help='The Log-Likelihood value should increase.')      
        self.addParam('doShowPmax', BooleanParam, label="Show maximum model probability?", default=False, 
                      help='Show the maximum probability for a model, this should tend to be a deltha function.')      
        self.addParam('doShowSignalChange', BooleanParam, label="Show plot for signal change?", default=False, 
                      help='Should approach to zero when convergence.')      
        self.addParam('doShowMirror', BooleanParam, label="Show mirror fraction for last iteration?", default=False, 
                      help='he the mirror fraction of each reference in last iteration.')      
        
class XmippML2DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtML2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _definition = XmippDefML2DViewer()
    _label = 'Xmipp Viewer ML2D'
    
    def _getVisualizeDict(self):
        return {'doShowClasses': self._viewIterRefs,
                'doShowLL': self._viewPlot,
                'doShowPmax': self._viewPlot,
                'doShowSignalChange': self._viewPlot,
                'doShowMirror': self._viewPlot}
        
    def _viewAll(self, *args):
        if self.doShowClasses:
            self._viewIterRefs()
        plots = [p for p in ['doShowLL', 'doShowPmax', 'doShowSignalChange', 'doShowMirror'] if self.getAttributeValue(p)]
        self.createPlots(self.protocol, plots).show()
        
    def _viewPlot(self, paramName):
        print "viewing param: ", paramName
        xplotter = self.createPlots(self.protocol, [paramName])
        if xplotter:
            xplotter.show()
        
    def _viewIterRefs(self, e=None):
        lastIter = self.protocol._lastIteration()
        runShowJ('classes@' + self.protocol._getIterClasses())
        print "lastIter: ", lastIter
        
    def createPlots(self, protML, selectedPlots):
        ''' Launch some plot for an ML2D protocol run '''
        protML._plot_count = 0
        lastIter = protML._lastIteration()
        if lastIter == 0:
            return
        refs = protML._getIterClasses(iter=lastIter, block='classes')
    #    if not exists(refs):
    #        return 
    #    blocks = getBlocksInMetaDataFile(refs)
    #    lastBlock = blocks[-1]
        
        def doPlot(plotName):
            return plotName in selectedPlots
    
        # Remove 'mirror' from list if DoMirror is false
        if doPlot('doShowMirror') and not protML.doMirror:
            selectedPlots.remove('doShowMirror')
            
        n = len(selectedPlots)
        if n == 0:
            #showWarning("ML2D plots", "Nothing to plot", protML.master)
            print "No plots"
            return 
        elif n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [2, 2]
            
        xplotter = XmippPlotter(*gridsize)
            
        # Create data to plot
        iters = range(1, lastIter+1)
        ll = []
        pmax = []
        for iter in iters:
            logs = protML._getIterClasses(iter=iter, block='info')
            md = xmipp.MetaData(logs)
            id = md.firstObject()
            ll.append(md.getValue(xmipp.MDL_LL, id))
            pmax.append(md.getValue(xmipp.MDL_PMAX, id))
                
        if doPlot('doShowLL'):
            a = xplotter.createSubPlot('Log-likelihood (should increase)', 'iterations', 'LL', yformat=True)
            a.plot(iters, ll)
    
        #Create plot of mirror for last iteration
        if doPlot('doShowMirror'):
            from numpy import arange
            from matplotlib.ticker import FormatStrFormatter
            md = xmipp.MetaData(refs)
            mirrors = [md.getValue(xmipp.MDL_MIRRORFRAC, id) for id in md]
            nrefs = len(mirrors)
            ind = arange(1, nrefs + 1)
            width = 0.85
            a = xplotter.createSubPlot('Mirror fractions on last iteration', 'references', 'mirror fraction')
            a.set_xticks(ind + 0.45)
            a.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
            a.bar(ind, mirrors, width, color='b')
            a.set_ylim([0, 1.])
            a.set_xlim([0.8, nrefs + 1])
            
        if doPlot('doShowPmax'):
            a = xplotter.createSubPlot('Probabilities distribution', 'iterations', 'Pmax/Psum') 
            a.plot(iters, pmax, color='green')
        
        if doPlot('doShowSignalChange'):
            md = xmipp.MetaData()
            for iter in iters:
                fn = protML._getIterClasses(iter=iter, block='classes')
                md2 = xmipp.MetaData(fn)
                md2.fillConstant(xmipp.MDL_ITER, str(iter))
                md.unionAll(md2)
            # 'iter(.*[1-9].*)@2D/ML2D/run_004/ml2d_iter_refs.xmd')
            #a = plt.subplot(gs[1, 1])
            #print "md:", md
            md2 = xmipp.MetaData()    
            md2.aggregate(md, xmipp.AGGR_MAX, xmipp.MDL_ITER, xmipp.MDL_SIGNALCHANGE, xmipp.MDL_MAX)
            signal_change = [md2.getValue(xmipp.MDL_MAX, id) for id in md2]
            xplotter.createSubPlot('Maximum signal change', 'iterations', 'signal change')
            xplotter.plot(iters, signal_change, color='green')
        
        return xplotter

