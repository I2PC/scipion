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
This module is mainly for the Viewer class, which 
serve as base for implementing visualization tools(Viewer sub-classes).
"""

from os.path import join
from protocol import Protocol

DESKTOP_TKINTER = 'tkinter'
WEB_DJANGO = 'django'


class Viewer(object):
    """ All visualization wrappers should user the Viewer class
    as base and provide the implementation to launch specific 
    command line tools in order to visualize objects.
    
    The _targets class property should contains a list of string
    with the class names that this viewer is able to visualize.
    For example: _targets = ['Image', 'SetOfImages']
    """
    _targets = []
    _environments = [DESKTOP_TKINTER]
    
    def __init__(self, tmpPath='./Tmp', **args):
        self._tmpPath = tmpPath
        
    def _getTmpPath(self, *paths):
        return join(self._tmpPath, *paths)
    
    def visualize(self, obj):
        """ This method should make the necessary convertions
        and call the command line utilities to visualize this
        particular object.
        """
        pass
    
    def getView(self):
        """ This method should return the string value of the view in web
        that will respond to this viewer. This method only should be implemented
        in those viewers that have WEB_DJANGO environment defined. 
        """
        return None


class Wizard(object):
    """ This is a special case of GUI to help the user
    selecting important parameters.
    The _targets will serve to define to which Definition and 
    parameters the Wizard is defined, it will be a list of tuples such as:
    _targets = [(DefImportMicrographs, ['voltage', sphericalAberration']),
                (DefCTFMicrographs, ['lowRes', 'highRes'])]
    The _environmets will serve to define when this wizard can be used.
    For example>
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    """
    _targets = []
    _environments = [DESKTOP_TKINTER]
    
    def show(self, form, *params):
        """ This will show up the wizard to select parametes.
        Params:
            form: the protocol form, given access to to all parameters.
                Some times the same wizard will modifify several elements
                in the form.
            *params: a list of params to modify, sometimes the wizard can 
                be generic and can be used for different parameters in the
                same form.
        """
        pass
    
    def getView(self):
        """ This method should return the string value of the view in web
        that will respond to this wizard. This method only should be implemented
        in those wizards that have WEB_DJANGO environment defined. 
        """
        return None
    
    
class ProtocolViewer(Protocol, Viewer):
    """ This class will serve as base for viwers that will have form and parameters. """
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        Viewer.__init__(self, **args)
        self.allowHeader.set(False)
    
    def visualize(self, obj, **args):
        """Open the Protocol GUI Form given a Protocol instance"""
        from gui.form import FormWindow
        self.protocol = obj
        self.windows = args.get('windows', None)
        self.formWindow = FormWindow("Protocol Viewer: " + self.getClassName(), self, 
                       self._viewAll, self.windows,
                       visualizeDict=self._getVisualizeDict(),
                       visualizeMode=True)
        self.formWindow.visualizeMode = True
        self.formWindow.show(center=True)     

    def _getVisualizeDict(self):
        """ Create the visualization dict for view individual params. """
        return {}
    
    def _viewAll(self, *args):
        """ Visualize all data give the parameters. """
        pass
    
    #TODO: This method should not be necessary, instead NumericListParam should return a list and not a String 
    def _getListFromRangeString(self, rangeStr):
        ''' Create a list of integer from a string with range definitions
        Examples:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]
        '''
        elements = rangeStr.split(',')
        values = []
        for e in elements:
            if '-' in e:
                limits = e.split('-')
                values += range(int(limits[0]), int(limits[1])+1)
            else:
                # If values are separated by comma also splitted 
                values += map(int, e.split())
        return values
    
    def getClassName(self):
        return self.__class__.__name__

def createPlots(protML, selectedPlots):
    ''' Launch some plot for an ML2D protocol run '''
    from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
    import xmipp
    
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
