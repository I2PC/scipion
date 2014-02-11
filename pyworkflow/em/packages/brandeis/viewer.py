# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This module implement the wrappers around xmipp_showj
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_frealign import ProtFrealign
from pyworkflow.em.packages.xmipp3.viewer import runShowJ
from pyworkflow.em.plotter import EmPlotter

import numpy as np

LAST_ITER = 0
ALL_ITER = 1
SELECTED_ITERS = 2

class FrealignViewer(ProtocolViewer):
    """ Visualization of Frealign."""
    
    _targets = [ProtFrealign]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer Frealign'
    
    def _defineParams(self, form):
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam, label="Which iteration do you want to visualize?", default=0, 
                      choices=['last','all','selection'], display=EnumParam.DISPLAY_LIST)
        form.addParam('selectedIters', StringParam, default='',
              label='Selected iterations', condition='iterToShow==2',
              help='Which iteration do you want to visualize.')     
        form.addParam('doShow3DRefsVolumes', BooleanParam, label="Visualize the 3D-references volumes?", default=True)
        form.addParam('doShow3DReconsVolumes', BooleanParam, label="Visualize the 3D-reconstructed volumes?", default=True)
        form.addParam('doShow3DMatchProj', BooleanParam, label="Visualize the matching projections of refinement?", default=True)
        form.addParam('doShowAngDist', BooleanParam, label="Plot angular distribution?", default=True)
        form.addParam('doShowDataDist', BooleanParam, label="Plot data distribution over 3d-references?", default=True)     
        
#         form.addSection(label='Overall Results')
#         form.addParam('doShowStatistics', BooleanParam, label="Plot overall convergence statistics?", default=True)
    
    def _getVisualizeDict(self):
        return {'doShow3DRefsVolumes': self._view3DRefsVolumes,
                'doShow3DReconsVolumes': self._view3DReconVolumes,
                'doShow3DMatchProj': self._viewMatchProj,
                'doShowAngDist': self._plotAngularDistribution,
#                'doShowDataDist': self_showDataDist
                }
    
    def _viewAll(self, *args):
        if self.doShow3DRefsVolumes:
            self._view3DRefsVolumes()
        if self.doShow3DReconsVolumes:
            self._view3DReconVolumes()
        if self._doShow3DMatchProj:
            self._viewMatchProj()
        if self.doShowAngDist:
            self._plotAngularDistribution()
    
    def _view3DRefsVolumes(self, e=None):
        files = self._getIterationFile("reference_volume_iter_%03d.mrc")
        self._doShowJ(files)
        
        #self._viewIterationFile("reference_volume_iter_%03d.mrc")
        
    def _view3DReconVolumes(self, e=None):
        files = self._getIterationFile("volume_iter_%03d.mrc")
        self._doShowJ(files)
    
    def _viewMatchProj(self, e=None):
        files = self._getIterationFile("particles_match_iter_%03d.mrc")
        self._doShowJ(files)
    
    def _plotAngularDistribution(self, e=None):
        plots, errors = self._createAngularDistributionPlots()
        self._showPlots(plots, errors)
        
            
    def _createAngularDistributionPlots(self):
        """ Plot the angular distributions for each reference and each iteration.
        Returns:
            plots: a list of all angular distribution plots.
            errors: a list of errors (if some iteration does not exists.
        """
        plots = []
        errors = self.setVisualizeIterations()
        
        if len(errors) == 0:
            for iteration in self.visualizeIters:
                pathFile = self.protocol._getExtraPath("iter_%03d" % iteration, "particles_iter_%03d.par" % iteration)
                print pathFile
                
                if not os.path.exists(pathFile):
                    errors.append('Iteration %s does not exist.' % iteration)
                else:
                    xplotter = self._createIterAngularDistributionPlot(iteration, pathFile)
                    xplotter.draw()
                    plots.append(xplotter)

        return plots, errors
    
    def _createIterAngularDistributionPlot(self, iteration, pathFile):
        # Create Angular plot for one iteration
        file = open(pathFile)
        phi = []
        theta = []
        
        for line in file:
            if not line.startswith('C'):
                lineList = line.split()
                phi.append(float(lineList[3]))
                theta.append(float(lineList[2]))
        gridsize = [1, 1]
        figsize = (4, 4)
        xplotter = EmPlotter(*gridsize, figsize=figsize, 
                                windowTitle="Angular distribution - iteration %d" % iteration)
        plot_title = 'iter %d' % iteration
        
        xplotter.plotAngularDistribution(plot_title, phi, theta)
        
        return xplotter
    
    
    def _doShowJ(self, files):
        if len(files) != 0:
            for f in files:
                runShowJ(f)
    
    def _getIterationFile(self, filePath):
        self.setVisualizeIterations()
        
        path = []
        for i, iter in enumerate(self.visualizeIters):
            pathDir = self.protocol._getExtraPath("iter_%03d" % iter)
            pathNew = join(pathDir, filePath % iter)
#             print "path=%s" % pathNew
            
            if os.path.exists(pathNew):
                path.append(pathNew)
#                runShowJ(path)
            else:
                self.formWindow.showError('Iteration %s does not exist.' % iter)
        
        return path
        
        
    def setVisualizeIterations(self):
        '''Validate and set the set of iterations to visualize.
        If not set is selected, only use the last iteration'''
        self.lastIter = self.protocol.numberOfIterations.get()
        self.visualizeIters = []
        
        if self.iterToShow.get() == LAST_ITER:
            self.visualizeIters = [self.lastIter]
            
        elif self.iterToShow.get() == ALL_ITER:
            self.visualizeIters = range(1, self.lastIter + 1)
        elif self.iterToShow.get() == SELECTED_ITERS:
            if self.selectedIters.empty():
                return ['Please select the iterations that you want to visualize.']
            else:
                try:
                    self.visualizeIters = self._getListFromRangeString(self.selectedIters.get())
                except Exception, ex:
                    return ['Invalid iterations range.']
        return [] # No errors resulted

    def getVisualizeDictWeb(self):
        return {'doShow3DRefsVolumes': 'doShow3DRefsVolumes',
                'doShow3DReconsVolumes': 'doShow3DReconsVolumes',
                'doShow3DMatchProj': 'doShow3DMatchProj',
                'doShowAngDist': 'doShowAngDist',
                'doShowDataDist': 'doShowDataDist'
                }

    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerFrealign"
    
