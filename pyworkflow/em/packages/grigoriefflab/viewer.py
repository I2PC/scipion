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
Visualization of the results of the Frealign protocol.
"""
from os.path import exists
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.em.plotter import EmPlotter

from protocol_refinement import ProtFrealign
from protocol_ml_classification import ProtFrealignClassify
from protocol_ctffind3 import ProtCTFFind, ProtRecalculateCTFFind


LAST_ITER = 0
ALL_ITER = 1
SELECTED_ITERS = 2


class FrealignViewer(ProtocolViewer):
    """ Visualization of Frealign."""    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Frealign'
    _targets = [ProtFrealign, ProtFrealignClassify]    
    
    def _defineParams(self, form):
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam, label="Which iteration do you want to visualize?", default=0, 
                      choices=['last','all','selection'], display=EnumParam.DISPLAY_LIST)
        form.addParam('selectedIters', StringParam, default='',
              label='Selected iterations', condition='iterToShow==2',
              help="""
*last*: only the last iteration will be visualized.
*all*: all iterations  will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                   """)
        form.addParam('doShow3DRefsVolumes', BooleanParam, label="Visualize the 3D-references volumes?", default=True)
        form.addParam('doShow3DReconsVolumes', BooleanParam, label="Visualize the 3D-reconstructed volumes?", default=True)
        form.addParam('doShow3DMatchProj', BooleanParam, label="Visualize the matching projections of refinement?", default=True)
        form.addParam('doShowAngDist', BooleanParam, label="Plot angular distribution?", default=True)
#         form.addParam('doShowDataDist', BooleanParam, label="Plot data distribution over 3d-references?", default=True)
        
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
        pass
#         if self.doShow3DRefsVolumes:
#             self._view3DRefsVolumes()
#         if self.doShow3DReconsVolumes:
#             self._view3DReconVolumes()
#         if self._doShow3DMatchProj:
#             self._viewMatchProj()
#         if self.doShowAngDist:
#             self._plotAngularDistribution()
    
    def _view3DRefsVolumes(self, e=None):
        files = self._getIterationFile("reference_volume_iter_%03d.mrc")
        path = self.protocol._getExtraPath('viewer_refvolumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project.getName(), self.protocol.strId(), path)]
        
    def _view3DReconVolumes(self, e=None):
        files = self._getIterationFile("volume_iter_%03d.mrc")
        path = self.protocol._getExtraPath('viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project.getName(), self.protocol.strId(), path)]
    
    def _viewMatchProj(self, e=None):
        files = self._getIterationFile("particles_match_iter_%03d.mrc")
        return [self.createDataView(files[0])]
    
    def _plotAngularDistribution(self, e=None):
        """ Plot the angular distributions for each reference and each iteration.
        Returns:
            views: a list of all angular distribution plots or some errors.
        """
        views = []
        errors = self.setVisualizeIterations()
        
        if len(errors) == 0:
            for iteration in self.visualizeIters:
                pathFile = self.protocol._getExtraPath("iter_%03d" % iteration, "particles_iter_%03d.par" % iteration)
                print pathFile
                
                if not os.path.exists(pathFile):
                    errors.append('Iteration %s does not exist.' % iteration)
                else:
                    xplotter = self._createIterAngularDistributionPlot(iteration, pathFile)
                    views.append(xplotter)
                    
        if errors:
            views.append(self.errorMessage('\n'.join(errors), "Visualization errors"))

        return views
    
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
    
    
    def createDataView(self, filename):
        return DataView(filename)
    
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


class ProtCTFFindViewer(Viewer):
    """ Visualization of Frealign."""    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer CtfFind'
    _targets = [ProtCTFFind, ProtRecalculateCTFFind]
    
    
    def __init__(self, **args):
        Viewer.__init__(self, **args)
        self._views = []
        
    def visualize(self, obj, **args):
        self._visualize(obj, **args)
        
        for v in self._views:
            v.show()
            
    def _visualize(self, obj, **args):
        cls = type(obj)
        
        def _getMicrographDir(mic):
            """ Return an unique dir name for results of the micrograph. """
            return obj._getExtraPath(removeBaseExt(mic.getFileName()))        
        
        def iterMicrographs(mics):
            """ Iterate over micrographs and yield
            micrograph name and a directory to process.
            """
            for mic in mics:
                micFn = mic.getFileName()
                micDir = _getMicrographDir(mic) 
                yield (micFn, micDir, mic)
        
        def visualizeObjs(obj, setOfMics):
                
                if exists(obj._getPath("ctfs_temporary.sqlite")):
                    os.remove(obj._getPath("ctfs_temporary.sqlite"))
                
                ctfSet = self.protocol._createSetOfCTF("_temporary")
                for fn, micDir, mic in iterMicrographs(setOfMics):
                    out = self.protocol._getCtfOutPath(micDir)
                    psdFile = self.protocol._getPsdPath(micDir)
                    
                    if exists(out) and exists(psdFile):
                        result = self.protocol._parseOutput(out)
                        defocusU, defocusV, defocusAngle = result
                        # save the values of defocus for each micrograph in a list
                        ctfModel = self.protocol._getCTFModel(defocusU, defocusV, defocusAngle, psdFile)
                        ctfModel.setMicrograph(mic)
                        ctfSet.append(ctfModel)
                
                if ctfSet.getSize() < 1:
                    raise Exception("Has not been completed the CTT estimation of any micrograph")
                else:
                    ctfSet.write()
                    ctfSet.close()
                    self._visualize(ctfSet)
        
        if issubclass(cls, ProtCTFFind) and not obj.hasAttribute("outputCTF"):
            mics = obj.inputMicrographs.get()
            visualizeObjs(obj, mics)

        elif issubclass(cls, ProtRecalculateCTFFind) and not obj.hasAttribute("outputCTF"):
            
            mics = obj.inputCtf.get().getMicrographs()
            visualizeObjs(obj, mics)
        
        elif obj.hasAttribute("outputCTF"):
            self._visualize(obj.outputCTF)
            
        else:
            fn = obj.getFileName()
            psdLabels = '_psdFile'
            labels = 'id enabled comment %s _defocusU _defocusV _defocusAngle _defocusRatio _micObj._filename' % psdLabels 
            self._views.append(ObjectView(self._project.getName(), obj.strId(), fn,
                                          viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels, ZOOM: 50, RENDER: psdLabels}))
            
        return self._views
