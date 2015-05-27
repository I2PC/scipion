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
import os
from os.path import exists, relpath
from pyworkflow.utils.path import cleanPath
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO, Viewer)
import pyworkflow.em as em
import pyworkflow.em.showj as showj
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (LabelParam, NumericRangeParam,IntParam,
                                        EnumParam, BooleanParam, FloatParam)
from protocol_refinement import ProtFrealign
from protocol_ml_classification import ProtFrealignClassify
from protocol_ctffind import ProtCTFFind
from convert import readCtfModel


LAST_ITER = 0
ALL_ITERS = 1
SELECTED_ITERS = 2

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CLASSES_ALL = 0
CLASSES_SEL = 1


class FrealignViewer(ProtocolViewer):
    """ Visualization of Frealign."""    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Frealign'
    _targets = [ProtFrealign, ProtFrealignClassify]    
    
    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        form.addParam('iterToShow', EnumParam, label="Iteration to visualize", default=0,
                      choices=['last','all','selection'], display=EnumParam.DISPLAY_LIST)
        form.addParam('selectedIters', NumericRangeParam, default='1',
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
        
        if self.protocol.IS_REFINE:
            group.addParam('showImagesAngularAssignment', LabelParam,
                           label='Particles angular assignment')
            group.addParam('show3DMatchProj', LabelParam, condition=self.protocol.writeMatchProjections,
                           label="Visualize the matching projections of refinement")
        else:
            group.addParam('showImagesInClasses', LabelParam,
                          label='Particles assigned to each Class', important=True,
                          help='Display the classes and the images associated.')
        
        group = form.addGroup('Volumes')
        
        if not self.protocol.IS_REFINE:
            group.addParam('showClasses3D', BooleanParam, default=CLASSES_ALL,
                           choices=['all', 'selection'], 
                           display=EnumParam.DISPLAY_HLIST,
                           label='CLASS 3D to visualize',
                           help='')
            group.addParam('class3DSelection', NumericRangeParam, default='1',
                          condition='showClasses3D == %d' % CLASSES_SEL,
                          label='Classes list',
                          help='')
        
        group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'], 
              default=VOLUME_SLICES, display=EnumParam.DISPLAY_HLIST, 
              label='Display volume with',
              help='*slices*: display volumes as 2D slices along z axis.\n'
                   '*chimera*: display volumes as surface with Chimera.')

        group.addParam('displayAngDist', EnumParam, choices=['2D plot', 'chimera'], 
              default=ANGDIST_2DPLOT, display=EnumParam.DISPLAY_HLIST, 
              label='Display angular distribution',
              help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                   '*chimera*: display angular distribution using Chimera with red spheres.')
        group.addParam('spheresScale', IntParam, default=100, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Spheres size',
                      help='')
        
        group = form.addGroup('Resolution')
        group.addParam('resolutionPlotsSSNR', LabelParam, default=True,
                      label='Display SSNR plots',condition=self.protocol.doAditionalStatisFSC,
                      help='Display signal to noise ratio plots (SSNR) ')
        group.addParam('resolutionPlotsFSC', LabelParam, default=True,
                      label='Display resolution plots (FSC)',
                      help='')
        group.addParam('resolutionThresholdFSC', FloatParam, default=0.143, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold in resolution plots',
                      help='')
    
    def _getVisualizeDict(self):
        self._load()
        return {'showImagesInClasses': self._showImagesInClasses,
                'showImagesAngularAssignment' : self._showImagesAngularAssignment,
                'show3DMatchProj': self._viewMatchProj,
                'displayVol': self._showVolumes,
                'displayAngDist': self._showAngularDistribution,
                'resolutionPlotsSSNR': self._showSSNR,
                'resolutionPlotsFSC': self._showFSC
                }

#===============================================================================
# showImagesInClasses     
#===============================================================================
            
    def _showImagesInClasses(self, paramName=None):
        """ Read frealign .par images file and 
        generate a new SetOfClasses3D.
        If the new metadata was already written, it is just shown.
        """
        views = []
        
        for it in self._iterations:
            fn = self._getIterClasses(it)
            v = self.createScipionView(fn)
            views.append(v)
        
        return views

#===============================================================================
# showImagesAngularAssignment     
#===============================================================================

    def _showImagesAngularAssignment(self, paramName=None):
        
        views = []
        
        for it in self._iterations:
            fn = self._getIterData(it)
            v = self.createScipionPartView(fn)
            views.append(v)
        
        return views

#===============================================================================
# show3DMatchProj     
#===============================================================================
    
    def _viewMatchProj(self, paramName=None):
        views = []
        for it in self._iterations:
            files = self.protocol._getFileName('match', iter=it)
            v = self.createDataView(files)
            views.append(v)
        return views
    
#===============================================================================
# ShowVolumes
#===============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()
    
    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """

        path = self.protocol._getExtraPath('frealign_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()

        files = []
        if self.protocol.IS_REFINE:
            for it in self._iterations:
                volFn = self.protocol._getFileName('iter_vol', iter=it)
                if exists(volFn):
                    files.append(volFn)
        else:
            for it in self._iterations:
                for ref3d in self._refsList:
                    volFn = self.protocol._getFileName('iter_vol_class', iter=it, ref3d=ref3d)
                if exists(volFn):
                    files.append(volFn)
        
        self.createVolumesSqlite(files, path, samplingRate)
        return [em.ObjectView(self._project, self.protocol.strId(), path)]
    
    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        
        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for vol in volumes:
                localVol = relpath(vol, self.protocol._getExtraPath())
                if exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = em.ChimeraView(cmdFile)
        else:
            from pyworkflow.em.viewer import ChimeraClientView
            #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volumes[0])
            view = ChimeraClientView(volumes[0], showProjection=False)
            
        return [view]
    
#===============================================================================
# showAngularDistribution
#===============================================================================
    def _showAngularDistribution(self, paramName=None):
        views = []
        
        if self.displayAngDist == ANGDIST_CHIMERA:
            for it in self._iterations:
                views.append(self._createAngDistChimera(it))
                        
        elif self.displayAngDist == ANGDIST_2DPLOT:
            for it in self._iterations:
                plot = self._createAngDist2D(it)
                if isinstance(plot, EmPlotter):
                    views.append(plot)
        return views
    
    def _createAngDistChimera(self, it):
        nparts = self.protocol.inputParticles.get().getSize()
        radius = self.spheresScale.get()

        volumes = self._getVolumeNames()
        
        if len(volumes) > 1:
            raise Exception("Please, select a single volume to show it's angular distribution")
        else:
            if self.protocol.IS_REFINE:
                data_angularDist = self.protocol._getFileName("output_par", iter=it)
                if exists(data_angularDist):
                    sqliteFn = self.protocol._getFileName('projections', iter=it)
                    self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(it, data_angularDist))
                    view = em.ChimeraClientView(volumes[0], showProjection=True, angularDistFile=sqliteFn, spheresDistance=radius)
            else:
                for ref3d in self._refsList:
                    data_angularDist = self.protocol._getFileName("output_par_class", iter=it, ref=ref3d)
                    if exists(data_angularDist):
                        sqliteFn = self.protocol._getFileName('projectionsClass', iter=it, ref=ref3d)
                        self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(it, data_angularDist))
                        view = em.ChimeraClientView(volumes[0], showProjection=True, angularDistFile=sqliteFn, spheresDistance=radius)
        return view
    
    def _createAngDist2D(self, it):
        nrefs = len(self._refsList)
        nparts = self.protocol.inputParticles.get().getSize()
        gridsize = self._getGridSize(nrefs)
        
        if self.protocol.IS_REFINE:
            data_angularDist = self.protocol._getFileName("output_par", iter=it)
            if exists(data_angularDist):
                plotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                                    mainTitle="Iteration %d" % it, windowTitle="Angular distribution")
                title = 'iter %d' % it
                sqliteFn = self.protocol._getFileName('projections', iter=it)
                self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(it, data_angularDist))
                plotter.plotAngularDistributionFromMd(sqliteFn, title)
                return plotter
            else:
                return
        else:
            for ref3d in self._refsList:
                data_angularDist = self.protocol._getFileName("output_par_class", iter=it, ref=ref3d)
                if exists(data_angularDist):
                    plotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                                        mainTitle="Iteration %d" % it, windowTitle="Angular distribution")
                    title = 'class %d' % ref3d
                    sqliteFn = self.protocol._getFileName('projectionsClass', iter=it, ref=ref3d)
                    self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(it, data_angularDist))
                    plotter.plotAngularDistributionFromMd(sqliteFn, title)
            return plotter
    
#===============================================================================
# plotFSC
#===============================================================================
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        xplotter = EmPlotter(x=gridsize[0], y=gridsize[1], windowTitle='Resolution FSC')
        
        if self.protocol.IS_REFINE:
            plot_title = 'FSC'
            a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'FSC', yformat=False)
            legends = []
            
            show = False
            for it in self._iterations:
                parFn = self.protocol._getFileName('output_vol_par', iter=it)
                if exists(parFn):
                    show = True
                    self._plotFSC(a, parFn)
                    legends.append('iter %d' % it)
            xplotter.showLegend(legends)
            
            if show:
                if threshold < self.maxFrc:
                    a.plot([self.minInv, self.maxInv],[threshold, threshold], color='black', linestyle='--')
                a.grid(True)
            else:
                raise Exception("Set a valid iteration to show its FSC")
        else:
            for ref3d in self._refsList:
                plot_title = 'class %s' % ref3d
                a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'FSC', yformat=False)
                legends = []
                
                for it in self._iterations:
                    parFn = self.protocol._getFileName('output_vol_par_class', iter=it, ref=ref3d)
                    if exists(parFn):
                        show = True
                        self._plotFSC(a, parFn)
                        legends.append('iter %d' % it)
                xplotter.showLegend(legends)
                if show:
                    if threshold < self.maxFrc:
                        a.plot([self.minInv, self.maxInv],[threshold, threshold], color='black', linestyle='--')
                    a.grid(True)
                else:
                    raise Exception("Set a valid iteration to show its FSC")
        
        return [xplotter]
    
    def _plotFSC(self, a, parFn):
        resolution_inv = self._getColunmFromFilePar(parFn, 2, invert=True)
        frc = self._getColunmFromFilePar(parFn, 5)
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])
        
# #===============================================================================
# # plotSSNR              
# #===============================================================================
    def _showSSNR(self, paramName=None):
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        xplotter = EmPlotter(x=gridsize[0], y=gridsize[1])
         
        for ref3d in self._refsList:
            plot_title = 'Resolution SSNR, for Class %s' % ref3d
            a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'sqrt(SSNR)', yformat=False)
            legendName = []
            for it in self._iterations:
                if self.protocol.IS_REFINE:
                    fn = self.protocol._getFileName('output_vol_par', iter=it)
                else:
                    fn = self.protocol._getFileName('output_vol_par_class', iter=it, ref=ref3d)
                if exists(fn):
                    self._plotSSNR(a, fn)
                legendName.append('iter %d' % it)
            xplotter.showLegend(legendName)
            a.grid(True)
         
        return [xplotter]
    
    def _plotSSNR(self, a, parFn):
        resolution_inv = self._getColunmFromFilePar(parFn, 2, invert=True)
        frc = self._getColunmFromFilePar(parFn, 8)
        
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)               
  
#===============================================================================
# Utils Functions
#===============================================================================
    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self._refsList = [1] 
        if not self.protocol.IS_REFINE:
            if self.showClasses3D == CLASSES_ALL:
                self._refsList = range(1, self.protocol.numberOfClasses.get()+1)
            else:
                self._refsList = self._getListFromRangeString(self.class3DSelection.get())
        
        self.protocol._createFilenameTemplates() # Load filename templates
        self.lastIter = self.protocol._getLastIter()
        self._iterations = []
        
        if self.iterToShow.get() == LAST_ITER:
            self._iterations = [self.lastIter]
        elif self.iterToShow.get() == ALL_ITERS:
            self._iterations = range(1, self.lastIter + 1)
        elif self.iterToShow.get() == SELECTED_ITERS:
            self._iterations = self._getListFromRangeString(self.selectedIters.get())
            
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
    
    def _validate(self):
        if self.lastIter is None:
            return ['There are not iterations completed.'] 
    
    def createDataView(self, filename, viewParams={}):
        return em.DataView(filename, env=self._env, viewParams=viewParams)
    
    def createScipionView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        return em.Classes3DView(self._project,
                                self.protocol.strId(), filename, other=inputParticlesId,
                                env=self._env, viewParams=viewParams)
        
    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        
        labels =  'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER:labels,
                      showj.VISIBLE: labels, showj.RENDER:'_filename',
                      'labels': 'id',
                      }
        return em.ObjectView(self._project,  self.protocol.strId(),
                             filename, other=inputParticlesId,
                             env=self._env, viewParams=viewParams)
    
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
        return "1/%0.2f" % inv
    
    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of references. """
        if n is None:
            n = len(self._refsList)
        
        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n+1)/2, 2]
            
        return gridsize
    
    def _getIterClasses(self, it, clean=False):
        """ Return the file with the classes for this iteration.
        If the file doesn't exists, it will be created. 
        """
        from convert import readSetOfClasses3D
        
        dataClasses = self.protocol._getFileName('classes_scipion', iter=it)
        
        if clean:
            cleanPath(dataClasses)
        
        if not exists(dataClasses):
            fileparList = []
            volumeList = []
            for ref in self._refsList:
                filepar = self.protocol._getFileName('output_par_class', iter=it, ref=ref)
                volFn = self.protocol._getFileName('iter_vol_class', iter=it, ref=ref)
                fileparList.append(filepar)
                volumeList.append(volFn)
            
            clsSet = em.SetOfClasses3D(filename=dataClasses)
            clsSet.setImages(self.inputParticles.get())
            readSetOfClasses3D(clsSet, fileparList, volumeList)
            clsSet.write()
            clsSet.close()

        return dataClasses
    
    def _getIterData(self, it, **kwargs):
        from convert import readSetOfParticles
        
        imgSqlite = self.protocol._getFileName('data_scipion', iter=it)
        
        if not exists(imgSqlite):
            imgPar = self.protocol._getFileName('output_par', iter=it)
            
            cleanPath(imgSqlite)
            imgSet = em.SetOfParticles(filename=imgSqlite)
            readSetOfParticles(self.protocol.inputParticles.get(), imgSet, imgPar)
            imgSet.write()
            
        return imgSqlite
    
    def _iterAngles(self, it, dataAngularDist):
        f = open(dataAngularDist)
        for line in f:
            if not line.startswith('C'):
                angles = map(float, line.split())
                rot = angles[1]
                tilt = angles[2]
                yield rot, tilt
        
        f.close()
    
    def _getColunmFromFilePar(self, parFn, col, invert=False):
        f1 = open(parFn)
        
        readLines = False
        value = []
        for l in f1:
            if "C  Average" in l:
                readLines = False
            if readLines:
                valList = l.split()
                val = float(valList[col])
                if invert:
                    val = 1/val
                value.append(val)
            if "NO.  RESOL" in l:
                readLines = True
        f1.close()
        return value

    def _getVolumeNames(self):
        volumes = []
        if self.protocol.IS_REFINE:
            for it in self._iterations:
                volFn = self.protocol._getFileName('iter_vol', iter=it)
                volumes.append(volFn)
        else:
            for it in self._iterations:
                for ref3d in self._refsList:
                    volFn = self.protocol._getFileName('iter_vol_class', iter=it, ref3d=ref3d)
                    volumes.append(volFn)
        return volumes


class ProtCTFFindViewer(Viewer):
    """ Visualization of Frealign."""
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer CtfFind'
    _targets = [ProtCTFFind]
     
     
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
            from pyworkflow.utils.path import removeBaseExt
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
                samplingRate = mic.getSamplingRate() * self.protocol.ctfDownFactor.get()
                mic.setSamplingRate(samplingRate)
                out = self.protocol._getCtfOutPath(micDir)
                psdFile = self.protocol._getPsdPath(micDir)
                
                if exists(out) and exists(psdFile):
                    ctfModel = em.CTFModel()
                
                    if not self.protocol.useCftfind4:
                        readCtfModel(ctfModel, out)
                    else:
                        readCtfModel(ctfModel, out, True)
                    
                    ctfModel.setPsdFile(psdFile)
                    ctfModel.setMicrograph(mic)
                    
                    ctfSet.append(ctfModel)
            
            if ctfSet.getSize() < 1:
                raise Exception("Has not been completed the CTF estimation of any micrograph")
            else:
                ctfSet.write()
                ctfSet.close()
                self._visualize(ctfSet)
         
        if issubclass(cls, ProtCTFFind) and not obj.hasAttribute("outputCTF"):
            mics = obj.inputMicrographs.get()
            visualizeObjs(obj, mics)
 
         
        elif obj.hasAttribute("outputCTF"):
            self._visualize(obj.outputCTF)
             
        else:
            fn = obj.getFileName()
            psdLabels = '_psdFile'
            labels = 'id enabled comment %s _defocusU _defocusV _defocusAngle _defocusRatio' % psdLabels
            if self.protocol.useCftfind4:
                labels = labels + ' _ctffind4_ctfResolution'
                self._views.append(em.ObjectView(self._project, obj.strId(), fn,
                                                 viewParams={showj.MODE: showj.MODE_MD,
                                                             showj.ORDER: labels,
                                                             showj.VISIBLE: labels,
                                                             showj.ZOOM: 50,
                                                             showj.RENDER: psdLabels,
                                                             showj.OBJCMDS: "'%s'" % showj.OBJCMD_CTFFIND4}))
            else:
                self._views.append(em.ObjectView(self._project, obj.strId(), fn,
                                                 viewParams={showj.MODE: showj.MODE_MD,
                                                             showj.ORDER: labels,
                                                             showj.VISIBLE: labels,
                                                             showj.ZOOM: 50,
                                                             showj.RENDER: psdLabels}))
        
        return self._views
    
def createCtfPlot(ctfSet, ctfId):
    from pyworkflow.utils.path import removeExt
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = removeExt(psdFn) + "_avrot.txt"

    gridsize = [1, 1]
    xplotter = EmPlotter(x=gridsize[0], y=gridsize[1], windowTitle='CTF Fitting')
    plot_title = "CTF Fitting"
    a = xplotter.createSubPlot(plot_title, 'pixels^-1', 'CTF', yformat=False)
    legendName = []
    for i in range(1, 5):
        _plotCurve(a, i, fn)
        if i == 1:
            legendName.append('rotational avg. No Astg')
        elif i == 2:
            legendName.append('rotational avg.')
        elif i == 3:
            legendName.append('CTF Fit')
        elif i == 4:
            legendName.append('Cross Correlation')
        elif i == 5:
            legendName.append('2sigma cross correlation of noise')
            
    xplotter.showLegend(legendName)
    a.grid(True)
    xplotter.show()

def _plotCurve(a, i, fn):
    freqs = _getValues(fn, 0)
    curv = _getValues(fn, i)
    a.plot(freqs, curv)

def _getValues(fn, row):
    f = open(fn)
    values = []
    i = 0
    for line in f:
        if not line.startswith("#"):
            if i == row:
                values = line.split()
                break
            i += 1
    f.close()
    return values
