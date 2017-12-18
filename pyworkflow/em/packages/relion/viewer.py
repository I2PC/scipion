# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import exists

import pyworkflow.em as em
import pyworkflow.em.showj as showj
import pyworkflow.em.metadata as md
from pyworkflow.em.data import SetOfParticles, SetOfImages
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.protocol.params as params
from pyworkflow.viewer import (Viewer, ProtocolViewer, DESKTOP_TKINTER,
                               WEB_DJANGO)

from protocol_classify2d import ProtRelionClassify2D
from protocol_classify3d import ProtRelionClassify3D
from protocol_refine3d import ProtRelionRefine3D
from protocol_polish import ProtRelionPolish
from protocol_postprocess import ProtRelionPostprocess
from protocol_autopick import ProtRelionAutopick, ProtRelionAutopickFom
from protocol_sort import ProtRelionSortParticles
from protocol_initialmodel import ProtRelionInitialModel
from protocol_localres import ProtRelionLocalRes
from convert import relionToLocation

ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CHIMERADATAVIEW = 0

CLASSES_ALL = 0
CLASSES_SEL = 1

FSC_CORRECTED = 0
FSC_UNMASKEDMAPS = 1
FSC_MASKEDMAPS = 2
FSC_RANDOMIZED = 3
FSC_ALL = 4



class RelionPlotter(EmPlotter):
    """ Class to create several plots with Xmipp utilities"""
    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        EmPlotter.__init__(self, x, y, mainTitle, **kwargs)
    
    def plotMdAngularDistribution(self, title, angularMd, color='blue'):
        """Create an special type of subplot, representing the angular
        distribution of weight projections. A metadata should be provided containing
        labels: RLN_ORIENT_ROT, RLN_ORIENT_TILT, MDL_WEIGHT """
        from math import radians
        
        rot = [radians(angularMd.getValue(md.RLN_ORIENT_ROT, objId)) for objId in angularMd]
        tilt = [angularMd.getValue(md.RLN_ORIENT_TILT, objId) for objId in angularMd]
        weight = [angularMd.getValue(md.MDL_WEIGHT, objId) for objId in angularMd]
        
        self.plotAngularDistribution(title, rot, tilt, weight)

    def plotMd(self, mdObj, mdLabelX, mdLabelY, color='g',**args):
        """ plot metadata columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        if mdLabelX:
            xx = []
        else:
            xx = range(1, mdObj.size() + 1)
        yy = []
        for objId in mdObj:
            if mdLabelX:
                xx.append(mdObj.getValue(mdLabelX, objId))
            yy.append(mdObj.getValue(mdLabelY, objId))
        
        nbins = args.pop('nbins', None)
        if nbins is None:
            self.plotData(xx, yy, color, **args)
        else:
            self.plotHist(yy, nbins, color, **args)
        
    def plotMdFile(self, mdFilename, mdLabelX, mdLabelY, color='g', **args):
        """ plot metadataFile columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        mdObj = md.MetaData(mdFilename)
        self.plotMd(mdObj, mdLabelX, mdLabelY, color='g',**args)
        

class RelionViewer(ProtocolViewer):
    """ This protocol serve to analyze the results of Relion runs.
    (for protocols classify 2d/3d, 3d auto-refine and initial model)
    The visualization tools follow the recommendations of Relion 1.3 tutorial:
    http://www2.mrc-lmb.cam.ac.uk/groups/scheres/relion13_tutorial.pdf
    """
    _targets = [ProtRelionClassify2D, ProtRelionClassify3D,
                ProtRelionRefine3D, ProtRelionInitialModel]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer relion'
    
    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize", 
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', params.NumericRangeParam, 
                      condition='viewIter==%d' % ITER_SELECTION, 
                      label="Iterations list", 
                      help="Write the iteration list to visualize.")

        changesLabel = 'Changes in Offset and Angles'
        
        group = form.addGroup('Particles')
        if self.protocol.IS_CLASSIFY:

            group.addParam('showImagesInClasses', params.LabelParam, 
                           label='Show classification in Scipion', important=True,
                           help='Display each class with the number of particles assigned. \n'
                                '*Note1*: The images of one class can be shown by \n'
                                'right-click on the class and select "Open images".\n'
                                '*Note2*: This option convert the Relion star file to\n'
                                'Scipion format and can take several minutes if the \n'
                                'number of particles is high.')
            group.addParam('showClassesOnly', params.LabelParam, 
                           label='Show classes only (*_model.star)',
                           help='Display the classes directly form the *_model.star file.')
            changesLabel = 'Changes in Offset, Angles and Classes'
        else:
            group.addParam('showImagesAngularAssignment', params.LabelParam, 
                           label='Particles angular assignment')
        group.addParam('showOptimiserFile', params.LabelParam, 
                       label='Show *_optimiser.star file')
        
        if self.protocol.IS_3D:
            group = form.addGroup('Volumes')
            
            if self.protocol.IS_CLASSIFY:
                group.addParam('showClasses3D', params.EnumParam, default=CLASSES_ALL,
                               choices=['all', 'selection'], 
                               display=params.EnumParam.DISPLAY_HLIST,
                               label='3D Class to visualize',
                               help='')
                group.addParam('class3DSelection', params.NumericRangeParam, default='1',
                               condition='showClasses3D == %d' % CLASSES_SEL,
                               label='Classes list',
                               help='')
            else:
                if self.protocol.IS_3D_INIT:
                    group.addHidden('showHalves', params.IntParam, default=3)
                else:
                    group.addParam('showHalves', params.EnumParam, default=0,
                                   choices=['half1', 'half2', 'both', 'final'],
                                   label='Volume to visualize',
                                   help='Select which half do you want to visualize.')
            
            group.addParam('displayVol', params.EnumParam, choices=['slices', 'chimera'], 
                           default=VOLUME_SLICES, display=params.EnumParam.DISPLAY_HLIST,
                           label='Display volume with',
                           help='*slices*: display volumes as 2D slices along z axis.\n'
                                '*chimera*: display volumes as surface with Chimera.')
            group.addParam('displayAngDist', params.EnumParam, choices=['2D plot', 'chimera'], 
                           default=ANGDIST_2DPLOT, display=params.EnumParam.DISPLAY_HLIST,
                           label='Display angular distribution',
                           help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                                '*chimera*: display angular distribution using Chimera with red spheres.')
            group.addParam('spheresScale', params.IntParam, default=100, 
                           expertLevel=LEVEL_ADVANCED,
                           label='Spheres size',
                           help='')

            group = form.addGroup('Resolution')
            group.addParam('figure', params.EnumParam, default=0,
                           choices=['new', 'active'],
                           label='Figure', display=params.EnumParam.DISPLAY_HLIST)
            group.addParam('resolutionPlotsSSNR', params.LabelParam, default=True,
                           label='Display SSNR plots',
                           help='Display signal to noise ratio plots (SSNR)')
            if not self.protocol.IS_CLASSIFY and not self.protocol.IS_3D_INIT:
                group.addParam('resolutionPlotsFSC', params.LabelParam, default=True,
                               label='Display resolution plots (FSC)',
                               help='')
                group.addParam('resolutionThresholdFSC', params.FloatParam, default=0.143, 
                               expertLevel=LEVEL_ADVANCED,
                               label='Threshold in resolution plots',
                               help='')
            
        form.addSection('Overall')      
        form.addParam('showPMax', params.LabelParam, default=True,
                      label="Show average PMax",
                      help='Average (per class) of the maximum value\n '
                           'of normalized probability function')
        form.addParam('showChanges', params.LabelParam, default=True,
                      label=changesLabel,
                      help='Visualize changes in orientation, offset and\n '
                           'number images assigned to each class')
        form.addParam('plotClassDistribution', params.LabelParam, default=True,
                      label='Prot class distribution over iterations',
                      help='Plot each class distribution over iterations as '
                           'bar plots.')
                                              
    def _getVisualizeDict(self):
        self._load()
        visualizeDict = {
                'showImagesInClasses': self._showImagesInClasses,
                'showClassesOnly': self._showClassesOnly,
                'showImagesAngularAssignment': self._showImagesAngularAssignment,
                'showOptimiserFile': self._showOptimiserFile,
                'showLL': self._showLL,
                'showPMax': self._showPMax,
                'showChanges': self._showChanges,
                'displayVol': self._showVolumes,
                'displayAngDist': self._showAngularDistribution,
                'resolutionPlotsSSNR': self._showSSNR,
                'resolutionPlotsFSC': self._showFSC,
                'plotClassDistribution': self._plotClassDistribution,
                }

        # If the is some error during the load, just show that instead
        # of any viewer
        if self._errors:
            for k in visualizeDict.keys():
                visualizeDict[k] = self._showErrors

        return visualizeDict

    def _showErrors(self, param=None):
        views = []
        self.errorList(self._errors, views)
        return views
        
    def _viewAll(self, *args):
        pass
    
#===============================================================================
# showImagesInClasses     
#===============================================================================
    def _getZoom(self):
        # Ensure that classes are shown at least at 128 px to 
        # properly see the rlnClassDistribution label. 
        dim = self.protocol.inputParticles.get().getDim()[0]
        if dim < 128:
            zoom = 128*100/dim
        else:
            zoom = 100
        return zoom
        
    def _showImagesInClasses(self, paramName=None):
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        views = []
        if (self.viewIter == ITER_LAST and
            getattr(self.protocol, 'outputClasses', None) is not None):
            fn = self.protocol.outputClasses.getFileName()
            v = self.createScipionView(fn)
            views.append(v)
        else:
            for it in self._iterations:
                fn = self.protocol._getIterClasses(it)
                v = self.createScipionView(fn)
                views.append(v)
        
        return views
    
    def _showClassesOnly(self, paramName=None):
        views = []
        viewParams = {showj.MODE: showj.MODE_GALLERY,
                      showj.RENDER: 'rlnReferenceImage',
                      showj.SORT_BY: 'rlnClassDistribution desc',
                      showj.LABELS: 'rlnClassDistribution',
                      showj.ZOOM: str(self._getZoom())
                      }
        
        for it in self._iterations:
            modelFile = self.protocol._getFileName('model', iter=it)
            v = self.createDataView('model_classes@' + modelFile, viewParams=viewParams)
            views.append(v)
        return views        

#===============================================================================
# showImagesAngularAssignment     
#===============================================================================
    def _showImagesAngularAssignment(self, paramName=None):
        views = []
        
        for it in self._iterations:
            fn = self.protocol._getIterData(it, alignType=em.ALIGN_PROJ)
            v = self.createScipionPartView(fn)
            views.append(v)
        
        return views
    
    def _showOptimiserFile(self,  paramName=None):
        views = []
        
        for it in self._iterations:
            optimiserFile = self.protocol._getFileName('optimiser', iter=it)
            v = self.createDataView(optimiserFile)
            views.append(v)
        return views
#=====================================================================
# showLLRelion
#=====================================================================
    def _showLL(self, paramName=None):
        views = []
        for it in self._iterations:
            fn = self.protocol._getIterData(it)
            views.append(self.createScipionView(fn))
            
        return views

#===============================================================================
# ShowPMax
#===============================================================================
    def _showPMax(self, paramName=None):
        labels = [md.RLN_MLMODEL_AVE_PMAX, md.RLN_PARTICLE_PMAX]
        
        mdIters = md.MetaData()
        iterations = range(self.firstIter, self.lastIter+1)
        
        for it in iterations:  # range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
            objId = mdIters.addObject()
            mdIters.setValue(md.MDL_ITER, it, objId)
            for i, prefix in enumerate(self.protocol.PREFIXES):
                fn = 'model_general@' + self.protocol._getFileName(prefix + 'model', iter=it)
                mdModel = md.RowMetaData(fn)
                pmax = mdModel.getValue(md.RLN_MLMODEL_AVE_PMAX)
                mdIters.setValue(labels[i], pmax, objId)
        fn = self.protocol._getFileName('all_avgPmax_xmipp')
        mdIters.write(fn)
            
        colors = ['g', 'b']

        xplotter = RelionPlotter()
        xplotter.createSubPlot("Avg PMax per Iterations", "Iterations", "Avg PMax")
        
        for label, color in zip(labels, colors):
            xplotter.plotMd(mdIters, md.MDL_ITER, label, color)
        
        if len(self.protocol.PREFIXES) > 1:
            xplotter.showLegend(self.protocol.PREFIXES)

        return [self.createDataView(fn), xplotter]

# ==============================================================================
# Get classes info per iteration
# ==============================================================================
    def _plotClassDistribution(self, paramName=None):
        labels = ["rlnClassDistribution", "rlnAccuracyRotations",
                  "rlnAccuracyTranslations"]
        iterations = range(self.firstIter, self.lastIter + 1)

        classInfo = {}

        for it in iterations:
            modelStar = self.protocol._getFileName('model', iter=it)

            for row in md.iterRows('%s@%s' % ('model_classes', modelStar)):
                i, fn = relionToLocation(row.getValue('rlnReferenceImage'))
                if i == em.NO_INDEX: # the case for 3D classes
                    # NOTE: Since there is not an proper ID value in
                    #  the clases metadata, we are assuming that class X
                    # has a filename *_classXXX.mrc (as it is in Relion)
                    # and we take the ID from there
                    index = int(fn[-7:-4])
                else:
                    index = i

                if index not in classInfo:
                    classInfo[index] = {}
                    for l in labels:
                        classInfo[index][l] = []

                for l in labels:
                    classInfo[index][l].append(row.getValue(l))

        xplotter = RelionPlotter()
        xplotter.createSubPlot("Classes distribution over iterations",
                               "Iterations", "Classes Distribution")

        # Empty list for each iteration
        iters = [[]] * len(iterations)

        l = labels[0]
        for index in sorted(classInfo.keys()):
            for it, value in enumerate(classInfo[index][l]):
                iters[it].append(value)

        ax = xplotter.getLastSubPlot()

        n = len(iterations)
        ind = range(n)
        bottomValues = [0] * n
        width = 0.45  # the width of the bars: can also be len(x) sequence

        def get_cmap(N):
            import matplotlib.cm as cmx
            import matplotlib.colors as colors
            """Returns a function that maps each index in 0, 1, ... N-1 to a distinct
            RGB color."""
            color_norm = colors.Normalize(vmin=0, vmax=N)#-1)
            scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')

            def map_index_to_rgb_color(index):
                return scalar_map.to_rgba(index)

            return map_index_to_rgb_color

        cmap = get_cmap(len(classInfo))

        for classId in sorted(classInfo.keys()):
            values = classInfo[classId][l]
            ax.bar(ind, values, width, bottom=bottomValues, color=cmap(classId))
            bottomValues = [a+b for a, b in zip(bottomValues, values)]

        ax.get_xaxis().set_ticks([i + 0.25 for i in ind])
        ax.get_xaxis().set_ticklabels([str(i) for i in ind])

        return [xplotter]
    
#===============================================================================
# ShowChanges    
#===============================================================================    
    def _showChanges(self, paramName=None):
        
        mdIters = md.MetaData()
        iterations = range(self.firstIter, self.lastIter+1)
        
        print " Computing average changes in offset, angles, and class membership"
        for it in iterations:
            print "Computing data for iteration; %03d" % it
            objId = mdIters.addObject()
            mdIters.setValue(md.MDL_ITER, it, objId)
            #agregar por ref3D
            fn = self.protocol._getFileName('optimiser', iter=it )
            mdOptimiser = md.RowMetaData(fn)
            for label in self.protocol.CHANGE_LABELS:
                mdIters.setValue(label, mdOptimiser.getValue(label), objId)
        fn = self.protocol._getFileName('all_changes_xmipp')
        mdIters.write(fn)
        
        return [self.createDataView(fn)]
    
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
        path = self.protocol._getExtraPath('relion_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if exists(volFn.replace(':mrc', '')):
                files.append(volFn)
        self.createVolumesSqlite(files, path, samplingRate)
        return [em.ObjectView(self._project, self.protocol.strId(), path)]
    
    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        
        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for volFn in volumes:
                # We assume that the chimera script will be generated
                # at the same folder than relion volumes
                vol = volFn.replace(':mrc', '')
                localVol = os.path.basename(vol)
                if exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = em.ChimeraView(cmdFile)
        else:
            #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volumes[0])
            view = em.ChimeraClientView(volumes[0])
            
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
                if isinstance(plot, RelionPlotter):
                    views.append(plot)
        return views
    
    def _createAngDistChimera(self, it):
        # Common variables to use
        nparts = self.protocol.inputParticles.get().getSize()
        radius = self.spheresScale.get()
        if radius < 0:
            radius = self.protocol.maskDiameterA.get()/2
            
        prefixes = self._getPrefixes()
        if len(self._refsList) == 1:
            # If just one reference we can show the angular distribution
            ref3d = self._refsList[0]
            volFn = self._getVolumeNames()[0]
            if exists(volFn.replace(":mrc", "")):
                for prefix in prefixes:
                    sqliteFn = self.protocol._getFileName('projections', iter=it, ref3d=ref3d, half=prefix)
                    if not exists(sqliteFn):
                        self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(self._getMdOut(it, prefix, ref3d)))
                    return em.ChimeraClientView(volFn, angularDistFile=sqliteFn, spheresDistance=radius)
            else:
                raise Exception("This class is Empty. Please try with other class")
        else:
            return self.infoMessage("Please select only one class to display angular distribution",
                                    "Input selection") 
    
    def _createAngDist2D(self, it):
        # Common variables to use
        nparts = self.protocol.inputParticles.get().getSize()
        prefixes = self._getPrefixes()
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        if prefixes[0] == "final":
            title = "Final"
        else:
            title = 'Iteration %d' % it
        
        plotter = RelionPlotter(x=gridsize[0], y=gridsize[1], 
                                 mainTitle=title, windowTitle="Angular Distribution")
        for prefix in prefixes:
            for ref3d in self._refsList:
                dataStar = self._getDataStar(prefix, it)
                randomSet = self._getRandomSet(prefix)
                if randomSet > 0:
                    title = '%s class %d' % (prefix, ref3d)
                else:
                    title = 'class %d' % ref3d
                sqliteFn = self.protocol._getFileName('projections', iter=it, ref3d=ref3d, half=prefix)
                if not exists(sqliteFn):
                    self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(self._getMdOut(it, prefix, ref3d)))
                plotter.plotAngularDistributionFromMd(sqliteFn, title)
        
        for prefix in prefixes:
            dataStar = self._getDataStar(prefix, it)
            if exists(dataStar):
                return plotter
            else:
                return
    
#===============================================================================
# plotSSNR              
#===============================================================================

    def _getFigure(self):
        return None if self.figure == 0 else 'active'

    def _showSSNR(self, paramName=None):
        prefixes = self._getPrefixes()        
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        gridsize = self._getGridSize(n)
        md.activateMathExtensions()
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1], figure=self._getFigure())
        
        for prefix in prefixes:
            for ref3d in self._refsList:
                plot_title = 'Resolution SSNR %s, for Class %s' % (prefix, ref3d)
                a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'log(SSNR)', yformat=False)
                blockName = 'model_class_%d@' % ref3d
                for it in self._iterations:
                    fn = self._getModelStar(prefix, it)
                    if exists(fn):
                        self._plotSSNR(a, blockName+fn, 'iter %d' % it)
                xplotter.legend()
                a.grid(True)
        
        return [xplotter]
    
    def _plotSSNR(self, a, fn, label):
        mdOut = md.MetaData(fn)
        mdSSNR = md.MetaData()
        # only cross by 1 is important
        mdSSNR.importObjects(mdOut, md.MDValueGT(md.RLN_MLMODEL_DATA_VS_PRIOR_REF, 0.9))
        mdSSNR.operate("rlnSsnrMap=log(rlnSsnrMap)")
        resolution_inv = [mdSSNR.getValue(md.RLN_RESOLUTION, id) for id in mdSSNR]
        frc = [mdSSNR.getValue(md.RLN_MLMODEL_DATA_VS_PRIOR_REF, id) for id in mdSSNR]
        a.plot(resolution_inv, frc, label=label)
        a.xaxis.set_major_formatter(self._plotFormatter)               
 
#===============================================================================
# plotFSC            
#===============================================================================
    def _showFSC(self, paramName=None):
        #self._iterations = self._getListFromRangeString(self.iterSelection.get())
        print("_showFSC_self._iterations",self._iterations)
        threshold = self.resolutionThresholdFSC.get()
        prefixes = self._getPrefixes()        
        nrefs = len(self._refsList)
        n = nrefs * len(prefixes)
        #gridsize = self._getGridSize(n)
        
        md.activateMathExtensions()
        
        fscViewer = em.FscViewer(project=self.protocol.getProject(),
                                 threshold=threshold,
                                 protocol=self.protocol,
                                 figure=self._getFigure(),
                                 addButton=True)
        fscSet = self.protocol._createSetOfFSCs()
        for prefix in prefixes:
            for ref3d in self._refsList:#ROB: I believe len(_refsList)==1
                #plot_title = prefix + 'class %s' % ref3d
                blockName = 'model_class_%d@' % ref3d
                for it in self._iterations:
                    model_star = self._getModelStar(prefix, it)

                    if exists(model_star):
                        #fnFSC = blockName + model_star
                        fsc = self._plotFSC(None, blockName + model_star, 'iter %d' % it)
                        fscSet.append(fsc)
        fscViewer.visualize(fscSet)
        return [fscViewer]
    
    def _plotFSC(self, a, model_star, label):
        mdStar = md.MetaData(model_star)
        resolution_inv = [mdStar.getValue(md.RLN_RESOLUTION, id) for id in mdStar]
        frc = [mdStar.getValue(md.RLN_MLMODEL_FSC_HALVES_REF, id) for id in mdStar]

        fsc = em.data.FSC(objLabel=label)
        fsc.setData(resolution_inv, frc)

        return fsc


    
#===============================================================================
# Utils Functions
#===============================================================================
    def _validate(self):
        if self.lastIter is None:
            return ['There are not iterations completed.'] 
    
    def createDataView(self, filename, viewParams={}):
        return em.DataView(filename, env=self._env, viewParams=viewParams)

    def createScipionView(self, filename):
        labels =  'enabled id _size _representative._filename '
        labels += '_rlnclassDistribution _rlnAccuracyRotations _rlnAccuracyTranslations'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, 
                      showj.RENDER:'_representative._filename',
                      showj.SORT_BY: '_size desc',
                      showj.ZOOM: str(self._getZoom())
                      }
        inputParticlesId = self.protocol.inputParticles.get().strId()
        ViewClass = em.ClassesView if self.protocol.IS_2D else em.Classes3DView
        view = ViewClass(self._project,
                          self.protocol.strId(), filename, other=inputParticlesId,
                          env=self._env,
                          viewParams=viewParams)

        return view

    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol._getInputParticles().strId()
        
        labels = 'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER:labels,
                      showj.VISIBLE: labels, showj.RENDER:'_filename',
                      'labels': 'id',
                      }
        return em.ObjectView(self._project, 
                             self.protocol.strId(), filename, other=inputParticlesId,
                             env=self._env,
                             viewParams=viewParams)

    def _getRange(self, var, label):
        """ Check if the range is not empty.
        :param var: The variable to retrieve the value
        :param label: the labe used for the message string
        :return: the list with the range of values, empty
        """
        value = var.get()
        if value is None or not value.strip():
            self._errors.append('Provide %s selection.' % label)
            result = []
        else:
            result = self._getListFromRangeString(value)

        return result

    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self._refsList = [1]
        self._errors = []

        if self.protocol.IS_3D and self.protocol.IS_CLASSIFY:
            if self.showClasses3D == CLASSES_ALL:
                self._refsList = range(1, self.protocol.numberOfClasses.get()+1)
            else:
                self._refsList = self._getRange(self.class3DSelection, 'classes 3d')
        self.protocol._initialize() # Load filename templates
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()
        
        halves = getattr(self, 'showHalves', None)
        if self.viewIter.get() == ITER_LAST or halves == 3:
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getRange(self.iterSelection, 'iterations')
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        
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
    
    def _getPrefixes(self):
        prefixes = self.protocol.PREFIXES
        halves = getattr(self, 'showHalves', None)
        if halves:
            if halves == 0:
                prefixes = ['half1_']
            elif halves == 1:
                prefixes = ['half2_']
            elif halves == 3:
                prefixes = ['final']
        return prefixes
    
    def _iterAngles(self, mdOut):
        
        for objId in mdOut:
            rot = mdOut.getValue(md.RLN_ORIENT_ROT, objId)
            tilt = mdOut.getValue(md.RLN_ORIENT_TILT, objId)
            yield rot, tilt
    
    def _getVolumePrefixes(self):
        prefixes = self._getPrefixes()
        if prefixes[0] == 'final':
            prefixes += ['final_half1_', 'final_half2_']
        if prefixes[0] == 'final' and self.protocol.IS_3D_INIT:
            prefixes = ['finalSGD']
        return prefixes

    def _getVolumeNames(self):
        vols = []
        prefixes = self._getVolumePrefixes()
        for it in self._iterations:
            for ref3d in self._refsList:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume', iter=it, ref3d=ref3d)
                    vols.append(volFn)
        return vols
    
    def _getMdOut(self, it, prefix, ref3d):
        randomSet = self._getRandomSet(prefix)
        dataStar = self._getDataStar(prefix, it)
        
        if randomSet > 0 and  randomSet < 3:
            mdAll = md.MetaData(dataStar)
            mdTmp = md.MetaData()
            mdTmp.importObjects(mdAll, md.MDValueEQ(md.RLN_PARTICLE_RANDOM_SUBSET, randomSet))
        else:
            mdTmp = md.MetaData(dataStar)
        
        mdOut = md.MetaData()
        mdOut.importObjects(mdTmp, md.MDValueEQ(md.RLN_PARTICLE_CLASS, ref3d))
        return mdOut
    
    def _getDataStar(self, prefix, it):
        randomSet = self._getRandomSet(prefix)
        if randomSet > 0 or self.protocol.IS_3D_INIT:
            return self.protocol._getFileName('data', iter=it)
        else:
            return self.protocol._getFileName('dataFinal')
    
    def _getModelStar(self, prefix, it):
        randomSet = self._getRandomSet(prefix)
        if self.protocol.IS_3D_INIT:
            return self.protocol._getFileName('model', iter=it)
        if randomSet > 0:
            return self.protocol._getFileName(prefix + 'model', iter=it)
        else:
            return self.protocol._getFileName('modelFinal')
    
    def _getRandomSet(self, prefix):
        if prefix == "final":
            return 0
        elif prefix == "half1_":
            return 1
        elif prefix == "half2_":
            return 2
        else:
            return 3


class PostprocessViewer(ProtocolViewer):
    """ Class to visualize Relion postprocess protocol """
    _targets = [ProtRelionPostprocess]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer postprocess relion'
    
    def setProtocol(self, protocol):
        ProtocolViewer.setProtocol(self, protocol)
        self.__defineParams(self._form)
        self._createVarsFromDefinition()
        from pyworkflow.em.packages.xmipp3 import getEnviron
        self._env = dict(getEnviron())
#        self._load()
        
    def _defineParams(self, form):
        self._form = form
        
    def __defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('3D analysis')
        
        group.addParam('displayVol', params.EnumParam, choices=['slices', 'chimera'], 
                      display=params.EnumParam.DISPLAY_HLIST, default=VOLUME_SLICES,
                      label='Display volume with',
                      help='*slices*: display volumes as 2D slices along z axis.\n'
                           '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayMaskedVol', params.EnumParam, choices=['slices', 'chimera'], 
                      display=params.EnumParam.DISPLAY_HLIST, default=VOLUME_SLICES,
                      label='Display masked volume with',
                      help='*slices*: display masked volume as 2D slices along z axis.\n'
                           '*chimera*: display masked volume as surface with Chimera.')
        group.addParam('figure', params.EnumParam, default=0,
                           choices=['new', 'active'],
                           label='Figure',
                       display=params.EnumParam.DISPLAY_HLIST,
                       help="plot in a new window vs the last open one")
        group.addParam('resolutionPlotsFSC', params.EnumParam,
                      choices=['Corrected', 'Unmasked Maps', 'Masked Maps', 'Phase Randomized Masked Maps', 'all'],
                      default=FSC_CORRECTED, display=params.EnumParam.DISPLAY_COMBO, 
                      label='Display resolution plots (FSC)',
                      help='') 
#         group.addParam('resolutionPlotsFSC', params.LabelParam, default=True,
#                       label='Display resolution plots (FSC) ?',
#                       help='')
        group.addParam('resolutionThresholdFSC', params.FloatParam, default=0.143, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold in resolution plots',
                      help='')
        group.addParam('guinierPlots', params.LabelParam,
              default=True, label='Display guinier plots',
              help='')

    def _getVisualizeDict(self):
        self._load()
        return {'displayVol': self._showVolume,
                'displayMaskedVol': self._showMaskedVolume,
                'guinierPlots': self._showGuinier,
                'resolutionPlotsFSC': self._showFSC
                }
#===============================================================================
# ShowVolumes
#===============================================================================
        
    def _showVolumeShowj(self, volPath):        
        return [em.DataView(volPath)]
    
    def _showVolumesChimera(self, volPath):
        """ Create a chimera script to visualize selected volumes. """
        #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volPath)
        view = em.ChimeraClientView(volPath)
        return [view]
            
    def _showVolume(self, paramName=None):
        volPath = self.protocol._getExtraPath('postprocess.mrc:mrc')
        
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera(volPath)
        
        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumeShowj(volPath)
                
    def _showMaskedVolume(self, paramName=None):
        volPath = self.protocol._getExtraPath('postprocess_masked.mrc:mrc')
        
        if self.displayMaskedVol == VOLUME_CHIMERA:
            return self._showVolumesChimera(volPath)
        
        elif self.displayMaskedVol == VOLUME_SLICES:
            return self._showVolumeShowj(volPath)
    
#===============================================================================
# plotFSC            
#===============================================================================
    def _getFigure(self):
        return None if self.figure == 0 else 'active'

    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()

        n = 1
        gridsize = [1, 1]

        md.activateMathExtensions()
        fscViewer = em.FscViewer(project=self.protocol.getProject(),
                                 threshold=threshold,
                                 protocol=self.protocol,
                                 figure=self._getFigure(),
                                 addButton=True)
        fscSet = self.protocol._createSetOfFSCs()

        modelStar = self.protocol._getExtraPath('postprocess.star')
        for label in self._getFSCLabels():
            if exists(modelStar):
                model = 'fsc@' + modelStar
                legend = self._getLegend(label)
                fsc = self._plotFSC(None, model, label, legend)
                fscSet.append(fsc)
        fscViewer.visualize(fscSet)
        return [fscViewer]

    #ROB this function is duplicated
    def _plotFSC(self, a, model_star, label, legend):
        mdStar = md.MetaData(model_star)
        resolution_inv = [mdStar.getValue(md.RLN_RESOLUTION, id) for id in mdStar]
        frc = [mdStar.getValue(label, id) for id in mdStar]

        fsc = em.data.FSC(objLabel=legend)
        fsc.setData(resolution_inv, frc)

        return fsc

#===============================================================================
# plotGuinier
#===============================================================================
    def _showGuinier(self, paramName=None):
        gridsize = [1, 1]
        md.activateMathExtensions()
        
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1], windowTitle='Guinier Plot')
        a = xplotter.createSubPlot("", 'Angstroms^-2', 'log(Amplitude)', yformat=False)
        legends = []
        modelStar = self.protocol._getExtraPath('postprocess.star')
        for label in self._getGuinerLabels():
            if exists(modelStar):
                model = 'guinier@' + modelStar
                self._plotGuinier(a, model, label)
                legends.append(self._getGuinerLegend(label))
            
        xplotter.showLegend(legends)
        a.grid(True)
        
        return [xplotter]
    
    def _plotGuinier(self, a, model, label):
        mdStar = md.MetaData(model)
        resolSqInv = [mdStar.getValue(md.RLN_POSTPROCESS_GUINIER_RESOL_SQUARED, id) for id in mdStar]
        logAmp = [mdStar.getValue(label, id) for id in mdStar]
        self.maxfsc = max(logAmp)
        self.minInv = min(resolSqInv)
        self.maxInv = max(resolSqInv)
        a.plot(resolSqInv, logAmp)
        a.xaxis.set_major_formatter(self._plotFormatter)
    
#===============================================================================
# Utils Functions
#===============================================================================
    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
        return "1/%0.2f" % inv
    
    def _getFSCLabels(self):
        if self.resolutionPlotsFSC.get() == 0:
            return [md.RLN_POSTPROCESS_FSC_TRUE]
        elif self.resolutionPlotsFSC.get() == 1:
            return [md.RLN_POSTPROCESS_FSC_UNMASKED]
        elif self.resolutionPlotsFSC.get() == 2:
            return [md.RLN_POSTPROCESS_FSC_MASKED]
        elif self.resolutionPlotsFSC.get() == 3:
            return [md.RLN_POSTPROCESS_FSC_RANDOM_MASKED]
        else:
            return [md.RLN_POSTPROCESS_FSC_TRUE, md.RLN_POSTPROCESS_FSC_UNMASKED,
                    md.RLN_POSTPROCESS_FSC_MASKED, md.RLN_POSTPROCESS_FSC_RANDOM_MASKED]
    
    def _getLegend(self, label):
        if label == md.RLN_POSTPROCESS_FSC_TRUE:
            return 'Corrected'
        elif label == md.RLN_POSTPROCESS_FSC_UNMASKED:
            return 'Unmasked Maps'
        elif label == md.RLN_POSTPROCESS_FSC_MASKED:
            return 'Masked Maps'
        else:
            return 'Phase Randomized Masked Maps'
    
    def _getGuinerLabels(self):
        return [md.RLN_POSTPROCESS_GUINIER_VALUE_IN,
                md.RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED,
                md.RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED,
                md.RLN_POSTPROCESS_GUINIER_VALUE_INVMTF,
                md.RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT]
    
    def _getGuinerLegend(self, label):
        if label == md.RLN_POSTPROCESS_GUINIER_VALUE_IN:
            return 'log(Amplitudes) Original'
        elif label == md.RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED:
            return 'log(Amplitudes) Weighted'
        elif label == md.RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED:
            return 'log(Amplitudes) Sharpened'
        elif label == md.RLN_POSTPROCESS_GUINIER_VALUE_INVMTF:
            return 'log(Amplitudes) MTF Corrected'
        else:
            return 'log(Amplitudes) Intercept'
        



class RelionAutopickViewerFOM(Viewer):
    """ Class to visualize Relion postprocess protocol """
    _targets = [ProtRelionAutopickFom]
    _environments = [DESKTOP_TKINTER]
    
    def _visualize(self, obj, **args):
        return [self.warnMessage('\n'.join(obj._summary()), 
                                 "No output expected")]
        

class RelionAutopickViewer(Viewer):
    """ Class to visualize Relion postprocess protocol """
    _targets = [ProtRelionAutopick]
    _environments = [DESKTOP_TKINTER]
    
    def visualize(self, obj, **args):
        micPath, coordPath = obj.writeXmippOutputCoords()
        import pyworkflow.em.packages.xmipp3 as xmipp3
        xmipp3.viewer.launchSupervisedPickerGUI(micPath, coordPath, self.protocol)


class RelionPolishViewer(ProtocolViewer):
    """ Viewer for Relion polishing protocol"""
    _targets = [ProtRelionPolish]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer relion polish'
    
    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        
        group = form.addGroup('Particles')
        group.addParam('displayShinyParticles', params.LabelParam, 
                       label='Display Shiny Particles')
        
        group = form.addGroup('Volumes')
        
        group.addParam('showHalves', params.EnumParam, choices=['half1', 'half2', 'both', 'final shiny'], default=0,
                       label='Volume to visualize',
                       help='Select which half do you want to visualize.')
        group.addParam('viewFrame', params.EnumParam, choices=['all', 'selection'], default=0, 
                       display=params.EnumParam.DISPLAY_HLIST, condition="showHalves<3",
                       label="Frame to visualize",
                       help="""
*all*: all frames volumes will be visualized.
*selection*: you may specify a range of frames.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        group.addParam('frameSelection', params.NumericRangeParam, 
                       condition='showHalves<3 and viewFrame==1' ,
                       label="Frames list", default = 1,
                       help="Write the frame list to visualize.")
        group.addParam('displayVol', params.EnumParam, choices=['slices', 'chimera'], 
                       default=VOLUME_SLICES, display=params.EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayAngDist', params.EnumParam, choices=['2D plot', 'chimera'], 
                       default=ANGDIST_2DPLOT, display=params.EnumParam.DISPLAY_HLIST,
                       label='Display angular distribution', condition="showHalves==3",
                       help='*2D plot*: display angular distribution as interactive 2D in matplotlib.\n'
                            '*chimera*: display angular distribution using Chimera with red spheres.')
        group.addParam('spheresScale', params.IntParam, default=100, 
                       expertLevel=LEVEL_ADVANCED, condition="showHalves==3 and displayAngDist==%d" % ANGDIST_CHIMERA,
                       label='Spheres size',
                       help='')
        group = form.addGroup('Plots')
        group.addParam('resolutionPlotsFSC', params.EnumParam,
                       choices=['Corrected', 'Unmasked Maps', 'Masked Maps', 'Phase Randomized Masked Maps', 'all'],
                       default=FSC_CORRECTED, display=params.EnumParam.DISPLAY_COMBO,
                       label='Display resolution plots (FSC)',
                       help='')
        group.addParam('resolutionThresholdFSC', params.FloatParam, default=0.143, 
                       expertLevel=LEVEL_ADVANCED,
                       label='Threshold in resolution plots',
                       help='')
        group.addParam('guinierPlots', params.LabelParam,
                       default=True, condition="showHalves<3",
                       label='Display guinier plots',
                       help='')
        group.addParam('bfactorsPlot', params.LabelParam,
                       default=True,
                       label='Display bfactors plot',
                       help='')
    
    def _getVisualizeDict(self):
        self._load()
        visualizeDict = {
                'displayShinyParticles': self._showShinyParticles,
                'displayVol': self._showVolumes,
                'displayAngDist': self._showAngularDistribution,
                'resolutionPlotsFSC': self._showFSC,
                'guinierPlots': self._showGuinier,
                'bfactorsPlot': self._showBfactors
                }

        # If the is some error during the load, just show that instead
        # of any viewer
        if self._errors:
            for k in visualizeDict.keys():
                visualizeDict[k] = self._showErrors

        return visualizeDict

    def _showErrors(self, e=None):
        views = []
        self.errorList(self._errors, views)
        return views

    def _viewAll(self, *args):
        pass
    
#===============================================================================
# displayShinyParticles     
#===============================================================================
    def _showShinyParticles(self, paramName=None):
        partsId = self.protocol.outputParticles.strId()
        fn = self.protocol.outputParticles.getFileName()
        view = self.createScipionPartView(fn, partsId)
        
        return [view]
    
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
        path = self.protocol._getExtraPath('relion_viewer_volumes.sqlite')
        samplingRate = self.protocol.outputParticles.getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if exists(volFn.replace(':mrc', '')):
                files.append(volFn)
        self.createVolumesSqlite(files, path, samplingRate)
        return [em.ObjectView(self._project, self.protocol.strId(), path)]
    
    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        
        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for volFn in volumes:
                # We assume that the chimera script will be generated
                # at the same folder than relion volumes
                vol = volFn.replace(':mrc', '')
                localVol = os.path.basename(vol)
                if exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = em.ChimeraView(cmdFile)
        else:
            view = em.ChimeraClientView(volumes[0])
            
        return [view]
    
#===============================================================================
# showAngularDistribution
#===============================================================================
    def _showAngularDistribution(self, paramName=None):
        if self.displayAngDist == ANGDIST_CHIMERA:
            view = self._createAngDistChimera()
        elif self.displayAngDist == ANGDIST_2DPLOT:
            view = self._createAngDist2D()
        return [view]
    
    def _createAngDistChimera(self):
        # Common variables to use
        nparts = self.protocol.outputParticles.getSize()
        radius = self.spheresScale.get()
        if radius < 0:
            radius = self.protocol.maskDiameterA.get()/2
            
        volFn = self._getVolumeNames()[0]
        if exists(volFn.replace(":mrc", "")):
            sqliteFn = self.protocol._getFileName('projections', iter=1, half="shiny_", ref3d=1)
            if not exists(sqliteFn):
                self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(self.protocol._getFileName('shiny')))
            return em.ChimeraClientView(volFn, angularDistFile=sqliteFn, spheresDistance=radius)
        else:
            raise Exception("They aren't shiny particles to show it's angular distribution.")
    
    def _createAngDist2D(self):
        # Common variables to use
        nparts = self.protocol.outputParticles.getSize()
        n = 1
        gridsize = self._getGridSize(n)
        
        plotter = RelionPlotter(x=gridsize[0], y=gridsize[1], 
                                 mainTitle="Shiny Particles", windowTitle="Angular Distribution")
        dataStar = self.protocol._getFileName('shiny')
        if exists(dataStar):
            sqliteFn = self.protocol._getFileName('projections', iter=1, half="shiny_", ref3d=1)
            if not exists(sqliteFn):
                self.createAngDistributionSqlite(sqliteFn, nparts, itemDataIterator=self._iterAngles(dataStar))
            plotter.plotAngularDistributionFromMd(sqliteFn, "")
            return plotter
        else:
            return
    
#===============================================================================
# plotFSC            
#===============================================================================
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        n = 1
        gridsize = [1, 1]
        
        md.activateMathExtensions()
        
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1], windowTitle='Resolution FSC')
        a = xplotter.createSubPlot("FSC Shiny Particles", 'Angstroms^-1', 'FSC', yformat=False)
        legends = []
        modelStar = self.protocol._getFileName('fsc_shiny')
        for label in self._getFSCLabels():
            if exists(modelStar):
                model = 'fsc@' + modelStar
                self._plotFSC(a, model, label)
                legends.append(self._getLegend(label))
        xplotter.showLegend(legends)
        if threshold < self.maxfsc:
            a.plot([self.minInv, self.maxInv], [threshold, threshold], color='black', linestyle='--')
        a.grid(True)
        
        return [xplotter]
    
    def _plotFSC(self, a, model, label):
        mdStar = md.MetaData(model)
        resolution_inv = [mdStar.getValue(md.RLN_RESOLUTION, id) for id in mdStar]
        fsc = [mdStar.getValue(label, id) for id in mdStar]
        self.maxfsc = max(fsc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, fsc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])
            
#===============================================================================
# plotGuinier
#===============================================================================
    def _showGuinier(self, paramName=None):
        gridsize = [1, 1]
        md.activateMathExtensions()
        
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1], windowTitle='Relative Guinier Plot')
        a = xplotter.createSubPlot("", 'Angstroms^-2', 'log(Amplitude)', yformat=False)
        legends = []
        for frame in self._frames:
            modelStar = self.protocol._getFileName('guinier_frame', frame=frame)
            if exists(modelStar):
                model = 'relative_guinier@' + modelStar
                self._plotGuinier(a, model)
                legends.append("frame %d" % frame)
        xplotter.showLegend(legends)
        a.grid(True)
        
        return [xplotter]
    
    def _plotGuinier(self, a, model):
        mdStar = md.MetaData(model)
        resolSqInv = [mdStar.getValue(md.RLN_POSTPROCESS_GUINIER_RESOL_SQUARED, id) for id in mdStar]
        logAmp = [mdStar.getValue(md.RLN_POSTPROCESS_GUINIER_VALUE_IN, id) for id in mdStar]
        self.maxfsc = max(logAmp)
        self.minInv = min(resolSqInv)
        self.maxInv = max(resolSqInv)
        a.plot(resolSqInv, logAmp)
        a.xaxis.set_major_formatter(self._plotFormatter)
    
#===============================================================================
# plotBfactors
#===============================================================================
    def _showBfactors(self, paramName=None):
        gridsize = self._getGridSize(2)
        md.activateMathExtensions()
        
        xplotter = RelionPlotter(x=gridsize[0], y=gridsize[1], windowTitle='Per Frame B-factors')
        modelStar = self.protocol._getFileName('bfactors')
        model = 'perframe_bfactors@' + modelStar
        for title in ['B-factor', 'Scale factor']:
            a = xplotter.createSubPlot(title, 'Frame', title, yformat=False)
            if title == "B-factor":
                label = md.RLN_POSTPROCESS_BFACTOR
            else:
                label = md.RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT
            self._plotBfactor(a, model, label)
        
        return [xplotter]
    
    def _plotBfactor(self, a, model, label):
        mdStar = md.MetaData(model)
        frames = [mdStar.getValue(md.RLN_IMAGE_FRAME_NR, id) for id in mdStar]
        bfactors = [mdStar.getValue(label, id) for id in mdStar]
        a.plot(frames, bfactors)
    
#===============================================================================
# Utils Functions
#===============================================================================
    def _load(self):
        self.protocol._initialize()  # Load filename templates
        self._errors = []
        self.lastIter = self.protocol._lastIter()
        halves = getattr(self, 'showHalves', None)
        if self.viewFrame.get() == 0 and halves < 3:
            self._frames = range(1, self.protocol._lastFrame())
        elif halves < 3:
            frameSelection = self.frameSelection.get()
            if frameSelection and frameSelection.strip():
                self._frames = self._getListFromRangeString(frameSelection)
            else:
                self._errors.append('Please provide FRAMES selection')
        else:
            self._frames = [1]
            
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        
    def createScipionPartView(self, filename, partsId, viewParams={}):
        labels = 'enabled id _filename _transform._matrix'
        viewParams = {showj.ORDER:labels,
                      showj.VISIBLE: labels, showj.RENDER:'_filename',
                      'labels': 'id',
                      }
        return em.ObjectView(self._project, 
                             self.protocol.strId(), filename, other=partsId,
                             env=self._env,
                             viewParams=viewParams)
    
    def _getPrefixes(self):
        prefixes = self.protocol.PREFIXES
        halves = getattr(self, 'showHalves', None)
        if halves:
            if halves == 0:
                prefixes = ['half1_']
            elif halves == 1:
                prefixes = ['half2_']
            elif halves == 3:
                prefixes = ['shiny']
        return prefixes
    
    def _getVolumeNames(self):
        vols = []
        prefixes = self._getPrefixes()
        if prefixes == ["shiny"]:
            volFn = self.protocol._getFileName('volume_shiny')
            if exists(volFn.replace(':mrc', '')):
                vols.append(volFn)
        else:
            for frm in self._frames:
                for prefix in prefixes:
                    volFn = self.protocol._getFileName('volume_frame', halve=prefix, frame=frm, ref3d=1)
                    if exists(volFn.replace(':mrc', '')):
                        vols.append(volFn)
        return vols
    
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        import math
        inv = 999.
        if value:
            inv = 1/math.sqrt(value)
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
    
    def _iterAngles(self, fn):
        mdOut = md.MetaData(fn)
        for objId in mdOut:
            rot = mdOut.getValue(md.RLN_ORIENT_ROT, objId)
            tilt = mdOut.getValue(md.RLN_ORIENT_TILT, objId)
            yield rot, tilt
    
    def _getFSCLabels(self):
        if self.resolutionPlotsFSC.get() == 0:
            return [md.RLN_POSTPROCESS_FSC_TRUE]
        elif self.resolutionPlotsFSC.get() == 1:
            return [md.RLN_POSTPROCESS_FSC_UNMASKED]
        elif self.resolutionPlotsFSC.get() == 2:
            return [md.RLN_POSTPROCESS_FSC_MASKED]
        elif self.resolutionPlotsFSC.get() == 3:
            return [md.RLN_POSTPROCESS_FSC_RANDOM_MASKED]
        else:
            return [md.RLN_POSTPROCESS_FSC_TRUE, md.RLN_POSTPROCESS_FSC_UNMASKED,
                    md.RLN_POSTPROCESS_FSC_MASKED, md.RLN_POSTPROCESS_FSC_RANDOM_MASKED]
    
    def _getLegend(self, label):
        if label == md.RLN_POSTPROCESS_FSC_TRUE:
            return 'Corrected'
        elif label == md.RLN_POSTPROCESS_FSC_UNMASKED:
            return 'Unmasked Maps'
        elif label == md.RLN_POSTPROCESS_FSC_MASKED:
            return 'Masked Maps'
        else:
            return 'Phase Randomized Masked Maps'


class RelionSortViewer(Viewer):
    """ Wrapper to visualize Relion sort protocol results
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionSortParticles]

    def _visualize(self, obj, **kwargs):
        views = []

        if obj.hasAttribute('outputParticles'):  # Protocol finished
            particles = obj.outputParticles
            labels = ('id enabled _index _filename _rlnSelectParticlesZscore '
                      '_coordinate._rlnAutopickFigureOfMerit _sampling '
                      '_ctfModel._defocusU _ctfModel._defocusV '
                      '_ctfModel._defocusAngle _transform._matrix')
            sortBy = '_rlnSelectParticlesZscore asc'
            strId = particles.strId()
            fn = particles.getFileName()
            views.append(em.ObjectView(self._project, strId, fn,
                                       viewParams={showj.ORDER: labels,
                                                   showj.VISIBLE: labels,
                                                   showj.SORT_BY: sortBy,
                                                   showj.RENDER: '_filename'}))

            fn = obj._getExtraPath('input_particles_sorted.star')
            mdFn = md.MetaData(fn)
            # If Zscore in output images plot Zscore particle sorting
            if mdFn.containsLabel(md.RLN_SELECT_PARTICLES_ZSCORE):
                # sort output by Z-score
                mdFn.sort(md.RLN_SELECT_PARTICLES_ZSCORE)
                xplotter = RelionPlotter(windowTitle="Zscore particles sorting")
                xplotter.createSubPlot("Particle sorting", "Particle number",
                                       "Zscore")
                xplotter.plotMd(mdFn, False,
                                mdLabelY=md.RLN_SELECT_PARTICLES_ZSCORE)
                views.append(xplotter)

        return views


# Color maps
COLOR_JET = 0
COLOR_TERRAIN = 1
COLOR_GIST_EARTH = 2
COLOR_GIST_NCAR = 3
COLOR_GNU_PLOT = 4
COLOR_GNU_PLOT2 = 5
COLOR_OTHER = 6

COLOR_CHOICES = em.OrderedDict()

COLOR_CHOICES[COLOR_JET] = 'jet'
COLOR_CHOICES[COLOR_TERRAIN] = 'terrain'
COLOR_CHOICES[COLOR_GIST_EARTH] = 'gist_earth'
COLOR_CHOICES[COLOR_GIST_NCAR] = 'gist_ncar'
COLOR_CHOICES[COLOR_GNU_PLOT] = 'gnuplot'
COLOR_CHOICES[COLOR_GNU_PLOT2] = 'gnuplot2'
COLOR_CHOICES[COLOR_OTHER] = 'other'

binaryCondition = ('(colorMap == %d) ' % (COLOR_OTHER))


# Axis code
AX_X = 0
AX_Y = 1
AX_Z = 2


class RelionLocalResViewer(ProtocolViewer):
    """
    Visualization tools for local resolution results.

    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionLocalRes]
    _label = 'viewer relion localres'

    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, *args, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('Colored resolution')
        group.addParam('colorMap', params.EnumParam,
                       choices=COLOR_CHOICES.values(),
                       default=COLOR_JET,
                       label='Color map',
                       help='Select the color map to apply to the resolution '
                            'map. http://matplotlib.org/1.3.0/examples/color/'
                            'colormaps_reference.html.')
        group.addParam('otherColorMap', params.StringParam, default='jet',
                       condition=binaryCondition,
                       label='Customized Color map',
                       help='Name of a color map to apply to the resolution '
                            'map. Valid names can be found at http://'
                            'matplotlib.org/1.3.0/examples/color/'
                            'colormaps_reference.html')

        group = form.addGroup('Slices')
        group.addParam('sliceAxis', params.EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Slice axis')
        group.addParam('doShowVolumeSlices', params.LabelParam,
                      label="Show volume slices")
        group.addParam('doShowChimera', params.LabelParam,
                       label="Show colored map in Chimera", default=True)

    def _getVisualizeDict(self):
        self.protocol._createFilenameTemplates()
        return {
               'doShowVolumeSlices': self._showVolumeSlices,
               'doShowChimera': self._showChimera,
               }

# ==============================================================================
# doShowVolumeSlices
# ==============================================================================
    def _showVolumeSlices(self, param=None):
        imageFile = self.protocol._getFileName('resolMap')
        imgData, minRes, maxRes = self._getImgData(imageFile)
        
        xplotter = RelionPlotter(x=2, y=2, mainTitle="Local Resolution Slices "
                                                     "along %s-axis."
                                                     %self._getAxis())
        for i in xrange(4):
            slice = self._getSlice(i+1, imgData)
            a = xplotter.createSubPlot("Slice %s" % (slice), '', '')
            matrix = self._getSliceImage(imgData, i+1, self._getAxis())
            plot = xplotter.plotMatrix(a, matrix, minRes, maxRes,
                                       cmap=self._getColorName(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]

# ==============================================================================
# showChimera
# ==============================================================================
    def _showChimera(self, param=None):
        cmdFile = self.protocol._getExtraPath('chimera_local_res.cmd')
        self._createChimeraScript(cmdFile)
        view = em.ChimeraView(cmdFile)
        return [view]

# ==============================================================================
# Utils Functions
# ==============================================================================
    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def _getImgData(self, imgFile):
        import numpy as np
        img = em.ImageHandler().read(imgFile+":mrc")
        imgData = img.getData()

        maxRes = np.amax(imgData)
        imgData2 = np.ma.masked_where(imgData < 0.1, imgData, copy=True)
        minRes = np.amin(imgData2)

        return imgData2, minRes, maxRes

    def _getSlice(self, index, volumeData):
        return int((index + 3) * volumeData.shape[0] / 9)

    def _getSliceImage(self, volumeData, index, dataAxis):
        slice = self._getSlice(index, volumeData)
        if dataAxis == 'y':
            imgSlice = volumeData[:, slice, :]
        elif dataAxis == 'x':
            imgSlice = volumeData[:, :, slice]
        else:
            imgSlice = volumeData[slice, :, :]
        return imgSlice
    
    def _getColorName(self):
        if self.colorMap.get() != COLOR_OTHER:
            return COLOR_CHOICES[self.colorMap.get()]
        else:
            return self.otherColorMap.get()

    def _createChimeraScript(self, scriptFile):
        import pyworkflow.gui.plotter as plotter
        from itertools import izip
        fhCmd = open(scriptFile, 'w')
        imageFile = os.path.abspath(self.protocol._getFileName('resolMap'))
        
        _, minRes, maxRes = self._getImgData(imageFile)
        
        stepColors = self._getStepColors(minRes, maxRes)
        colorList = plotter.getHexColorList(stepColors, self._getColorName())
        
        fnVol = os.path.abspath(self.protocol._getFileName('finalMap'))

        fhCmd.write("background solid white\n")
        
        fhCmd.write("open %s\n" % fnVol)
        fhCmd.write("open %s\n" % (imageFile))
        
        sampRate = self.protocol.outputVolume.getSamplingRate()
        fhCmd.write("volume #0 voxelSize %s\n" % (str(sampRate)))
        fhCmd.write("volume #1 voxelSize %s\n" % (str(sampRate)))
        fhCmd.write("volume #1 hide\n")

        scolorStr = ''
        for step, color in izip(stepColors, colorList):
            scolorStr += '%s,%s:' % (step, color)
        scolorStr = scolorStr[:-1]
        line = ("scolor #0 volume #1 perPixel false cmap " + scolorStr + "\n")
        fhCmd.write(line)

        scolorStr2 = ''
        for step, color in izip(stepColors, colorList):
            indx = stepColors.index(step)
            if ((indx % 4) != 0):
                scolorStr2 += '" " %s ' % color
            else:
                scolorStr2 += '%s %s ' % (step, color)
        line = ("colorkey 0.01,0.05 0.02,0.95 labelColor None "
                + scolorStr2 + " \n")
        fhCmd.write(line)
        fhCmd.close()

    def _getStepColors(self, minRes, maxRes, numberOfColors=13):
        inter = (maxRes - minRes) / (numberOfColors - 1)
        rangeList = []
        for step in range(0, numberOfColors):
            rangeList.append(round(minRes + step * inter, 2))
        return rangeList
