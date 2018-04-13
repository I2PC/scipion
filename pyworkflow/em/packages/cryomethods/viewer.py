# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Javier Vargas Balbuena (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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
from pyworkflow.em.plotter import EmPlotter
import pyworkflow.protocol.params as params
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
from initial_volume_selector import ProtInitialVolumeSelector
from convert import relionToLocation

ITER_LAST = 0
ITER_SELECTION = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

CHIMERADATAVIEW = 0

CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS,
                 md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
                 md.RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES]

class VolSelPlotter(EmPlotter):
    """ Class to create several plots with Xmipp utilities"""
    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        EmPlotter.__init__(self, x, y, mainTitle, **kwargs)

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


class VolumeSelectorViewer(ProtocolViewer):
    """ This protocol serve to analyze the results of Initial
    Volume Selector protocol.
    """
    _targets = [ProtInitialVolumeSelector]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer volume selector'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help='*last*: only the last iteration will be '
                           'visualized.\n'
                           '*selection*: you may specify a range of '
                           'iterations.\n'
                           'Examples:\n'
                           '"1,5-8,10" -> [1,5,6,7,8,10]\n'
                           '"2,6,9-11" -> [2,6,9,10,11]\n'
                           '"2 5, 6-8" -> [2,5,6,7,8] ')

        form.addParam('iterSelection', params.NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Volumes')
        group.addParam('displayVol', params.EnumParam,
                       choices=['slices', 'chimera'], default=VOLUME_SLICES,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z '
                            'axis.\n'
                            '*chimera*: display volumes as surface with '
                            'Chimera.')

        group = form.addGroup('Resolution')
        group.addParam('figure', params.EnumParam, default=0,
                       choices=['new', 'active'],
                       label='Figure', display=params.EnumParam.DISPLAY_HLIST)
        group.addParam('resolutionPlotsSSNR', params.LabelParam, default=True,
                       label='Display SSNR plots',
                       help='Display signal to noise ratio plots (SSNR)')
        group = form.addGroup('Overall')
        group.addParam('showPMax', params.LabelParam, default=True,
                      label="Show average PMax",
                      help='Average (per class) of the maximum value of '
                           'normalized probability function')
        group.addParam('showChanges', params.LabelParam, default=True,
                      label='Changes in Offset, Angles and Classes',
                      help='Visualize changes in orientation, offset and '
                           'number images assigned to each class')
        group.addParam('plotClassDistribution', params.LabelParam, default=True,
                      label='Plot class distribution over iterations',
                      help='Plot each class distribution over iterations as '
                           'bar plots.')

    def _getVisualizeDict(self):
        visualizeDict = {'displayVol': self._showVolumes,
                         'showPMax': self._showPMax,
                         'showChanges': self._showChanges,
                         'resolutionPlotsSSNR': self._showSSNR,
                         'plotClassDistribution': self._plotClassDistribution
                         }
        self._load()

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

# ==============================================================================
# ShowVolumes
# ==============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumesSqlite()

    def _showVolumesSqlite(self):
        """ Write (if it is needed) an sqlite with all volumes selected for
        visualization. """

        view = []
        if (self.viewIter == ITER_LAST and
                getattr(self.protocol, 'outputVolumes', None) is not None):
            fn = self.protocol.outputVolumes.getFileName()

            view.append(self.createView(filename=fn,
                                        viewParams=self._getViewParams()))
        else:
            for it in self._iterations:
                volSqlte = self.protocol._getIterVolumes(it)
                view.append(self.createView(filename=volSqlte,
                                            viewParams=self._getViewParams()))
        return view

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
# ShowPMax
#===============================================================================
    def _showPMax(self, paramName=None):
        labels = [md.RLN_MLMODEL_AVE_PMAX, md.RLN_PARTICLE_PMAX]

        mdIters = md.MetaData()
        iterations = range(self.firstIter, self.lastIter+1)

        # range (firstIter, # self._visualizeLastIteration+1):
        # #always list all iteration
        for it in iterations:
            objId = mdIters.addObject()
            mdIters.setValue(md.MDL_ITER, it, objId)
            for i, prefix in enumerate(self.protocol.PREFIXES):
                filename = self.protocol._getFileName(prefix + 'model', iter=it)
                fn = 'model_general@' + filename
                mdModel = md.RowMetaData(fn)
                pmax = mdModel.getValue(md.RLN_MLMODEL_AVE_PMAX)
                mdIters.setValue(labels[i], pmax, objId)
        fn = self.protocol._getFileName('all_avgPmax_xmipp')
        mdIters.write(fn)

        colors = ['g', 'b']

        xplotter = VolSelPlotter()
        xplotter.createSubPlot("Avg PMax per Iterations", "Iterations", "Avg PMax")

        for label, color in zip(labels, colors):
            xplotter.plotMd(mdIters, md.MDL_ITER, label, color)

        if len(self.protocol.PREFIXES) > 1:
            xplotter.showLegend(self.protocol.PREFIXES)

        return [self.createView(fn), xplotter]

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

        xplotter = VolSelPlotter()
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
            """Returns a function that maps each index in 0, 1, ... N-1 
            to a distinct RGB color."""
            color_norm = colors.Normalize(vmin=0, vmax=N)#-1)
            scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')

            def map_index_to_rgb_color(index):
                return scalar_map.to_rgba(index)

            return map_index_to_rgb_color

        cmap = get_cmap(len(classInfo))

        for classId in sorted(classInfo.keys()):
            values = classInfo[classId][l]
            ax.bar(ind, values, width, label='class %s' % classId,
                   bottom=bottomValues, color=cmap(classId))
            bottomValues = [a+b for a, b in zip(bottomValues, values)]

        ax.get_xaxis().set_ticks([i + 0.25 for i in ind])
        ax.get_xaxis().set_ticklabels([str(i) for i in ind])
        ax.legend(loc='upper left', fontsize='xx-small')

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
            for label in CHANGE_LABELS:
                mdIters.setValue(label, mdOptimiser.getValue(label), objId)
        fn = self.protocol._getFileName('all_changes_xmipp')
        mdIters.write(fn)

        return [self.createView(fn)]

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
        xplotter = VolSelPlotter(x=gridsize[0], y=gridsize[1],
                                 figure=self._getFigure())

        for prefix in prefixes:
            for ref3d in self._refsList:
                plot_title = 'Resolution SSNR %s, for Class %s' % (prefix, ref3d)
                a = xplotter.createSubPlot(plot_title, 'Angstroms^-1',
                                           'log(SSNR)', yformat=False)
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
# Utils Functions
#===============================================================================
    def _getZoom(self):
        # Ensure that classes are shown at least at 128 px to
        # properly see the rlnClassDistribution label.[[
        dim = self.protocol.inputParticles.get().getDim()[0]
        if dim < 128:
            zoom = 128*100/dim
        else:
            zoom = 100
        return zoom

    def _validate(self):
        if self.lastIter is None:
            return ['There are not iterations completed.']

    def _getViewParams(self):
        labels = ('enabled id _filename _rlnClassDistribution '
                 '_rlnAccuracyRotations _rlnAccuracyTranslations '
                  '_rlnEstimatedResolution')
        viewParams = {showj.ORDER: labels,
                      showj.MODE: showj.MODE_MD,
                      showj.VISIBLE: labels,
                      showj.RENDER: '_filename',
                      showj.SORT_BY: '_rlnClassDistribution desc',
                      showj.ZOOM: str(self._getZoom())
                      }
        return viewParams


    def createView(self, filename, viewParams={}):
        return em.ObjectView(self._project, self.protocol.strId(),
                             filename, viewParams=viewParams)

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

        volSize = self.protocol.inputVolumes.get().getSize()
        self._refsList = range(1, volSize+1)

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

    def _getVolumeNames(self):
        vols = []
        for it in self._iterations:
            for ref3d in self._refsList:
                volFn = self.protocol._getFileName('volume', iter=it, ref3d=ref3d)
                vols.append(volFn)
        return vols

    def _getModelStar(self, prefix, it):
        return self.protocol._getFileName(prefix + 'model', iter=it)