# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
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
"""
This module implement the first version of viewers using 
around xmipp_showj visualization program.
"""
import os, math

from pyworkflow.gui.project import ProjectWindow
import pyworkflow.gui.text as text
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER,
                               WEB_DJANGO)
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
import pyworkflow.em.showj as showj
from pyworkflow.em.viewer import (ObjectView, DataView,
                                  ChimeraView, ChimeraClientView, ClassesView)
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.executor import StepExecutor
from pyworkflow.protocol.params import (LabelParam, NumericRangeParam,
                                        EnumParam, FloatParam, IntParam,
                                        BooleanParam)
import pyworkflow.utils as pwutils

from protocol_boxing import EmanProtBoxing
from protocol_ctf import EmanProtCTFAuto
from protocol_initialmodel import EmanProtInitModel
from protocol_refine2d import EmanProtRefine2D
from protocol_refine2d_bispec import EmanProtRefine2DBispec
from protocol_refineasy import EmanProtRefine
from protocol_tiltvalidate import EmanProtTiltValidate
from convert import loadJson
from eman2 import getEmanProgram

LAST_ITER = 0
ALL_ITERS = 1
SELECTED_ITERS = 2

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

TILT_SCATTER = 0
TILT_CONTOUR = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

FSC_UNMASK = 0
FSC_MASK = 1
FSC_MASKTIGHT = 2
FSC_ALL = 3

HALF_EVEN = 0
HALF_ODD = 1
FULL_MAP = 2
ALL_MAPS = 3

OBJCMD_CLASSAVG_PROJS = 'Show class-averages/projections'
OBJCMD_PROJS = 'Show only projections'
OBJCMD_INITVOL = 'Show initial random volume'


class EmanViewer(XmippViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER]
    _targets = [EmanProtBoxing, EmanProtInitModel]

    def _visualize(self, obj, **args):

        if isinstance(obj, EmanProtBoxing):
            coords = obj.getCoords()
            if coords:
                XmippViewer._visualize(self, obj.outputCoordinates)

        elif isinstance(obj, EmanProtInitModel):
            obj = obj.outputVolumes
            fn = obj.getFileName()
            labels = 'id enabled comment _filename '
            objCommands = "'%s' '%s' '%s'" % (OBJCMD_CLASSAVG_PROJS,
                                              OBJCMD_PROJS,
                                              OBJCMD_INITVOL)

            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={showj.MODE: showj.MODE_MD,
                                                      showj.VISIBLE: labels,
                                                      showj.RENDER: '_filename',
                                                      showj.OBJCMDS: objCommands},
                                          ))


def showExtraFile(volumeSet, volId, suffix):
    vol = volumeSet[volId]
    volFn = vol.getFileName().replace('.hdf', suffix)
    DataView(volFn).show()


def showClassAvgProjs(volumeSet, volId):
    showExtraFile(volumeSet, volId, '_aptcl.hdf')


def showProjs(volumeSet, volId):
    showExtraFile(volumeSet, volId, '_proj.hdf')


def showInitialRandomVolume(volumeSet, volId):
    showExtraFile(volumeSet, volId, '_init.hdf')


ProjectWindow.registerObjectCommand(OBJCMD_CLASSAVG_PROJS, showClassAvgProjs)
ProjectWindow.registerObjectCommand(OBJCMD_PROJS, showProjs)
ProjectWindow.registerObjectCommand(OBJCMD_INITVOL, showInitialRandomVolume)

class Refine2DViewer(ProtocolViewer):
    """ Visualization of e2refine2d results. """

    _targets = [EmanProtRefine2D, EmanProtRefine2DBispec]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer refine 2d'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam,
                      label="Which iteration do you want to visualize?",
                      default=0, choices=['last', 'all', 'selection'],
                      display=EnumParam.DISPLAY_LIST)
        form.addParam('iterSelection', NumericRangeParam, default='1',
                      label='Selected iterations',
                      condition='iterToShow==%d' % SELECTED_ITERS,
                      help="""
*last*: only the last iteration will be visualized.
*all*: all iterations  will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]
                   """)
        form.addParam('showClasses', LabelParam,
                      label='Show classification in Scipion',
                      help='Display final unaligned class averages')
        form.addParam('showAllRef', LabelParam,
                      condition='_protocolIsNotBispec',
                      label='Show sorted and aligned class averages',
                      help='')
        form.addParam('showBasis', LabelParam,
                      label='Show MSA vectors',
                      help='Display MSA basis vectors, which may be useful '
                           'if looking for signs of specific symmetries, etc.')
        form.addParam('showAliRef', LabelParam,
                      condition='_protocolIsNotBispec',
                      label='Show alignment references',
                      help='May be useful to look at which averages were '
                           'used as 2-D alignment references')

    def _getVisualizeDict(self):
        self._load()
        return {'showClasses': self._showClasses,
                'showAllRef': lambda paramName: self._showMisc(key='allrefs'),
                'showBasis': lambda paramName: self._showMisc(key='basis'),
                'showAliRef': lambda paramName: self._showMisc(key='alirefs')
                }

    def _showClasses(self, paramName=None):
        views = []
        if (self.iterToShow == LAST_ITER and
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

    def createScipionView(self, filename):
        labels =  'enabled id _size _representative._filename '
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.RENDER:'_representative._filename',
                      showj.SORT_BY: '_size desc'
                      }
        inputParticlesId = self.protocol.inputParticles.get().strId()
        view = ClassesView(self._project,
                           self.protocol.strId(), filename, other=inputParticlesId,
                           env=self._env,
                           viewParams=viewParams)

        return view

    def _showMisc(self, key):
        views = []
        if self.iterToShow.get() == LAST_ITER:
            last = self.protocol._lastIter()
            fn = self.protocol._getFileName(key, iter=last)
            v = self.createScipionView(fn)
            views.append(v)
        else:
            for it in self._iterations:
                fn = self.protocol._getFileName(key, iter=it)
                v = self.createScipionView(fn)
                views.append(v)

        return views

    def _load(self):
        """ Load selected iterations and classes 2D for visualization mode. """
        self.protocol._createFilenameTemplates()
        self.protocol._createIterTemplates()
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()

        if self.iterToShow.get() == LAST_ITER:
            self._iterations = [self.lastIter]
        elif self.iterToShow.get() == ALL_ITERS:
            self._iterations = range(1, self.lastIter + 1)
        elif self.iterToShow.get() == SELECTED_ITERS:
            self._iterations = self._getListFromRangeString(
                self.iterSelection.get())

    def _protocolIsNotBispec(self):
        return self.protocol.getClassName() == 'EmanProtRefine2D'


class RefineEasyViewer(ProtocolViewer):
    """ Visualization of e2refine_easy results. """

    _targets = [EmanProtRefine]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Refine Easy'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam,
                      label="Which iteration do you want to visualize?",
                      default=0, choices=['last', 'all', 'selection'],
                      display=EnumParam.DISPLAY_LIST)
        form.addParam('iterSelection', NumericRangeParam, default='1',
                      label='Selected iterations',
                      condition='iterToShow==%d' % SELECTED_ITERS,
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

        group.addParam('showHalves', EnumParam,
                       choices=['half even', 'half odd', 'full map',
                                'all maps'],
                       default=HALF_EVEN,
                       label='Map to visualize',
                       help='Select which map do you want to visualize.')
        group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'],
                       default=VOLUME_SLICES, display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayAngDist', EnumParam,
                       choices=['2D plot', 'chimera'],
                       default=ANGDIST_2DPLOT, display=EnumParam.DISPLAY_HLIST,
                       label='Display angular distribution',
                       help='*2D plot*: display angular distribution as interative 2D in'
                            ' matplotlib.\n*chimera*: display angular distribution using'
                            ' Chimera with red spheres.')
        group.addParam('spheresScale', IntParam, default=100,
                       expertLevel=LEVEL_ADVANCED,
                       label='Spheres size',
                       help='')

        group = form.addGroup('Resolution')

        group.addParam('resolutionPlotsFSC', EnumParam,
                       choices=['unmasked', 'masked', 'masked tight', 'all'],
                       default=FSC_UNMASK, display=EnumParam.DISPLAY_COMBO,
                       label='Display resolution plots (FSC)',
                       help='*unmasked*: display FSC of unmasked maps.\n'
                            '*masked*: display FSC of masked maps.\n'
                            '*masked tight*: display FSC of masked tight maps.')

        group.addParam('resolutionThresholdFSC', FloatParam, default=0.143,
                       expertLevel=LEVEL_ADVANCED,
                       label='Threshold in resolution plots',
                       help='')

        form.addParam('showHtmlReport', LabelParam,
                      label='Show HTML report',
                      help='An HTML report file will be generated as this '
                           'program runs, telling you exactly what it decided '
                           'to do and why, as well as giving information about '
                           'runtime, etc while the job is still running.')

    def _getVisualizeDict(self):
        self._load()
        return {'showImagesAngularAssignment': self._showImagesAngularAssignment,
                'displayVol': self._showVolumes,
                'displayAngDist': self._showAngularDistribution,
                'resolutionPlotsFSC': self._showFSC,
                'showHtmlReport': self._showHtmlReport
                }

    # ===============================================================================
    # showImagesAngularAssignment
    # ===============================================================================

    def _showImagesAngularAssignment(self, paramName=None):
        views = []

        for it in self._iterations:
            fn = self.protocol._getIterData(it)
            v = self.createScipionPartView(fn)
            views.append(v)

        return views

    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        labels = 'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, showj.RENDER: '_filename',
                      'labels': 'id',
                      }
        return ObjectView(self._project,
                          self.protocol.strId(), filename,
                          other=inputParticlesId,
                          env=self._env, viewParams=viewParams)

    # ===============================================================================
    # ShowVolumes
    # ===============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()

        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for vol in volumes:
                # We assume that the chimera script will be generated
                # at the same folder than eman volumes
                if os.path.exists(vol):
                    localVol = os.path.relpath(vol,
                                               self.protocol._getExtraPath())
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = ChimeraView(cmdFile)
        else:
            view = ChimeraClientView(volumes[0])

        return [view]

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """

        path = self.protocol._getExtraPath('eman2_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()

        vols = self._getVolumeNames()
        files = []
        for vol in vols:
            if os.path.exists(vol):
                files.append(vol)
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    # ===============================================================================
    # showAngularDistribution
    # ===============================================================================
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
        radius = self.spheresScale.get()
        volumes = self._getVolumeNames()
        self.protocol._execEmanProcess(self.protocol._getRun(), it)
        angularDist = self.protocol._getFileName("angles", iter=it)

        def createSqlite(prefix):
            nparts = self._getNumberOfParticles(it, prefix)
            sqliteFn = self.protocol._getFileName('projections', iter=it,
                                                  half=prefix)
            self.createAngDistributionSqlite(sqliteFn, nparts,
                                             itemDataIterator=self._iterAngles(
                                                 it, prefix))
            return sqliteFn

        if len(volumes) > 1:
            raise Exception(
                "Please, select a single volume to show it's angular "
                "distribution")
        elif not os.path.exists(angularDist):
            raise Exception(
                "Please, select a valid iteration to show the angular "
                "distribution")
        else:
            if self.showHalves.get() == HALF_EVEN:
                sqliteFn = createSqlite('even')
            elif self.showHalves.get() == HALF_ODD:
                sqliteFn = createSqlite('odd')
            elif self.showHalves.get() == FULL_MAP:
                sqliteFn = createSqlite('full')
            else:
                raise Exception(
                    "Please, select a single volume to show it's angular "
                    "distribution")

            view = ChimeraClientView(volumes[0], showProjection=True,
                                     angularDistFile=sqliteFn,
                                     spheresDistance=radius)
        return view

    def _createAngDist2D(self, it):
        nrefs = self._getNumberOfRefs()
        gridsize = self._getGridSize(nrefs)
        self.protocol._execEmanProcess(self.protocol._getRun(), it)
        angularDist = self.protocol._getFileName("angles", iter=it)

        if os.path.exists(angularDist):
            xplotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                                 mainTitle="Iteration %d" % it,
                                 windowTitle="Angular distribution")

            def plot(prefix):
                nparts = self._getNumberOfParticles(it, prefix)
                title = '%s particles' % prefix
                sqliteFn = self.protocol._getFileName('projections', iter=it,
                                                      half=prefix)
                self.createAngDistributionSqlite(sqliteFn, nparts,
                                                 itemDataIterator=self._iterAngles(
                                                     it, prefix))
                xplotter.plotAngularDistributionFromMd(sqliteFn, title)

            if self.showHalves.get() == HALF_EVEN:
                plot('even')
            elif self.showHalves.get() == HALF_ODD:
                plot('odd')
            elif self.showHalves.get() == FULL_MAP:
                plot('full')
            else:
                for prefix in ['even', 'odd', 'full']:
                    plot(prefix)
            return xplotter
        else:
            return

    # ===============================================================================
    # plotFSC
    # ===============================================================================
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        gridsize = self._getGridSize(1)
        xplotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                             windowTitle='Resolution FSC')

        plot_title = 'FSC'
        a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'FSC',
                                   yformat=False)
        legends = []

        show = False
        for it in self._iterations:
            if self.resolutionPlotsFSC.get() == FSC_UNMASK:
                fscUnmask = self.protocol._getFileName('fscUnmasked',
                                                       run=self.protocol._getRun(),
                                                       iter=it)
                if os.path.exists(fscUnmask):
                    show = True
                    self._plotFSC(a, fscUnmask)
                    legends.append('unmasked map it %d' % it)
                xplotter.showLegend(legends)

            elif self.resolutionPlotsFSC.get() == FSC_MASK:
                fscMask = self.protocol._getFileName('fscMasked',
                                                     run=self.protocol._getRun(),
                                                     iter=it)
                if os.path.exists(fscMask):
                    show = True
                    self._plotFSC(a, fscMask)
                    legends.append('masked map it %d' % it)
                xplotter.showLegend(legends)

            elif self.resolutionPlotsFSC.get() == FSC_MASKTIGHT:
                fscMaskTight = self.protocol._getFileName('fscMaskedTight',
                                                          run=self.protocol._getRun(),
                                                          iter=it)
                if os.path.exists(fscMaskTight):
                    show = True
                    self._plotFSC(a, fscMaskTight)
                    legends.append('masked tight map it %d' % it)
                xplotter.showLegend(legends)
            elif self.resolutionPlotsFSC.get() == FSC_ALL:
                fscUnmask = self.protocol._getFileName('fscUnmasked',
                                                       run=self.protocol._getRun(),
                                                       iter=it)
                fscMask = self.protocol._getFileName('fscMasked',
                                                     run=self.protocol._getRun(),
                                                     iter=it)
                fscMaskTight = self.protocol._getFileName('fscMaskedTight',
                                                          run=self.protocol._getRun(),
                                                          iter=it)
                if os.path.exists(fscUnmask):
                    show = True
                    self._plotFSC(a, fscUnmask)
                    legends.append('unmasked map it %d' % it)
                    self._plotFSC(a, fscMask)
                    legends.append('masked map it %d' % it)
                    self._plotFSC(a, fscMaskTight)
                    legends.append('masked tight map it %d' % it)
                xplotter.showLegend(legends)

        if show:
            if threshold < self.maxFrc:
                a.plot([self.minInv, self.maxInv], [threshold, threshold],
                       color='black', linestyle='--')
            a.grid(True)
        else:
            raise Exception("Set a valid iteration to show its FSC")

        return [xplotter]

    def _plotFSC(self, a, fscFn):
        resolution_inv = self._getColunmFromFilePar(fscFn, 0)
        frc = self._getColunmFromFilePar(fscFn, 1)
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])

    def _showHtmlReport(self, paramName=None):
        reportPath = self.protocol._getFileName('reportHtml',
                                                run=self.protocol._getRun())
        if pwutils.exists(reportPath):
            text._open_cmd(reportPath, self.getTkRoot())
        else:
            self.showInfo('Your html report is not ready yet. Please try again in a minute.')

    # ===============================================================================
    # Utils Functions
    # ===============================================================================
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
            self._iterations = self._getListFromRangeString(
                self.iterSelection.get())

        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq)

    @staticmethod
    def _formatFreq(value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1 / value
        return "1/%0.2f" % inv

    def _getVolumeNames(self):
        vols = []
        runType = self.protocol._getRun()
        for it in self._iterations:
            if self.showHalves.get() == HALF_EVEN:
                volFn = self.protocol._getFileName('mapEven', run=runType,
                                                   iter=it)
                vols.append(volFn)
            elif self.showHalves.get() == HALF_ODD:
                volFn = self.protocol._getFileName('mapOdd', run=runType,
                                                   iter=it)
                vols.append(volFn)
            elif self.showHalves.get() == FULL_MAP:
                volFn = self.protocol._getFileName('mapFull', run=runType,
                                                   iter=it)
                vols.append(volFn)
            else:
                volFn = self.protocol._getFileName('mapEven', run=runType,
                                                   iter=it)
                vols.append(volFn)
                volFn = self.protocol._getFileName('mapOdd', run=runType,
                                                   iter=it)
                vols.append(volFn)
                volFn = self.protocol._getFileName('mapFull', run=runType,
                                                   iter=it)
                vols.append(volFn)
        return vols

    def _getNumberOfRefs(self):
        if self.showHalves.get() == ALL_MAPS:
            refs = 3
        else:
            refs = 1

        return refs

    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of
        references. """
        if n is None:
            n = len(self._refsList)

        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n + 1) / 2, 2]

        return gridsize

    def _getColunmFromFilePar(self, fscFn, col):
        f1 = open(fscFn)
        value = []
        for l in f1:
            valList = l.split()
            val = float(valList[col])
            value.append(val)
        f1.close()
        return value

    def _getNumberOfParticles(self, it, prefix='full'):
        f = open(self.protocol._getFileName('angles', iter=it))
        nLines = int(f.readlines()[-1].split()[0]) + 1

        if prefix == 'full':
            return nLines
        else:
            return nLines / 2

    def _iterAngles(self, it, half="full"):
        f = open(self.protocol._getFileName('angles', iter=it))
        rest = 0 if half == 'even' else 1

        for i, line in enumerate(f):
            if '#' not in line:
                angles = map(float, line.split())
                if angles[1] != 0:  # skip disabled images
                    rot = float("{0:.2f}".format(angles[2]))
                    tilt = float("{0:.2f}".format(angles[3]))
                    if half == 'full' or i % 2 == rest:
                        yield rot, tilt

        f.close()


class TiltValidateViewer(ProtocolViewer):
    """ Visualization of Eman2 tilt validate protocol."""

    _targets = [EmanProtTiltValidate]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer tilt validate'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('radcut', FloatParam,
                      default=45.0,
                      label='Max radius (deg.)',
                      help='Truncate the polar plt at this radius value. '
                           '-1 means no limit.')
        form.addParam('planethres', FloatParam,
                      default=360.0,
                      label='Max out of plane threshold (deg.)',
                      help='Maximum out of plane threshold for the tilt axis. '
                           '0 = perfectly in plane, 360 = normal to plane.')
        form.addParam('displayPlot', EnumParam,
                      choices=['scatter plot', 'contour plot'],
                      default=TILT_SCATTER, display=EnumParam.DISPLAY_HLIST,
                      label='Display tilt geometry plot',
                      help='*scatter plot*: display polar plot of computed '
                           ' tilt geometry per particle pair.\n'
                           '*contour plot*: display contour plot similar '
                           'to fig.6 in Henderson paper.')
        form.addParam('displayEmanPlot', LabelParam,
                      label='Display scatter plot with EMAN2 GUI')
        form.addParam('colozaxis', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition='displayPlot == 0',  # scatter plot
                      label='Color Z-axis',
                      help='Color scatter dots by Z axis')

    def _getVisualizeDict(self):
        self._load()
        return {'displayPlot': self._showPlot,
                'displayEmanPlot': self._showEmanPlot}

    def _showPlot(self, paramName=None):
        views = []
        color = self.colozaxis
        rmax = self.radcut.get()

        if self.displayPlot == TILT_SCATTER:
            views.append(self._createScatterPlot(rmax, colorzaxis=color))
        elif self.displayPlot == TILT_CONTOUR:
            plotFn = self.protocol._getFileName('outputContourPlot')
            if pwutils.exists(plotFn):
                views.append(DataView(plotFn))
            else:
                raise Exception("Contour plot file not found: %s" % plotFn)

        return views

    def _createScatterPlot(self, rmax, colorzaxis=False):
        gridsize = self._getGridSize(1)
        xplotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                             windowTitle='Tilt geometry plot')
        plot_title = 'Tilt pair parameter plot'
        a = xplotter.createSubPlot(plot_title, 'Tilt axis', 'Tilt angle',
                                   projection='polar')

        datap, r, theta, zaxis = self._getValues()

        if colorzaxis:
            a.scatter(theta, r, c=zaxis)
        else:
            a.scatter(theta, r)
        a.set_rmax(rmax)

        return xplotter

    def _showEmanPlot(self, paramName=None):
        program = getEmanProgram('e2tiltvalidate.py')
        args = "--path=TiltValidate_01 --radcut=%0.2f --gui --planethres=%0.2f" % (
            self.radcut.get(), self.planethres.get())
        if self.colozaxis:
            args += " --colorzaxis"

        hostConfig = self.protocol.getHostConfig()
        # Create the steps executor
        executor = StepExecutor(hostConfig)
        self.protocol.setStepsExecutor(executor)
        # Finally run the protocol

        self.protocol.runJob(program, args, cwd=self.protocol._getExtraPath(),
                             numberOfMpi=1, numberOfThreads=1)

        return []

    # ===============================================================================
    # Utils Functions
    # ===============================================================================
    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self.protocol._createFilenameTemplates()

    def _getValues(self):
        fileName = self.protocol._getFileName('outputAngles')
        planethres = self.planethres.get()

        r = []
        theta = []
        datap = []
        zaxis = []

        if pwutils.exists(fileName):
            jsonPosDict = loadJson(fileName)
            if jsonPosDict.has_key("particletilt_list"):
                tiltpairs = jsonPosDict["particletilt_list"]
                maxcolorval = max(tiltpairs, key=lambda x: x[3])[3]

                for tp in tiltpairs:
                    if tp[3] > planethres:
                        continue
                    datap.append(tp[0])
                    r.append(tp[1])
                    theta.append(math.radians(tp[2]))
                    # Color the Z axis out of planeness
                    zaxis.append(self._computeRGBcolor(tp[3], 0, maxcolorval))

                return datap, r, theta, zaxis

    def _computeRGBcolor(self, value, minval, maxval):
        """ From e2tiltvalidate.py:
    Author: John Flanagan (jfflanag@bcm.edu)
    Copyright (c) 2000-2011 Baylor College of Medicine
    Modified by Stephen Murray (scmurray@bcm.edu) 3/28/13
    Compute a RGB value to represent a data range.
    Basically convert Hue to GSB with I=0.33 and S=1.0
   """
        # Normalize from 0 to 1
        normval = (value - minval) / (maxval - minval)
        radval = normval * 2 * math.pi
        if radval < 2 * math.pi / 3:
            B = 0.0
            R = 0.33 * (1 + math.cos(radval) / math.cos(math.pi / 3 - radval))
            G = 1.0 - R
            return "#%02x%02x%02x" % (255 * R, 255 * G, 255 * B)
        if 2 * math.pi / 3 < radval < 4 * math.pi / 3:
            hue = radval - 2 * math.pi / 3
            R = 0.0
            G = 0.33 * (1 + math.cos(hue) / math.cos(math.pi / 3 - hue))
            B = 1.0 - G
            return "#%02x%02x%02x" % (255 * R, 255 * G, 255 * B)
        if radval > 4 * math.pi / 3:
            hue = radval - 4 * math.pi / 3
            G = 0
            B = 0.33 * (1 + math.cos(hue) / math.cos(math.pi / 3 - hue))
            R = 1.0 - B
            return "#%02x%02x%02x" % (255 * R, 255 * G, 255 * B)

    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of references."""
        if n is None or n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n + 1) / 2, 2]

        return gridsize


class CtfViewer(ProtocolViewer):
    """ Visualization of Eman2 ctf auto protocol."""

    _targets = [EmanProtCTFAuto]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer ctf'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('outputType', EnumParam,
                      default=0,
                      choices=self._getOutputs(),
                      label='Select output set',
                      help='Choose output particle set to display')
        form.addParam('displayCtf', LabelParam,
                      label='Display particle set with Scipion')
        form.addParam('displayEmanCtf', LabelParam,
                      label='Display all results with EMAN2 GUI')

    def _getVisualizeDict(self):
        self._load()
        return {'displayCtf': self._showCtf,
                'displayEmanCtf': self._showEmanCtf}

    def _showCtf(self, paramName=None):
        views = []
        obj = "obj = self.protocol." + self.getEnumText('outputType')
        exec (obj)
        strId = obj.strId()
        fn = obj.getFileName()
        particlesView = ObjectView(self._project, strId, fn)
        views.append(particlesView)
        return views

    def _showEmanCtf(self, paramName=None):
        program = getEmanProgram('e2ctf.py')
        args = '--allparticles --minptcl=0 --minqual=0'
        args += ' --gui --constbfactor=-1.0 --sf=auto'

        hostConfig = self.protocol.getHostConfig()
        # Create the steps executor
        executor = StepExecutor(hostConfig)
        self.protocol.setStepsExecutor(executor)
        # Finally run the protocol

        self.protocol.runJob(program, args, cwd=self.protocol._getExtraPath(),
                             numberOfMpi=1, numberOfThreads=1)

        # Open dialog to request confirmation to overwrite output
        saveChanges = askYesNo("Save output changes?",
                               "Do you want to overwrite output particles with new CTF values?",
                               self.getTkRoot())
        if saveChanges:
            self.protocol.createOutputStep()

    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self.protocol._createFilenameTemplates()

    def _getOutputs(self):
        outputList = []
        for attrName, _ in self.protocol.iterOutputEM():
            outputList.append(attrName)
        return outputList
