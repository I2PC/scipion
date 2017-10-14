# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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
This module implements the viewer for Xmipp mltomo protocol
"""

import os
from os.path import exists

from pyworkflow.em import (DataView, Classes3DView, FscViewer, FSC,
                           ObjectView, ChimeraView, ChimeraClientView)
import pyworkflow.em.showj as showj
import pyworkflow.em.metadata as md
from pyworkflow.protocol.params import EnumParam, NumericRangeParam, LabelParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from protocol_mltomo import XmippProtMLTomo

ITER_LAST = 0
ITER_SELECTION = 1

CLASSES_ALL = 0
CLASSES_SEL = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1


class XmippMLTomoViewer(ProtocolViewer):
    """ Visualization of the MLTomo results. """
    _label = 'viewer mltomo'
    _targets = [XmippProtMLTomo]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        form.addParam('viewIter', EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help="""
        *last*: only the last iteration will be visualized.
        *selection*: you may specify a range of iterations.
        Examples:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]                      
                                   """)
        form.addParam('iterSelection', NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Classification')
        group.addParam('showImagesInClasses', LabelParam,
                      label='Show classification in Scipion', important=True,
                      help='Display each class with the number of volumes assigned. \n'
                           '*Note1*: The volumes of one class can be shown by \n'
                           'right-click on the class and select "Open images".\n'
                           '*Note2*: This option can take several minutes if the \n'
                           'number of volumes is high.')
        group.addParam('showLogFile', LabelParam,
                       label='Show *_log.xmd file')

        group = form.addGroup('Volumes')
        group.addParam('showClasses3D', EnumParam, default=CLASSES_ALL,
                       choices=['all', 'selection'],
                       display=EnumParam.DISPLAY_HLIST,
                       label='3D Class to visualize')
        group.addParam('class3DSelection', NumericRangeParam, default='1',
                       condition='showClasses3D == %d' % CLASSES_SEL,
                       label='Classes list')
        group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'],
                       default=VOLUME_SLICES, display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')

        group = form.addGroup('Resolution')
        group.addParam('figure', EnumParam, default=0,
                       choices=['new', 'active'],
                       label='Figure', display=EnumParam.DISPLAY_HLIST)
        group.addParam('resolutionPlotsFSC', LabelParam, default=True,
                       label='Display resolution plots (FSC)')

    def _getVisualizeDict(self):
        self._load()
        return {'showImagesInClasses': self._showImagesInClasses,
                'showLogFile': self._showLogFile,
                'displayVol': self._showVolumes,
                'resolutionPlotsFSC': self._showFSC
                }

    def _showImagesInClasses(self, paramName=None):
        views = []
        if (self.viewIter == ITER_LAST and
                    getattr(self.protocol, 'outputClasses', None) is not None):
            fn = self.protocol.outputClasses.getFileName()
            v = self.createScipionView(fn)
            views.append(v)
        elif self.viewIter == ITER_LAST:  # run not finished
            for it in self._iterations:
                fn = self.protocol._getIterClasses(it)
                v = self.createScipionView(fn)
                views.append(v)
        else:
            self._iterations = self._getRange(self.iterSelection, 'iterations')
            for it in self._iterations:
                fn = self.protocol._getIterClasses(it)
                v = self.createScipionView(fn)
                views.append(v)
        return views

    def _showLogFile(self, paramName=None):
        views = []

        for it in self._iterations:
            logFile = self.protocol._getFileName('log_it', iter=it)
            v = self.createDataView(logFile)
            views.append(v)
        return views

# =============================================================================
# ShowVolumes
# =============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('xmipp_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputVols.get().getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if exists(volFn):
                files.append(volFn)
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()

        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for volFn in volumes:
                # We assume that the chimera script will be generated
                # in the same folder as xmipp volumes
                localVol = os.path.basename(volFn)
                if exists(volFn):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = ChimeraView(cmdFile)
        else:
            view = ChimeraClientView(volumes[0])

        return [view]

# =============================================================================
# plotFSC
# =============================================================================
    def _showFSC(self, paramName=None):
        md.activateMathExtensions()
        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=0.5,
                              protocol=self.protocol,
                              figure=self._getFigure(),
                              addButton=True)
        fscSet = self.protocol._createSetOfFSCs()
        for ref3d in self._refsList:
            blockName = 'class_%06d@' % ref3d
            for it in self._iterations:
                fscFn = self.protocol._getFileName('fsc_it', iter=it)
                if exists(fscFn):
                    fsc = self._plotFSC(None, blockName + fscFn, 'iter %d, class %d' % (it, ref3d))
                    fscSet.append(fsc)
        fscViewer.visualize(fscSet)

        return [fscViewer]

    def _plotFSC(self, a, fscFn, label):
        mdStar = md.MetaData(fscFn)
        resolution_inv = [mdStar.getValue(md.MDL_RESOLUTION_FREQ, id) for id in mdStar]
        frc = [mdStar.getValue(md.MDL_RESOLUTION_FRC, id) for id in mdStar]
        fsc = FSC(objLabel=label)
        fsc.setData(resolution_inv, frc)

        return fsc

# =============================================================================
# Utils Functions
# =============================================================================

    def createDataView(self, filename, viewParams={}):
        return DataView(filename, env=self._env, viewParams=viewParams)

    def createScipionView(self, filename):
        labels = 'enabled id _size _representative._filename '
        labels += '_xmippWeight _xmippSignalChange'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.RENDER: '_representative._filename',
                      showj.SORT_BY: '_size desc',
                      showj.ZOOM: str(self._getZoom())
                      }
        inputVolumesId = self.protocol.inputVols.get().strId()
        view = Classes3DView(self._project,
                             self.protocol.strId(), filename, other=inputVolumesId,
                             env=self._env,
                             viewParams=viewParams)
        return view

    def _getZoom(self):
        # Ensure that classes are shown at least at 128 px to
        # properly see the other xmipp labels.
        dim = self.protocol.inputVols.get().getDim()[0]
        if dim < 128:
            zoom = 128*100/dim
        else:
            zoom = 100
        return zoom

    def _getRange(self, var, label):
        """ Check if the range is not empty.
        var: The variable to retrieve the value
        label: the label used for the message string
        return: the list with the range of values, empty
        """
        value = var.get()
        if value is None or not value.strip():
            self.formWindow.showError('Provide %s selection.' % label)
            result = []
        else:
            result = self._getListFromRangeString(value)
        return result

    def _getVolumeNames(self):
        vols = []
        for it in self._iterations:
            for ref3d in self._refsList:
                volFn = self.protocol._getFileName('volume', iter=it, ref3d=ref3d)
                vols.append(volFn)
        return vols

    def _load(self):
        """ Load selected iterations and classes 3D
        for visualization mode. """
        self._refsList = [1]
        if self.showClasses3D == CLASSES_ALL:
            self._refsList = range(1, self.protocol.numberOfReferences.get() + 1)
        else:
            self._refsList = self._getRange(self.class3DSelection, 'classes 3d')

        self.protocol._initialize()  # Load filename templates
        self.firstIter = self.protocol._firstIter()
        self.lastIter = self.protocol._lastIter()
        if self.viewIter.get() == ITER_LAST:
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getRange(self.iterSelection, 'iterations')

    def _getFigure(self):
        return None if self.figure == 0 else 'active'
