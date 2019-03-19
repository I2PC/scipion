# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.convert import ImageHandler

import views
import showj


class DataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]
    _targets = [
        em.Image,
        em.SetOfClasses2D,
        em.SetOfClasses3D,
        em.SetOfCoordinates,
        em.SetOfCTF,
        em.SetOfImages,
        em.SetOfMovies,
        em.SetOfNormalModes,
        em.SetOfPDBs,
        em.ProtParticlePicking,
        em.ProtImportMovies,
        # TiltPairs related data
        em.CoordinatesTiltPair,
        em.MicrographsTiltPair,
        em.ParticlesTiltPair
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _addObjView(self, obj, fn, viewParams={}):
        objView = views.ObjectView(self._project, obj.strId(), fn,
                                   viewParams=viewParams)
        self._views.append(objView)
        return objView

    def _visualize(self, obj, **kwargs):
        cls = type(obj)

        if issubclass(cls, em.Volume):
            fn = ImageHandler.locationToXmipp(obj)
            self._addObjView(obj, fn,
                             {showj.RENDER: 'image',
                              showj.SAMPLINGRATE: obj.getSamplingRate()})

        elif issubclass(cls, em.Image):
            fn = ImageHandler.locationToXmipp(obj)
            self._addObjView(obj, fn)

        elif issubclass(cls, em.SetOfPDBs):
            fn = obj.getFileName()
            labels = 'id _filename '
            self._addObjView(obj, fn, {showj.ORDER: labels,
                                       showj.VISIBLE: labels,
                                       showj.MODE: showj.MODE_MD,
                                       showj.RENDER: "no"})

        elif issubclass(cls, em.SetOfMovies):
            fn = obj.getFileName()
            # Enabled for the future has to be available
            labels = ('id _filename _samplingRate _acquisition._dosePerFrame '
                      '_acquisition._doseInitial ')
            moviesView = self._addObjView(obj, fn, {showj.ORDER: labels,
                                                    showj.VISIBLE: labels,
                                                    showj.MODE: showj.MODE_MD,
                                                    showj.RENDER: "no"})
            # For movies increase the JVM memory by 1 GB, just in case
            moviesView.setMemory(showj.getJvmMaxMemory() + 1)

        elif issubclass(cls, em.SetOfMicrographs):
            self._views.append(views.MicrographsView(self._project, obj, **kwargs))

        elif issubclass(cls, em.MicrographsTiltPair):
            labels = 'id enabled _untilted._filename _tilted._filename'
            renderLabels = '_untilted._filename _tilted._filename'
            self._addObjView(obj, obj.getFileName(), {showj.ORDER: labels,
                                                      showj.VISIBLE: labels,
                                                      showj.MODE: showj.MODE_MD,
                                                      showj.RENDER: renderLabels})

        elif issubclass(cls, em.ParticlesTiltPair):
            labels = 'id enabled _untilted._filename _tilted._filename'
            renderLabels = '_untilted._filename _tilted._filename'
            self._addObjView(obj, obj.getFileName(), {showj.ORDER: labels,
                                                      showj.VISIBLE: labels,
                                                      showj.RENDER: renderLabels,
                                                      showj.MODE: showj.MODE_MD})

        elif issubclass(cls, em.SetOfCoordinates):
            # FIXME: Remove dependency on xmipp3 plugin to visualize coordinates
            xmipp3 = pwutils.importFromPlugin('xmipp3',
                                              errorMsg="xmipp3 plugin is required "
                                                       "now to visualize coordinates.")
            micSet = obj.getMicrographs()  # accessing mics to provide metadata file
            if micSet is None:
                raise Exception('visualize: SetOfCoordinates has no micrographs set.')

            mdFn = getattr(micSet, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:  # happens if protocol is not an xmipp one
                fn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                xmipp3.convert.writeSetOfMicrographs(micSet, fn)
            tmpDir = self._getTmpPath(obj.getName())
            pwutils.cleanPath(tmpDir)
            pwutils.makePath(tmpDir)
            # FIXME: (JMRT) We are always writing the SetOfCoordinates and removing
            # the tmpDir, we need to take into account if the user have pick
            # some particles in the tmpDir and have not save them, that now
            # will loose all picked particles.
            # A possible solution could be to alert that changes have not been
            # written during modification of tmpDir or create a new Xmipp picking
            # protocol to continue picking later without loosing the coordinates.
            xmipp3.convert.writeSetOfCoordinates(tmpDir, obj)
            self._views.append(views.CoordinatesObjectView(self._project, fn,
                                                           tmpDir, self.protocol,
                                                           inTmpFolder=True))

        elif issubclass(cls, em.SetOfParticles):
            fn = obj.getFileName()
            labels = ('id enabled _index _filename _xmipp_zScore _xmipp_cumulativeSSNR '
                      '_sampling _xmipp_scoreByVariance _xmipp_scoreEmptiness '
                      '_ctfModel._defocusU _ctfModel._defocusV _ctfModel._defocusAngle '
                      '_transform._matrix')
            self._addObjView(obj, fn, {showj.ORDER: labels,
                                       showj.VISIBLE: labels,
                                       showj.SORT_BY: '_xmipp_zScore asc',
                                       showj.RENDER: '_filename'})

        elif issubclass(cls, em.SetOfVolumes):
            fn = obj.getFileName()
            labels = 'id enabled comment _filename '
            self._addObjView(obj, fn, {showj.MODE: showj.MODE_MD,
                                       showj.ORDER: labels,
                                       showj.VISIBLE: labels,
                                       showj.RENDER: '_filename'})

        elif issubclass(cls, em.SetOfClasses2D):
            self._views.append(views.ClassesView(self._project, obj.strId(),
                                                 obj.getFileName(), **kwargs))

        elif issubclass(cls, em.SetOfClasses3D):
            self._views.append(views.Classes3DView(self._project, obj.strId(),
                                                   obj.getFileName()))

        elif issubclass(cls, em.SetOfImages):
            self._views.append(views.ObjectView(self._project, obj.strId(),
                                                obj.getFileName(), **kwargs))

        elif issubclass(cls, em.SetOfCTF):
            self._views.append(views.CtfView(self._project, obj))

        elif issubclass(cls, em.CoordinatesTiltPair):
            # FIXME: Remove dependency on xmipp3 plugin to visualize coordinates
            xmipp3 = pwutils.importFromPlugin('xmipp3',
                                              errorMsg="xmipp3 plugin is required "
                                                       "now to visualize coordinates.")
            tmpDir = self._getTmpPath(obj.getName())
            pwutils.makePath(tmpDir)

            mdFn = os.path.join(tmpDir, 'input_micrographs.xmd')
            xmipp3.convert.writeSetOfMicrographsPairs(
                obj.getUntilted().getMicrographs(),
                obj.getTilted().getMicrographs(), mdFn)
            parentProtId = obj.getObjParentId()
            parentProt = self.getProject().mapper.selectById(parentProtId)
            extraDir = parentProt._getExtraPath()

            #           extraDir = parentProt._getExtraPath()
            # TODO: Review this if ever a non Xmipp CoordinatesTiltPair is available
            xmipp3.convert.writeSetOfCoordinates(tmpDir, obj.getUntilted())
            xmipp3.convert.writeSetOfCoordinates(tmpDir, obj.getTilted())
            showj.launchTiltPairPickerGUI(mdFn, tmpDir, self.protocol)

        elif issubclass(cls, em.ProtParticlePicking):
            if obj.getOutputsSize() >= 1:
                self._visualize(obj.getCoords())

        elif issubclass(cls, em.ProtImportMovies):
            movs = obj.outputMovies
            self._visualize(movs)
            gainFn = movs.getGain()
            if gainFn is not None:
                if os.path.exists(gainFn):
                    self._views.append(views.DataView(gainFn))

        return self._views

