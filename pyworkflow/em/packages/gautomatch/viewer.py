# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This module implements the viewer for Gautomatch program
"""

from pyworkflow.protocol.params import *
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.packages.xmipp3.convert import *
from pyworkflow.em.packages.xmipp3.viewer import launchSupervisedPickerGUI
from protocol_gautomatch import ProtGautomatch


class GautomatchViewer(ProtocolViewer):
    """ Visualization of Gautomatch protocol. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtGautomatch]
    _label = 'viewer Gautomatch'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowAutopick', LabelParam,
                      label="Show auto-picked particles?", default=True,
                      help="Show the auto-picked particles")
        form.addParam('doShowRejected', LabelParam,
                      label="Show rejected particles?", default=True,
                      help="Rejected particles may be useful for adjusting "
                           "picking parameters.")

        form.addSection(label='Debug output')
        form.addParam('doShowCC', LabelParam,
                      label="Show cross-correlation files?", default=True,
                      condition='_writeCC')
        form.addParam('doShowFilt', LabelParam,
                      label="Show pre-filtered micrographs?", default=True,
                      condition='_writeFilt')
        form.addParam('doShowBgEst', LabelParam,
                      label="Show background estimation?", default=True,
                      condition='_writeBg')
        form.addParam('doShowBgSub', LabelParam,
                      label="Show background-subtracted micrographs?", default=True,
                      condition='_writeBgSub')
        form.addParam('doShowSigma', LabelParam,
                      label="Show local sigma micrographs?", default=True,
                      condition='_writeSigma')
        form.addParam('doShowMask', LabelParam,
                      label="Show auto-detected mask?", default=True,
                      condition='_writeMsk')

    def _writeCC(self):
        return True if self.protocol.writeCC.get() else False

    def _writeFilt(self):
        return True if self.protocol.writeFilt.get() else False

    def _writeBg(self):
        return True if self.protocol.writeBg.get() else False

    def _writeBgSub(self):
        return True if self.protocol.writeBgSub.get() else False

    def _writeSigma(self):
        return True if self.protocol.writeSigma.get() else False

    def _writeMsk(self):
        return True if self.protocol.writeMsk.get() else False

    def _getVisualizeDict(self):
        return {'doShowAutopick': self._viewParam,
                'doShowRejected': self._viewParam,
                'doShowCC': self._viewParam,
                'doShowFilt': self._viewParam,
                'doShowBgEst': self._viewParam,
                'doShowBgSub': self._viewParam,
                'doShowSigma': self._viewParam,
                'doShowMask': self._viewParam
                }

    def _viewParam(self, param=None):
        micSet = self.protocol.getInputMicrographs()
        tmpDir = self.protocol._getTmpPath()
        pwutils.cleanPath(tmpDir)
        pwutils.makePath(tmpDir)
        # FIXME: (JMRT) We are always writing the SetOfCoordinates and removing
        # the tmpDir, we need to take into account if the user has picked
        # some particles in the tmpDir and has not saved them, that now he
        # will lose all picked particles.
        # A possible solution could be to alert that changes have not been
        # written during modification of tmpDir or create a new Xmipp picking
        # protocol to continue picking later without losing the coordinates.

        if micSet is None:
            raise Exception('visualize: SetOfCoordinates has no micrographs set.')

        micsFn = pwutils.join(tmpDir, micSet.getName() + '_micrographs.xmd')
        writeSetOfMicrographs(micSet, micsFn)
        inTmpFolder = True
        view = []

        if param == 'doShowAutopick':
            self._convertCoords(micSet, tmpDir, coordsType='autopick')
            launchSupervisedPickerGUI(micsFn, tmpDir, self.protocol, mode='review',
                                      inTmpFolder=inTmpFolder)
        elif param == 'doShowRejected':
            self._convertCoords(micSet, tmpDir, coordsType='rejected')
            launchSupervisedPickerGUI(micsFn, tmpDir, self.protocol, mode='review',
                                      inTmpFolder=inTmpFolder)
        elif param == 'doShowCC':
            fn = self.protocol._getPath('micrographs_ccmax.sqlite')
            view.append(ObjectView(self._project, self.protocol.strId(), fn))
            return view
        elif param == 'doShowFilt':
            fn = self.protocol._getPath('micrographs_pref.sqlite')
            view.append(ObjectView(self._project, self.protocol.strId(), fn))
            return view
        elif param == 'doShowBgEst':
            fn = self.protocol._getPath('micrographs_bg.sqlite')
            view.append(ObjectView(self._project, self.protocol.strId(), fn))
            return view
        elif param == 'doShowBgSub':
            fn = self.protocol._getPath('micrographs_bgfree.sqlite')
            view.append(ObjectView(self._project, self.protocol.strId(), fn))
            return view
        elif param == 'doShowSigma':
            fn = self.protocol._getPath('micrographs_lsigma.sqlite')
            view.append(ObjectView(self._project, self.protocol.strId(), fn))
            return view
        elif param == 'doShowMask':
            fn = self.protocol._getPath('micrographs_mask.sqlite')
            view.append(ObjectView(self._project, self.protocol.strId(), fn))
            return view

    def _convertCoords(self, micSet, tmpDir, coordsType):
        """ Link specified coord set to tmpDir folder and convert it to .pos files"""
        coordTypes = {'autopick': 'coordinates.sqlite',
                      'rejected': 'coordinates_rejected.sqlite'}
        coordsFnIn = self.protocol._getPath(coordTypes[coordsType])
        coordsFnOut = pwutils.join(tmpDir, 'coordinates.sqlite')
        pwutils.createLink(coordsFnIn, coordsFnOut)
        coordSet = SetOfCoordinates(filename=coordsFnOut)
        coordSet.setMicrographs(micSet)
        writeSetOfCoordinates(tmpDir, coordSet, ismanual=False)
