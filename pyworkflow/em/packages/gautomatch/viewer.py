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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module implements the viewer for Gautomatch program
"""

import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import *
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import CoordinatesObjectView, ObjectView
from pyworkflow.em.data import SetOfCoordinates
from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates, writeSetOfMicrographs
import pyworkflow.em.showj as showj
from protocol_gautomatch import ProtGautomatch


class GautomatchViewer(ProtocolViewer):
    """ Visualization of Gautomatch protocol. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtGautomatch]
    _label = 'viewer Gautomatch'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowAutopick', LabelParam,
                      label="Show auto-picked particles?", default=True)
        form.addParam('doShowRejected', LabelParam,
                      label="Show rejected particles?", default=True,
                      help="Rejected particles may be useful for adjusting "
                           "picking parameters.")

        form.addSection(label='Debug output')
        form.addParam('doShowCC', LabelParam,
                      label="Show cross-correlation files?", default=True)
        form.addParam('doShowFilt', LabelParam,
                      label="Show pre-filtered micrographs?", default=True,
                      condition='_doPrefilt')
        form.addParam('doShowBgEst', LabelParam,
                      label="Show background estimation?", default=True)
        form.addParam('doShowBgSub', LabelParam,
                      label="Show background-subtracted micrographs?", default=True)
        form.addParam('doShowSigma', LabelParam,
                      label="Show local sigma micrographs?", default=True)
        form.addParam('doShowMask', LabelParam,
                      label="Show auto-detected mask?", default=True)

    def _doPrefilt(self):
        return True if self.protocol.preFilt.get() else False

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
        project = self.protocol.getProject()
        micSet = self.protocol.getInputMicrographs()
        coordsDir = project.getTmpPath(micSet.getName())

        if micSet is None:
            raise Exception('visualize: SetOfCoordinates has no micrographs set.')

        micfn = project.getTmpPath(micSet.getName() + '_micrographs.xmd')
        writeSetOfMicrographs(micSet, micfn)
        labels = 'enabled id _size _filename'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, showj.RENDER: '_filename',
                      'labels': 'id',
                      }
        fn = ''
        if param == 'doShowAutopick':
            self._convertCoords(micSet, coordsDir, coordsType='autopick')
            view = CoordinatesObjectView(project, micfn, coordsDir, self.protocol)
            return [view]
        elif param == 'doShowRejected':
            self._convertCoords(micSet, coordsDir, coordsType='rejected')
            view = CoordinatesObjectView(project, micfn, coordsDir, self.protocol)
            return [view]
        elif param == 'doShowCC':
            fn = self.protocol.createDebugOutput(suffix='_ccmax')
        elif param == 'doShowFilt':
            fn = self.protocol.createDebugOutput(suffix='_pref')
        elif param == 'doShowBgEst':
            fn = self.protocol.createDebugOutput(suffix='_bg')
        elif param == 'doShowBgSub':
            fn = self.protocol.createDebugOutput(suffix='_bgfree')
        elif param == 'doShowSigma':
            fn = self.protocol.createDebugOutput(suffix='_lsigma')
        elif param == 'doShowMask':
            fn = self.protocol.createDebugOutput(suffix='_mask')

        view = ObjectView(project, self.protocol.strId(), fn, viewParams=viewParams)

        return [view]

    def _convertCoords(self, micSet, coordsDir, coordsType):
        """ Link specified coord set to Tmp folder and convert it to .pos files"""
        pwutils.cleanPath(coordsDir)
        pwutils.makePath(coordsDir)
        coordTypes = {'autopick': 'coordinates.sqlite',
                      'rejected': 'coordinates_rejected.sqlite'
                      }
        coordsFnIn = self.protocol._getPath(coordTypes[coordsType])
        coordsFnOut = pwutils.join(coordsDir, 'coordinates.sqlite')
        pwutils.createLink(coordsFnIn, coordsFnOut)
        coordSet = SetOfCoordinates(filename=coordsFnOut)
        coordSet.setMicrographs(micSet)
        writeSetOfCoordinates(coordsDir, coordSet, ismanual=False)
