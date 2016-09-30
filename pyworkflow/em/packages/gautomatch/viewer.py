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
from pyworkflow.em.viewer import CoordinatesObjectView
from pyworkflow.em.data import SetOfCoordinates
from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates, writeSetOfMicrographs
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
        form.addParam('doShowNonunique', LabelParam,
                      label="Show non-unique particles?", default=True)
        form.addParam('doShowRejected', LabelParam,
                      label="Show rejected particles?", default=True)

        form.addSection(label='Debug output')
        form.addParam('doShowCC', LabelParam, expertLevel=LEVEL_ADVANCED,
                      label="Show cross-correlation files?", default=True)
        form.addParam('doShowBgEst', LabelParam, expertLevel=LEVEL_ADVANCED,
                      label="Show background estimation?", default=True)
        form.addParam('doShowBgSub', LabelParam, expertLevel=LEVEL_ADVANCED,
                      label="Show background-subtracted micrographs?", default=True)
        form.addParam('doShowSigma', LabelParam, expertLevel=LEVEL_ADVANCED,
                      label="Show local sigma micrographs?", default=True)
        form.addParam('doShowMask', LabelParam, expertLevel=LEVEL_ADVANCED,
                      label="Show auto-detected mask?", default=True)

    def _getVisualizeDict(self):
        return {'doShowAutopick': self._viewParam,
                'doShowNonunique': self._viewParam,
                'doShowRejected': self._viewParam,
                'doShowCC': self._viewParam,
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

        if param == 'doShowAutopick':
            self._convertCoords(micSet, coordsDir, coordsType='autopick')
            view = CoordinatesObjectView(project, micfn, coordsDir, self.protocol)
        elif param == 'doShowNonunique':
            self._convertCoords(micSet, coordsDir, coordsType='nonunique')
            view = CoordinatesObjectView(project, micfn, coordsDir, self.protocol)
        elif param == 'doShowRejected':
            self._convertCoords(micSet, coordsDir, coordsType='rejected')
            view = CoordinatesObjectView(project, micfn, coordsDir, self.protocol)
        elif param == 'doShowCC':
            pass
            #view = DataView(self.protocol._getFileName('mic_CC'))
        elif param == 'doShowBgEst':
            pass
            #view = DataView(self.protocol._getFileName('mic_BgEst'))
        elif param == 'doShowBgSub':
            pass
            #view = DataView(self.protocol._getFileName('mic_BgSub'))
        elif param == 'doShowSigma':
            pass
            #view = DataView(self.protocol._getFileName('mic_Sigma'))
        elif param == 'doShowMask':
            pass
            #view = DataView(self.protocol._getFileName('mic_Mask'))

        return [view]


    def _convertCoords(self, micSet, coordsDir, coordsType):
        """ Link specified coord set to Tmp folder and convert it to .pos files"""
        pwutils.cleanPath(coordsDir)
        pwutils.makePath(coordsDir)
        coordTypes = {'autopick': 'coordinates.sqlite',
                      'nonunique': 'coordinates_nonunique.sqlite',
                      'rejected': 'coordinates_rejected.sqlite'
                      }

        coordsFnIn = self.protocol._getPath(coordTypes[coordsType])
        coordsFnOut = pwutils.join(coordsDir, 'coordinates.sqlite')
        pwutils.createLink(coordsFnIn, coordsFnOut)
        coordSet = SetOfCoordinates(filename=coordsFnOut)
        coordSet.setMicrographs(micSet)
        writeSetOfCoordinates(coordsDir, coordSet, ismanual=False)
