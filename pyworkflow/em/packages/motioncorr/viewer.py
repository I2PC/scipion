# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import MicrographsView
import pyworkflow.em.showj as showj

from protocol_motioncorr import ProtMotionCorr


class ProtMotioncorrViewer(Viewer):
    _targets = [ProtMotionCorr]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer motioncorr'

    def _visualize(self, obj, **kwargs):
        views = []

        labelsDef = 'enabled id _filename'
        viewParamsDef = {showj.MODE: showj.MODE_MD,
                         showj.ORDER: labelsDef,
                         showj.VISIBLE: labelsDef,
                         showj.RENDER: None
                         }

        outputLabels = ['outputMicrographs', 'outputMicrographsDoseWeighted',
                        'outputMovies']

        if not any(obj.hasAttribute(l) for l in outputLabels):
            return [self.infoMessage("Output (micrographs or movies) have "
                                     "not been produced yet.")]

        # Display only the first available output, showing all of them
        # can be confusing and not useful.
        # The user can still double-click in the specific output
        for l in outputLabels:
            if obj.hasAttribute(l):
                output = getattr(obj, l)
                if 'Micrographs' in l:
                    return [MicrographsView(self.getProject(), output)]
                else:  # Movies case
                    return [self.objectView(output, viewParams=viewParamsDef)]

        return views
