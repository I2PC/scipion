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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.protocol import ProtClassesConsensus

from .views import DataView


class ViewerClassesConsensus(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtClassesConsensus]

    def _visualize(self, obj, **kwargs):
        labels = ('class1.id class1._representative._filename class2.id '
                  'class2._representative._filename jaccard intersection union')
        return [DataView(obj.outputConsensus.getFileName(),
                         viewParams={'order': labels, 'mode': 'metadata',
                                     'visible': labels,
                                     'render': 'class1._representative._filename class2._representative._filename'
                                     })
                ]

    def visualize(self, obj, **kwargs):
        self._visualize(obj, **kwargs)[0].show()
