# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
This module implement the viewer for Xmipp mltomo protocol
"""

import glob
from os.path import exists

from pyworkflow.em import DataView
from pyworkflow.protocol.params import LabelParam, StringParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from protocol_mltomo import XmippProtMLTomo


class XmippMLTomoViewer(ProtocolViewer):
    """ Visualization of the MLTomo results. """
    _label = 'viewer mltomo'
    _targets = [XmippProtMLTomo]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowLastIter', LabelParam, label="Visualize last iteration.")
        form.addParam('showSeveralIters', StringParam, default='',
                      label='Visualize several iterations', condition='not doShowLastIter',
                      help='Create a  list of iterations like: 0,1,3 or 0-3 ')

    def _getVisualizeDict(self):
        return {'doShowLastIter': self._viewLastIter,
                'showSeveralIters': self._viewIterFiles}

    def _viewLastIter(self, e=None):
        views = []
        lastIterFile = self.protocol._getExtraPath("results/mltomo_ref.xmd")
        views.append(DataView("classes@" + lastIterFile))

        return views

    def _viewIterFiles(self, e=None):
        views = []
        errors = []
        iterFiles = glob.glob(self.protocol._getExtraPath("results/mltomo_it*_ref.xmd"))

        if iterFiles:
            iterFiles.sort()
            if self.showSeveralIters.empty():
                self.formWindow.showError('Please select the iterations that you want to visualize.')
            else:
                listOfIters = []
                try:
                    listOfIters = self._getListFromRangeString(self.showSeveralIters.get())
                except Exception, ex:
                    errors.append('Invalid iterations range.')

                for iterNum in listOfIters:
                    fn = self.protocol._getExtraPath("results/mltomo_it%06d_ref.xmd" % iterNum)
                    if exists(fn):
                        views.append(DataView("classes@" + fn))
                    else:
                        self.formWindow.showError('Iteration %s does not exist.' % iterNum)
            self.errorList(errors, views)

        return views 
