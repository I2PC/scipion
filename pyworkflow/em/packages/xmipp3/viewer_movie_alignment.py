# **************************************************************************
# *
# * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
# *              J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Airen Zaldivar  (airenzp@gmail.com)
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
This module implement the wrappers aroung Xmipp ML2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.protocol.params import BooleanParam, IntParam, EnumParam

from protocol_movie_alignment import ProtMovieAlignment

from pyworkflow.em.viewer import ObjectView
from pyworkflow.em.showj import OBJCMDS
from pyworkflow.em.packages.xmipp3.protocol_movie_alignment import createPlots

class XmippMovieAlignViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    #_targets = [ProtMovieAlignment]
    _targets = []
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer movie alignment'
    #_plotVars = ['doShowLL', 'doShowPmax', 'doShowSignalChange', 'doShowMirror']

    def _defineParams(self, form):
        form.addSection('Visualization')
        form.addParam('showMovies', BooleanParam, default=False,
                      important=True,
                      label="Display movies",
                      help='Check to see the list of available movies.')
        form.addParam('movieId', IntParam,
                      default=1,
                      label='Movie number',
                      help='The movie number for which the plots are requested.')
        form.addParam('plotToShow', EnumParam, choices=PLOT_CHOICES,
                      default=PLOT_POLAR, display=EnumParam.DISPLAY_COMBO,
                      label="Plot type",
                      help="Select which type of plot to use")

    def _getVisualizeDict(self):
        return {'showMovies': self._showMics,
                'plotToShow': self.showPlots}

    def _showMics(self, paramName):
        #return XmippViewer(project=self._project)._visualize(self.protocol.outputMicrographs)
        outputMics = self.protocol.outputMicrographs
        #objCommands = '%s' % (OBJCMD_MOVIE_ALIGNPOLAR)
        return [ObjectView(self._project, outputMics.strId(), outputMics.getFileName(), self.protocol.strId())]

    def showPlots(self, paramName):
        createPlots(self.plotToShow.get(), self.protocol, self.movieId.get())