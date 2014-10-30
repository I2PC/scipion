# **************************************************************************
# *
# * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
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
This module implement the wrappers aroung Xmipp ML2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.protocol.params import BooleanParam, IntParam, EnumParam
from pyworkflow.gui.plotter import Plotter
from protocol_movie_alignment import ProtMovieAlignment
import matplotlib.pyplot as plt
from viewer import XmippViewer
import numpy as np

PLOT_CART = 0
PLOT_POLAR = 1
PLOT_POLARCART = 2
PLOT_CHOICES = ['cartesian', 'polar', 'both']

class XmippMovieAlignViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [ProtMovieAlignment]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer movie alignment'
    #_plotVars = ['doShowLL', 'doShowPmax', 'doShowSignalChange', 'doShowMirror']

    def _defineParams(self, form):
        form.addSection('Visualization')
        form.addParam('showMovies', BooleanParam,
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
        return {'showMovies': self._showMovies,
                'plotToShow': self._createPlots}

    def _showMovies(self, paramName):
        return XmippViewer(project=self._project)._visualize(self.protocol.outputMovies)

    def _createPlots(self, paramName):

        import xmipp
        meanX = []
        meanY = []
        stdX = []
        stdY = []
        colors = []
        gr = 255.0
        plotType = self.plotToShow.get()
        movieId = self.movieId.get()
        movie = self.protocol.outputMovies[movieId]

        if movie is None:
            return [self.errorMessage("Invalid movie id *%d*" % movieId,
                                      title="Invalid input")]

        alignedMovie = self.protocol.outputMovies[movieId].alignMetaData
        md = xmipp.MetaData(alignedMovie)
        colorDist = 255 / md.size()

        cartPosition = None
        polarPosition = None
        colorBarPosition = 122
        figureSize = (8, 6)

        if plotType == PLOT_CART:
            cartPosition = 121
        elif plotType == PLOT_POLAR:
            polarPosition = 121
        elif plotType == PLOT_POLARCART:
            cartPosition = 132
            polarPosition = 131
            colorBarPosition = 133

        plotter = Plotter(*figureSize)
        figure = plotter.getFigure()

        # Plot the color bar
        ax = figure.add_subplot(colorBarPosition, aspect='equal', xlim=[0, 6])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        colorBarX = np.array([1, 1])
        colorBarY = np.array([2, 4])

        # Read the shifts information from the Metadata
        for objId in md:
            meanX.append(md.getValue(xmipp.MDL_OPTICALFLOW_MEANX, objId))
            meanY.append(md.getValue(xmipp.MDL_OPTICALFLOW_MEANY, objId))
            stdX.append(md.getValue(xmipp.MDL_OPTICALFLOW_STDY, objId))
            stdY.append(md.getValue(xmipp.MDL_OPTICALFLOW_STDY, objId))
            colors.append((1, gr / 255.0, 0))
            ax.plot(colorBarX, colorBarY, c=(1, gr / 255.0, 0), linewidth=10)
            ax.text(2, np.mean(colorBarY), str(objId)+'-'+str(objId+1))
            colorBarY += 2
            gr -= colorDist
        area = (np.sqrt(np.power(np.asarray(stdX), 2) + np.power(np.asarray(stdY), 2)))*700

        # Plot in polar if needed
        if polarPosition:
            r = np.sqrt(np.power(np.asarray(meanX), 2) + np.power(np.asarray(meanY), 2))
            theta = np.arctan2(meanY, meanX) * 180 / np.pi
            ax = figure.add_subplot(polarPosition, projection='polar')
            ax.set_title('Polar representation')
            c = ax.scatter(theta, r, c=colors, s=area, cmap=plt.cm.hsv)
            c.set_alpha(0.75)
            ax.plot(theta, r, '-^')

        # Plot in cartesian if needed
        if cartPosition:
            ax = figure.add_subplot(cartPosition)
            ax.grid()
            ax.set_title('Cartesian representation')
            c = ax.scatter(np.asarray(meanX), np.asarray(meanY), c=colors, s=area, cmap=plt.cm.hsv)
            c.set_alpha(0.75)
            ax.plot(np.asarray(meanX), np.asarray(meanY), '-^')

        #plt.tight_layout()
        #plt.show()
        return [plotter]