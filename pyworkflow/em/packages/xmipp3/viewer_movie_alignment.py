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
import pyworkflow.em.showj as showj

from protocol_movie_opticalflow import (XmippProtOFAlignment,
                                        OBJCMD_MOVIE_ALIGNCARTESIAN)


class XmippMovieAlignViewer(Viewer):
    _targets = [XmippProtOFAlignment]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    _label = 'viewer optical alignment'

    def _visualize(self, obj, **kwargs):
        views = []

        plotLabels = ('psdCorr._filename plotPolar._filename '
                      'plotCart._filename')
        labels = plotLabels + ' _filename '
        viewParams = {showj.MODE: showj.MODE_MD,
                      showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.RENDER: plotLabels,
                      showj.ZOOM: 50,
                      showj.OBJCMDS: "'%s'" % OBJCMD_MOVIE_ALIGNCARTESIAN
                      }

        if obj.hasAttribute('outputMicrographs'):
            views.append(self.objectView(obj.outputMicrographs,
                                         viewParams=viewParams))
        elif obj.hasAttribute('outputMovies'):
            views.append(self.objectView(obj.outputMovies,
                                         viewParams=viewParams))
        else:
            views.append(self.infoMessage("Output (micrographs or movies) has "
                                          "not been produced yet."))

        return views

