# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import matplotlib.pyplot as plt

from pyworkflow.protocol.params import LabelParam, EnumParam, IntParam
from pyworkflow.utils import exists
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.em.viewer import DataView, ChimeraClientView
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from protocol_cryoef import ProtCryoEF
from convert import iterAngles

VOL_RS_PSF = 0
VOL_FS_PSF = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1


class CryoEFViewer(ProtocolViewer):
    """
    Visualization tools for cryoEF results.
    
    cryoEF is software tool for analysing the orientation
    distribution of single-particle EM data.

    Find more information at http://www.mrc-lmb.cam.ac.uk/crusso/cryoEF/
    """
           
    _environments = [DESKTOP_TKINTER]
    _targets = [ProtCryoEF]
    _label = 'viewer cryoEF'

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('Volumes')
        group.addParam('displayVol', EnumParam, choices=['slices', 'chimera'],
                       default=VOLUME_SLICES, display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('doShowOutVol', EnumParam, default=VOL_RS_PSF,
                       choices=['real space PSF', 'fourier space PSF'],
                       display=EnumParam.DISPLAY_COMBO,
                       label='PSF volume to display',
                       help='Display output volumes:\n'
                            '1) First one contains the shape of the point spread '
                            'function (PSF) corresponding to the geometry of '
                            'the orientation distribution.\n'
                            '2) Second one contains the Fourier space (k-space) '
                            'information coverage of the orientation '
                            'distribution. Ideally, it should be spherically '
                            'symmetric.')
        form.addParam('displayAngDist', LabelParam,
                      label='Display angular distribution',
                      help='Display angular distribution as '
                           'interactive 2D in matplotlib.')
        form.addParam('spheresScale', IntParam, default=100,
                      expertLevel=LEVEL_ADVANCED,
                      label='Spheres size')
        form.addParam('doShowHistogram', LabelParam,
                      label="Show histogram")
        form.addParam('doShowLog', LabelParam,
                      label="Show output log")

    def _getVisualizeDict(self):
        self.protocol._initialize()  # Load filename templates
        return {'doShowOutVol': self._showVolumes,
                'displayAngDist': self._showAngularDistribution,
                'doShowHistogram': self._showHistogram,
                'doShowLog': self._showLogFile
                }

# =============================================================================
# ShowVolumes
# =============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumeShowj()

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        view = ChimeraClientView(volumes)

        return [view]

    def _showVolumeShowj(self):
        vols = str(self._getVolumeNames())
        return [DataView(vols)]

# ===============================================================================
# showAngularDistribution
# ===============================================================================
    def _showAngularDistribution(self, paramName=None):
        views = []
        plot = self._createAngDist2D()
        views.append(plot)

        return views

    def _createAngDist2D(self):
        # Common variables to use
        nparts = self.protocol._getInputParticles().getSize()
        title = "Angular Distribution"
        plotter = EmPlotter(x=1, y=1,
                            windowTitle=title)
        sqliteFn = self.protocol._getFileName('projections')
        if not exists(sqliteFn):
            self.createAngDistributionSqlite(sqliteFn, nparts,
                                             itemDataIterator=iterAngles(
                                                 self.protocol._getFileName('anglesFn')))
        plotter.plotAngularDistributionFromMd(sqliteFn, title)

        return plotter

# ===============================================================================

    def _showHistogram(self, param=None):
        fn = self.protocol._getFileName('output_hist')
        f = open(fn)

        views = []
        numberOfBins = 10
        plotter = EmPlotter()
        plotter.createSubPlot("PSF Resolution histogram",
                              "Resolution (A)", "%")
        resolution = [float(line.strip()) for line in f]
        f.close()
        plotter.plotHist(resolution, nbins=numberOfBins)
        plotter.show()

        return views.append(plotter)

    def _showLogFile(self, param=None):
        view = self.textView([self.protocol._getFileName('output_log')],
                             "Output log file")
        return [view]

    def _getVolumeNames(self):
        if self.doShowOutVol.get() == VOL_RS_PSF:
            vols = self.protocol._getFileName('real space PSF')
        else:  # VOL_FS_PSF
            vols = self.protocol._getFileName('fourier space PSF')

        return vols
