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

import os

from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.utils import exists
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.em.viewer import (ImageView, ChimeraView,
                                  ChimeraClientView, ObjectView)
from protocol_3dfsc import Prot3DFSC

VOL_ORIG = 0
VOL_TH = 1
VOL_THBIN = 2
VOL_ALL = 3

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1


class ThreedFscViewer(ProtocolViewer):
    """
    Visualization tools for 3D FSC results.
    
    3D FSC is software tool for quantifying directional
    resolution using 3D Fourier shell correlation volumes.

    Find more information at https://github.com/nysbc/Anisotropy
    """
           
    _environments = [DESKTOP_TKINTER]
    _targets = [Prot3DFSC]
    _label = 'viewer 3D FSC'

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
        group.addParam('doShowOutVol', EnumParam, default=VOL_ORIG,
                       choices=['original', 'thresholded',
                                'thresholded and binarized', 'all'],
                       display=EnumParam.DISPLAY_COMBO,
                       label='3D FSC volume to display')

        form.addParam('doShowHistogram', LabelParam,
                      label="Show histogram and directional FSC plot")
        form.addParam('doShowPlotFT', LabelParam,
                      label="Show Fourier Transform Power plot")
        form.addParam('doShowPlot3DFSC', LabelParam,
                      label="Show 3D FSC plots")
        form.addParam('doShowChimera', LabelParam,
                      label="Show Chimera animation", default=True,
                      help="Display 3D FSC and coloring original map by "
                           "angular resolution.")
        
    def _getVisualizeDict(self):
        self.protocol._initialize()  # Load filename templates
        return {'doShowOutVol': self._showVolumes,
                'doShowHistogram': self._showHistogram,
                'doShowPlotFT': self._showPlotFT,
                'doShowPlot3DFSC': self._showPlot3DFSC,
                'doShowChimera': self._showChimera
                }

# =============================================================================
# ShowVolumes
# =============================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()

        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for volFn in volumes:
                # We assume that the chimera script will be generated
                # at the same folder as 3DFSC volumes
                vol = volFn.replace(':mrc', '')
                localVol = os.path.basename(vol)
                if exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = ChimeraView(cmdFile)
        else:
            # view = CommandView('xmipp_chimera_client --input "%s"
            # --mode projector 256 &' % volumes[0])
            view = ChimeraClientView(volumes[0])

        return [view]

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('3DFSC_viewer_volumes.sqlite')
        samplingRate = self.protocol.inputVolume.get().getSamplingRate()

        vols = self._getVolumeNames()
        files = []
        for vol in vols:
            if os.path.exists(vol):
                files.append(vol)
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

# =============================================================================

    def _showHistogram(self, param=None):
        return [ImageView(self.protocol._getFileName('out_histogram'))]

    def _showPlotFT(self, param=None):
        return [ImageView(self.protocol._getFileName('out_plotFT'))]

    def _showPlot3DFSC(self, param=None):
        return [ImageView(self.protocol._getFileName('out_plot3DFSC'))]

    def _showChimera(self, param=None):
        #os.system('chimera "%s" &' % self.protocol._getFileName('out_cmdChimera'))
        cmdFile = self.protocol._getFileName('out_cmdChimera')
        view = ChimeraView(cmdFile)
        return [view]

    def _getVolumeNames(self):
        volsFn = ['out_vol3DFSC', 'out_vol3DFSC-th', 'out_vol3DFSC-thbin']

        if self.doShowOutVol.get() == VOL_ORIG:
            vols = [self.protocol._getFileName(volsFn[0])]
        elif self.doShowOutVol.get() == VOL_TH:
            vols = [self.protocol._getFileName(volsFn[1])]
        elif self.doShowOutVol.get() == VOL_THBIN:
            vols = [self.protocol._getFileName(volsFn[2])]
        else:
            vols = [self.protocol._getFileName(f) for f in volsFn]

        return vols
