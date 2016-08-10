# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              
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

import os
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, Viewer
from protocol_import import ProtImportTiltSeries
from pyworkflow.em.viewer import DataView
from pyworkflow.em.showj import RENDER, ORDER, VISIBLE, MODE, MODE_MD


class ProtImportTiltSeriesViewer(Viewer):
    
    _label = 'viewer import tiltSeries'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtImportTiltSeries]
        
    def _visualize(self, e=None):
        viewLabels = 'id enabled image angleTilt _ctfModel._defocusU'
        fnMd = self.protocol._defineOutputMdName()
        
        if not os.path.exists(fnMd):
            return [self.errorMessage('The necessary metadata was not produced.\n'
                                      'Please wait if protocol is still running, '
                                      'otherwise execute it again.',
                                      title='Missing result file')]
        
        cm = DataView(fnMd, viewParams={MODE : MODE_MD,
                                        ORDER : viewLabels,
                                        RENDER :'image',
                                        VISIBLE : viewLabels})    
        return [cm]
    
    
    