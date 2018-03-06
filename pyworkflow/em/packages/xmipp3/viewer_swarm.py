# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Nov 2014
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

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer, ObjectView, DataView
import pyworkflow.em.showj as showj
import xmipp

from protocol_reconstruct_swarm import XmippProtReconstructSwarm

class XmippReconstructSwarmViewer(XmippViewer):
    """ Visualize the output of protocol reconstruct swarm """
    _label = 'viewer reconstruct swarm'
    _targets = [XmippProtReconstructSwarm]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

    def _visualize(self, obj, **args):
        import os
        fnVolume = self.protocol._getExtraPath("volumeAvg.vol")
        if os.path.exists(fnVolume):
            fnDir = self.protocol._getExtraPath()
            samplingRate=self.protocol.readInfoField(fnDir,"sampling",xmipp.MDL_SAMPLINGRATE)
            self._views.append(ObjectView(self._project, None, fnVolume, viewParams={showj.RENDER: 'image', showj.SAMPLINGRATE: samplingRate}))
        
        fnSwarm = self.protocol._getExtraPath("swarm.xmd")
        if os.path.exists(fnSwarm):
            self._views.append(DataView('bestByVolume@' + fnSwarm, viewParams = {showj.MODE: showj.MODE_MD}))
        
