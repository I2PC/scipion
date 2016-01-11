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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView, ChimeraView
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer

from protocol_volume_strain import XmippProtVolumeStrain

class XmippVolumeStrainViewer(XmippViewer):
    """ Visualize the output of protocol volume strain """
    _label = 'viewer split volume'
    _targets = [XmippProtVolumeStrain]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

    def _visualize(self, obj, **args):
        import os
        fnCmd = self.protocol._getPath('result_strain_chimera.cmd')
        if os.path.exists(fnCmd):
            self._views.append(ChimeraView(fnCmd))
        fnCmd = self.protocol._getPath('result_localrot_chimera.cmd')
        if os.path.exists(fnCmd):
            self._views.append(ChimeraView(fnCmd))
        fnCmd = self.protocol._getPath('result_morph_chimera.cmd')
        if os.path.exists(fnCmd):
            self._views.append(ChimeraView(fnCmd))

