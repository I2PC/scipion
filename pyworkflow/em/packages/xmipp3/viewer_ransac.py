# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
"""
This module implement the wrappers aroung Xmipp Ransac Protocol
visualization program.
"""
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_ransac import XmippProtRansac
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer

class XmippViewerRansac(XmippViewer):
    """ Wrapper to visualize Ransac results """
    _label = 'viewer ransac'
    _targets = [XmippProtRansac]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

    def _visualize(self, obj, **args):
        if hasattr(self.protocol, "outputVolumes"):
            obj = self.protocol.outputVolumes
            fn = obj.getFileName()
            labels = 'id enabled comment _filename _xmipp_volScoreSum _xmipp_volScoreSumTh _xmipp_volScoreMean _xmipp_volScoreMin'
            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={showj.ORDER: labels, 
                                                      showj.VISIBLE: labels, 
                                                      showj.MODE: showj.MODE_MD,
                                                      showj.RENDER:'_filename'}))
