# **************************************************************************
# *
# * Authors:     C.O.S. Sorzano (coss@cnb.csic.es)
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
This module implement the wrappers around xmipp_showj
visualization program.
"""


from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import ObjectView
import pyworkflow.em.showj as showj 

from protocol_solid_angles import XmippProtSolidAngles



class SolidAnglesViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtSolidAngles]
    
    def _visualize(self, obj, **kwargs):
        views = []
        outputClasses = getattr(obj, 'outputClasses', None)

        if outputClasses is not None:
            renderLabels = '_representative._filename _reprojection._filename'
            labels = 'id enabled %s _representative._xmipp_angleRot _representative._xmipp_angleTilt _representative._xmipp_classCount' % renderLabels
    
            views.append(ObjectView(self._project, outputClasses.strId(),
                                    outputClasses.getFileName(),
                                    viewParams={showj.ORDER: labels, 
                                                showj.VISIBLE: labels,
                                                showj.RENDER: renderLabels,
                                                showj.MODE: showj.MODE_MD}))

            outputAverages = getattr(obj, 'outputAverages', None)
            if outputAverages is not None:
                renderLabels = '_filename'
                labels = 'id enabled %s _xmipp_angleRot _xmipp_angleTilt _xmipp_classCount' % renderLabels
                views.append(ObjectView(self._project, outputAverages.strId(),
                                        outputAverages.getFileName(),
                                        viewParams={showj.ORDER: labels, 
                                                    showj.VISIBLE: labels,
                                                    showj.RENDER: renderLabels,
                                                    showj.MODE: showj.MODE_MD}))            
        else:
            views.append(self.infoMessage("No output has been generate yet"))
        
        return views

        
    
    

        
