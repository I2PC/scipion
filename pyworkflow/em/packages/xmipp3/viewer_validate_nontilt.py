# **************************************************************************
# *
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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
from pyworkflow.em.showj import RENDER, ORDER, VISIBLE, MODE, MODE_MD
from pyworkflow.em.viewer import DataView
import pyworkflow.em.metadata as md
from pyworkflow.protocol.params import LabelParam, StringParam
from protocol_validate_nontilt import XmippProtValidateNonTilt
from plotter import XmippPlotter


IMAGE_INDEX = 'Image index'
P_INDEX = 'Cluster tendency parameter'

class XmippValidateNonTiltViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer validate_nontilt'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtValidateNonTilt]
    
    def _defineParams(self, form):
        form.addSection(label='Show Results Validate NonTilt')
        form.addParam('volForCurve', StringParam, default=1, 
                      label="Show info for volume")
        form.addParam('doShowVolume', LabelParam,
                      label="Display the volume with quality parameter?")
        form.addParam('doShowP', LabelParam,
                      label="Display the clustering tendency curve?")
        
    def _getVisualizeDict(self):
        return {'doShowVolume': self._viewVolume,
                'doShowP': self._viewP
                }        
        
    def _viewVolume(self, e=None):
        import pyworkflow.em.metadata as md 
        viewLabels = 'id enabled _index _filename _xmipp_weight'
        cm = DataView(self.protocol.outputVolumes.getFileName(), viewParams={MODE : MODE_MD,
                                                                             ORDER : viewLabels,
                                                                             RENDER :'_filename',
                                                                             VISIBLE : viewLabels})
        return [cm]
        
    def _viewP(self, e=None):
        volId = int(self.volForCurve.get())
        vol = self.protocol.outputVolumes[volId]
        
        if vol is None: # Wrong volume selection
            return [self.errorMessage("Invalid volume id *%d*" % volId)]

        pFn = vol.clusterMd.get()
        mdCls = md.MetaData(vol.clusterMd.get())
        return [self._viewPlot("Cluster tendency parameter for each image", IMAGE_INDEX, P_INDEX, 
                       mdCls, md.MDL_IMAGE_IDX, md.MDL_WEIGHT, color='b'),
                                DataView(pFn)]
        
    def _viewPlot(self, title, xTitle, yTitle, mdFn, mdLabelX, mdLabelY, color='g'):        
        #xplotter = XmippPlotter(1, 1, figsize=(4,4), windowTitle="Plot")
        #xplotter.createSubPlot(title, xTitle, yTitle)
        #xplotter.plotMdFile(md, mdLabelX, mdLabelY, color)
        mdFn.sort(md.MDL_WEIGHT)
        plotter = XmippPlotter()
        plotter.createSubPlot("", IMAGE_INDEX, P_INDEX)
        plotter.plotMdFile(mdFn, None, md.MDL_WEIGHT)
        return plotter
 
