# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import ObjectView, MODE, MODE_MD, ORDER, VISIBLE
from pyworkflow.em.plotter import EmPlotter

from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy



class XmippCTFDiscrepancyViewer(Viewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer CTF Discrepancy'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtCTFDiscrepancy]
    
    def _visualize(self, obj):
        views = []
        fn = obj.outputCTF.getFileName()
    #            self._views.append(DataView(fn, viewParams={MODE: 'metadata'}))
        labels = 'id enabled _micObj._filename discrepancy _defocusU _defocusV _defocusAngle' 
        views.append(ObjectView(self._project.getName(), obj.strId(), fn, 
                                      viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        
        views.append(self._createPlots(obj))
        
        return views 
    
    def _createPlots(self, obj):
        plotter = EmPlotter()
        plotter.createSubPlot("Discrepancies histogram", 
                      "Resolution", "# of Micrographs")
        discrepancies = [ctf.discrepancy.get() for ctf in obj.outputCTF]
        n = len(discrepancies)
        print "n/20", n/20
        plotter.plotHist(discrepancies, nbins=min(n, 1+n/20))
        return plotter

    
    def _viewEstructureFactor(self, e=None):
        strFactFn = self.protocol._defineStructFactorName()
        md = MetaData(strFactFn)
        return [self._viewPlot("Structure Factor", FREQ_LABEL, 'Structure Factor', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_STRUCTURE_FACTOR),
                self._viewPlot("Structure Factor", FREQ_LABEL, 'log(Structure Factor)', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_LOG_STRUCTURE_FACTOR),
                DataView(strFactFn)]        
    
    def _viewSsnr(self, e=None):
        pass
    
    def _viewVssnr(self, e=None):
        pass
    
