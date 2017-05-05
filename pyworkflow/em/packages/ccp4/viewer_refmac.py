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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from protocol_refmac import CCP4ProtRunRefmac
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer, Viewer
from pyworkflow.gui.text import _open_cmd
from pyworkflow.em.data import EMSet, EMObject
from pyworkflow.object import Float, String
from pyworkflow.em.viewer import ObjectView, TableView

from pyworkflow.em.viewer import ImageView, ChimeraView

class CCP4ProtRunRefmacViewer(ProtocolViewer):
    """ Viewer for CCP4 program refmac
    """
    _label = 'Refmac Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CCP4ProtRunRefmac]

    def __init(self):
        ProtocolViewer.__init__()

    # ROB: do we need this memory for something?
    # _memory = False
    # temporary metadata file with ctf that has some resolution greathan than X
    # tmpMetadataFile = 'viewersTmp.sqlite'

    def _defineParams(self, form):
        form.addSection(label='Visualization of Refmac results')
        # group = form.addGroup('Overall results')
        form.addParam('displayMask', LabelParam,
                      label="PDB based Mask",
                      help="Display Masked map")
        form.addParam('showFinalResults', LabelParam,
                      label="Final Results Table",
                      help="Table of Final Results from refine.log file")
        form.addParam('showLogFile', LabelParam,
                      label="Show log file",
                      help="open refmac log file in a text editor")
        form.addParam('showLastIteration', LabelParam,
                      label="Results Table (last iteration)",
                      help="Table stored in log file summarizing the last iteration")
        form.addParam('displayRFactorPlot', LabelParam,
                      label="R-factor vs. iteration",
                      help="Plot R-factor as a function of the iteration")
        form.addParam('displayFOMPlot', LabelParam,
                      label="FOM vs. iteration",
                      help="Plot Figure Of Merit as a function of the iteration")
        form.addParam('displayLLPlot', LabelParam,
                      label="-LL vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayLLfreePlot', LabelParam,
                      label="-LLfree vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayGeometryPlot', LabelParam,
                      label="Geometry vs. iteration",
                      help="""Plot Geometry as a function of the iteration:
Geometry includes rmsBOND (root mean square bond lengths)
zBOND (zscore of the deviation of bond lengths)
rmsANGL (root mean square bond angles)
zANGL (zscore of the deviation of bond angles)
and rmsCHIRAL (root mean square of chiral index""")

    def _getVisualizeDict(self):
        return {
            'displayMask': self._visualizeMask,
            'showFinalResults': self._visualizeFinalResults,
            'showLogFile': self._visualizeLogFile,
            'showLastIteration': self._visualizeLastIteration,
            'displayRFactorPlot': self._visualizeRFactorPlot,
            'displayFOMPlot': self._visualizeFOMPlot,
            'displayLLPlot': self._visualizeLLPlot,
            'displayLLfreePlot': self._visualizeLLfreePlot,
            'displayGeometryPlot': self._visualizeGeometryPlot
        }

    def _visualizeMask(self, e=None):
        pass

    def _visualizeFinalResults(self, e=None):
        """
        views = []
        labels = '_1 _2'
        emSet = EMSet(filename="/tmp/kk.sqlite")
        emObject = EMObject()
        emObject._1 = String('first parameter')
        emObject._2 = Float(12.)
        emSet.append(emObject)
        emObject = EMObject()
        emObject._1 = String('second parameter')
        emObject._2 = Float(22.)
        emSet.append(emObject)
        emSet.write()
        views.append(ObjectView(self._project,
                                self.protocol.strId(),
                                "/tmp/kk.sqlite",
                                viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        return views
"""
        headerList = ['property', 'value']
        dataList = [
        ('1John', 'Smith1') ,
        ('2Larry', 'Black') ,
        ('3Walter', 'White') ,
        ('4Fred', 'Becker'),
        ('5John', 'Smith') ,
        ('6Larry', 'Black') ,
        ('7Walter', 'White') ,
        ('8Fred', 'Becker'),
        ('9John', 'Smith') ,
        ('10Larry', 'Black') ,
        ('11Walter', 'White') ,
        ('12Fred', 'Becker'),
        ('13John', 'Smith') ,
        ('14Larry', 'Black') ,
        ('15Walter', 'White') ,
        ('16Fred', 'Becker')
        ]

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg="This is a list with many important persons. This comment will be a bit longer",
                  title= "Refmac: Final Results Summary",
                  height=len(dataList))

    def _visualizeLogFile(self, e=None):
        """Show refmac log file."""
        refineLogFileName = self.protocol._getExtraPath(self.protocol.refineLogFileName)
        _open_cmd(refineLogFileName, self.getTkRoot())

    def _visualizeLastIteration(self, e=None):
        pass

    def _visualizeRFactorPlot(self, e=None):
        pass

    def _visualizeFOMPlot(self, e=None):
        pass

    def _visualizeLLPlot(self, e=None):
        pass

    def _visualizeLLfreePlot(self, e=None):
        pass

    def _visualizeGeometryPlot(self, e=None):
        pass
