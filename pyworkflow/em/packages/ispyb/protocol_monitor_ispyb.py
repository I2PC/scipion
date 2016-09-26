# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtMonitor, Monitor, MonitorCTF
from pyworkflow.em.protocol import ProtCTFMicrographs

from report_ispyb import MonitorISPyB


class ProtMonitorISPyB(ProtMonitor):
    """ Provide some summary of the basic steps of the Scipion-Box:
    - Import movies
    - Align movies (global and/or local)
    - CTF estimation.
    """
    _label = 'monitor to ISPyB'

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        group = form.addGroup('Experiment')
        group.addParam('groupid', params.StringParam,
                      label="Group Id",
                      help="Group Id")
        group.addParam('visit', params.StringParam,
                      label="Visit",
                      help="Visit")
        group.addParam('sampleid', params.StringParam,
                      label="Sample Id",
                      help="Sample Id")
        group.addParam('detectorid', params.StringParam,
                      label="Detector Id",
                      help="Detector Id")

        form.addParam('db', params.EnumParam,
                      choices=["production", "devel", "test"],
                      label="Database",
                      help="Select which ISPyB database you want to use.")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        monitor = MonitorISPyB(self, refreshSecs=self.samplingInterval.get(),
                               workingDir=self._getPath(),
                               samplingInterval=self.samplingInterval.get(),
                               monitorTime=100
                               )
        monitor.initLoop()
        monitor.loop()

