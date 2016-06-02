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

import time

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtCTFMicrographs

from protocol_monitor import ProtMonitor
from protocol_monitor_ctf import MonitorCTF



class ProtMonitorSummary(ProtMonitor):
    """ Provide some summary of the basic steps of the Scipion-Box:
    - Import movies
    - Align movies (global and/or local)
    - CTF estimation.
    """
    _label = 'monitor summary'

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        form.addSection('CTF Monitor')
        form.addParam('maxDefocus', params.FloatParam,default=40000,
              label="Raise Alarm if maximum defocus (A) >",
              help="Raise alarm if defocus is greater than given value")
        form.addParam('minDefocus', params.FloatParam,default=1000,
              label="Raise Alarm if minimum defocus (A) <",
              help="Raise alarm if defocus is smaller than given value")
        form.addParam('astigmatism', params.FloatParam,default=0.2,
              label="Raise Alarm if astigmatism >",
              help="Raise alarm if astigmatism is greater than given value")

        form.addParam('monitorTime', params.FloatParam, default=300,
              label="Total Logging time (min)",
              help="Log during this interval")

        form.addSection('Mail settings')
        ProtMonitor._sendMailParams(self, form)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        # finished = False
        # completedDict = {}
        # f = open(self._getPath('log.txt'), 'w')
        #
        # def log(msg):
        #     print >> f, msg
        #     f.flush()
        #
        # log("Inside monitorStep......")
        #
        # while not finished:
        #     log("Iterating over protocols...")
        #     for protPointer in self.inputProtocols:
        #         prot = protPointer.get()
        #         log(" prot.getObjId(): %s" % prot.getObjId())
        #         for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
        #             outSet.load()
        #             outSet.loadAllProperties()
        #             completedDict[(prot.getObjId(), outName)] = outSet.getSize()
        #             outSet.close()
        #     log("=" * 80)
        #     log(pwutils.prettyTime())
        #     log(completedDict)
        #     # Wait some seconds before checking for new data
        #     time.sleep(self.samplingInterval.get())
        #
        # f.close()

        self.createMonitor().loop()


    def _getCtfProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtCTFMicrographs):
                return prot
        return None

    def createMonitor(self):
        ctfProt = self._getCtfProtocol()
        ctfProt.setProject(self.getProject())

        ctfMonitor = MonitorCTF(ctfProt,
                                workingDir=self.workingDir.get(),
                                samplingInterval=self.samplingInterval.get(),
                                monitorTime=self.monitorTime.get(),
                                email=self.createEmailNotifier(),
                                stdout=True,
                                minDefocus=self.minDefocus.get(),
                                maxDefocus=self.maxDefocus.get(),
                                astigmatism=self.astigmatism.get())
        return ctfMonitor
