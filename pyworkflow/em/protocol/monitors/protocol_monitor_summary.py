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
from protocol_monitor import ProtMonitor, Monitor
from protocol_monitor_ctf import MonitorCTF
from protocol_monitor_system import MonitorSystem
from pyworkflow import VERSION_1_1
from pyworkflow.em.protocol import ProtCTFMicrographs
from pyworkflow.em.protocol.monitors.report_html import ReportHtml
import getnifs

class ProtMonitorSummary(ProtMonitor):
    """ Provide some summary of the basic steps of the Scipion-Box:
    - Import movies
    - Align movies (global and/or local)
    - CTF estimation.
    """
    _label = 'monitor summary'

    _lastUpdateVersion = VERSION_1_1
    nifs = getnifs.get_network_interfaces()
    nifsNameList = [nif.getName() for nif in nifs]

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

        form.addParam('monitorTime', params.FloatParam, default=30000,
              label="Total Logging time (min)",
              help="Log during this interval")

        form.addSection('System Monitor')
        form.addParam('cpuAlert', params.FloatParam, default=101,
                      label="Raise Alarm if CPU > XX%",
                      help="Raise alarm if memory allocated is greater "
                           "than given percentage")

        form.addParam('memAlert', params.FloatParam, default=101,
                      label="Raise Alarm if Memory > XX%",
                      help="Raise alarm if cpu allocated is greater "
                           "than given percentage")
        form.addParam('swapAlert', params.FloatParam, default=101,
                      label="Raise Alarm if Swap > XX%",
                      help="Raise alarm if swap allocated is greater "
                           "than given percentage")
        group = form.addGroup('GPU')
        group.addParam('doGpu', params.BooleanParam, default=False,
                       label="Check GPU",
                       help="Set to true if you want to monitor the GPU")
        group.addParam('gpusToUse', params.StringParam, default='0',
                          label='Which GPUs to use:', condition='doGpu',
                          help='Providing a list of which GPUs '
                               '(0,1,2,3, etc). Default is monitor GPU 0 only')
        group = form.addGroup('NETWORK')
        group.addParam('doNetwork', params.BooleanParam, default=False,
                       label="Check Network",
                       help="Set to true if you want to monitor the Network")
        group.addParam('netInterfaces', params.EnumParam, choices=self.nifsNameList,
                      default=1,#usually 0 is the loopback
                      label="Interface", condition='doNetwork',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry"
                           " for a description of the symmetry groups format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )

        group = form.addGroup('Disk')
        group.addParam('doDiskIO', params.BooleanParam, default=False,
                       label="Check Disk IO",
                       help="Set to true if you want to monitor the Disk Acces")


        form.addSection('Mail settings')
        ProtMonitor._sendMailParams(self, form)

        form.addSection('HTML Report')
        form.addParam('publishCmd', params.StringParam, default='',
                      label="Publish command",
                      help="Specify a command to publish the template. "
                           "You can use the special token %(REPORT_FOLDER)s "
                           "that will be replaced with the report folder. "
                           "For example: \n"
                           "rsync -av %(REPORT_FOLDER)s scipion@webserver:public_html/")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        ctfMonitor = self.createCtfMonitor()
        sysMonitor = self.createSystemMonitor()
        reportHtml = self.createHtmlReport(ctfMonitor, sysMonitor)

        monitor = Monitor(workingDir=self.workingDir.get(),
                          samplingInterval=self.samplingInterval.get(),
                          monitorTime=self.monitorTime.get())

        def initAll():
            if ctfMonitor is not None:
                ctfMonitor.initLoop()
            sysMonitor.initLoop()

        def stepAll():
            finished = False
            try:
                if ctfMonitor is not None:
                    # FIXME: finished should be True if all monitored protocols
                    # are finished, failed or aborted
                    finished = ctfMonitor.step()
                sysMonitor.step()
                reportHtml.generate(finished)
            except Exception as ex:
                print "An error happened: %s" % ex
            # Stop when the CTF monitor is done
            return finished

        monitor.initLoop = initAll
        monitor.step = stepAll

        monitor.loop()

    def _getCtfProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtCTFMicrographs):
                return prot
        return None

    def createCtfMonitor(self):
        ctfProt = self._getCtfProtocol()

        if ctfProt is None:
            return None

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

    def createSystemMonitor(self):
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)

        sysMonitor = MonitorSystem(protocols,
                                   workingDir=self.workingDir.get(),
                                   samplingInterval=self.samplingInterval.get(),
                                   monitorTime=self.monitorTime.get(),
                                   email=self.createEmailNotifier(),
                                   stdout=True,
                                   cpuAlert=self.cpuAlert.get(),
                                   memAlert=self.memAlert.get(),
                                   swapAlert=self.swapAlert.get(),
                                   doGpu=self.doGpu.get(),
                                   gpusToUse = self.gpusToUse.get(),
                                   doNetwork=self.doNetwork.get(),
                                   doDiskIO=self.doDiskIO.get(),
                                   nif=self.nifsNameList[self.netInterfaces.get()],
                                   )

        return sysMonitor

    def createHtmlReport(self, ctfMonitor=None, sysMonitor=None):
        ctfMonitor = ctfMonitor or self.createCtfMonitor()
        sysMonitor = sysMonitor or self.createSystemMonitor()

        return ReportHtml(self, ctfMonitor, sysMonitor,  self.publishCmd.get(),
                          refreshSecs=self.samplingInterval.get())
