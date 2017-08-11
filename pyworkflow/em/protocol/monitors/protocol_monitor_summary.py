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
from protocol_monitor_movie_gain import MonitorMovieGain
from protocol_monitor_system import MonitorSystem
from pyworkflow import VERSION_1_1
from pyworkflow.em.protocol import ProtCTFMicrographs, ProtProcessMovies
from pyworkflow.em.protocol.monitors.report_html import ReportHtml


class ProtMonitorSummary(ProtMonitor):
    """ Provide some summary of the basic steps of the Scipion-Box:
    - Import movies
    - Align movies (global and/or local)
    - CTF estimation
    - Movie gain estimation.
    """
    _label = 'monitor summary'
    _lastUpdateVersion = VERSION_1_1

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        form.addSection('MovieGain Monitor')
        form.addParam('stddevValue', params.FloatParam, default=0.04,
                      label="Raise Alarm if residual gain standard "
                            "deviation >",
                      help="Raise alarm if residual gain standard deviation "
                           "is greater than given value")
        form.addParam('ratio1Value', params.FloatParam, default=1.15,
                      label="Raise Alarm if the ratio between the 97.5 "
                            "and 2.5 percentiles >",
                      help="Raise alarm if the ratio between the 97.5 "
                           "and 2.5 percentiles is greater than given value")
        form.addParam('ratio2Value', params.FloatParam, default=4.5,
                      label="Raise Alarm if the ratio between the maximum "
                            "gain value and the 97.5 percentile >",
                      help="Raise alarm if the ratio between the maximum "
                           "gain value and the 97.5 percentile is greater "
                           "than given value")

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
        movieGainMonitor = self.createMovieGainMonitor()
        ctfMonitor = self.createCtfMonitor()
        sysMonitor = self.createSystemMonitor()
        reportHtml = self.createHtmlReport(ctfMonitor, sysMonitor, movieGainMonitor)

        monitor = Monitor(workingDir=self.workingDir.get(),
                          samplingInterval=self.samplingInterval.get(),
                          monitorTime=self.monitorTime.get())

        def initAll():
            if ctfMonitor is not None:
                ctfMonitor.initLoop()
            if movieGainMonitor is not None:
                movieGainMonitor.initLoop()
            sysMonitor.initLoop()

        def stepAll():
            finished = False
            try:
                if ctfMonitor is not None:
                    # FIXME: finished should be True if all monitored protocols
                    # are finished, failed or aborted
                    finished = ctfMonitor.step()

                # sysmonitor watches all input protocols so
                # when sysmonitor done we done
                finished = sysMonitor.step()

                if movieGainMonitor is not None:
                    finished = movieGainMonitor.step()
                    
                reportHtml.generate(finished)
            except Exception as ex:
                print("An error happened: %s" % ex)
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

    def _getMovieGainProtocol(self):
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            if isinstance(prot, ProtProcessMovies):
                return prot
        return None

    def createMovieGainMonitor(self):
        movieGainProt = self._getMovieGainProtocol()

        if movieGainProt is None:
            return None

        movieGainProt.setProject(self.getProject())

        movieGainMonitor = MonitorMovieGain(movieGainProt,
                                            workingDir=self.workingDir.get(),
                                            samplingInterval=self.samplingInterval.get(),
                                            monitorTime=self.monitorTime.get(),
                                            email=self.createEmailNotifier(),
                                            stdout=True,
                                            stddevValue=self.stddevValue.get(),
                                            ratio1Value = self.ratio1Value.get(),
                                            ratio2Value = self.ratio2Value.get())
        return movieGainMonitor

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
        protocols = self.getInputProtocols()

        sysMonitor = MonitorSystem(protocols,
                                   workingDir=self.workingDir.get(),
                                   samplingInterval=self.samplingInterval.get(),
                                   monitorTime=self.monitorTime.get(),
                                   email=self.createEmailNotifier(),
                                   stdout=True,
                                   cpuAlert=self.cpuAlert.get(),
                                   memAlert=self.memAlert.get(),
                                   swapAlert=self.swapAlert.get())
        return sysMonitor

    def createHtmlReport(self, ctfMonitor=None, sysMonitor=None, movieGainMonitor=None):
        ctfMonitor = ctfMonitor or self.createCtfMonitor()
        sysMonitor = sysMonitor or self.createSystemMonitor()
        movieGainMonitor = movieGainMonitor or self.createMovieGainMonitor()

        return ReportHtml(self, ctfMonitor, sysMonitor, movieGainMonitor, self.publishCmd.get(),
                          refreshSecs=self.samplingInterval.get())