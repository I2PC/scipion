# **************************************************************************
# *
# * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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

import os
import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor, Monitor
from pyworkflow.protocol.constants import STATUS_RUNNING
from pyworkflow import VERSION_1_1


class ProtMonitorMovieGain(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'movie gain monitor'
    _version = VERSION_1_1

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        # ProtMonitor._defineParams(self, form)
        form.addSection(label='Input')

        form.addParam('inputProtocol', params.PointerParam,
                      label="Input protocol", important=True,
                      pointerClass='XmippProtMovieGain',
                      help="This protocol will be monitorized")

        form.addParam('stddevValue', params.FloatParam, default=0.04,
                      label="Raise Alarm if residual gain standard "
                            "deviation >",
                      help="Raise alarm if residual gain standard deviation "
                           "is greater than given value")
        form.addParam('ratio1Value', params.FloatParam, default=1.15,
                      label="Raise Alarm if the ratio between the 97.5 and "
                            "2.5 percentiles >",
                      help="Raise alarm if the ratio between the 97.5 and "
                           "2.5 percentiles is greater than given value")
        form.addParam('ratio2Value', params.FloatParam, default=4.5,
                      label="Raise Alarm if the ratio between the maximum "
                            "gain value and the 97.5 percentile >",
                      help="Raise alarm if the ratio between the maximum "
                           "gain value and the 97.5 percentile is greater "
                           "than given value")

        form.addParam('monitorTime', params.FloatParam, default=300,
                      label="Total Logging time (min)",
                      help="Log during this interval")

        form.addParam('samplingInterval', params.IntParam, default=60,
                      label="Sampling Interval (sec)",
                      pointerClass='EMProtocol',
                      help="Take one sample each SamplinInteval seconds")
        ProtMonitor._sendMailParams(self, form)

    # -------------------------- STEPS functions ------------------------------
    def monitorStep(self):
        self.createMonitor().loop()


    def createMonitor(self):

        movieGainProt = self.inputProtocol.get()
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

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        prot = self.inputProtocol.get()
        fnWarning = prot._getPath("warningsMonitor.txt")
        if not os.path.exists(fnWarning):
            summary = ["Monitor movie gain, no warnings yet."]
        else:
            fhWarning = open(fnWarning, "r")
            summary = []
            for line in fhWarning.readlines():
                summary.append(line.rstrip())
            fhWarning.close()
        return summary


class MonitorMovieGain(Monitor):
    """ This will be monitoring a movie gain estimation protocol.
    It will internally handle a database to store produced
    movie gain values.
    """
    def __init__(self, protocol, **kwargs):
        Monitor.__init__(self, **kwargs)

        # The movieGain protocol to monitor
        self.protocol = protocol

        self.stddevValue = kwargs['stddevValue']
        self.ratio1Value = kwargs['ratio1Value']
        self.ratio2Value = kwargs['ratio2Value']

    def warning(self, msg):
        self.notify("Scipion Movie Gain Monitor WARNING", msg)

    def initLoop(self):
        pass

    def step(self):
        prot = self.protocol
        fnSummary = prot._getPath("summaryForMonitor.txt")
        fnWarning = prot._getPath("warningsMonitor.txt")
        if not os.path.exists(fnSummary) or os.path.getsize(fnSummary) < 1:
            return False
        fhSummary = open(fnSummary, "r")
        if not os.path.exists(fnWarning):
            fhWarning = open(fnWarning, "w")
        else:
            fhWarning = open(fnWarning, "a")

        line = fhSummary.readlines()[-1]
        stddev, perc25, perc975, maxVal = map(float, line.split()[1:])
        movie_name = map(str, line.split()[0])
        values = line.split()

        if float(values[1]) > self.stddevValue:
            self.warning("Residual gain standard deviation is %f."
                         % stddev)
            fhWarning.write("%s: Residual gain standard deviation is %f.\n"
                         % (movie_name, stddev))

        if (perc975 / perc25) > self.ratio1Value:
            self.warning("The ratio between the 97.5 and 2.5 "
                         "percentiles is %f."
                         % (perc975 / perc25))
            fhWarning.write("%s: The ratio between the 97.5 and 2.5 "
                            "percentiles is %f.\n"
                            % (movie_name, (perc975 / perc25)))

        if (maxVal / perc975) > self.ratio2Value:
            self.warning("The ratio between the maximum gain value "
                         "and the 97.5 percentile is %f."
                         % (maxVal / perc975))
            fhWarning.write("%s: The ratio between the maximum gain value "
                            "and the 97.5 percentile is %f.\n"
                            % (movie_name, (maxVal / perc975)))
        fhSummary.close()
        fhWarning.close()
        return prot.getStatus() != STATUS_RUNNING


    def getData(self):
        idValues = []
        stddevValues = []
        ratio1Values = []
        ratio2Values = []

        prot = self.protocol
        fnSummary = prot._getPath("summaryForMonitor.txt")
        if not os.path.exists(fnSummary) < 1:
            fhSummary = open(fnSummary, "r")
            for idx, line in enumerate(fhSummary.readlines()):
                stddev, perc25, perc975, maxVal = map(float, line.split()[1:])
                idValues.append(idx)
                stddevValues.append(stddev)
                ratio1Values.append(perc975 / perc25)
                ratio2Values.append(maxVal / perc975)
            fhSummary.close()

        data = {
            'idValues': idValues,
            'standard_deviation': stddevValues,
            'ratio1': ratio1Values,
            'ratio2': ratio2Values
        }
        return data

