# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Pablo Conesa (pconesa@cnb.csic.es)
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

import json
import os
from os.path import join, exists

import pyworkflow.utils as pwutils
from summary_provider import SummaryProvider


class ReportHtml():
    """ Create an html report with a summary of the processing.
    The report will be updated with a given frequency.
    """
    def __init__(self, protocol, ctfMonitor, sysMonitor, **kwargs):
        # The CTF protocol to monitor
        self.protocol = protocol
        self.provider = SummaryProvider(protocol)
        self.ctfMonitor = ctfMonitor
        self.sysMonitor = sysMonitor
        # Get the html template to be used, by default use the one
        # in scipion/config/templates
        self.template = kwargs.get('template',
                                   join(pwutils.getTemplatesFolder(),
                                        'execution.summary.template.html'))

    def getHTMLReportText(self):
        if exists(self.template):
            return open(self.template).read()
        else:
            return ""

    def info(self, msg):
        if self.protocol._log is not None:
            self.protocol.info(msg)
        else:
            print msg

    def generate(self):
        reportTemplate = self.getHTMLReportText()

        if not reportTemplate:
            raise Exception("HTML template file '%s' not found. "
                            % self.template)

        reportName = 'index.html'
        reportDir = self.protocol._getExtraPath('Report')

        self.info("Creating report directory: %s" % reportDir)
        pwutils.cleanPath(reportDir)
        pwutils.makePath(reportDir)

        reportPath = join(reportDir, reportName)

        acquisitionLines = ''

        for item in self.provider.acquisition:
            if not acquisitionLines == '':
                acquisitionLines += ','

            acquisitionLines += '{propertyName:"%s", propertyValue:"%s"}' % item

        runLines = ''
        wasProtocol = None
        for obj in self.provider.getObjects():

            # If it's a protocol
            isProtocol = True if obj.name else False

            if isProtocol:
                if runLines != '': runLines += ']},'
                runLines += '{protocolName: "%s", output:[' % obj.name
            else:
                if not wasProtocol: runLines += ','
                runLines += '{name: "%s",  size:"%s"}' % (obj.output, obj.outSize)

            wasProtocol = isProtocol

        # End the runLines JSON object
        runLines += ']}'

        # Ctf monitor chart data
        data = self.ctfMonitor.getData()
        ctfData = json.dumps(data)

        # system monitor chart data
        data = self.sysMonitor.getData()
        systemData = json.dumps(data)

        args = {'acquisitionLines': acquisitionLines,
                'runLines': runLines,
                'dateStr': pwutils.prettyTime(secs=True),
                'projectName': self.protocol.getProject().getShortName(),
                'scipionVersion': os.environ['SCIPION_VERSION'],
                'ctfData': ctfData,
                'systemData': systemData
                }

        self.info("Writing report html to: %s" % reportPath)
        reportFile = open(reportPath, 'w')
        reportFile.write(reportTemplate % args)
        reportFile.close()

        return reportPath
