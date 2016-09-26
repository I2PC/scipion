# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es) [1]
# *              Kevin Savage () [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2]
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

from os.path import abspath

import pyworkflow.utils as pwutils
from pyworkflow.em.protocol.monitors import SummaryProvider


class ReportISPyB:
    """ Create an
    """
    def __init__(self, protocol, ctfMonitor, publishCmd=None,
                 **kwargs):
        # The CTF protocol to monitor
        self.protocol = protocol
        self.provider = SummaryProvider(protocol)
        self.ctfMonitor = ctfMonitor

        # Get the html template to be used, by default use the one
        # in scipion/config/templates
        self.publishCmd = publishCmd

    def info(self, msg):
        if self.protocol._log is not None:
            self.protocol.info(msg)
        else:
            print msg

    def generate(self, finished):
        project = self.protocol.getProject()

        projName = project.getShortName()
        reportDir = abspath(self.protocol._getExtraPath(projName))

        self.info("Creating report directory: %s" % reportDir)
        pwutils.cleanPath(reportDir)
        pwutils.makePath(reportDir)

        acquisitionLines = ''
        self.provider.refreshObjects()

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
                runLines += '{name: "%s",  size:"%s"}' % (obj.output,
                                                          obj.outSize)

            wasProtocol = isProtocol

        print "\n\n====================================="
        print runLines

        # Ctf monitor chart data
        data = [] if self.ctfMonitor is None else self.ctfMonitor.getData()
