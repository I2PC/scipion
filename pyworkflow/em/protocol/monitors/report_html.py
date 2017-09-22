# -*- coding: utf-8 -*-
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
from os.path import join, exists, abspath, basename
from numpy import histogram, array, insert
from operator import itemgetter

from pyworkflow import getTemplatePath
import pyworkflow.utils as pwutils
from summary_provider import SummaryProvider
from pyworkflow.em.convert import ImageHandler



class ReportHtml:
    """ Create an html report with a summary of the processing.
    The report will be updated with a given frequency.
    """
    def __init__(self, protocol, ctfMonitor, sysMonitor, movieGainMonitor, publishCmd=None,
                 **kwargs):
        # The CTF protocol to monitor
        self.protocol = protocol
        self.provider = SummaryProvider(protocol)
        self.ctfMonitor = ctfMonitor
        self.sysMonitor = sysMonitor
        self.movieGainMonitor = movieGainMonitor

        # Get the html template to be used, by default use the one
        # in scipion/config/templates
        defaultTemplate = getTemplatePath('execution.summary.template.html')
        self.template = kwargs.get('template', defaultTemplate)

        self.publishCmd = publishCmd
        self.refreshSecs = kwargs.get('refreshSecs', 60)

    def getHTMLReportText(self):
        if exists(self.template):
            return open(self.template, 'rb').read().decode('utf-8')
        else:
            return ""

    def info(self, msg):
        if self.protocol._log is not None:
            self.protocol.info(msg)
        else:
            print(msg)

    @staticmethod
    def generateReportImages(imgPathList, targetDir, imageType, relPath=None, ext="png", **kwargs):
        outputPaths = []
        ih = ImageHandler()
        if not exists(targetDir):
            pwutils.makePath(targetDir)

        for imgPath in imgPathList:
            if imageType == "psd":
                movie = basename(os.path.dirname(imgPath))
                imgOutPath = join(targetDir,
                                  "%s_%s" % (movie, pwutils.replaceExt(basename(imgPath), ext)))
                if not exists(imgOutPath):
                    ih.computeThumbnail(imgPath, imgOutPath, scaleFactor=1)

            elif imageType == "mic":
                imgOutPath = join(targetDir, pwutils.replaceExt(basename(imgPath), ext))
                if not exists(imgOutPath):
                    scaleFactor = kwargs.get('scaleFactor', 6)
                    ih.computeThumbnail(imgPath, imgOutPath, scaleFactor=scaleFactor)

            elif imageType == "shift":
                imgOutPath = join(targetDir,  pwutils.replaceExt(basename(imgPath), ext))
                if not exists(imgOutPath):
                    pwutils.copyFile(imgPath, targetDir)

            if relPath is not None:
                # append path relative to relPath
                outputPaths.append(os.path.relpath(imgOutPath, relPath))
            else:
                # append absolute path
                outputPaths.append(imgOutPath)

        return outputPaths

    def processDefocusValues(self, defocusList):
        maxDefocus = self.protocol.maxDefocus.get()
        minDefocus = self.protocol.minDefocus.get()
        edges = array(range(0, int(maxDefocus)+1, 2500))
        edges = insert(edges[edges > minDefocus], 0, minDefocus)
        values, binEdges = histogram(defocusList, bins=edges, range=(minDefocus, maxDefocus))
        belowThresh = 0
        aboveThresh = 0
        labels = ["%d-%d" % (x[0], x[1]) for x in zip(binEdges, binEdges[1:])]
        for v in defocusList:
            if v < minDefocus:
                belowThresh += 1
            elif v > maxDefocus:
                aboveThresh += 1
        zipped = zip(values, labels)
        zipped[:0] = [(-belowThresh, "0-%d" % minDefocus)]  # Negative values for highcharts heatmap color scale
        zipped.append((-aboveThresh, "> %d" % maxDefocus))

        return zipped

    def generate(self, finished):
        reportTemplate = self.getHTMLReportText()

        if not reportTemplate:
            raise Exception("HTML template file '%s' not found. "
                            % self.template)

        project = self.protocol.getProject()

        reportName = 'index.html'
        projName = project.getShortName()
        reportDir = abspath(self.protocol._getExtraPath(projName))

        if not exists(reportDir):
            self.info("Creating report directory: %s" % reportDir)
            pwutils.makePath(reportDir)

        reportPath = join(reportDir, reportName)
        pwutils.cleanPath(reportPath)

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

        # End the runLines JSON object
        runLines += ']}'

        # Ctf monitor chart data
        data = [] if self.ctfMonitor is None else self.ctfMonitor.getData()
        if data:
            # generate mic thumbnails
            micThumbsDir = join(reportDir, 'imgMicThumbs')
            data['imgMicThumbs'] = self.generateReportImages(data['imgMicPath'], micThumbsDir,
                                                             'mic', relPath=reportDir)
            # generate psd thumbnails
            psdThumbsDir = join(reportDir, 'imgPsdThumbs')
            data['imgPsdThumbs'] = self.generateReportImages(data['imgPsdPath'], psdThumbsDir,
                                                             'psd', relPath=reportDir)

            # copy shift plots
            if data['imgShiftPath'] and data['imgShiftPath'][0] != "":
                shiftPlotsDir = join(reportDir, 'imgShiftPlots')
                data['imgShiftCopyPath'] = self.generateReportImages(data['imgShiftPath'], shiftPlotsDir,
                                                                     'shift', relPath=reportDir)
            if len(data['defocusU']) < 100:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'])
            else:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'][:-50])
                data['defocusCoverageLast50'] = self.processDefocusValues(data['defocusU'][-50:])

        ctfData = json.dumps(data)

        # Movie gain monitor chart data
        data = [] if self.movieGainMonitor is None else self.movieGainMonitor.getData()

        movieGainData = json.dumps(data)

        # system monitor chart data
        data = self.sysMonitor.getData()
        systemData = json.dumps(data)

        args = {'projectName': projName,
                'startTime': pwutils.dateStr(project.getCreationTime(), secs=True),
                'dateStr': pwutils.prettyTime(secs=True),
                'projectDuration': pwutils.prettyDelta(project.getElapsedTime()),
                'projectStatus': "FINISHED" if finished else "RUNNING",
                'scipionVersion': os.environ['SCIPION_VERSION'],
                'acquisitionLines': acquisitionLines,
                'runLines': runLines,
                'ctfData': ctfData,
                'movieGainData': movieGainData,
                'systemData': systemData,
                'refresh': '<META http-equiv="refresh" content="%s" >' % self.refreshSecs if not finished else '',
                }

        self.info("Writing report html to: %s" % abspath(reportPath))
        reportFile = open(reportPath, 'w')
        reportTemplate = reportTemplate % args
        reportFile.write(reportTemplate.encode('utf-8'))
        reportFile.close()

        if self.publishCmd:
            self.info("Publishing the report:")
            cmd = self.publishCmd % {'REPORT_FOLDER': reportDir}
            self.info(cmd)
            os.system(cmd)

        return reportPath
