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
import multiprocessing

from pyworkflow import getTemplatePath
import pyworkflow.utils as pwutils
from summary_provider import SummaryProvider
from pyworkflow.em.convert import ImageHandler

####################### CONSTANTS ##########################
# These constants are the keys used in the ctfMonitor function
# getData() to store paths to micrographs, psd files and shift plots
MIC_PATH = 'imgMicPath'
PSD_PATH = 'imgPsdPath'
SHIFT_PATH = 'imgShiftPath'
# These constants are the name of the folders where thumbnails
# for the html report will be stored. They are also the keys to
# the attribute thumbPaths.
MIC_THUMBS = 'imgMicThumbs'
PSD_THUMBS = 'imgPsdThumbs'
SHIFT_THUMBS = 'imgShiftThumbs'


class ReportHtml:
    """ Create an html report with a summary of the processing.
    The report will be updated with a given frequency.
    """
    def __init__(self, protocol, ctfMonitor, sysMonitor, movieGainMonitor, publishCmd=None, **kwargs):
        # The CTF protocol to monitor
        self.protocol = protocol
        self.reportPath = protocol.reportPath
        self.reportDir = protocol.reportDir
        self.provider = SummaryProvider(protocol)
        self.ctfMonitor = ctfMonitor
        self.sysMonitor = sysMonitor
        self.movieGainMonitor = movieGainMonitor
        self.lastThumbIndex = 0
        self.thumbsReady = 0
        self.thumbPaths = {MIC_THUMBS: [],
                           PSD_THUMBS: [],
                           SHIFT_THUMBS: []}

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

    def checkNewThumbsReady(self):
        for i in range(self.thumbsReady, self.lastThumbIndex):
            pathsReady = [exists(join(self.reportDir, self.thumbPaths[k][i])) for k in self.thumbPaths.keys()]
            if all(pathsReady):
                self.thumbsReady += 1
        return self.thumbsReady

    def getThumbnailPaths(self, ctfData, ext="png"):
        """Adds to self.thumbPaths the paths to the report thumbnails,
        and creates their folders in the reportDir if they don't exist.

        ===== Params =====
        - ctfData: dict resulting from calling ctfMonitor.getData() containing
                   paths to the original images.
        - ext: extension of the thumbnail images. Defaults to png.

        """
        copyShiftPlots = SHIFT_PATH in ctfData

        # If we're in the first round, create thumbnail folders
        if self.lastThumbIndex == 0:
            if not copyShiftPlots:
                self.thumbPaths.pop(SHIFT_PATH)
            folderPaths = [join(self.reportDir, f) for f in self.thumbPaths.keys()]
            pwutils.makePath(*folderPaths)

        for i in range(self.lastThumbIndex, len(ctfData[MIC_PATH])):
            micPath = ctfData[MIC_PATH][i]
            micThumb = join(MIC_THUMBS, pwutils.replaceExt(basename(micPath), ext))
            self.thumbPaths[MIC_THUMBS].append(micThumb)

            psdPath = ctfData[PSD_PATH][i]
            movie = basename(os.path.dirname(psdPath))
            psdThumb = join(PSD_THUMBS, "%s_%s" % (movie, pwutils.replaceExt(basename(psdPath), ext)))
            self.thumbPaths[PSD_THUMBS].append(psdThumb)

            if copyShiftPlots:
                shiftPath = ctfData[SHIFT_PATH][i]
                shiftCopy = join(SHIFT_THUMBS, pwutils.replaceExt(basename(shiftPath), ext))
                self.thumbPaths[SHIFT_THUMBS].append(shiftCopy)

        return


    def generateReportImages(self, ctfData, firstThumbIndex=0, micScaleFactor=6):
        """ Function to generate thumbnails for the report.

        ===== Params =====
        - ctfData: dict resulting from calling ctfMonitor.getData()
        - firstThumbIndex: index from which we start generating thumbnails
        - micScaleFactor: how much to reduce in size the micrographs.
        """
        ih = ImageHandler()

        numMics = len(ctfData[MIC_PATH])
        numShiftPlots = len(ctfData.get(SHIFT_PATH, []))

        for i in range(firstThumbIndex, numMics):
            print('Generating images for mic %d' % i)
            # mic thumbnails
            dstImgPath = join(self.reportDir, self.thumbPaths[MIC_THUMBS][i])
            if not exists(dstImgPath):
                ih.computeThumbnail(ctfData[MIC_PATH][i], dstImgPath, scaleFactor=micScaleFactor)

            # psd thumbnails
            dstImgPath = join(self.reportDir, self.thumbPaths[PSD_THUMBS][i])
            if not exists(dstImgPath):
                ih.computeThumbnail(ctfData[PSD_PATH][i], dstImgPath, scaleFactor=1)

            # shift plots
            if numShiftPlots:
                dstImgPath = join(self.reportDir, self.thumbPaths[SHIFT_THUMBS][i])
                if not exists(dstImgPath):
                    pwutils.createAbsLink(ctfData[SHIFT_PATH][i], dstImgPath)

        return

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
        projName = project.getShortName()
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
        reportFinished = True

        if data:
            numImages = len(data[MIC_PATH])

            if self.lastThumbIndex < numImages:  # if there are new images to generate in this round
                # get the thumbnail paths
                self.getThumbnailPaths(data)
                newImages = numImages - self.lastThumbIndex
                if newImages < 10:
                    # we have few new images, eg streaming mode, generate thumbnails now
                    self.generateReportImages(data, self.lastThumbIndex)
                else:
                    # we have many images, generate thumbs in a separate process
                    processName = 'Images %d to %d' % (self.lastThumbIndex, numImages)
                    process = multiprocessing.Process(name=processName, target=self.generateReportImages,
                                                      args=(data, self.lastThumbIndex))
                    process.start()
                # update number of thumbnails we generated
                self.lastThumbIndex = len(data[MIC_PATH])

            self.thumbsReady = self.checkNewThumbsReady()

            # check if we generated any new images in this round
            if self.thumbsReady < self.lastThumbIndex:
                reportFinished = False

            if len(data['defocusU']) < 100:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'])
            else:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'][:-50])
                data['defocusCoverageLast50'] = self.processDefocusValues(data['defocusU'][-50:])

            # remove data of original pictures, report doesn't need it
            for k in [MIC_PATH, PSD_PATH, SHIFT_PATH]:
                data.pop(k)
            # send over only thumbnails of the mics that have been fully processed
            for k in self.thumbPaths:
                data[k] = self.thumbPaths[k][:self.thumbsReady]

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

        self.info("Writing report html to: %s" % abspath(self.reportPath))
        pwutils.cleanPath(self.reportPath)
        reportFile = open(self.reportPath, 'w')
        reportTemplate = reportTemplate % args
        reportFile.write(reportTemplate.encode('utf-8'))
        reportFile.close()

        if self.publishCmd:
            self.info("Publishing the report:")
            cmd = self.publishCmd % {'REPORT_FOLDER': self.reportDir}
            self.info(cmd)
            os.system(cmd)

        return reportFinished
