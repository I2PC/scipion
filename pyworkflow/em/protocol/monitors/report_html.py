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
        self.thumbsDone = 0

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

    def getThumbnailPaths(self, ctfData, ext="png"):
        """Adds to the dict ctfData the paths to the report thumbnails,
        and creates their folders in the reportDir if they don't exist.
        Also checks which thumbnails are ready. This function doesn't
        actually generate thumbnails, only adds their paths to ctfData.

        ===== Params =====
        - ctfData: dict resulting from calling ctfMonitor.getData()
        - ext: extension of the thumbnail images. Defaults to png

        ===== Returns =====
        - micsReady: int with the number of micrographs that are ready
                     to be displayed in the report i.e. how many rows
                     in the table have all their thumbnails generated.
        """
        micsReady = 0
        micsFolder = 'imgMicThumbs'
        psdFolder = 'imgPsdThumbs'
        shiftPlotsFolder = 'imgShiftThumbs'
        ctfData[micsFolder] = []
        ctfData[psdFolder] = []
        micThumbDir = join(self.reportDir, micsFolder)
        psdThumbDir = join(self.reportDir, psdFolder)
        # Create folders if they don't exist
        if 'imgShiftPath' in ctfData:
            shiftPlotDir = join(self.reportDir, shiftPlotsFolder)
            pwutils.makePath(micThumbDir, psdThumbDir, shiftPlotDir)
            ctfData[shiftPlotsFolder] = []
            copyShiftPlots = True
        else:
            pwutils.makePath(micThumbDir, psdThumbDir)
            copyShiftPlots = False

        for i in range(len(ctfData.get('imgMicPath', []))):
            micPath = ctfData['imgMicPath'][i]
            micThumb = join(micsFolder, pwutils.replaceExt(basename(micPath), ext))
            ctfData[micsFolder].append(micThumb)

            psdPath = ctfData['imgPsdPath'][i]
            movie = basename(os.path.dirname(psdPath))
            psdThumb = join(psdFolder, "%s_%s" % (movie, pwutils.replaceExt(basename(psdPath), ext)))
            ctfData[psdFolder].append(psdThumb)

            if copyShiftPlots:
                shiftPath = ctfData['imgShiftPath'][i]
                shiftCopy = join(shiftPlotsFolder, pwutils.replaceExt(basename(shiftPath), ext))
                ctfData[shiftPlotsFolder].append(shiftCopy)
                shiftPlotReady = exists(join(self.reportDir, shiftCopy))
            else:
                shiftPlotReady = True

            # Check if all thumbnails are available
            if exists(join(self.reportDir, micThumb)) and exists(join(self.reportDir, psdThumb)) and shiftPlotReady:
                micsReady += 1

        return micsReady

    def generateReportImages(self, ctfData, firstThumbIndex=0, micScaleFactor=6):
        """ Function to generate thumbnails for the report.

        ===== Params =====
        - ctfData: dict resulting from calling ctfMonitor.getData() after being
                   modified in getThumbnailPaths.
        - firstThumbIndex: index from which we start generating thumbnails
        - micScaleFactor: how much to reduce in size the micrographs.
        """
        ih = ImageHandler()
        micsFolder = 'imgMicThumbs'
        psdFolder = 'imgPsdThumbs'
        shiftPlotsFolder = 'imgShiftThumbs'

        numMics = len(ctfData['imgMicPath'])
        numShiftPlots = len(ctfData.get('imgShiftPath', []))
        print('Generating thumbnails for micrographs %s to %s' % (firstThumbIndex, numMics))

        for i in range(firstThumbIndex, numMics):
            print('Generating images for mic %d' % i)
            # mic thumbnails
            dstImgPath = join(self.reportDir, ctfData[micsFolder][i])
            if not exists(dstImgPath):
                ih.computeThumbnail(ctfData['imgMicPath'][i], dstImgPath, scaleFactor=micScaleFactor)

            # psd thumbnails
            dstImgPath = join(self.reportDir, ctfData[psdFolder][i])
            if not exists(dstImgPath):
                ih.computeThumbnail(ctfData['imgPsdPath'][i], dstImgPath, scaleFactor=1)

            # shift plots
            if numShiftPlots:
                dstImgPath = join(self.reportDir, ctfData[shiftPlotsFolder][i])
                if not exists(dstImgPath):
                    pwutils.createAbsLink(ctfData['imgShiftPath'][i], dstImgPath)

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
            # get the thumbnail paths
            micsReady = self.getThumbnailPaths(data)
            # generate actual images in a separate thread
            processName = 'Images %d to %d' % (self.thumbsDone, len(data['imgMicPath']))
            process = multiprocessing.Process(name=processName, target=self.generateReportImages,
                                              args=(data, self.thumbsDone))
            process.start()
            # update number of processed thumbnails
            self.thumbsDone = len(data['imgMicPath'])

            # check if we generated any new images in this round
            if micsReady < self.thumbsDone:
                reportFinished = False

            if len(data['defocusU']) < 100:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'])
            else:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'][:-50])
                data['defocusCoverageLast50'] = self.processDefocusValues(data['defocusU'][-50:])

            # remove data of original pictures, report doesn't need it
            for k in ['imgMicPath', 'imgPsdPath', 'imgShiftPath']:
                data.pop(k)
            # send over only thumbnails of the mics that have been fully processed
            for k in ['imgMicThumbs', 'imgPsdThumbs', 'imgShiftThumbs']:
                data[k] = data[k][:micsReady]

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
