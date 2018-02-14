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
import numpy as np
import multiprocessing
from datetime import datetime

from pyworkflow.protocol import getUpdatedProtocol
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
# used in the execution.summary.template.html to read data (where
# they will need to be changed if they're changed here)
MIC_THUMBS = 'imgMicThumbs'
PSD_THUMBS = 'imgPsdThumbs'
SHIFT_THUMBS = 'imgShiftThumbs'
MIC_ID = 'micId'
DEFOCUS_HIST_BIN_WIDTH = 0.5
RESOLUTION_HIST_BIN_WIDTH = 0.5


class ReportHtml:
    """ Create an html report with a summary of the processing.
    The report will be updated with a given frequency.
    """
    def __init__(self, protocol, ctfMonitor, sysMonitor, movieGainMonitor, publishCmd=None, **kwargs):
        # The CTF protocol to monitor
        self.protocol = protocol
        self.ctfProtocol = protocol._getCtfProtocol()
        self.alignProtocol = protocol._getAlignProtocol()
        self.micThumbSymlinks = False
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
                           SHIFT_THUMBS: [],
                           MIC_PATH: [],
                           SHIFT_PATH: [],
                           PSD_PATH: [],
                           MIC_ID: []}

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
        thumbKeys = [MIC_THUMBS]
        if PSD_THUMBS in self.thumbPaths.keys():
            lastThumb = min(len(self.thumbPaths[PSD_THUMBS]), len(self.thumbPaths[MIC_THUMBS]))
            thumbKeys.append(PSD_THUMBS)
        else:
            lastThumb = len(self.thumbPaths[MIC_THUMBS])
        if SHIFT_THUMBS in self.thumbPaths.keys():
            thumbKeys.append(SHIFT_THUMBS)
        for i in range(self.thumbsReady, lastThumb):
            pathsReady = [exists(join(self.reportDir, self.thumbPaths[k][i])) for k in thumbKeys]
            if all(pathsReady):
                self.thumbsReady += 1
        return self.thumbsReady

    def setUp(self):
        """Setup actions for the html report: create directories to store
        thumbnails, check if thumbnails are already generated in the
        alignment protocol.
        """
        # make report folders
        pwutils.makePath(join(self.reportDir, MIC_THUMBS),
                         join(self.reportDir, PSD_THUMBS),
                         join(self.reportDir, SHIFT_THUMBS))
        # check if align protocol already has thumbnails
        if (hasattr(self.alignProtocol, 'doComputeMicThumbnail')
            and self.alignProtocol._doComputeMicThumbnail()):
            self.micThumbSymlinks = True

    def getThumbPaths(self, thumbsDone=0, ctfData=None, ext='png', micIdSet=None):
        """Adds to self.thumbPaths the paths to the report thumbnails
           that come from the alignment and/or ctf protocol.

            ===== Params =====
            - thumbsDone: how many thumbnails have already been generated.
                          we will get paths starting from this index
            - ctfData: dict resulting from ctfMonitor.getData()
            - ext: extension of the thumbnail images. Defaults to png.
            - micIdSet: mic indexes to use
        """
        # get psd thumbs from ctfData
        if ctfData is not None:
            for i in range(thumbsDone, len(ctfData[PSD_PATH])):
                psdPath = ctfData[PSD_PATH][i]
                movie = basename(os.path.dirname(psdPath))
                psdThumb = join(PSD_THUMBS, "%s_%s" % (movie, pwutils.replaceExt(basename(psdPath), ext)))
                self.thumbPaths[PSD_THUMBS].append(psdThumb)
                self.thumbPaths[PSD_PATH].append(psdPath)

        # get alignment and mic thumbs
        if self.alignProtocol is not None:
            getMicFromCTF = False
            updatedProt = getUpdatedProtocol(self.alignProtocol)
            if hasattr(updatedProt, 'outputMicrographs'):
                outputSet = updatedProt.outputMicrographs
                if micIdSet is None:
                    micIdSet = list(outputSet.getIdSet())
            else:
                return

        elif self.ctfProtocol is not None:
            getMicFromCTF = True
            updatedProt = getUpdatedProtocol(self.ctfProtocol)
            if hasattr(updatedProt, 'outputCTF'):
                outputSet = updatedProt.outputCTF
                if micIdSet is None:
                    micIdSet = list(outputSet.getIdSet())
            else:
                return
        else:
            return

        for micId in micIdSet[thumbsDone:]:
            mic = outputSet[micId]
            if getMicFromCTF:
                mic = mic.getMicrograph()
            if hasattr(mic, 'thumbnail'):
                srcMicFn = abspath(mic.thumbnail.getFileName())
            else:
                srcMicFn = abspath(mic.getFileName())
            micThumbFn = join(MIC_THUMBS, pwutils.replaceExt(basename(srcMicFn), ext))
            self.thumbPaths[MIC_PATH].append(srcMicFn)
            self.thumbPaths[MIC_THUMBS].append(micThumbFn)

            shiftPlot = (getattr(mic, 'plotCart', None) or getattr(mic, 'plotGlobal', None))
            if shiftPlot is not None:
                shiftPath = "" if shiftPlot is None else abspath(shiftPlot.getFileName())
                shiftCopy = "" if shiftPlot is None else join(SHIFT_THUMBS, pwutils.replaceExt(basename(shiftPath), ext))
                self.thumbPaths[SHIFT_PATH].append(shiftPath)
                self.thumbPaths[SHIFT_THUMBS].append(shiftCopy)
            else:
                if SHIFT_PATH in self.thumbPaths:
                    self.thumbPaths.pop(SHIFT_PATH, None)
                if SHIFT_THUMBS in self.thumbPaths:
                    self.thumbPaths.pop(SHIFT_THUMBS, None)

            self.thumbPaths[MIC_ID].append(micId)

            if self.ctfProtocol is None:
                psdPath = mic.psdJpeg.getFileName() if hasattr(mic, 'psdJpeg') else mic.psdCorr.getFileName()
                psdThumb = join(PSD_THUMBS, pwutils.replaceExt(basename(psdPath), ext))
                self.thumbPaths[PSD_THUMBS].append(psdThumb)
                self.thumbPaths[PSD_PATH].append(psdPath)

    def generateReportImages(self, firstThumbIndex=0, micScaleFactor=6):
        """ Function to generate thumbnails for the report. Uses data from
        self.thumbPaths.

        ===== Params =====
        - firstThumbIndex: index from which we start generating thumbnails
        - micScaleFactor: how much to reduce in size the micrographs.

        """
        ih = ImageHandler()

        numMics = len(self.thumbPaths[MIC_PATH])

        for i in range(firstThumbIndex, numMics):
            print('Generating images for mic %d' % (i+1))
            # mic thumbnails
            dstImgPath = join(self.reportDir, self.thumbPaths[MIC_THUMBS][i])
            if not exists(dstImgPath):
                if self.micThumbSymlinks:
                    pwutils.createAbsLink(self.thumbPaths[MIC_PATH][i], dstImgPath)
                else:
                    ih.computeThumbnail(self.thumbPaths[MIC_PATH][i], dstImgPath, scaleFactor=micScaleFactor)

            # shift plots
            if SHIFT_THUMBS in self.thumbPaths:
                dstImgPath = join(self.reportDir, self.thumbPaths[SHIFT_THUMBS][i])
                if not exists(dstImgPath):
                    pwutils.createAbsLink(self.thumbPaths[SHIFT_PATH][i], dstImgPath)

            # Psd thumbnails
            if self.ctfProtocol is None:
                srcImgPath = self.thumbPaths[PSD_PATH][i]
                dstImgPath = join(self.reportDir, self.thumbPaths[PSD_THUMBS][i])
                if not exists(dstImgPath):
                    if srcImgPath.endswith('psd'):
                        psdImg1 = ih.read(srcImgPath)
                        psdImg1.convertPSD()
                        psdImg1.write(dstImgPath)
                        ih.computeThumbnail(dstImgPath, dstImgPath, scaleFactor=1)
                    else:
                        pwutils.createAbsLink(srcImgPath, dstImgPath)
            else:
                dstImgPath = join(self.reportDir, self.thumbPaths[PSD_THUMBS][i])
                if not exists(dstImgPath):
                    ih.computeThumbnail(self.thumbPaths[PSD_PATH][i], dstImgPath, scaleFactor=1)

        return

    def processDefocusValues(self, defocusList):
        maxDefocus = self.protocol.maxDefocus.get()*1e-4
        minDefocus = self.protocol.minDefocus.get()*1e-4
        # Convert defocus values to microns
        defocusList = [i*1e-4 for i in defocusList]
        edges = np.arange(0, maxDefocus+DEFOCUS_HIST_BIN_WIDTH, DEFOCUS_HIST_BIN_WIDTH)
        edges = np.insert(edges[edges > minDefocus], 0, minDefocus)
        values, binEdges = np.histogram(defocusList, bins=edges, range=(minDefocus, maxDefocus))
        belowThresh = 0
        aboveThresh = 0
        labels = ["%0.2f-%0.2f" % (x[0], x[1]) for x in zip(binEdges, binEdges[1:])]
        for v in defocusList:
            if v < minDefocus:
                belowThresh += 1
            elif v > maxDefocus:
                aboveThresh += 1
        zipped = zip(values, labels)
        zipped[:0] = [(belowThresh, "0-%0.2f" % (minDefocus))]
        zipped.append((aboveThresh, "> %0.2f" % (maxDefocus)))

        return zipped

    def getResolutionHistogram(self, resolutionValues):
        if len(resolutionValues) == 0:
            return []
        maxValue = int(np.ceil(max(resolutionValues)))
        edges = np.append(np.arange(0, maxValue, RESOLUTION_HIST_BIN_WIDTH), maxValue)
        values, binEdges = np.histogram(resolutionValues, bins=edges, range=(0, maxValue))
        return zip(values, binEdges)

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
        data = {} if self.ctfMonitor is None else self.ctfMonitor.getData()

        if data:
            numMicsDone = len(self.thumbPaths[PSD_THUMBS])
            numMics = len(data[PSD_PATH])
            numMicsToDo = numMics - numMicsDone
            self.getThumbPaths(ctfData=data, thumbsDone=numMicsDone, micIdSet=data['idValues'])

            if len(data['defocusU']) < 100:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'])
            else:
                data['defocusCoverage'] = self.processDefocusValues(data['defocusU'][:-50])
                data['defocusCoverageLast50'] = self.processDefocusValues(data['defocusU'][-50:])

            data['resolutionHistogram'] = self.getResolutionHistogram(data['resolution'])

        else:
            # Thumbnails for Micrograph Table
            numMicsDone = len(self.thumbPaths[MIC_THUMBS])
            self.getThumbPaths(thumbsDone=numMicsDone)
            numMics = len(self.thumbPaths[MIC_PATH])
            numMicsToDo = numMics - numMicsDone

        if numMicsToDo <= 10:
            # we have few new images, eg streaming mode, generate thumbnails now
            self.generateReportImages(firstThumbIndex=numMicsDone)
        else:
            # we have many images, generate thumbs in a separate process
            process = multiprocessing.Process(target=self.generateReportImages,
                                              args=(numMicsDone,))
            process.start()

        # send over only thumbnails of the mics that have been fully processed
        self.thumbsReady = self.checkNewThumbsReady()
        thumbsLoading = numMics - self.thumbsReady
        for k in [MIC_THUMBS, SHIFT_THUMBS, PSD_THUMBS]:
            if k in self.thumbPaths:
                data[k] = self.thumbPaths[k][:self.thumbsReady] + ['']*thumbsLoading

        data[MIC_ID] = self.thumbPaths[MIC_ID]

        reportFinished = self.thumbsReady == numMics

        ctfData = json.dumps(data)

        # Movie gain monitor chart data
        data = [] if self.movieGainMonitor is None else self.movieGainMonitor.getData()

        movieGainData = json.dumps(data)

        # system monitor chart data
        data = self.sysMonitor.getData()
        systemData = json.dumps(data)
        tnow = datetime.now()
        args = {'projectName': projName,
                'startTime': pwutils.dateStr(project.getCreationTime(), secs=True),
                'dateStr': pwutils.prettyTime(dt=tnow, secs=True),
                'projectDuration': pwutils.prettyDelta(tnow-project.getCreationTime()),
                'projectStatus': "FINISHED" if finished else "RUNNING",
                'scipionVersion': os.environ['SCIPION_VERSION'],
                'acquisitionLines': acquisitionLines,
                'runLines': runLines,
                'ctfData': ctfData,
                'movieGainData': movieGainData,
                'systemData': systemData,
                'refresh': self.refreshSecs
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
