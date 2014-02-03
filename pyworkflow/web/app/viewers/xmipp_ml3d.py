# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
from pyworkflow.web.app.views_util import *
from pyworkflow.em.packages.xmipp3.protocol_ml3d import *


def viewerML3D(request, protocolViewer):
    ioDict = {}
    
    if protocolViewer.doShowGreyScaleRefVol:
        typeUrl, url = viewCorrectedVols(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowFilterRefVol:
        typeUrl, url = viewFilteredVolsProtocol(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowSeedsVols:
        typeUrl, url = viewGeneratedVols(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShow2DAvgs:
        typeUrl, url = view2DAvgs(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShow3DRefsVolumes:
        typeUrl, url = view3DRefsVolumes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowAngDist:
        typeUrl, url = doPlotAngularDistribution(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowDataDist:
        typeUrl, url = doPlotClassDistribution(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowStatistics: 
        typeUrl, url = doPlotStatistics(request, protocolViewer)
        ioDict[typeUrl]= url
        
    return ioDict

def viewCorrectedVols(request, protocolViewer):
    path = str(protocolViewer.protocol._getExtraPath("corrected_volumes.stk"))
    return "showj", "/visualize_object/?path="+ path

def viewFilteredVolsProtocol(request, protocolViewer):
    path = str(protocolViewer.protocol._getExtraPath("filtered_volumes.stk"))
    return "showj", "/visualize_object/?path="+ path

def viewGeneratedVols(request, protocolViewer):
    path = str(protocolViewer.protocol._getExtraPath("generated_volumes.stk"))
    return "showj", "/visualize_object/?path="+ path

def view2DAvgs(request, protocolViewer):
    files = protocolViewer._getIterationFile("iter_classes.xmd", extra="mlf2dextra")
    urls = buildShowjPath(files)
    return "showjs", urls

def view3DRefsVolumes(request, protocolViewer):
    files = protocolViewer._getIterationFile("vol000001.vol", extra="extra")
    urls = buildShowjPath(files)
    return "showjs", urls

def doPlotAngularDistribution(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    protocolViewer.setVisualizeIterations()
    plots = []
    for iter in protocolViewer.visualizeIters:
        extra = '%s2d' % protocolViewer.protocol.getProgramId() + 'extra'
        path = protocolViewer.protocol._getPath(extra + "/iter%03d/iter_classes.xmd" % iter)
        url = "/view_plots/?function=plotAngularDistribution&protViewerClass="+ protViewerClass + "&path="+ path + "&iter="+ str(iter) + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
        plots.append(str(url))
    return "plots", plots

def plotAngularDistribution(request, protocolViewer):
    iter = request.GET.get('iter', None)
    path = request.GET.get('path', None)
    xplotter = protocolViewer._createIterAngularDistributionPlot(int(iter), path)
    return xplotter
    
def doPlotClassDistribution(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    protocolViewer.setVisualizeIterations()
    plots = []
    for iter in protocolViewer.visualizeIters:
        extra = '%s2d' % protocolViewer.protocol.getProgramId() + 'extra'
        path = protocolViewer.protocol._getPath(extra + "/iter%03d/iter_classes.xmd" % iter)
        url = "/view_plots/?function=plotClassDistribution&protViewerClass="+ protViewerClass + "&path="+ path + "&iter="+ str(iter) + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
        plots.append(str(url))
    return "plots", plots

def plotClassDistribution(request, protocolViewer):
    iter = request.GET.get('iter', None)
    path = request.GET.get('path', None)
    xplotter = protocolViewer._createIterClassDistributionPlots(int(iter), path)
    return xplotter

def doPlotStatistics(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(2)
    return "plot","/view_plots/?function=plotStatistics&protViewerClass="+ protViewerClass + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)

def plotStatistics(request, protocolViewer):
    xplotter = protocolViewer._plotStatistics()
    return xplotter

            