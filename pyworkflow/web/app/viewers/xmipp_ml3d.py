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

LAST_ITER = 0
ALL_ITER = 1
SELECTED_ITERS = 2

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
    extra = '%s2d' % protocolViewer.protocol.getProgramId() + 'extra'
    file = extra + "/iter%03d/iter_classes.xmd"
#    viewIterationFile("ml2dextra/iter%03d/iter_classes.xmd")
    pass

def view3DRefsVolumes(request, protocolViewer):
    file = "extra/iter%03d/vol000001.vol"
#    viewIterationFile("extra/iter%03d/vol000001.vol")
    pass



def doPlotAngularDistribution(request, protocolViewer):
    iterToShow = str(protocolViewer.iterToShow.get())
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    return "plot","/view_plots/?function=plotAngularDistribution&protViewerClass="+ protViewerClass + "&protId="+ protId +"&iterToShow="+ iterToShow + "&width=" + str(width) + "&height="+ str(height)

def plotAngularDistribution(request, protocolViewer):
    protocolViewer.iterToShow.set(request.GET.get('iterToShow', None))
    plots, errors = protocolViewer._createAngularDistributionPlots()
    if len(errors) != 0:
        pass
    else:
        return plots
    
def doPlotClassDistribution(request, protocolViewer):
    iterToShow = str(protocolViewer.iterToShow.get())
    protocolViewer.iterToShow.set(request.GET.get('iterToShow', None))
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    return "plot","/view_plots/?function=plotClassDistribution&protViewerClass="+ protViewerClass + "&protId="+ protId +"&iterToShow="+ iterToShow + "&width=" + str(width) + "&height="+ str(height)

def plotClassDistribution(request, protocolViewer):
    protocolViewer.iterToShow.set(request.GET.get('iterToShow', None))
    plots, errors = protocolViewer._createClassDistributionPlots()
    xplotter = plots[0]
    if len(errors) != 0:
        pass
    else:
        return plots

def doPlotStatistics(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    return "plot","/view_plots/?function=plotStatistics&protViewerClass="+ protViewerClass + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)

def plotStatistics(request, protocolViewer):
    xplotter = protocolViewer._plotStatistics()
    return xplotter

            