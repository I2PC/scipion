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

from pyworkflow.web.app.views_util import *

def viewerCAPCA(request, protocolViewer):
    ioDict = {}
    # SHOWJS
    if protocolViewer.doShowEigenImages and protocolViewer.doShowReconsImages:
        typeUrl, url = doShowImagesCAPCA(request, protocolViewer)
        ioDict[typeUrl]= url
    else:   
        if protocolViewer.doShowEigenImages:
            typeUrl, url = doShowEigenImages(request, protocolViewer)
            ioDict[typeUrl]= url
        elif protocolViewer.doShowReconsImages:
            typeUrl, url = doShowReconsImages(request, protocolViewer)
            ioDict[typeUrl]= url
    # PLOTS
    if protocolViewer.doShowHistogram and protocolViewer.doShowFactorMaps:
        typeUrl, url = doPlotsCAPCA(request, protocolViewer)
        ioDict[typeUrl]= url
    else:
        if protocolViewer.doShowHistogram:
            typeUrl, url = doPlotsCAPCA(request, protocolViewer)
            ioDict[typeUrl]= url
        elif protocolViewer.doShowFactorMaps:
            typeUrl, url = doPlotFactorMaps(request, protocolViewer)
            ioDict[typeUrl]= url
    # FILE VIEWER
    if protocolViewer.doShowPcaFile:
        typeUrl, url = doShowPcaFile(request, protocolViewer)
        ioDict[typeUrl]= url
        
    return ioDict

def doShowImagesCAPCA(request, protocolViewer):
    _, eigenUrl = doShowEigenImages(request, protocolViewer)
    _, reconsUrl = doShowReconsImages(request, protocolViewer)
    return "showjs", [str(eigenUrl) , str(reconsUrl)]

def doShowEigenImages(request, protocolViewer):
    path = protocolViewer.protocol._getFileName('eigenimages')
    return "showj", "/visualize_object/?path="+ path

def doShowReconsImages(request, protocolViewer):
    path = protocolViewer.protocol._getFileName('reconstituted')
    return "showj", "/visualize_object/?path="+ path

def doPlotsCAPCA(request, protocolViewer):
    _, histogram = doPlotHistogram(request, protocolViewer)
    _, factorMaps = doPlotFactorMaps(request, protocolViewer)
    return "plots", [str(histogram) , str(factorMaps)]

def doPlotHistogram(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=plotHistogram&protViewerClass="+ protViewerClass + "&protId="+ protId

def plotHistogram(request, protocolViewer):
    xplotter = protocolViewer._plotHistogram()
    return xplotter

def doPlotFactorMaps(request, protocolViewer):
    first = str(protocolViewer.firstFactor.get())
    second = str(protocolViewer.secondFactor.get())
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=plotFactorMaps&protViewerClass="+ protViewerClass + "&protId="+ protId +"&first="+ first +"&second="+second 

def plotFactorMaps(request, protocolViewer):
    protocolViewer.firstFactor.set(request.GET.get('first', None))
    protocolViewer.secondFactor.set(request.GET.get('second', None))
    xplotter = protocolViewer._plotFactorMaps()
    return xplotter

def doShowPcaFile(request, protocolViewer):
    html = textfileViewer("PCA files", [protocolViewer.protocol.imcFile.filename.get()])
    return "html", html
