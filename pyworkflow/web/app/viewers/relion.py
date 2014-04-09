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
from pyworkflow.em.packages.relion.viewer import *

#===============================================================================
# doShowImagesInClasses
#===============================================================================

def doShowImagesInClasses(request, protocolViewer):
    paths = protocolViewer._createImagesInClasses()
    urls = buildShowjPath(paths)
    return "showjs", urls

#===============================================================================
# doShowLLRelion
#===============================================================================

def doShowLLRelion(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    
    # Need to set iterations
    protocolViewer._load()
    
    urls = []

    for it in protocolViewer._iterations:
        
        # FILES
        path = protocolViewer._getFilesPerIterationLL(it)
        url_showj = "/visualize_object/?path="+ path
        # PLOTTERS
        url_plot = "/view_plots/?function=plotShowLLRelion&protViewerClass="+ protViewerClass + "&path="+ path + "&iter="+ str(it) + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
        
        # Add to the url list
        urls.append(str(url_showj))
        urls.append(str(url_plot))
        
    return "urls", urls

def plotShowLLRelion(request, protocolViewer):
    iter = request.GET.get('iter', None)
    path = request.GET.get('path', None)
    xplotter = protocolViewer._createPlotPerIterationLL(int(iter), path)
    return xplotter

#===============================================================================
# doShowPMax
#===============================================================================

def doShowPMax(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    
    # Need to set iterations
    protocolViewer._load()

    # FILE
    path = protocolViewer._createFilePMax()
    url_showj = "/visualize_object/?path="+ path
    # PLOT
    url_plot = "/view_plots/?function=plotShowPMax&protViewerClass="+ protViewerClass + "&path="+ path + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
    
    return "urls", [url_showj, url_plot]

def plotShowPMax(request, protocolViewer):
    path = request.GET.get('path', None)
    xplotter = protocolViewer._createPlotPMax(path)
    return xplotter


#===============================================================================
# doShowChanges        
#===============================================================================

def doShowChanges(request, protocolViewer):
    filename = protocolViewer._createChanges()
    return "showj","/visualize_object/?path="+ filename
    
#===============================================================================
# doShowVolumes
#===============================================================================
    
def doShowVolumes(request, protocolViewer):
    files = protocolViewer._createVolumes()
    urls = buildShowjPath(files)
    return "showjs", urls

#===============================================================================
# doShowAngularDistributionRelion
#===============================================================================

def doShowAngularDistributionRelion(request, protocolViewer):
    plotters, arguments = protocolViewer._createAngularDistribution()
    return "",""

#===============================================================================
# doPlotsSSNR
#===============================================================================

def doPlotsSSNR(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    
    # PLOT
    url_plot = "/view_plots/?function=plotSSNR&protViewerClass="+ protViewerClass + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
    
    return "plot", url_plot

def plotSSNR(request, protocolViewer):
    xplotter = protocolViewer._createSSNR()
    return xplotter

#===============================================================================
# doPlotsFSC
#===============================================================================

def doPlotsFSC(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    
    # PLOT
    url_plot = "/view_plots/?function=plotFSC&protViewerClass="+ protViewerClass + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
    
    return "plot", url_plot

def plotFSC(request, protocolViewer):
    xplotter = protocolViewer._createFSC()
    return xplotter




