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
from pyworkflow.em.packages.brandeis.viewer_frealign import *


def viewerFrealign(request, protocolViewer):
    ioDict = {}
    
    if protocolViewer.iterToShow:
        pass
    if protocolViewer.selectedIters:
        pass
    if protocolViewer.doShow3DRefsVolumes:
        typeUrl, url = doShow3DRefsVolumes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShow3DReconsVolumes:
        typeUrl, url = doShow3DReconsVolumes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShow3DMatchProj:
        typeUrl, url = doShow3DMatchProj(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowAngDist:
        typeUrl, url = doShowAngDist(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowDataDist:
        pass
    
    return ioDict

def doShow3DRefsVolumes(request, protocolViewer):
    files = protocolViewer._getIterationFile("reference_volume_iter_%03d.mrc")
    urls = buildPath(files)
    return "showjs", urls


def doShow3DReconsVolumes(request, protocolViewer):
    files = protocolViewer._getIterationFile("volume_iter_%03d.mrc")
    urls = buildPath(files)
    return "showjs", urls


def doShow3DMatchProj(request, protocolViewer):
    files = protocolViewer._getIterationFile("particles_match_iter_%03d.mrc")
    urls = buildPath(files)
    return "showjs", urls

def doShowAngDist(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    iter = str(protocolViewer.iterToShow.get())
    return "plot","/view_plots/?function=plotShowAngDist&protViewerClass="+ protViewerClass + "&iter="+ iter + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)

#    return "showj","/visualize_object/?path="+ path

def plotShowAngDist(request, protocolViewer):
    iter = request.GET.get('iter', None)
    protocolViewer.iterToShow.set(iter)
    xplotter, errors = protocolViewer._createAngularDistributionPlots()
    print xplotter
    
    return xplotter

def buildPath(files):
    urls = []
    for f in files:
        url = "/visualize_object/?path="+ f
        urls.append(url)
        
    return urls


