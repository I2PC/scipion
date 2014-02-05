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
        typeUrl, url = doShowDataDist(request, protocolViewer)
        ioDict[typeUrl]= url
    
    return ioDict

def doShow3DRefsVolumes(request, protocolViewer):
    files = protocolViewer._getIterationFile("reference_volume_iter_%03d.mrc")
    urls = buildShowjPath(files)
    return "showjs", urls


def doShow3DReconsVolumes(request, protocolViewer):
    files = protocolViewer._getIterationFile("volume_iter_%03d.mrc")
    urls = buildShowjPath(files)
    return "showjs", urls


def doShow3DMatchProj(request, protocolViewer):
    files = protocolViewer._getIterationFile("particles_match_iter_%03d.mrc")
    urls = buildShowjPath(files)
    return "showjs", urls

def doShowAngDist(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
#    iter = str(protocolViewer.iterToShow.get())
    protocolViewer.setVisualizeIterations()
    plots = []
    for iter in protocolViewer.visualizeIters:
        path = protocolViewer.protocol._getExtraPath("iter_%03d" % iter, "particles_iter_%03d.par" % iter)
        url = "/view_plots/?function=plotShowAngDist&protViewerClass="+ protViewerClass + "&path="+ path + "&iter="+ str(iter) + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)
        plots.append(str(url))
    return "plots", plots

def plotShowAngDist(request, protocolViewer):
    iter = request.GET.get('iter', None)
    path = request.GET.get('path', None)
    xplotter = protocolViewer._createIterAngularDistributionPlot(int(iter), path)
    return xplotter

def doShowDataDist(request, protocolViewer):
    pass





