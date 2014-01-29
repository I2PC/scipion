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
from pyworkflow.em.packages.brandeis.viewer_frealign import setVisualizeIterations

LAST_ITER = 0
ALL_ITER = 1
SELECTED_ITERS = 2

def viewerFrealign(request, protocolViewer):
    ioDict = {}
    
    if protocolViewer.doShow3DRefsVolumes:
        typeUrl, url = view3DRefsVolumes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShow3DReconVolumes:
        typeUrl, url = view3DReconVolumes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowAngDist:
        typeUrl, url = doPlotAngularDistribution(request, protocolViewer)
        ioDict[typeUrl]= url
        
    return ioDict

def doShow3DRefsVolumes(request, protocolViewer):
    path = str("reference_volume_iter_%03d.mrc")
    sourcePath = self._viewIterationFile(path)
    return "showj", "/visualize_object/?path="+ sourcePath

def doShow3DReconVolumes(request, protocolViewer):
    path = str("volume_iter_%03d.mrc")
    sourcePath = self._viewIterationFile(path)
    return "showj", "/visualize_object/?path="+ sourcePath

def doPlotAngularDistribution(request, protocolViewer):
    iterToShow = str(protocolViewer.iterToShow.get())
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(-1)
    return "plot","/view_plots/?function=plotAngularDistribution&protViewerClass="+ protViewerClass + "&protId="+ protId +"&iterToShow="+ iterToShow + "&width=" + str(width) + "&height="+ str(height)

def plotAngularDistribution(request, protocolViewer):
    protocolViewer.iterToShow.set(request.GET.get('iterToShow', None))
    plots, errors = protocolViewer._createAngularDistributionPlots()
    
#    xplotter = plots[0]

    if len(errors) != 0:
        pass
    else:
        return plots

def _viewIterationFile(self, filePath):
    self.setVisualizeIterations()
    for iter in self.visualizeIters:
        pathDir = self.protocol._getExtraPath("iter_%03d" % iter)
        path = join(pathDir, filePath % iter)
        if os.path.exists(path):
            return path       



