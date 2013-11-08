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
# *  e-mail address 'jose.gutierrez@cnb.csic.es'
# *
# **************************************************************************

from views_util import *

def viewerWARD(request, protocolViewer):
    ioDict = {}
    
    if protocolViewer.doShowClasses:
        typeUrl, url = doVisualizeClasses(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowDendrogram:
        typeUrl, url = doVisualizeDendrogram(request, protocolViewer)
        ioDict[typeUrl]= url

    return ioDict

def doVisualizeClasses(request, protocolViewer):
    pass

def doVisualizeDendrogram(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=visualizeDendrogram&protViewerClass="+ protViewerClass + "&protId="+ protId

def visualizeDendrogram(request, protocolViewer):
    xplotter = protocolViewer._plotDendrogram()
    return xplotter
