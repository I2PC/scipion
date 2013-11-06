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

from pyworkflow.viewer import createPlots

############## VIEWER XMIPP ML2D ##############
def viewerML2D(request, protocolViewer):
    ioDict = {}

    if protocolViewer.doShowClasses:
        typeUrl, url = doShowClasses(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowPlots:
        typeUrl, url = doAllPlotsML2D(request, protocolViewer)
        ioDict[typeUrl]= url
    else:
        typeUrl, url = doSomePlotsML2D(protocolViewer)
        if url != None:
            ioDict[typeUrl]= url
        
    return ioDict

def doShowClasses(request, protocolViewer):
    objId = str(protocolViewer.protocol.outputClasses.getObjId())
    return "showj","/visualize_object/?objectId="+ objId
            
def doAllPlotsML2D(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plotsComposite","/view_plots/?function=allPlotsML2D&protViewerClass="+ protViewerClass + "&protId="+ protId

def doSomePlotsML2D(protocolViewer):
    plots=""
    for p in protocolViewer._plotVars:
        if protocolViewer.getAttributeValue(p):
            plots = plots + p + "-"
    
    if plots != "":
        return "plotsComposite","/view_plots/?function=somePlotsML2D&plots="+str(plots)+"&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocolViewer.protocol.getObjId())
    else:
        return "", None

def doShowLL(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowLL&protViewerClass="+ protViewerClass + "&protId="+ protId

def doShowPmax(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowPmax&protViewerClass="+ protViewerClass + "&protId="+ protId

def doShowSignalChange(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowSignalChange&protViewerClass="+ protViewerClass + "&protId="+ protId
    
def doShowMirror(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowMirror&protViewerClass="+ protViewerClass + "&protId="+ protId

def allPlotsML2D(request, protocolViewer):
    xplotter = createPlots(protocolViewer.protocol, protocolViewer._plotVars)
    return xplotter

def somePlotsML2D(request, protocolViewer):
    plots = request.GET.get('plots', None)
    plots = str(plots).split("-")
    if len(plots) > 1:
        plots.remove(plots[len(plots)-1])
    
    xplotter = createPlots(protocolViewer.protocol, plots)
    return xplotter