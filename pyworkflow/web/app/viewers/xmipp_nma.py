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

from pyworkflow.em.packages.xmipp3.viewer_nma import *
from pyworkflow.web.app.views_util import *

def viewerNMA(request, protocolViewer):
    ioDict = {}
    if protocolViewer.isEm:
        if protocolViewer.displayPseudoAtom:
            #CHIMERA
#            typeUrl, url = doDisplayPseudoAtom(request, protocolViewer)
#            ioDict[typeUrl]= url
            pass
        if protocolViewer.displayPseudoAtomAproximation:
            typeUrl, url = doDisplayPseudoAtomAproximation(request, protocolViewer)
            ioDict[typeUrl]= url
    if protocolViewer.displayModes:
        typeUrl, url = doDisplayModes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.displayMaxDistanceProfile:
        typeUrl, url = doDisplayMaxDistanceProfile(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.displayDistanceProfile:
        typeUrl, url = doDisplayDistanceProfile(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.singleMode.hasValue():
        #VMD
#       typeUrl, url = doSingleMode(request, protocolViewer)
#       ioDict[typeUrl]= url
        pass
        
    return ioDict

def doDisplayPseudoAtom(request, protocolViewer):
#    chimera(self.protocol._getPath("chimera.cmd"))
    pass

def doDisplayPseudoAtomAproximation(request, protocolViewer):
    file = protocolViewer.protocol.inputStructure.get().getFirstItem().getFileName()
    extra = protocolViewer.protocol._getExtraPath('pseudoatoms_approximation.vol')

    url1 = "/visualize_object/?path="+ file
    url2 = "/visualize_object/?path="+ extra
    
    return "showjs", [str(url1) , str(url2)]

def doDisplayModes(request, protocolViewer):
    path = str(protocolViewer.protocol._getPath('modes.xmd'))
    return "showj","/visualize_object/?path="+ path

def doDisplayMaxDistanceProfile(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    return "plot","/view_plots/?function=plotMaxDistanceProfile&protViewerClass="+ protViewerClass + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)

def plotMaxDistanceProfile(request, protocolViewer):
    fn = str(protocolViewer.protocol._getExtraPath("maxAtomShifts.xmd"))
    xplotter = protocolViewer._createShiftPlot(fn, "Maximum atom shifts", "maximum shift")
    return xplotter

def doDisplayDistanceProfile(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    width, height = getSizePlotter(1)
    mode = str(protocolViewer.singleMode.get())
    return "plot","/view_plots/?function=plotDistanceProfile&protViewerClass="+ protViewerClass + "&mode="+ mode + "&protId="+ protId + "&width=" + str(width) + "&height="+ str(height)

def plotDistanceProfile(request, protocolViewer):
    mode = int(request.GET.get('mode', None))
    fn = str(protocolViewer.protocol._getExtraPath("distanceProfiles","vec%d.xmd" % mode))
    xplotter = protocolViewer._createShiftPlot(fn, "Atom shifts for mode %d" % mode, "shift")
    return xplotter

def doSingleMode(request, protocolViewer):
#    vmdFile = self.protocol._getExtraPath("animations", "animated_mode_%03d.vmd" % self.singleMode.get())
#    os.system("vmd -e %s" % vmdFile)
    pass


