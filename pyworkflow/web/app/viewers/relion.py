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

def viewerRelion(request, protocolViewer):
    
    ioDict = {}
    
    if protocolViewer.showImagesInClasses:
#        typeUrl, url = doShowImagesInClasses(request, protocolViewer)
#        ioDict[typeUrl]= url
        pass
    if protocolViewer.showLL:
        typeUrl, url = doShowLL(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.showPMax:
        typeUrl, url = doShowPMax(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.showChanges:
        typeUrl, url = doShowChanges(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.displayVol:
        typeUrl, url = doShowVolumes(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.displayAngDist:
        typeUrl, url = doShowAngularDistribution(request, protocolViewer)
        ioDict[typeUrl]= url
    
    return ioDict

def doShowImagesInClasses(request, protocolViewer):
    pass

def doShowLL(request, protocolViewer):
    plotters, files = protocolViewer._createLL()
    pass

def doShowPMax(request, protocolViewer):
    xplotter, fn = protocolViewer._createPMax()
    pass

def doShowChanges(request, protocolViewer):
    fn = protocolViewer._createChanges()
    pass
#    urls = buildShowjPath(fn)
#    return "showjs", urls
    
def doShowVolumes(request, protocolViewer):
    files = protocolViewer._createVolumes()
    pass
#    return xplotter

def doShowAngularDistribution(request, protocolViewer):
    plotters, arguments = protocolViewer._createAngularDistribution()
    pass


