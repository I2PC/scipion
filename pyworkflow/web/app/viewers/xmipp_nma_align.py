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

from pyworkflow.em.packages.xmipp3.viewer_nma_alignment import *
from pyworkflow.web.app.views_util import *

def viewerNMAAlign(request, protocolViewer):
    ioDict = {}
    if protocolViewer.displayRawDeformation:
        typeUrl, url = doDisplayRawDeformation(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.analyzeMatlab:
        typeUrl, url = doAnalyzeMatlab(request, protocolViewer)
        ioDict[typeUrl]= url

    return ioDict

def doDisplayRawDeformation(request, protocolViewer):
    pass

def doAnalyzeMatlab(request, protocolViewer):
    pass
