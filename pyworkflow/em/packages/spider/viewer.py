# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
This module implement the wrappers around xmipp_showj
visualization program.
"""
import Tkinter as tk
from pyworkflow.protocol.params import *
from pyworkflow.viewer import Viewer, ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.utils.graph import Graph
from pyworkflow.gui.graph import LevelTree
from pyworkflow.gui.canvas import Canvas, ImageBox
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer, runShowJ
from pyworkflow.gui.text import showTextfileViewer

from spider import PcaFile
from protocol_filters import SpiderProtFilter
from protocol_custommask import SpiderProtCustomMask
from protocol_align_apsr import SpiderProtAlignAPSR

    
class SpiderViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """
    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PcaFile, SpiderProtCustomMask, SpiderProtFilter]
    _label = 'viewer'

    def visualize(self, obj, **args):
        
        if isinstance(obj, PcaFile):
            showTextfileViewer("PCA file", [obj.filename.get()])
        
        elif isinstance(obj, SpiderProtFilter):
            XmippViewer().visualize(obj.outputParticles)
            
        elif isinstance(obj, SpiderProtCustomMask):
            mask = obj.outputMask
            XmippViewer().visualize(mask)
            # Remove location to visualize the whole stack
            runShowJ(mask.getFileName())
            
    @classmethod
    def getView(self):
        return "viewerSpider"
            
