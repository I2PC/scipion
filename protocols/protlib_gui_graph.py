'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

import Tkinter as tk
import ttk
from protlib_utils import reportError
import numpy as np
from pyworkflow.drawing import *


def showDependencyTree(runsDict):
    ''' This function will create a Canvas to display
    the protocols dependency tree''' 

    XSIZE, YSIZE = 8, 6
    DPI = 100
    XDIM, YDIM = XSIZE*DPI, YSIZE*DPI
    DY = 56
    DX = 50
    FONT = "sans-serif"
    FONTSIZE = 9
    colors = ['#D9F1FA', '#D9F1FA', '#FCCE62', '#D2F5CB', '#F5CCCB', '#F3F5CB', '#416FF0']
    # Store the right horizontal x position for better packaging of
    # the graph, assuming max of 100 y-levels
    from numpy import zeros
    hLimits = zeros(100) + DX/2
    
    root = tk.Toplevel()
    root.withdraw()
    root.title("Dependency tree")
    root.columnconfigure(0, weight=1, minsize=800)
    root.rowconfigure(0, weight=1, minsize=500)
    
    parent = root
    canvas = Canvas(parent)
    canvas.grid(row=0, column=0, sticky='nsew')
    parent.grid_columnconfigure(0, weight=1)
    parent.grid_rowconfigure(0, weight=1)
    
    def showNode(dd, x, y):
        if dd.prot is None:
            nodeText = dd.extRunName
        else:
            nodeText = "%s\n%s" % (dd.protName, dd.runName)
        
        t = canvas.createTextbox(nodeText, x, y, bgColor=colors[dd.state])
        dd.width, dd.height = t.getDimensions()
        dd.start = x
        dd.t = t
        
        return t
        
    def showLevel(dd, level):
        y = level * DY
        
        if len(dd.deps):
            #width = (xmax - xmin) / n
            childs = [runsDict[rn] for rn in dd.deps]
            for c in childs:
                showLevel(c, level + 1)
                
            firstChild = childs[0]
            lastChild = childs[-1]
            
            t = showNode(dd, 0, y)
            dd.start = (lastChild.start + lastChild.width + firstChild.start - dd.width) / 2
            dd.start = max(dd.start, hLimits[level] + DX)
            t.moveTo(dd.start, y)            
            hLimits[level] = dd.start + dd.width     
            for c in childs:
                canvas.createEdge(t, c.t)
                #canvas.createEdge(src, dst)
                #a.add_line(mlines.Line2D(xx, yy, lw=2., alpha=0.4))
        else:
            t = showNode(dd, hLimits[level] + DX, y)
            half = dd.width / 2
            #dd.hLimits = [(-half, half)]
            hLimits[level] = dd.start + dd.width
    rootNode = runsDict['PROJECT']
    showLevel(rootNode, 1)  
    
    from protlib_gui_ext import ToolTip, centerWindows
    centerWindows(root)
    root.deiconify()
    root.mainloop()  
    #xplotter.show()



