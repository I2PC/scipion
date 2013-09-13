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

# All these drawing classes should be used
# in the pyworkflow.drawing new module
class Canvas(tk.Frame):
    '''Canvas to draw some objects.
    It will really contains a Frame, a Canvas and scrollbars'''
    def __init__(self, parent, **args):
        tk.Frame.__init__(self, parent)        
        h = tk.Scrollbar(self, orient=tk.HORIZONTAL)
        v = tk.Scrollbar(self, orient=tk.VERTICAL)
        self.canvas = tk.Canvas(self, scrollregion=(0, 0, 1000, 1000), bg='white',
                                yscrollcommand=v.set, xscrollcommand=h.set)
        h['command'] = self.canvas.xview
        v['command'] = self.canvas.yview
        self.canvas.grid(column=0, row=0, sticky='nsew')
        h.grid(column=0, row=1, sticky='we')
        v.grid(column=1, row=0, sticky='ns')
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.lastItem = None # Track last item selected
        self.lastPos = (0, 0) # Track last clicked position
        self.items = {} # Keep a dictionary with high-level items
        self.onClickCallback = None
        self.onDoubleClickCallback = None
        self.onRightClickCallback = None
        
        # Add bindings
        self.canvas.bind("<Button-1>", self.onClick)
        self.canvas.bind("<Button-3>", self.onRightClick)
        self.canvas.bind("<Double-Button-1>", self.onDoubleClick)
        self.canvas.bind("<B1-Motion>", self.onDrag)
        #self.canvas.bind("<MouseWheel>", self.onScroll)
    
    def getCoordinates(self, event):
        '''Converts the events coordinates to canvas coordinates'''
        # Convert screen coordinates to canvas coordinates
        xc = self.canvas.canvasx(event.x)
        yc = self.canvas.canvasy(event.y)
        return (xc, yc)
    
    def onClick(self, event):
        xc, yc = self.getCoordinates(event)
        items = self.canvas.find_overlapping(xc - 1, yc - 1,  xc + 1, yc + 1)
        if self.lastItem:
            self.lastItem.setSelected(False)
            self.lastItem = None
        self.lastPos = (0, 0)
        for i in items:
            if i in self.items:
                self.lastItem = self.items[i]
                self.lastItem.setSelected(True)
                if self.onClickCallback:
                    self.onClickCallback(self.lastItem.text)
                self.lastPos = (xc, yc)
                break
            
    def onRightClick(self, event):
        xc, yc = self.getCoordinates(event)
        items = self.canvas.find_overlapping(xc - 1, yc - 1,  xc + 1, yc + 1)
        if self.lastItem:
            self.lastItem.setSelected(False)
            self.lastItem = None
        self.lastPos = (0, 0)
        for i in items:
            if i in self.items:
                self.lastItem = self.items[i]
                self.lastItem.setSelected(True)
                if self.onRightClickCallback:
                    event.text = self.lastItem.text
                    self.onRightClickCallback(event)
                self.lastPos = (xc, yc)
                break
    
    def onDoubleClick(self, event):
        xc, yc = self.getCoordinates(event)
        items = self.canvas.find_overlapping(xc - 1, yc - 1,  xc + 1, yc + 1)
        self.lastItem = None
        self.lastPos = (0, 0)
        for i in items:
            if i in self.items:
                self.lastItem = self.items[i]
                if self.onDoubleClickCallback:
                    self.onDoubleClickCallback(self.lastItem.text)
                self.lastPos = (xc, yc)
                break

    def onDrag(self, event):
        if self.lastItem:
            xc, yc = self.getCoordinates(event)
            self.lastItem.move(xc-self.lastPos[0], yc-self.lastPos[1])
            self.lastPos = (xc, yc)  
            
    def onScroll(self, event):
        if event.num == 5 or event.delta < 0:
            count = 1
        if event.num == 4 or event.delta > 0:
            count = -1
        self.canvas.yview("scroll", count, "units")          
        
    def createTextbox(self, text, x, y, bgColor="#99DAE8", textColor='black', font=None):
        tb = TextBox(self.canvas, text, x, y, bgColor, textColor, font)
        self.items[tb.id] = tb
        return tb
    
    def createEdge(self, src, dst):
        edge = Edge(self.canvas, src, dst)
        #self.items[edge.id] = edge
        return edge
    
    def clear(self):
        '''Clear all items from the canvas'''
        self.canvas.delete(tk.ALL)
       
        
class TextBox():
    '''This class will serve to paint and store
    rectange boxes with some text'''
    def __init__(self, canvas, text, x, y, bgColor, textColor='black', font=None):
        self.bgColor = bgColor
        self.textColor = textColor
        self.text = text
        self.font = font
        self.margin = 3
        self.canvas = canvas
        self.x = x
        self.y = y
        self.paint()
        self.listeners = []
        
    def paint(self):
        '''Paint the object in a specific position.'''
        self.id_text = self.canvas.create_text(self.x, self.y, text=self.text, 
                                               justify=tk.CENTER, fill=self.textColor, font=self.font)
        xr, yr, w, h = self.canvas.bbox(self.id_text)
        m = self.margin
        xr -= m
        yr -= m

        self.id = self.canvas.create_rectangle(xr, yr, w+m, h+m, 
                                               fill=self.bgColor)
        self.canvas.tag_raise(self.id_text)
        
    def getDimensions(self):
        x, y, x2, y2 = self.canvas.bbox(self.id)
        return (x2-x, y2-y)
        
    def move(self, dx, dy):
        '''Move TextBox to new position
        dx and y are differences to current position'''
        self.canvas.move(self.id_text, dx, dy)
        self.canvas.move(self.id, dx, dy)
        self.x += dx
        self.y += dy
        for listenerFunc in self.listeners:
            listenerFunc(dx, dy)
            
    def moveTo(self, x, y):
        '''Move TextBox to new position
        dx and dy are the new position'''
        self.move(x-self.x, y-self.y)
        
    def addPositionListener(self, listenerFunc):
        self.listeners.append(listenerFunc)
        
    def setSelected(self, value):
        bw = 1
        if value:
            bw = 2
        self.canvas.itemconfig(self.id, width=bw)
        

class Edge():
    '''Edge between two objects'''
    def __init__(self, canvas, source, dest):
        self.srcX, self.srcY = source.x, source.y
        self.dstX, self.dstY = dest.x, dest.y
        self.canvas = canvas
        source.addPositionListener(self.updateSrc)
        dest.addPositionListener(self.updateDst)
        self.id = None
        self.paint()
        
    def paint(self):
        if self.id:
            self.canvas.delete(self.id)
        self.id = self.canvas.create_line(self.srcX, self.srcY, 
                                          self.dstX, self.dstY,
                                          width=2)
        self.canvas.tag_lower(self.id)
        
    def updateSrc(self, dx, dy):
        self.srcX += dx
        self.srcY += dy
        self.paint()
        
    def updateDst(self, dx, dy):
        self.dstX += dx
        self.dstY += dy
        self.paint()

DY = 15
DX = 15
FONT = "sans-serif"
FONTSIZE = 9
colors = ['#D9F1FA', '#D9F1FA', '#FCCE62', '#D2F5CB', '#F5CCCB', '#F3F5CB', '#124EB0']


def showDependencyTree(canvas, runsDict, rootName):
    ''' This function will create a Canvas to display
    the protocols dependency tree''' 
    
#    root = tk.Toplevel()
#    root.withdraw()
#    root.title("Dependency tree")
#    root.columnconfigure(0, weight=1, minsize=800)
#    root.rowconfigure(0, weight=1, minsize=500)
#    
#    parent = root
    canvas.clear()
    #canvas.grid(row=0, column=0, sticky='nsew')
    
    def showNode(dd, y):
        if dd.prot is None:
            nodeText = dd.extRunName
        else:
            nodeText = "%s\n%s" % (dd.protName, dd.runName)
        
        textColor = 'black'
        if nodeText.startswith('Project'):
            textColor='white'
        
        from protlib_gui_ext import Fonts
        t = canvas.createTextbox(nodeText, 100, y, bgColor=colors[dd.state], textColor=textColor, font=Fonts['normal'])
        dd.width, dd.height = t.getDimensions()
        dd.half = dd.width / 2
        dd.hLimits = [[-dd.half, dd.half]]
        dd.t = t
        dd.y = y + dd.height / 2
        dd.offset = 0
        
        return t
        
    def printHLimits(dd, msg):
        print "\n=====%s========" % msg
        print " dd: %s" % dd.t.text.replace('\n', '_')
        print "  offset: %d, width: %d" % (dd.offset, dd.width)
        print "  hlimits:"
        for l, r in dd.hLimits:
            print "   [%d, %d]" % (l, r)
            
    def getHLimits(dd):
        '''This function will traverse the tree
        from dd to build the left and right profiles(hLimits)
        for each level of the tree'''
        dd.hLimits = [[-dd.half, dd.half]]
        #print "getHLimits, parent: ", dd.t.text
        for c in [runsDict[rn] for rn in dd.deps]:
            count = 1
            #printHLimits(c, " child")
            
            for l, r in c.hLimits:
                l += c.offset
                r += c.offset
                #print "  level ", count
                #print "   l, r", l, r
                if count < len(dd.hLimits):
                    #print "   dd.hLimits(A)", dd.hLimits[count]
                    if l < dd.hLimits[count][0]:
                        dd.hLimits[count][0] = l
                    if r > dd.hLimits[count][1]:
                        dd.hLimits[count][1] = r
                else:
                    dd.hLimits.append([l, r])
                    #print "   dd.hLimits(A)", dd.hLimits[count]
                #print "   dd.hLimits(B)", dd.hLimits[count]
                
                count += 1

        
    def getSeparation(child1, child2):
        '''Calcualte separation between siblings
        at each height level'''
        sep = 0 
        hL1 = child1.hLimits
        hL2 = child2.hLimits
        n1 = len(hL1)
        n2 = len(hL2)
        h = min(n1, n2)
            
        for i in range(h):
            right = hL1[i][1]
            left = hL2[i][0]            
            if left + sep < right:
                sep = right - left                
  
        return sep + DX
        
    def showLevel(dd, level, y):
        n = len(dd.deps)
        
        showNode(dd, y)
        ny = y + dd.height + DY
        
        if n > 0:
            #width = (xmax - xmin) / n
            childs = [runsDict[rn] for rn in dd.deps]
            for c in childs:
                showLevel(c, level + 1, ny)
                
            if n > 1:
                offset = 0
                for i in range(n-1):
                    sep = getSeparation(childs[i], childs[i+1])
                    offset += sep
                    c = childs[i+1]
                    c.offset = offset
                
#                print "\n\n----A------"
#                print "Parent: %s" % dd.t.text.replace('\n', '_')
#                for c in childs:
#                    print "  child: %s, width: %d, offset: %d" % (c.t.text.replace('\n', '_'), c.width, c.offset)
                
                total = childs[0].half + offset + childs[-1].half
                half = total/2
                for c in childs:
                    c.offset -= half - childs[0].half
                
#                print "\n----B------"
#                print "childs[0].half: ", childs[0].half
#                print "Parent: %s" % dd.t.text.replace('\n', '_')
#                for c in childs:
#                    print "  child: %s, width: %d, offset: %d" % (c.t.text.replace('\n', '_'), c.width, c.offset)
            else:
                childs[0].offset = 0
            getHLimits(dd)
#            print "\n=====C========"
#            print "Parent: %s" % dd.t.text.replace('\n', '_')
#            print " offset: %d, width: %d" % (dd.offset, dd.width)
#            for l, r in dd.hLimits:
#                print "[%d, %d]" % (l, r)
        return dd

    def showDD(dd, x):
        nx = x + dd.offset
        dd.t.moveTo(nx, dd.y)
        #print "dd: ", dd.t.text, " x:", nx
        childs = [runsDict[rn] for rn in dd.deps]
        for c in childs:
            showDD(c, nx)  
            canvas.createEdge(dd.t, c.t)
    
    rootNode = runsDict[rootName]
    dd = showLevel(rootNode, 1, DY)
    m = 9999
    for left, right in dd.hLimits:
        m = min(m, left)
    
    showDD(dd, -m + DY)
    
    return canvas
    
#    from protlib_gui_ext import ToolTip, centerWindows
#    centerWindows(root)
#    root.deiconify()
#    root.mainloop()  
    #xplotter.show()



