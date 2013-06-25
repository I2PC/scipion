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
This module extends the functionalities of a normal Tkinter Canvas.
The new Canvas class allows to easily display Texboxes and Edges
that can be interactively dragged and clicked.
"""

import Tkinter as tk


class Canvas(tk.Frame):
    """Canvas to draw some objects.
    It will really contains a Frame, a Canvas and scrollbars"""
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
        """Converts the events coordinates to canvas coordinates"""
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
        print "scrolling"
        if event.num == 5 or event.delta < 0:
            count = 1
        if event.num == 4 or event.delta > 0:
            count = -1
        self.canvas.yview("scroll", count, "units")          
        
    def createTextbox(self, text, x, y, bgColor="#99DAE8", textColor='black'):
        tb = TextBox(self.canvas, text, x, y, bgColor, textColor)
        self.items[tb.id] = tb
        return tb
    
    def createEdge(self, src, dst):
        edge = Edge(self.canvas, src, dst)
        #self.items[edge.id] = edge
        return edge
    
    def clear(self):
        """Clear all items from the canvas"""
        self.canvas.delete(tk.ALL)
       
        
class TextBox():
    """This class will serve to paint and store
    rectange boxes with some text"""
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        self.bgColor = bgColor
        self.textColor = textColor
        self.text = text
        self.margin = 3
        self.canvas = canvas
        self.x = x
        self.y = y
        self.paint()
        self.listeners = []
        
    def paint(self):
        """Paint the object in a specific position."""
        self.id_text = self.canvas.create_text(self.x, self.y, text=self.text, 
                                               justify=tk.CENTER, fill=self.textColor)
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
        """Move TextBox to new position
        dx and y are differences to current position"""
        self.canvas.move(self.id_text, dx, dy)
        self.canvas.move(self.id, dx, dy)
        self.x += dx
        self.y += dy
        for listenerFunc in self.listeners:
            listenerFunc(dx, dy)
            
    def moveTo(self, x, y):
        """Move TextBox to new position
        dx and dy are the new position"""
        self.move(x-self.x, y-self.y)
        
    def addPositionListener(self, listenerFunc):
        self.listeners.append(listenerFunc)
        
    def setSelected(self, value):
        bw = 1
        if value:
            bw = 2
        self.canvas.itemconfig(self.id, width=bw)
        

class Edge():
    """Edge between two objects"""
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
        
#class MouseMover():
#
#  def __init__(self):
#    self.unselect()
#    #self.item = 0; self.previous = (0, 0)
#
#  def unselect(self):
#    self.item = None
#    self.previous = None
#
#  def select(self, event):
#    widget = event.widget                       # Get handle to canvas 
#    # Convert screen coordinates to canvas coordinates
#    xc = widget.canvasx(event.x)
#    yc = widget.canvasx(event.y)
#    rect = (xc - 1, yc - 1,  xc + 1, yc + 1)
#    items = widget.find_overlapping(xc - 1, yc - 1,  xc + 1, yc + 1)
#    if len(items):
#      self.item = items[0]        # ID for closest
#      self.previous = (xc, yc)
#      print((xc, yc, self.item))
#    else:
#      self.unselect()
#
#  def drag(self, event):
#    if self.item:
#      widget = event.widget
#      xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
#      canvas.move(self.item, xc-self.previous[0], yc-self.previous[1])
#      self.previous = (xc, yc)
 
if __name__ == '__main__':
    root = tk.Tk()
    canvas = Canvas(root, width=400, height=400)
    canvas.grid(row=0, column=0, sticky='nsew')
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)
    
    tb1 = canvas.createTextbox("Project", 100, 100, "blue")
    tb2 = canvas.createTextbox("aqui estoy yo\ny tu tb", 200, 200)
    tb3 = canvas.createTextbox("otro mas\n", 100, 200, "red")
    e1 = canvas.createEdge(tb1, tb2)
    e2 = canvas.createEdge(tb1, tb3)
    
    tb3.moveTo(100, 300)

    root.mainloop()
