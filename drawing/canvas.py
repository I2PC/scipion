#!/usr/bin/env xmipp_python

import Tkinter as tk


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
        
        # Add bindings
        self.canvas.bind("<Button-1>", self.onClick)
        self.canvas.bind("<Double-Button-1>", self.onDoubleClick)
        self.canvas.bind("<B1-Motion>", self.onDrag)
    
    def getCoordinates(self, event):
        '''Converts the events coordinates to canvas coordinates'''
        # Convert screen coordinates to canvas coordinates
        xc = self.canvas.canvasx(event.x)
        yc = self.canvas.canvasx(event.y)
        return (xc, yc)
    
    def onClick(self, event):
        xc, yc = self.getCoordinates(event)
        items = self.canvas.find_overlapping(xc - 1, yc - 1,  xc + 1, yc + 1)
        self.lastItem = None
        self.lastPos = (0, 0)
        for i in items:
            if i in self.items:
                self.lastItem = self.items[i]
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
                print self.lastItem.text
                self.lastPos = (xc, yc)
                break

    def onDrag(self, event):
        if self.lastItem:
            xc, yc = self.getCoordinates(event)
            self.lastItem.move(xc-self.lastPos[0], yc-self.lastPos[1])
            self.lastPos = (xc, yc)
            
        
    def createTextbox(self, text, x, y, bgColor="#99DAE8"):
        tb = TextBox(self.canvas, text, x, y, bgColor)
        self.items[tb.id] = tb
        return tb
    
    def createEdge(self, src, dst):
        edge = Edge(self.canvas, src, dst)
        #self.items[edge.id] = edge
        return edge
       
        
class TextBox():
    '''This class will serve to paint and store
    rectange boxes with some text'''
    def __init__(self, canvas, text, x, y, bgColor):
        self.bgColor = bgColor
        self.text = text
        self.margin = 3
        self.canvas = canvas
        self.pos = (x, y)
        self.paint()
        self.listeners = []
        
    def paint(self):
        '''Paint the object in a specific position.'''
        self.id_text = self.canvas.create_text(self.pos[0], self.pos[1], text=self.text)
        xr, yr, w, h = self.canvas.bbox(self.id_text)
        m = self.margin
        xr -= m
        yr -= m
        self.id = self.canvas.create_rectangle(xr, yr, w+2*m, h+2*m, fill=self.bgColor)
        self.canvas.tag_raise(self.id_text)
        
    def move(self, x, y):
        '''Move TextBox to new position'''
        self.canvas.move(self.id_text, x, y)
        self.canvas.move(self.id, x, y)
        for listenerFunc in self.listeners:
            listenerFunc(x, y)
        
    def addPositionListener(self, listenerFunc):
        self.listeners.append(listenerFunc)
        

class Edge():
    '''Edge between two objects'''
    def __init__(self, canvas, source, dest):
        self.posSrc = source.pos
        self.posDst = dest.pos
        self.canvas = canvas
        source.addPositionListener(self.updateSrc)
        dest.addPositionListener(self.updateDst)
        self.id = None
        self.paint()
        
    def paint(self):
        if self.id:
            self.canvas.delete(self.id)
        self.id = self.canvas.create_line(self.posSrc[0], self.posSrc[1], 
                                          self.posDst[0], self.posDst[1],
                                          width=2, splinesteps=25)
        self.canvas.tag_lower(self.id)
        
    def updateSrc(self, x, y):
        self.posSrc = (self.posSrc[0]+x, self.posSrc[1]+y)
        self.paint()
        
    def updateDst(self, x, y):
        self.posDst = (self.posDst[0]+x, self.posDst[1]+y)
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
 

root = tk.Tk()
canvas = Canvas(root, width=400, height=200)
canvas.grid(row=0, column=0, sticky='nsew')
root.grid_columnconfigure(0, weight=1)
root.grid_rowconfigure(0, weight=1)

tb1 = canvas.createTextbox("aqui", 100, 100, "blue")
tb2 = canvas.createTextbox("aqui estoy yo\ny tu tb", 200, 200)

e = canvas.createEdge(tb1, tb2)

#canvas.create_oval(10, 10, 110, 60, fill="grey")
#canvas.create_text(60, 35, text="Oval")
#canvas.create_rectangle(10, 100, 110, 150, outline="blue", fill="green")
#canvas.create_text(60, 125, text="Rectangle")
#canvas.create_line(60, 60, 60, 100, width=3)
## Get an instance of the MouseMover object
#mm = MouseMover()
# 
## Bind mouse events to methods (could also be in the constructor)
#canvas.bind("<Button-1>", mm.select)
#canvas.bind("<B1-Motion>", mm.drag)
root.mainloop()
