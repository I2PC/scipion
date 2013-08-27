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
import math
import operator

# !!!! A connector can be a Socket (female) or a Plug (male)
# !!!! Different connectors paint socket and plug differently
# !!!! A cable connects 2 connectors. It only needs to know the coordinates of the connectors to paint itself


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
    
    def _handleMouseEvent(self, event, callback):
        xc, yc = self.getCoordinates(event)
        items = self.canvas.find_overlapping(xc-1, yc-1,  xc+1, yc+1)
        if self.lastItem:
            self.lastItem.setSelected(False)
            self.lastItem = None
        self.lastPos = (0, 0)
        for i in items:
            if i in self.items:
                self.lastItem = self.items[i]
                self.lastItem.setSelected(True)
                if callback:
                    callback(self.lastItem)
                self.lastPos = (xc, yc)
                break
        
    def onClick(self, event):
        self._handleMouseEvent(event, self.onClickCallback)
            
    def onRightClick(self, event):
        # RightClick callback will not work not, as it need
        # the event information to know the coordinates
        self._handleMouseEvent(event, self.onRightClickCallback)
    
    def onDoubleClick(self, event):
        self._handleMouseEvent(event, self.onDoubleClickCallback)
#        xc, yc = self.getCoordinates(event)
#        items = self.canvas.find_overlapping(xc - 1, yc - 1,  xc + 1, yc + 1)
#        self.lastItem = None
#        self.lastPos = (0, 0)
#        for i in items:
#            if i in self.items:
#                self.lastItem = self.items[i]
#                if self.onDoubleClickCallback:
#                    self.onDoubleClickCallback(self.lastItem)
#                self.lastPos = (xc, yc)
#                break

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

    def createTextCircle(self, text, x, y, bgColor="#99DAE8", textColor='black'):
        tb = TextCircle(self.canvas, text, x, y, bgColor, textColor)
        self.items[tb.id] = tb
        return tb

    def createRoundedTextbox(self, text, x, y, bgColor="#99DAE8", textColor='black'):
        tb = RoundedTextBox(self.canvas, text, x, y, bgColor, textColor)
        self.items[tb.id] = tb
        return tb

    
    def createEdge(self, srcItem, dstItem):
        edge = Edge(self.canvas, srcItem, dstItem)
        #self.items[edge.id] = edge
        return edge
    

    def clear(self):
        """Clear all items from the canvas"""
        self.canvas.delete(tk.ALL)
       

def findClosestPoints(list1,list2):
    candidates=[]
    for c1 in list1:
        for c2 in list2:
            candidates.append([c1,c2, math.hypot(c2[0] - c1[0], c2[1] - c1[1])])
    closestTuple=min(candidates,key=operator.itemgetter(2))
    return closestTuple[0],closestTuple[1]



# !!!! sometimes, in a hierarchichal graph, it's better the upper (or lower) connector, than the closest...
def findClosestConnectors(item1,item2):
    srcConnectors=item1.getConnectorsCoordinates()
    dstConnectors=item2.getConnectorsCoordinates()
    c1Coords,c2Coords=findClosestPoints(srcConnectors,dstConnectors)
    item1.setActiveConnector(c1Coords)
    item2.setActiveConnector(c2Coords)
    return c1Coords,c2Coords



# !!!! refactor code repeated in subclasses, probably into this Item class
# !!!! An item can have many Sockets, in the top or bottom, in different positions (order)
class Item(object):      
    def __init__(self,canvas):
        self.activeConnector=None
        self.canvas=canvas
        self.sockets=[]

    def getConnectorsCoordinates(self):
        x1,y1,x2,y2=self.getCorners()
        xc=(x2+x1)/2.0
        yc=(y2+y1)/2.0
        return [(xc,y1), (x2,yc), (xc,y2), (x1,yc)]

    def getCorners(self):
        return self.canvas.bbox(self.id)

    def setActiveConnector(self,coords):
        if self.activeConnector:
            self.activeConnector.erase()
        self.activeConnector=Connector(self.canvas,coords[0],coords[1])
        # add the connector to the items list in order to track it for events

    def addSocket(self,name,socketClass,verticalLocation,position=None):
        self.sockets.append({"name": name, "verticalLocation": verticalLocation, "class":socketClass})
        self.paint()

    def paintSocket(self,socket):
        x1,y1,x2,y2=self.getCorners()
        xc=(x2+x1)/2.0
        if socket["verticalLocation"] == "top":
            y=y1
        else:
            y=y2
        print "paint socket at %d %d" %(xc,y)
        #class_=globals()[socket["class"]]
        socket["class"].paintSocket(self.canvas,xc,y)

    def paintSockets(self):
        print self.sockets
        for socket in self.sockets:
            self.paintSocket(socket)


        
class TextItem(Item):
    """This class will serve to paint and store rectange boxes with some text.
       x and y are the coordinates of the center of this item"""
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        super(TextItem,self).__init__(canvas)
        self.bgColor = bgColor
        self.textColor = textColor
        self.text = text
        self.margin = 8
        self.canvas = canvas
        self.x = x
        self.y = y
        self.paint()
        self.listeners = []
        
    def _paintBounds(self, x, y, w, h, fillColor):
        """ Subclasses should implement this method 
        to paint the bounds to the text.
        Normally the bound are: rectangle or circles.
        Params:
            x, y: top left corner of the bounding box
            w, h: width and height of the box
            fillColor: color to fill the background
        Returns:
            should return the id of the created shape
        """
        pass
        
    def paint(self):
        """Paint the object in a specific position."""
        self.id_text = self.canvas.create_text(self.x, self.y, text=self.text, 
                                               justify=tk.CENTER, fill=self.textColor)
        xr, yr, w, h = self.canvas.bbox(self.id_text)
        m = self.margin

        self.id = self._paintBounds(xr-m, yr-m, w+m, h+m, fillColor=self.bgColor)
        self.canvas.tag_raise(self.id_text)
    
        self.paintSockets()


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
   
        
class TextBox(TextItem):
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        super(TextBox,self).__init__(canvas, text, x, y, bgColor, textColor)

    def _paintBounds(self, x, y, w, h, fillColor):
        return self.canvas.create_rectangle(x, y, w, h, fill=fillColor) 

class RoundedTextBox(TextItem):
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        super(RoundedTextBox,self).__init__(canvas, text, x, y, bgColor, textColor)


    def _paintBounds(self, upperLeftX, upperLeftY , bottomRightX, bottomRightY, fillColor):
        d=5
        # When smooth=1, you define a straight segment by including its ends twice
        return self.canvas.create_polygon(upperLeftX+d+1,upperLeftY, #1
                                          upperLeftX+d,upperLeftY, #1
                                          bottomRightX-d,upperLeftY, #2
                                          bottomRightX-d,upperLeftY, #2
                                          # bottomRightX-d+1,upperLeftY, #2b
                                          bottomRightX,upperLeftY+d-1, #3b
                                          bottomRightX,upperLeftY+d, #3
                                          bottomRightX,upperLeftY+d, #3
                                          bottomRightX, bottomRightY-d, #4
                                          bottomRightX, bottomRightY-d, #4
                                          bottomRightX-d,bottomRightY, #5
                                          bottomRightX-d,bottomRightY, #5
                                          upperLeftX+d,bottomRightY, #6
                                          upperLeftX+d,bottomRightY, #6
                                          # upperLeftX+d-1,bottomRightY, #6b
                                          upperLeftX, bottomRightY-d+1, #7b
                                          upperLeftX, bottomRightY-d, #7
                                          upperLeftX, bottomRightY-d, #7
                                          upperLeftX, upperLeftY+d, #8
                                          upperLeftX, upperLeftY+d, #8
                                          # upperLeftX, upperLeftY+d-1, #8b
                                          upperLeftX+d-1,upperLeftY, #1b
                                          fill=fillColor, outline='black',smooth=1) 

        def getDimensions(self):
            return self.canvas.bbox(self.id)
    
class TextCircle(TextItem):
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        super(TextCircle,self).__init__(canvas, text, x, y, bgColor, textColor)

    def _paintBounds(self, x, y, w, h, fillColor):
        return self.canvas.create_oval(x, y, w, h, fill=fillColor) 
        

class Connector():
    def __init__(self,canvas,x,y,paintFn=None):
        self.id=None
        self.x=x
        self.y=y
        self.canvas=canvas
        if paintFn == None:
            self.paintFn=self.defaultPaint
        self.listeners=[]
        self.paint()
        
    def addPositionListener(self,listenerFunc):
        self.listeners.append(listenerFunc)        

    def defaultPaint(self):
        """ Right now, the default is not painting de connector"""
        return None

    def ovalPaint(self):
        return self.canvas.create_oval(self.x,self.y,self.x+5,self.y+5)

    @classmethod
    def paintPlug(cls,canvas,x,y):
        """Should be implemented by the subclasses"""
        pass

    # !!!! improve connector paint (sometimes is not well centered against its item)
    def paint(self):
        self.id=self.paintFn()

    def erase(self):
        if self.id:
            self.canvas.delete(self.id)

    def moveTo(self, x, y):
        self.move(x-self.x, y-self.y)

    def move(self, dx, dy):
        self.canvas.move(self.id, dx, dy)
        self.x += dx
        self.y += dy
        for listenerFunc in self.listeners:
            listenerFunc(dx, dy)

    def setSelected(self, value):
        bw = 1
        if value:
            bw = 2
        self.canvas.itemconfig(self.id, width=bw)
        print self.x,self.y

class RoundConnector(Connector):
    @classmethod
    def paintSocket(cls,canvas,x,y):
        return canvas.create_oval(x,y,x+5,y+5)

class Edge():
    """Edge between two objects"""
    def __init__(self, canvas, source, dest):
        self.source=source
        self.dest=dest
        self.srcX, self.srcY = source.x, source.y
        self.dstX, self.dstY = dest.x, dest.y
        self.canvas = canvas
        source.addPositionListener(self.updateSrc)
        dest.addPositionListener(self.updateDst)
        self.id = None
        self.paint()
        
    def paint(self):
        c1Coords,c2Coords=findClosestConnectors(self.source,self.dest)

        if self.id:
           self.canvas.delete(self.id)

        self.id = self.canvas.create_line(c1Coords[0], c1Coords[1], 
                                          c2Coords[0], c2Coords[1],
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
        
 
if __name__ == '__main__':
    root = tk.Tk()
    canvas = Canvas(root, width=800, height=600)
    canvas.grid(row=0, column=0, sticky='nsew')
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)
    
    tb1 = canvas.createTextCircle("Project", 100, 100, "blue")
    tb2 = canvas.createTextbox("This is an intentionally quite big, big box,\nas you may appreciate looking carefully\nat it,\nas many times\nas you might need", 300, 200)
    tb2.addSocket("output1",RoundConnector, "bottom")
    tb3 = canvas.createRoundedTextbox("otro mas\n", 100, 200, "red")
    tb4 = canvas.createRoundedTextbox("tb4", 300, 300, "yellow")
    tb4.addSocket("input1",RoundConnector, "top")
    e1 = canvas.createEdge(tb1, tb2)
    e2 = canvas.createEdge(tb1, tb3)
    
    tb3.moveTo(100, 300)

    root.mainloop()
