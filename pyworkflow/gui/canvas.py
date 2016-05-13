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
import math
import Tkinter as tk

import gui
import operator
from widgets import Scrollable


DEFAULT_CONNECTOR_FILL = "blue"
DEFAULT_CONNECTOR_OUTLINE = "black"


class Canvas(tk.Canvas, Scrollable):
    """Canvas to draw some objects.
    It will really contains a Frame, a Canvas and scrollbars"""
    _images = {}
    
    def __init__(self, parent, tooltipCallback=None, tooltipDelay=1500, **kwargs):
        defaults = {'bg': 'white'}
        defaults.update(kwargs)
        Scrollable.__init__(self, parent, tk.Canvas, **defaults)
        
        self.lastItem = None # Track last item selected
        self.lastPos = (0, 0) # Track last clicked position
        self.firstPos = None  # Track first clicked position (for a drag action)
        self.items = {} # Keep a dictionary with high-level items
        self.cleanSelected = True
        
        self.onClickCallback = None
        self.onDoubleClickCallback = None
        self.onRightClickCallback = None
        self.onControlClickCallback = None
        self.onAreaSelected = None
        
        # Add bindings
        self.bind("<Button-1>", self.onClick)
        self.bind("<ButtonRelease-1>", self.onButton1Release)
        self.bind("<Button-3>", self.onRightClick)
        self.bind("<Double-Button-1>", self.onDoubleClick)
        self.bind("<B1-Motion>", self.onDrag)
                # Hide the right-click menu
        self.bind('<FocusOut>', self._unpostMenu)
        self.bind("<Key>", self._unpostMenu)
        self.bind("<Control-1>", self.onControlClick)
        #self.bind("<MouseWheel>", self.onScroll)
        
        self._tooltipId = None
        self._tooltipOn = False # True if the tooltip is displayed
        self._tooltipCallback = tooltipCallback
        self._tooltipDelay = tooltipDelay
        
        if tooltipCallback:
            self.bind('<Motion>', self.onMotion)
            #self.bind('<Leave>', self.onLeave)
            self._createTooltip() # This should set
        
        self._menu = tk.Menu(self, tearoff=0)
        
    def _createTooltip(self):
        """ Create a Tooltip window to display tooltips in 
        the canvas.
        """
        tw = tk.Toplevel(self)
        tw.withdraw() # hidden by default
        tw.wm_overrideredirect(1) # Remove window decorations
        tw.bind("<Leave>", self.hideTooltip)
        
        self._tooltip = tw
        
    def _showTooltip(self, x, y, item):
        # check that the mouse is still in the position
        nx = self.winfo_pointerx()
        ny = self.winfo_pointery()
        if x == nx and y == ny:
            self._tooltipOn = True
            tw = self._tooltip # short notation
            self._tooltipCallback(tw, item)
            tw.update_idletasks()
            tw.wm_geometry("+%d+%d" % (x, y))
            tw.deiconify()
        
    def hideTooltip(self, e=None):
        if self._tooltipOn:
            self._tooltipOn = False
            tw = self._tooltip # short notation
            tw.withdraw()
        
    def getImage(self, img):
        return gui.getImage(img, self._images)
    
    def _unpostMenu(self, e=None):
        self._menu.unpost()
      
    def getCoordinates(self, event):
        """Converts the events coordinates to canvas coordinates"""
        # Convert screen coordinates to canvas coordinates
        xc = self.canvasx(event.x)
        yc = self.canvasy(event.y)
        return (xc, yc)
    
    def selectItem(self, item):
        if self.lastItem:
            self.lastItem.setSelected(False)
        self.lastItem = item
        item.setSelected(True)

    def _findItem(self, xc, yc):
        """ Find if there is any item in the canvas
        in the mouse event coordinates.
        Return None if not Found
        """
        items = self.find_overlapping(xc-1, yc-1,  xc+1, yc+1)
        for i in items:
            if i in self.items:
                return self.items[i]
        return None
             
    def _handleMouseEvent(self, event, callback):
        xc, yc = self.getCoordinates(event)
        self.lastItem = self._findItem(xc, yc)
        self.callbackResults = None
        self.lastPos = (0, 0)
        if self.lastItem is not None:
            if callback:
                self.callbackResults = callback(self.lastItem)
            self.lastPos = (xc, yc)

    def onClick(self, event):
        self.cleanSelected = True
        self._unpostMenu()
        self._handleMouseEvent(event, self.onClickCallback)
        
    def onControlClick(self, event):
        self.cleanSelected = False
        self._unpostMenu()
        self._handleMouseEvent(event, self.onControlClickCallback)
            
    def onRightClick(self, e=None):
        # RightClick callback will not work not, as it need
        # the event information to know the coordinates
        self._handleMouseEvent(e, self.onRightClickCallback)
        unpost = True
        # If the callback return a list of actions
        # we will show up a menu with them
        actions = self.callbackResults
        
        if actions:
            self._menu.delete(0, tk.END)
            for a in actions:
                if a is None: 
                    self._menu.add_separator()
                else:
                    img = ''
                    if len(a) > 2: # image for the action
                        img = self.getImage(a[2])
                    self._menu.add_command(label=a[0], command=a[1], 
                                          image=img, compound=tk.LEFT)
            self._menu.post(e.x_root, e.y_root)
            unpost = False
        if unpost:
            self._menu.unpost()
    
    def onDoubleClick(self, event):
        self._handleMouseEvent(event, self.onDoubleClickCallback)

    def onDrag(self, event):
        try:
            if self.lastItem:
                xc, yc = self.getCoordinates(event)
                self.lastItem.move(xc-self.lastPos[0], yc-self.lastPos[1])
                self.lastPos = (xc, yc)

            elif self.firstPos is None:
                self.firstPos = (event.x, event.y)
                # print "onDrag position captured."
        except Exception, ex:
            # JMRT: We are having a weird exception here.
            # Presumably because there is concurrency between the onDrag
            # event and the refresh one. For now, just ignore it.
            pass

    def onButton1Release(self, event):

        if self.firstPos is not None:

            self.onAreaSelected(self.firstPos[0], self.firstPos[1], event.x, event.y)

            self.firstPos = None

    def onMotion(self, event):
        self.onLeave(event) # Hide tooltip and cancel schedule
            
        xc, yc = self.getCoordinates(event)
        item = self._findItem(xc, yc)
        if item is not None:
            self._tooltipId = self.after(self._tooltipDelay, 
                                         lambda: self._showTooltip(event.x_root,
                                                                   event.y_root,
                                                                   item))  
        
    def onLeave(self, event):
        if self._tooltipId:
            self.after_cancel(self._tooltipId)
            self.hideTooltip()
            
    def createTextbox(self, text, x, y, bgColor="#99DAE8", textColor='black'):
        tb = TextBox(self, text, x, y, bgColor, textColor)
        self.items[tb.id] = tb
        return tb

    def createTextCircle(self, text, x, y, bgColor="#99DAE8", textColor='black'):
        tb = TextCircle(self, text, x, y, bgColor, textColor)
        self.items[tb.id] = tb
        return tb

    def createRoundedTextbox(self, text, x, y, bgColor="#99DAE8", textColor='black'):
        tb = RoundedTextBox(self, text, x, y, bgColor, textColor)
        self.items[tb.id] = tb
        return tb
    
    def addItem(self, item):
        self.items[item.id] = item
    
    def createEdge(self, srcItem, dstItem):
        edge = Edge(self, srcItem, dstItem)
        #self.items[edge.id] = edge
        return edge
    
    def createCable(self,src,srcSocket,dst,dstSocket):
        return Cable(self,src,srcSocket,dst,dstSocket)

    def clear(self):
        """ Clear all items from the canvas """
        self.delete(tk.ALL)

    def updateScrollRegion(self):
        self.update_idletasks()
        self.config(scrollregion=self.bbox("all"))
        
    def drawGraph(self, graph, layout=None, drawNode=None):
        """ Draw a graph in the canvas.
        nodes in the graph should have x and y.
        If layout is not None, it will be used to 
        reorganize the node positions.
        Provide drawNode if you want to customize how
        to create the boxes for each graph node.
        """
        
        if drawNode is None:
            self.drawNode = self._drawNode
        else:
            self.drawNode = drawNode
        
        self._drawNodes(graph.getRoot(), {})
        
        if layout is not None:
            layout.draw(graph)
            
        # Update node positions
        self._updatePositions(graph.getRoot(), {})
        self.updateScrollRegion()
        
    def _drawNode(self, canvas, node):
        """ Default implementation to draw nodes as textboxes. """
        return TextBox(self, node.getLabel(), 0, 0,
                       bgColor="#99DAE8", textColor='black')
        
    def _drawNodes(self, node, visitedDict={}):
        nodeName = node.getName()
        
        if nodeName not in visitedDict:
            visitedDict[nodeName] = True
            item = self.drawNode(self, node)
            node.width, node.height = item.getDimensions()
            node.item = item
            item.node = node
            self.addItem(item)

            if getattr(node, 'expanded', True):
                for child in node.getChilds():
                    self._drawNodes(child, visitedDict)
            else:
                self._setupParentProperties(node, visitedDict)
                
    def _setupParentProperties(self, node, visitedDict):
        """ This methods is used for collapsed nodes, in which 
        the properties (width, height, x and y) is propagated
        to the hidden childs.
        """
        for child in node.getChilds():
            if child.getName() not in visitedDict:
                child.width = node.width
                child.height = node.height
                child.x = node.x
                child.y = node.y
                self._setupParentProperties(child, visitedDict)
        
    def _updatePositions(self, node, visitedDict={}):
        """ Update position of nodes and create the edges. """
        nodeName = node.getName()
        
        if nodeName not in visitedDict:
            visitedDict[nodeName] = True
            item = node.item
            item.moveTo(node.x, node.y)
            
            if getattr(node, 'expanded', True):
                for child in node.getChilds():
                    self.createEdge(item, child.item)
                    self._updatePositions(child, visitedDict)
            
        
        

def findClosestPoints(list1, list2):
    candidates=[]
    for c1 in list1:
        for c2 in list2:
            candidates.append([c1,c2, math.hypot(c2[0] - c1[0], c2[1] - c1[1])])
    closestTuple=min(candidates,key=operator.itemgetter(2))
    return closestTuple[0],closestTuple[1]


def findClosestConnectors(item1, item2):
    return findUpDownClosestConnectors(item1,item2)


def findUpDownClosestConnectors(item1, item2):
    srcConnectors = item1.getUpDownConnectorsCoordinates()
    dstConnectors = item2.getUpDownConnectorsCoordinates()
    if srcConnectors and dstConnectors:
        c1Coords, c2Coords = findClosestPoints(srcConnectors,dstConnectors)
        return c1Coords, c2Coords
    return None


def findStrictClosestConnectors(item1, item2):
    srcConnectors = item1.getConnectorsCoordinates()
    dstConnectors = item2.getConnectorsCoordinates()
    c1Coords,c2Coords = findClosestPoints(srcConnectors,dstConnectors)
    return c1Coords, c2Coords


def getConnectors(itemSource, itemDest):

    srcConnector = itemSource.getOutputConnectorCoordinates()
    dstConnector = itemDest.getInputConnectorCoordinates()

    return srcConnector, dstConnector


class Item(object):
    socketSeparation = 12
    verticalFlow = True

    TOP = 0
    RIGHT = 1
    BOTTOM = 2
    LEFT = 3

    def __init__(self, canvas, x, y):
        self.activeConnector = None
        self.canvas = canvas
        self.x = x
        self.y = y
        self.sockets = {}
        self.listeners = []
        self.selectionListeners = []

    def getCenter(self,x1,y1,x2,y2):
        xc=(x2+x1)/2.0
        yc=(y2+y1)/2.0
        return (xc,yc)

    def getConnectorsCoordinates(self):
        x1,y1,x2,y2=self.getCorners()
        xc,yc=self.getCenter(x1,y1,x2,y2)
        return [(xc,y1), (x2,yc), (xc,y2), (x1,yc)]

    def getTopConnectorCoordinates(self):

        return self._getConnectorCoordinates(self.TOP)

    def getBottomConnectorCoordinates(self):

        return self._getConnectorCoordinates(self.BOTTOM)

    def getLeftConnectorCoordinates(self):

        return self._getConnectorCoordinates(self.LEFT)

    def getRightConnectorCoordinates(self):

        return self._getConnectorCoordinates(self.RIGHT)

    def _getConnectorCoordinates(self, side):
        fourConnectors = self.getConnectorsCoordinates()
        return fourConnectors[side]

    def getInputConnectorCoordinates(self):

        if self.verticalFlow:
            return self.getTopConnectorCoordinates()
        else:
            return self.getLeftConnectorCoordinates()

    def getOutputConnectorCoordinates(self):

        if self.verticalFlow:
            return self.getBottomConnectorCoordinates()
        else:
            return self.getRightConnectorCoordinates()

    def getUpDownConnectorsCoordinates(self):
        corners = self.getCorners()
        if corners:
            x1,y1,x2,y2 = self.getCorners()
            xc, yc = self.getCenter(x1,y1,x2,y2)
            return [(xc,y1), (xc,y2)]
        return None

    def getCorners(self): 
        return self.canvas.bbox(self.id)

    def countSockets(self,verticalLocation):
        return len(self.getSocketsAt(verticalLocation))

    def addSocket(self,name,socketClass,verticalLocation,fillColor=DEFAULT_CONNECTOR_FILL,outline=DEFAULT_CONNECTOR_OUTLINE,position=None):
        count=self.countSockets(verticalLocation) + 1
        if position==None:
            position=count
        self.relocateSockets(verticalLocation, count)
        x,y=self.getSocketCoordsAt(verticalLocation, count, count)
        self.sockets[name]={"object": socketClass(self.canvas,x,y,name,fillColor=fillColor,outline=outline), "verticalLocation": verticalLocation, "position":position}
        self.paintSocket(self.getSocket(name))

    def getSocket(self,name):
        return self.sockets[name]["object"]

    def getSocketsAt(self,verticalLocation):
        return filter(lambda s: s["verticalLocation"] == verticalLocation, self.sockets.values())

    def getSocketCoords(self,name):
        socket=self.sockets[name]
        return self.getSocketCoordsAt(socket["verticalLocation"], socket["position"], self.countSockets(socket["verticalLocation"]))

    def getSocketCoordsAt(self,verticalLocation,position=1,socketsCount=1):
        x1,y1,x2,y2=self.getCorners()
        xc=(x2+x1)/2.0
        socketsGroupSize=(socketsCount-1)*self.socketSeparation
        socketsGroupStart=xc - (socketsGroupSize / 2)
        x=socketsGroupStart+(position-1)*self.socketSeparation
        if verticalLocation == "top":
            y=y1
        else:
            y=y2
        return (x,y)

    def relocateSockets(self,verticalLocation,count):
        sockets=self.getSocketsAt(verticalLocation)
        for socket in sockets:
            o=socket["object"]
            x,y=self.getSocketCoordsAt(verticalLocation,socket["position"],count)
            o.moveTo(x,y)


    def paintSocket(self,socket):
        # x,y=self.getSocketCoords(socket["name"])
        socket.paintSocket()

    def paintSockets(self):
        for name in self.sockets.keys():
            self.paintSocket(self.getSocket(name))

    def getDimensions(self):
        x, y, x2, y2 = self.canvas.bbox(self.id)
        return (x2-x, y2-y)
        
    def move(self, dx, dy):
        if hasattr(self,"id"):
            self.canvas.move(self.id, dx, dy)
        self.x += dx
        self.y += dy
        for name in self.sockets.keys():
            socket=self.sockets[name]
            socket["object"].move(dx,dy)
        for listenerFunc in self.listeners:
            listenerFunc(dx, dy)
            
    def moveTo(self, x, y):
        """Move TextBox to a new position (x,y)"""
        self.move(x-self.x, y-self.y)
        
    def addPositionListener(self, listenerFunc):
        self.listeners.append(listenerFunc)

    def addSelectionListener(self, listenerFunc):
        self.selectionListeners.append(listenerFunc)

    def _notifySelectionListeners(self, value):

        for listenerFunc in self.selectionListeners:
            listenerFunc(value)

    def setSelected(self, value):
        bw = 0
        bc = 'black'
        if value:
            bw = 2
            self.lift()
            # bc = 'Firebrick'
        self.canvas.itemconfig(self.id, width=bw, outline=bc)
        self._notifySelectionListeners(value)


    def lift(self):
        self.canvas.lift(self.id)
        
class TextItem(Item):
    """This class will serve to paint and store rectangle boxes with some text.
       x and y are the coordinates of the center of this item"""
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        super(TextItem,self).__init__(canvas,x,y)
        self.bgColor = bgColor
        self.textColor = textColor
        self.text = text
        self.margin = 8
        self.paint()
        
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
        self.canvas.lift(self.id_text)

    def move(self, dx, dy):
        super(TextItem, self).move(dx,dy)
        self.canvas.move(self.id_text, dx, dy)

    def lift(self):
        super(TextItem, self).lift()
        self.canvas.lift(self.id_text)
   
        
class TextBox(TextItem):
    def __init__(self, canvas, text, x, y, bgColor, textColor='black'):
        super(TextBox,self).__init__(canvas, text, x, y, bgColor, textColor)

    def _paintBounds(self, x, y, w, h, fillColor):
        return self.canvas.create_rectangle(x, y, w, h, fill=fillColor, outline=fillColor)

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
    
    
class ImageBox(Item):
    def __init__(self, canvas, imgPath, x=0, y=0, text=None):
        Item.__init__(self, canvas, x, y)
        # Create the image
        from pyworkflow.gui import getImage, getImageFromPath
        
        if imgPath is None:
            self.image = getImage('no-image.png')
        else:
            self.image = getImageFromPath(imgPath)

        if text is not None:
            self.label = tk.Label(canvas, image=self.image, text=text, 
                                  compound=tk.TOP, bg='gray')
            self.id = self.canvas.create_window(x, y, window=self.label)
            self.label.bind('<Button-1>', self._onClick)
        else:
            self.id = self.canvas.create_image(x, y, image=self.image)
         
        
    def setSelected(self, value): #Ignore selection highlight
        pass
    
    def _onClick(self, e=None):
        pass

class Connector(Item):
    """ Default connector has no graphical representation (hence, it'ss invisible). Subclasses offer different looks"""
    def __init__(self,canvas,x,y,name):
        super(Connector,self).__init__(canvas, x, y)
        self.name= name

    def paintSocket(self):
        """Should be implemented by the subclasses"""
        pass

    def paintPlug(self, canvas, x, y):
        """Should be implemented by the subclasses"""
        pass

    def move(self,dx,dy):
        super(Connector,self).move(dx,dy)
        if hasattr(self,"socketId"):
            self.canvas.move(self.socketId, dx, dy)
        if hasattr(self,"plugId"):
            self.canvas.move(self.plugId, dx, dy)

class ColoredConnector(Connector):
    def __init__(self,canvas,x,y,name,fillColor=DEFAULT_CONNECTOR_FILL,outline=DEFAULT_CONNECTOR_OUTLINE):
        super(ColoredConnector,self).__init__(canvas,x,y,name)
        self.fillColor=fillColor
        self.outline=outline

class RoundConnector(ColoredConnector):
    radius=3
    def paintSocket(self):
        self.socketId= self.canvas.create_oval(self.x-self.radius, self.y-self.radius, self.x+self.radius, self.y+self.radius, outline=self.outline)

    def paintPlug(self):
        self.plugId= self.canvas.create_oval(self.x-self.radius, self.y-self.radius, self.x+self.radius, self.y+self.radius, fill=self.fillColor,width=0)


class SquareConnector(ColoredConnector):
    halfside=3
    def paintSocket(self):
        self.socketId= self.canvas.create_rectangle(self.x-self.halfside, self.y-self.halfside, self.x+self.halfside, self.y+self.halfside, outline=self.outline)

    def paintPlug(self):
        self.plugId= self.canvas.create_rectangle(self.x-self.halfside, self.y-self.halfside, self.x+self.halfside, self.y+self.halfside, fill=self.fillColor,width=0)


# !!!! other figures: half circle, diamond...
class Oval():

    """Oval or circle"""
    def __init__(self, canvas, x, y, radio, color='green', anchor=None):

        self.anchor = anchor
        self.X, self.Y = x, y
        self.radio = radio
        self.color = color
        self.canvas = canvas
        anchor.addPositionListener(self.updateSrc)
        anchor.addSelectionListener(self.selectionListener)
        self.id = None
        self.paint()

    def paint(self):

        if self.id:
            self.canvas.delete(self.id)

        self.id = self.canvas.create_oval(self.X, self.Y, self.X + self.radio, self.Y + self.radio, fill=self.color, outline='black')
        # self.canvas.tag_raise(self.id)

    def updateSrc(self, dx, dy):
        self.X += dx
        self.Y += dy
        self.paint()

    def selectionListener(self, value):
        if value:
            self.canvas.lift(self.id)


class Rectangle():

    def __init__(self, canvas, x, y, width, height=None, color='green', anchor=None):

        self.anchor = anchor
        self.X, self.Y = x, y
        self.width = width
        self.height = height or width
        self.color = color
        self.canvas = canvas
        anchor.addPositionListener(self.updateSrc)
        anchor.addSelectionListener(self.selectionListener)
        self.id = None
        self.paint()

    def paint(self):

        if self.id:
            self.canvas.delete(self.id)

        self.id = self.canvas.create_rectangle(self.X, self.Y, self.X + self.width, self.Y + self.height, fill=self.color, outline=self.color)
        # self.canvas.tag_raise(self.id)

    def updateSrc(self, dx, dy):
        self.X += dx
        self.Y += dy
        self.paint()

    def selectionListener(self, value):
        if value:
            self.canvas.lift(self.id)


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
        # coords = findClosestConnectors(self.source,self.dest)
        coords = getConnectors(self.source, self.dest)

        if coords:
            c1Coords, c2Coords = coords
    
            if self.id:
                self.canvas.delete(self.id)
    
            self.id = self.canvas.create_line(c1Coords[0], c1Coords[1], 
                                              c2Coords[0], c2Coords[1],
                                              width=1., fill='#ccc')
            self.canvas.tag_lower(self.id)
        
    def updateSrc(self, dx, dy):
        self.srcX += dx
        self.srcY += dy
        self.paint()
        
    def updateDst(self, dx, dy):
        self.dstX += dx
        self.dstY += dy
        self.paint()


# !!!! Interaction: allow to reconnect cables dynamically

# !!!! Antialiasing for the line - it seems TkInter does not support antialiasing...
# Although Tk 8.5 supports anti-aliasing if the Cairo library is installed:
# @see http://wiki.tcl.tk/10630

class Cable():
    def __init__(self,canvas,src,srcConnector,dst,dstConnector):
        self.id=None
        self.canvas=canvas
        self.srcPlug=src.getSocket(srcConnector)
        self.dstPlug=dst.getSocket(dstConnector) 
        self.srcX=self.srcPlug.x
        self.srcY=self.srcPlug.y
        self.dstX,self.dstY=dst.getSocketCoords(dstConnector)
        src.addPositionListener(self.srcMoved)
        dst.addPositionListener(self.dstMoved)
        self.paint()

    def srcMoved(self,dx,dy):
        self.srcX=self.srcX+dx
        self.srcY=self.srcY+dy
        self.updateCoords()

    def updateCoords(self):
        self.canvas.coords(self.id, self.srcX, self.srcY,self.dstX,self.dstY)


    def dstMoved(self,dx,dy):
        self.dstX=self.dstX+dx
        self.dstY=self.dstY+dy
        self.updateCoords()

    def paint(self):
        self.id = self.canvas.create_line(self.srcX,self.srcY,self.dstX,self.dstY, width=2)
        self.canvas.tag_lower(self.id)
        self.paintPlugs()
        
    def paintPlugs(self):
        self.srcPlug.paintPlug()
        self.dstPlug.paintPlug()

 
if __name__ == '__main__':
    root = tk.Tk()
    canvas = Canvas(root, width=800, height=600)
    canvas.frame.grid(row=0, column=0, sticky='nsew')
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)
    
    def canvasExample1():    
        tb1 = canvas.createTextCircle("Project", 100, 100, "blue")
        tb2 = canvas.createTextbox("This is an intentionally quite big, big box,\nas you may appreciate looking carefully\nat it,\nas many times\nas you might need", 300, 200)
        tb2.addSocket("output1",RoundConnector, "bottom",fillColor="green")
        tb2.addSocket("output2",SquareConnector, "bottom",fillColor="yellow")
        tb2.addSocket("output3",SquareConnector, "bottom",fillColor="blue")
        tb3 = canvas.createRoundedTextbox("otro mas\n", 100, 200, "red")
        tb4 = canvas.createRoundedTextbox("tb4", 300, 300, "yellow")
        tb4.addSocket("input1",SquareConnector, "top",outline="red")
        tb5 = canvas.createTextCircle("tb5", 400, 300, "grey")
        tb5.addSocket("input1",SquareConnector, "top")
        e1 = canvas.createEdge(tb1, tb2)
        e2 = canvas.createEdge(tb1, tb3)
        c1= canvas.createCable(tb2,"output2",tb4,"input1")
        c2= canvas.createCable(tb2,"output3",tb5,"input1")
        tb3.moveTo(100, 300)
        
       
    canvasExample1()

    root.mainloop()
