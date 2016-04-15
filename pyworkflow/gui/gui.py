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
This module contains all configuration
settings and gui general functions 
"""
import os
import Tkinter as tk
import tkFont
import Queue
from pyworkflow.object import OrderedObject
from pyworkflow.utils.path import findResource
from pyworkflow.utils.properties import Message, Color, Icon
from widgets import Button

"""
Some GUI CONFIGURATION parameters
"""
#TODO: read font size and name from config file
cfgFontName = os.environ.get('SCIPION_FONT_NAME', "Helvetica")
cfgFontSize = int(os.environ.get('SCIPION_FONT_SIZE', 10))  

#TextColor
cfgCitationTextColor = "dark olive green"
cfgLabelTextColor = "black"
cfgSectionTextColor = "blue4"
#Background Color
cfgBgColor = "light grey"
cfgLabelBgColor = "white"
cfgHighlightBgColor = cfgBgColor
cfgButtonFgColor = "white"
cfgButtonActiveFgColor = "white"
cfgButtonBgColor = Color.RED_COLOR
cfgButtonActiveBgColor = "#A60C0C"
cfgEntryBgColor = "lemon chiffon" 
cfgExpertLabelBgColor = "light salmon"
cfgSectionBgColor = cfgButtonBgColor
#Color
cfgListSelectColor = "DeepSkyBlue4"
cfgBooleanSelectColor = "white"
cfgButtonSelectColor = "DeepSkyBlue2"
#Dimensions limits
cfgMaxHeight = 650
cfgMaxWidth = 800
cfgMaxFontSize = 14
cfgMinFontSize = 6
cfgWrapLenght = cfgMaxWidth - 50


class Config(OrderedObject):
    pass

def saveConfig(filename):
    from pyworkflow.mapper import SqliteMapper
    from pyworkflow.object import String, Integer
    
    mapper = SqliteMapper(filename)
    o = Config()
    for k, v in globals().iteritems():
        if k.startswith('cfg'):
            if type(v) is str:
                value = String(v)
            else: 
                value = Integer(v)
            setattr(o, k, value)
    mapper.insert(o)
    mapper.commit()
            
   
"""
FONT related variables and functions 
"""
def setFont(fontKey, update=False, **opts):
    """Register a tkFont and store it in a globals of this module
    this method should be called only after a tk.Tk() windows has been
    created."""
    if not hasFont(fontKey) or update:
        globals()[fontKey] = tkFont.Font(**opts)
        
    return globals()[fontKey] 
    
def hasFont(fontKey):
    return fontKey in globals()

def aliasFont(fontAlias, fontKey):
    """Set a fontAlias as another alias name of fontKey"""
    g = globals()
    g[fontAlias] = g[fontKey] 
    
def setCommonFonts(windows=None):
    """Set some predifined common fonts.
    Same conditions of setFont applies here."""
    f = setFont('fontNormal', family=cfgFontName, size=cfgFontSize)
    aliasFont('fontButton', 'fontNormal')
    fb = setFont('fontBold', family=cfgFontName, size=cfgFontSize, weight='bold')
    fi = setFont('fontItalic', family=cfgFontName, size=cfgFontSize, slant='italic')
    setFont('fontLabel', family=cfgFontName, size=cfgFontSize+1, weight='bold')
    if windows:
        windows.fontBig = tkFont.Font(size=cfgFontSize+2, family=cfgFontName, weight='bold')
        windows.font = f
        windows.fontBold = fb
        windows.fontItalic = fi 

def changeFontSizeByDeltha(font, deltha, minSize=-999, maxSize=999):
    size = font['size']
    new_size = size + deltha
    if new_size >= minSize and new_size <= maxSize:
        font.configure(size=new_size)
            
def changeFontSize(font, event, minSize=-999, maxSize=999):
    deltha = 2
    if event.char == '-':
        deltha = -2
    changeFontSizeByDeltha(font, deltha, minSize, maxSize)

"""
IMAGE related variables and functions 
"""

def getImage(imageName, imgDict=None, tkImage=True, percent=100, maxheight=None):
    """ Search for the image in the RESOURCES path list. """
    if imageName is None:
        return None
    if imgDict is not None and imageName in imgDict:
        return imgDict[imageName]
    if not os.path.isabs(imageName):
        imagePath = findResource(imageName)
    else:
        imagePath = imageName
    image = None
    if imagePath:
        from PIL import Image
        image = Image.open(imagePath)
        w, h = image.size
        newSize = None
        if percent != 100: # Display image with other dimensions
            fp = float(percent)/100.0
            newSize = int(fp * w), int(fp * h)
        elif maxheight and h > maxheight:
            newSize = int(w * float(maxheight)/h), maxheight
        if newSize:
            image.thumbnail(newSize, Image.ANTIALIAS)
        if tkImage:
            from PIL import ImageTk
            image = ImageTk.PhotoImage(image)
        if imgDict is not None:
            imgDict[imageName] = image
    return image

def getPILImage(imageXmipp, dim=None, normalize=True):
    """ Given an image read by Xmipp, convert it to PIL. """
    from PIL import Image
    import xmipp
    
    if normalize:
        imageXmipp.convert2DataType(xmipp.DT_UCHAR, xmipp.CW_ADJUST)
        
    imageData = imageXmipp.getData()
    image = Image.fromarray(imageData)
    if dim:
        size = int(dim), int(dim)
        image.thumbnail(size, Image.ANTIALIAS)
    return image

def getTkImage(imageXmipp, filename, dim):
    from PIL import ImageTk
    imageXmipp.readPreview(filename, dim)
    return ImageTk.PhotoImage(getPILImage(imageXmipp))

def getImageFromPath(imagePath):
    """ Read an image using Xmipp, convert to PIL
    and then return as expected by Tk.
    """
    import xmipp
    img = xmipp.Image(imagePath)
    imgPIL = getPILImage(img)
    from PIL import ImageTk
    imgTk = ImageTk.PhotoImage(imgPIL)
    
    return imgTk

"""
Windows geometry utilities
"""
def getGeometry(win):
    ''' Return the geometry information of the windows
    It will be a tuple (width, height, x, y)
    '''
    return win.winfo_reqwidth(), win.winfo_reqheight(), win.winfo_x(), win.winfo_y()

def centerWindows(root, dim=None, refWindows=None):
    """Center a windows in the middle of the screen 
    or in the middle of other windows(refWindows param)"""
    root.update_idletasks()
    if dim is None:
        gw, gh, gx, gy = getGeometry(root)
    else:
        gw, gh = dim
    if refWindows:
        rw, rh, rx, ry = getGeometry(refWindows)
        x = rx + (rw - gw) / 2
        y = ry + (rh - gh) / 2 
    else:
        w = root.winfo_screenwidth()
        h = root.winfo_screenheight()
        x = (w - gw) / 2
        y = (h - gh) / 2
        
    root.geometry("%dx%d+%d+%d" % (gw, gh, x, y))
    
def configureWeigths(widget, row=0, column=0):
    """This function is a shortcut to a common
    used pair of calls: rowconfigure and columnconfigure
    for making childs widgets take the space available"""
    widget.columnconfigure(column, weight=1)
    widget.rowconfigure(row, weight=1)
    
    
class Window():
    """Class to manage a Tk windows.
    It will encapsulates some basic creation and 
    setup functions. """
    
    def __init__(self, title='', masterWindow=None, weight=True, minsize=(500, 300),
                 icon=None, **kwargs):
        """Create a Tk window.
        title: string to use as title for the windows.
        master: if not provided, the windows create will be the principal one
        weight: if true, the first col and row will be configured with weight=1
        minsize: a minimum size for height and width
        icon: if not None, set the windows icon
        """
        if masterWindow is None:
            Window._root = self
            self.root = tk.Tk()
            self._images = {}
        else:
            self.root = tk.Toplevel(masterWindow.root)
            self._images = masterWindow._images
            
        self.root.withdraw()
        self.root.title(title)
        
        if weight:
            configureWeigths(self.root)
        if not minsize is None:
            self.root.minsize(minsize[0], minsize[1])
        if not icon is None:
            path = findResource(icon)          
            abspath = os.path.abspath(path)
            self.root.iconbitmap("@" + abspath)
            
        self.root.protocol("WM_DELETE_WINDOW", self._onClosing)
        self._w, self._h, self._x, self._y = 0, 0, 0, 0
        self.root.bind("<Configure>", self._configure)
        self.master = masterWindow
        setCommonFonts(self)
        
        if kwargs.get('enableQueue', False):
            self.queue = Queue.Queue(maxsize=0)
        else:
            self.queue = None
            
    def __processQueue(self):#called from main frame
        if not self.queue.empty():
            func = self.queue.get(block=False)
            # executes graphic interface function
            func()
        self._queueTimer = self.root.after(500, self.__processQueue)
        
    def enqueue(self, func):
        """ Put some function to be executed in the GUI main thread. """
        self.queue.put(func)
        
    def getRoot(self):
        return self.root
    
    def desiredDimensions(self):
        """Override this method to calculate desired dimensions."""
        return None
    
    def _configure(self, e):
        """ Filter event and call appropriate handler. """
        if self.root != e.widget:
            return
        
        _, _, x, y = getGeometry(self.root)
        w, h = e.width, e.height
        
        if w != self._w or h != self._h:
            self._w, self._h = w, h
            self.handleResize() 
        
        if x != self._x or y != self._y:
            self._x, self._y = x, y
            self.handleMove()    
        
    def handleResize(self):
        """Override this method to respond to resize events."""
        pass

    def handleMove(self):
        """Override this method to respond to move events."""
        pass

    def show(self, center=True):
        """This function will enter in the Tk mainloop"""
        if center:
            if self.master is None:
                refw = None
            else:
                refw = self.master.root
            centerWindows(self.root, dim=self.desiredDimensions(),
                          refWindows=refw)
        self.root.deiconify()
        self.root.focus_set()
        if self.queue is not None:
            self._queueTimer = self.root.after(1000, self.__processQueue)
        self.root.mainloop()
        
    def close(self, e=None):
        self.root.destroy()
        # JMRT: For some reason when Tkinter has an exception
        # it does not exit the application as expected and
        # remains in the mainloop, so here we are forcing
        # to exit the whole system (only applies for the main window)
        if self.master is None:
            import sys 
            sys.exit()
        
    def _onClosing(self):
        """Do some cleaning before closing."""
        if self.master is None: 
            pass
        else:
            self.master.root.focus_set()
        if self.queue is not None:
            self.root.after_cancel(self._queueTimer)
        self.close()
        
    def getImage(self, imgName, percent=100, maxheight=None):
        return getImage(imgName, self._images, percent=percent, maxheight=maxheight)
    
    def createMainMenu(self, menuConfig):
        """Create Main menu from the given MenuConfig object."""
        menu = tk.Menu(self.root)
        self._addMenuChilds(menu, menuConfig)
        self.root.config(menu=menu)
        return menu
        
    def _addMenuChilds(self, menu, menuConfig):
        """Add entries of menuConfig in menu
        (using add_cascade or add_command for sub-menus and final options)."""
        # Helper function to create the main menu.
        for sub in menuConfig:
            menuLabel = sub.text.get()
            if not menuLabel: # empty or None label means a separator
                menu.add_separator()
            elif len(sub) > 0:  # sub-menu
                submenu = tk.Menu(self.root, tearoff=0)
                menu.add_cascade(label=menuLabel, menu=submenu)
                self._addMenuChilds(submenu, sub)  # recursive filling
            else:  # menu option
                # If there is an entry called "Browse files", when clicked it
                # will call the method onBrowseFiles() (it has to be defined!)
                def callback(name):
                    """Return a callback function named "on<Name>"."""
                    f = "on%s" % "".join(x.capitalize() for x in name.split())
                    return lambda: getattr(self, f)()
                menu.add_command(label=menuLabel, compound=tk.LEFT,
                                 image=self.getImage(sub.icon.get()),
                                 command=callback(name=sub.text.get()))

    def showError(self, msg, header="Error"):
        from dialog import showError
        showError(header, msg, self.root)
        
    def showInfo(self, msg, header="Info"):
        from dialog import showInfo
        showInfo(header, msg, self.root)
        
    def showWarning(self, msg, header='Warning'):
        from dialog import showWarning
        showWarning(header, msg, self.root)
        
    def askYesNo(self, title, msg):
        from dialog import askYesNo
        return askYesNo(title, msg, self.root)
        
    def createCloseButton(self, parent):
        """ Create a button for closing the window, setting
        the proper label and icon. 
        """
        return Button(parent, Message.LABEL_BUTTON_CLOSE, Icon.ACTION_CLOSE, 
                          command=self.close)
        
    def configureWeights(self, row=0, column=0):
        configureWeigths(self.root, row, column)

