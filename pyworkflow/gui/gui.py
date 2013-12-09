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
import os, sys
import Tkinter as tk
import tkFont

from pyworkflow.object import OrderedObject
from pyworkflow.utils.path import findResource

from os.path import join, exists, basename

"""
Some GUI CONFIGURATION parameters
"""
cfgFontName = "Helvetica"
cfgFontSize = 10  

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
cfgButtonBgColor = "#7D0709"
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
    
def hasFont(fontKey):
    return fontKey in globals()

def aliasFont(fontAlias, fontKey):
    """Set a fontAlias as another alias name of fontKey"""
    g = globals()
    g[fontAlias] = g[fontKey] 
    
def setCommonFonts():
    """Set some predifined common fonts.
    Same conditions of setFont applies here."""
    setFont('fontNormal', family=cfgFontName, size=cfgFontSize)
    aliasFont('fontButton', 'fontNormal')
    setFont('fontBold', family=cfgFontName, size=cfgFontSize, weight='bold')
    setFont('fontItalic', family=cfgFontName, size=cfgFontSize, slant='italic')
    setFont('fontLabel', family=cfgFontName, size=cfgFontSize+1, weight='bold')

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

def getImage(imageName, imgDict=None, tk=True, percent=100):
    """ Search for the image in the RESOURCES path list. """
    if imageName is None:
        return None
    if imgDict is not None and imageName in imgDict:
        return imgDict[imageName]
    imagePath = findResource(imageName)
    image = None
    if imagePath:
        from PIL import Image, ImageTk
        image = Image.open(imagePath)
        if percent != 100: # Display image with other dimensions
            (w, h) = image.size
            fp = float(percent)/100.0
            newSize = int(fp * w), int(fp * h)
            image.thumbnail(newSize, Image.ANTIALIAS)
        if tk:
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
    
def configureWeigths(widget):
    """This function is a shortcut to a common
    used pair of calls: rowconfigure and columnconfigure
    for making childs widgets take the space available"""
    widget.columnconfigure(0, weight=1)
    widget.rowconfigure(0, weight=1)
    
class Window():
    """Class to manage a Tk windows.
    It will encapsulates some basic creation and 
    setup functions. """
    
    def __init__(self, title, masterWindow=None, weight=True, minsize=(500, 300),
                 icon=None, **args):
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
        setCommonFonts()
        
    def desiredDimensions(self):
        """This method should be used by subclasses
        to calculate desired dimensions"""
        return None
    
    def _configure(self, e):
        """ Filter event and call appropiate handler. """
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
        """ This method should be overriden by subclasses
        in order to response to resize event. """
        pass
        
    def handleMove(self):
        """ This method should be overriden by subclasses
        in order to response to move envet. """
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
        self.root.mainloop()
        
    def close(self):
        self.root.destroy()
        
    def _onClosing(self):
        """Do some cleanning before closing"""
        if self.master is None: 
            pass
        else:
            self.master.root.focus_set()
        self.close()
        
    def getImage(self, imgName, percent=100):
        return getImage(imgName, self._images, percent=percent)
    
    def createMainMenu(self, menuConfig):
        """Create Main menu from a given configuration"""
        menu = tk.Menu(self.root)
        self._addMenuChilds(menu, menuConfig)
        self.root.config(menu=menu)
        return menu
        
    def _addMenuChilds(self, menu, menuConfig):
        """Helper function for creating main menu"""
        for sub in menuConfig:
            print "sub",sub
            if len(sub):
                submenu = tk.Menu(self.root, tearoff=0)
                menu.add_cascade(label=sub.text.get(), menu=submenu)
                self._addMenuChilds(submenu, sub)
            else:
                menu.add_command(label=sub.text.get(), compound=tk.LEFT,
                                 image=self.getImage(sub.icon.get()))
                
    def showError(self, msg, header="Error"):
        from dialog import showError
        showError(header, msg, self.root)
        
    def showInfo(self, msg, header="Info"):
        from dialog import showInfo
        showInfo(header, msg, self.root)


VIEW_PROJECTS = 'Projects'
VIEW_PROTOCOLS = 'Protocols'
VIEW_DATA = 'Data'
VIEW_HOSTS = 'Hosts'
VIEW_LIST = [VIEW_PROJECTS, VIEW_PROTOCOLS, VIEW_DATA, VIEW_HOSTS]   

     
class WindowBase(Window):
    """Base Template Window
    It extends from Window and add some layout functions (header and footer)
    """
    def __init__(self, title, masterWindow=None, weight=True, minsize=(500, 300),
                 icon="scipion_bn.xbm", **args):
        Window.__init__(self, title, masterWindow, weight=weight, icon=icon, minsize=(900,500))
        
        content = tk.Frame(self.root)
        content.columnconfigure(0, weight=1)
        content.rowconfigure(1, weight=1)
        content.grid(row=0, column=0, sticky='news')
        self.content = content
        
        Window.createMainMenu(self, self.menuCfg)
        
        self.header = self.createHeaderFrame(content)
        self.header.grid(row=0, column=0, sticky='new')
        
        self.footer = tk.Frame(content, bg='white')
        self.footer.grid(row=1, column=0, sticky='news') 
        
        self.view, self.viewWidget = None, None
        
        
        from pyworkflow.apps.pw_manager import ProjectsView
        from pyworkflow.apps.pw_project import DataView
        from pyworkflow.apps.pw_project_viewprotocols import ProtocolsView
        from pyworkflow.apps.pw_project_viewhosts import HostsView
        
#        self.viewFuncs = {VIEW_PROJECTS: createProjectsView,
#                          VIEW_PROTOCOLS: createProtocolsView,
#                          VIEW_DATA: createDataView,
#                          VIEW_HOSTS: createHostsView
#                          }
        self.viewFuncs = {VIEW_PROJECTS: ProjectsView,
                          VIEW_PROTOCOLS: ProtocolsView,
                          VIEW_DATA: DataView,
                          VIEW_HOSTS: HostsView
                          }
#        self.switchView(VIEW_PROJECTS)
        
#    def __init__(self, path, master=None):   
#        # Load global configuration
#        self.projName = 'Project ' + basename(path)
#        self.projPath = path
#        self.loadProject()
#        self.icon = self.generalCfg.icon.get()
#        self.selectedProtocol = None
#        self.showGraph = False
#        
#        Window.__init__(self, self.projName, master, icon=self.icon, minsize=(900,500))
#        
#        content = tk.Frame(self.root)
#        content.columnconfigure(0, weight=1)
#        content.rowconfigure(1, weight=1)
#        content.grid(row=0, column=0, sticky='news')
#        self.content = content
#        
#        self.createMainMenu()
#        
#        header = self.createHeaderFrame(content)
#        header.grid(row=0, column=0, sticky='new')
#        
#        self.view, self.viewWidget = None, None
#        self.viewFuncs = {VIEW_PROTOCOLS: self.createProtocolsView,
#                          VIEW_DATA: self.createDataView,
#                          VIEW_HOSTS: self.createHostsView
#                          }
#        self.switchView(VIEW_PROTOCOLS)
        
    def createHeaderFrame(self, parent):
        
        """ Create the Header frame at the top of the windows.
        It has (from left to right):
            - Main application Logo
            - Project Name
            - View selection combobox
        """
        header = tk.Frame(parent, bg='white')        
        header.columnconfigure(1, weight=1)
        header.columnconfigure(2, weight=1)
        # Create the SCIPION logo label
        logoImg = self.getImage(self.generalCfg.logo.get(), percent=50)
        logoLabel = tk.Label(header, image=logoImg, 
                             borderwidth=0, anchor='nw', bg='white')
        logoLabel.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        
        # Create the Project Name label
        self.projNameFont = tkFont.Font(size=-28, family='helvetica')
        projLabel = tk.Label(header, text=self.projName if 'projName' in locals() else "", font=self.projNameFont,
                             borderwidth=0, anchor='nw', bg='white', fg='#707070')
        projLabel.grid(row=0, column=1, sticky='sw', padx=(20, 5), pady=10)
        
        # Create view selection frame
        viewFrame = tk.Frame(header, bg='white')
        viewFrame.grid(row=0, column=2, sticky='se', padx=5, pady=10)
        
        # Create gradient
        from widgets import GradientFrame
        GradientFrame(header, height=8, borderwidth=0).grid(row=1, column=0, columnspan=3, sticky='new')
#        viewLabel = tk.Label(viewFrame, text='View:', bg='white')
#        viewLabel.grid(row=0, column=0, padx=5)
#        self.viewVar = tk.StringVar()
#        self.viewVar.set(VIEW_PROTOCOLS)
        
#        viewCombo = ttk.Combobox(viewFrame, textvariable=self.viewVar, state='readonly')
#        viewCombo['values'] = [VIEW_PROTOCOLS, VIEW_DATA, VIEW_HOSTS]
#        viewCombo.grid(row=0, column=1)
#        viewCombo.bind('<<ComboboxSelected>>', self._viewComboSelected)
        
        def addLink(elementText):
            btn = tk.Label(viewFrame, text=elementText, cursor='hand2', fg="#6F3232", bg="white")
            btn.bind('<Button-1>', lambda e:self._viewComboSelected(elementText))
            return btn
        
        def addTube():        
            tube = tk.Label(viewFrame, text="|", fg="#6F3232", bg="white", padx=5)
            return tube
        
        for i, elementText in enumerate(VIEW_LIST):
            btn = addLink(elementText)
            btn.grid(row=0, column=i*2)
            
            if i < len(VIEW_LIST)-1:
                tube = addTube()
                tube.grid(row=0, column=(i*2)+1)
        
        # Create header line
#        headerLine = tk.Frame(header,bg='red', height=10, width=100%)
#        headerLine = Canvas(header, height=-1, bg='red')
#        headerLine.grid(row=1, column=0, columnspan=3, sticky='nw')
        
        return header
    
    
    def _viewComboSelected(self, elementText):
        if elementText != self.view:
            self.switchView(elementText)

    def switchView(self, newView):
        # Destroy the previous view if existing:
        if self.viewWidget:
            self.viewWidget.grid_forget()
            self.viewWidget.destroy()
        # Create the new view
        self.viewWidget = self.viewFuncs[newView](self.footer, self)
        # Grid in the second row (1)
        self.viewWidget.grid(row=0, column=0, columnspan=10, sticky='news')
        self.footer.rowconfigure(0, weight=1)
        self.footer.columnconfigure(0, weight=1)
        #header.columnconfigure(2, weight=1)
        self.view = newView    
#    def _viewComboSelected(self, e=None):
#        if self.viewVar.get() != self.view:
#            self.switchView(self.viewVar.get())

        
