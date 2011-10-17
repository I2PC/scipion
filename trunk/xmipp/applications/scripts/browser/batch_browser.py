#!/usr/bin/env xmipp_python

from Tkinter import *
import ttk
import os, sys, stat
from os.path import join, getsize, basename
import matplotlib.cm as cm
from numpy import zeros
import xmipp
from protlib_filesystem import getXmippPath
from protlib_utils import runImageJPlugin, pretty_date
from protlib_gui_ext import MyButton

RESOURCES = getXmippPath('resources')

# Some helper functions
def showj(filename, mode):
    runImageJPlugin("512m", "XmippBrowser.txt", "-i %s --mode %s" % (filename, mode), batchMode=True)
    
def chimera(filename):
    os.system('chimera spider:%s &' % filename)

def fileInfo(browser):
    msg = "<size:> %d bytes\n" % browser.stat.st_size
    msg += "<modified:> %s\n" % pretty_date(int(browser.stat.st_mtime))
    return msg
    
# Event handlers for each action and each type of file
def defaultOnClick(filename, browser):
    fn = join(RESOURCES, 'png', 'folder.png')
    print "calling preview with fn: ", fn
    browser.updatePreview(fn)
    return ''
        
def defaultFillMenu( filename, browser):
    return False

def defaultOnDoubleClick(filename, browser):
    os.system('kde-open %s &' % filename)

def getMdString(filename, browser):
    md = xmipp.MetaData(filename)
    labels = md.getActiveLabels()
    msg = "  <labels:>" + ''.join(["\n   - %s" % xmipp.label2Str(l) for l in labels])
    if xmipp.MDL_IMAGE in labels:
        browser.updatePreview(md.getValue(xmipp.MDL_IMAGE, md.firstObject()))
    return msg
    
def mdOnClick(filename, browser):
    if '@' not in filename:
        msg = "<Metadata File>\n"
        blocks = xmipp.getBlocksInMetaDataFile(filename)
        nblocks = len(blocks)
        if nblocks <= 1:
            msg += "  <single block>\n" + getMdString(filename, browser)
        else:
            msg += "  <%d blocks:>" % nblocks + ''.join(["\n  - %s" % b for b in blocks])
            # Insert blocks in metadata as independent items
            if len(browser.tree.get_children(filename)) == 0:
                fm = browser.managers['md']
                for b in blocks:
                    bname = "%s@%s" % (b, filename)                    
                    btext = "%s@%s" % (b, basename(filename))
                    browser.tree.insert(filename, 'end', bname, text=btext, image=fm.image)
    else:
        msg = "<Metadata Block>\n" + getMdString(filename, browser)
    return msg
        
def mdFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'metadata'))
    menu.add_command(label="Open as Images table", command=lambda:showj(filename, 'gallery'))
    menu.add_command(label="Open as ImageJ gallery", command=lambda: showj(filename, 'image'))
    menu.add_separator()
    menu.add_command(label="Delete file", command=None) 
    return True

def mdOnDoubleClick(filename, browser):
    showj(filename, 'metadata')

def imgOnClick(filename, browser):
    x, y, z, n = xmipp.SingleImgSize(filename)
    dimMsg = "<Image>\n  <dimensions:> %(x)d x %(y)d" 
    expMsg = "Columns x Rows "
    if z > 1: 
        dimMsg += " x %(z)d"
        expMsg += " x Slices"
    if n > 1:
        dimMsg += " x %(n)d" 
        expMsg += " x Objects"
    browser.updatePreview(filename)
    return (dimMsg + "\n" + expMsg) % locals()
        
def imgFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'image'))
    menu.add_separator()
    menu.add_command(label="Delete file", command=None) 
    return True

def imgOnDoubleClick(filename, browser):
    showj(filename, 'image')

def stackFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'gallery'))
    menu.add_command(label="Open in ImageJ", command=lambda:showj(filename, 'image'))
    menu.add_command(label="Open as MetaData", command=lambda: showj(filename, 'metadata'))
    menu.add_separator()
    menu.add_command(label="Delete file", command=None) 
    return True

def stackOnDoubleClick(filename, browser):
    showj(filename, 'gallery')
        
def volFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'gallery'))
    menu.add_command(label="Open as ImageJ gallery", command=lambda:showj(filename, 'image'))
    menu.add_command(label="Open with Chimera", command=lambda:chimera(filename))
    menu.add_separator()
    menu.add_command(label="Delete file", command=None) 
    return True

def volOnDoubleClick(filename, browser):
    showj(filename, 'gallery')
    
class FileManager():
    ''' Class to handle different types of files '''
    def __init__(self, **attributes):
        for k, v in attributes.iteritems():
            setattr(self, k, v)
            
class XmippBrowser():
    def __init__(self, initialDir='.', parent=None):
        self.dir = initialDir
        self.parent = parent
        if parent:
            self.root = Toplevel(parent)
        else:
            self.root = Tk()
        
        self.createFileManagers()
        self.createGUI()
        
    def addFileManager(self, key, icon, extensions, 
                       fillMenu=defaultFillMenu, 
                       onClick=defaultOnClick,
                       onDoubleClick=defaultOnDoubleClick):
        img = PhotoImage(file=join(RESOURCES, icon))
        fm = FileManager(key=key, icon=icon, image=img, 
                         fillMenu=fillMenu,
                         onClick=onClick,
                         onDoubleClick=onDoubleClick)
        self.managers[key] = fm
        for ext in extensions:
            self.extSet[ext] = fm
        
    def createFileManagers(self):
        self.managers = {}
        self.extSet = {}
        addFm = self.addFileManager
        self.addFileManager('md', 'md.gif', ['.xmd', '.sel', '.doc', '.ctfparam', '.ctfdat'], 
                            mdFillMenu, mdOnClick, mdOnDoubleClick)
        self.addFileManager('stk', 'stack.gif', ['.stk', '.mrcs'],
                            stackFillMenu, imgOnClick, stackOnDoubleClick)
        self.addFileManager('img', 'image.gif', ['.xmp', '.tif', '.spi'],
                            imgFillMenu, imgOnClick, imgOnDoubleClick)
        self.addFileManager('vol', 'vol.gif', ['.vol'], 
                            volFillMenu, imgOnClick, volOnDoubleClick)
        self.addFileManager('pyfile', 'python_file.gif', ['.py'])
        self.addFileManager('folder', 'folderopen.gif', [])
            
    def createGUI(self):
        self.root.withdraw()
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.minsize(600, 400)
        self.root.title("Xmipp Browser")

        #Create frame for tree
        frame = ttk.Frame(self.root, padding="3 3 12 12")
        frame.grid(column=0, row=0, sticky=(N, W, E, S))
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        tree = ttk.Treeview(frame)#, columns=('items'))
        tree.grid(column=0, row=0, sticky=(N,W, E, S))
        self.tree = tree
        #tree.column('items', width=60, anchor='e')
        #tree.heading('items', text='Items')

        #Create frame for details
        details = ttk.Frame(self.root, padding="3 3 12 12")
        details.grid(column=1, row=0, sticky=(N, W, E, S))
        details.columnconfigure(0, weight=1)
        details.rowconfigure(0, weight=1, minsize=100)  
        details.rowconfigure(1, weight=1)
        # Add details preview
        self.canvas = None
        #self.createPreviewCanvas(details)
        #
        # Add details text
        from protlib_gui_ext import TaggedText
        self.text = TaggedText(details, width=20, height=5)
        self.text.grid(column=0, row=1, sticky=(N, W, E, S))
        self.details = details
        
        #Create bottom frame
        filterFrame = ttk.Frame(self.root, padding="3 3 0 0")
        filterFrame.grid(column=0, row=1, columnspan=2)
        ttk.Label(filterFrame, text="Filter").grid(column=0, row=0, sticky=W, padx=2)
        self.filterVar = StringVar()
        filter = ttk.Entry(filterFrame, width=25, textvariable=self.filterVar)
        filter.grid(column=1, row=0, sticky=W, padx=2)
        btnFilter = MyButton(filterFrame, "Search", 'search.gif')
        btnFilter.grid(column=2, row=0, sticky=W, padx=2)
        btnOk = MyButton(filterFrame, "Select")
        btnOk.grid(column=3, row=0, sticky=E, padx=2)
        btnOk = MyButton(filterFrame, "Cancel")
        btnOk.grid(column=4, row=0, sticky=E, padx=2)
        
        #create menu
        # create a popup menu
        self.menu = Menu(self.root, tearoff=0)
        self.menu.add_command(label="Open", command=self.onDoubleClick)
        self.menu.add_command(label="Redo", command=self.onDoubleClick) 
             
        # add bindings
        tree.bind("<Double-1>", self.onDoubleClick)
        tree.bind("<Return>", self.onDoubleClick)
        tree.bind("<Button-3>", self.onRightClick)
        self.root.bind("<Button-1>", self.onClick)
        self.root.bind("<Key>", self.onKeyPress)
        
        #Create a dictionary with extensions and icon type
        self.insertFiles(self.dir)
        from protlib_gui_ext import centerWindows
        centerWindows(self.root, refWindows=self.parent)
        self.root.deiconify()
        self.root.mainloop()

    def insertElement(self, root, elem, isFolder=False):
        if root == dir: parent = ''
        else: parent = root
        fm = None
        if isFolder:
            fm = self.managers['folder']
        else:
            ext = os.path.splitext(elem)[1]
            if self.extSet.has_key(ext):
                fm = self.extSet[ext]
        if fm:
            self.tree.insert(parent, 'end', join(root, elem), text=elem, image=fm.image)
        else:
            self.tree.insert(parent, 'end', join(root, elem), text=elem)

    def insertFiles(self, dir):
      for root, dirs, files in os.walk(dir, followlinks=True    ):
        #if root != dir:
        #    self.tree.set(root, 'items', len(dirs) + len(files))
        for d in dirs: 
            self.insertElement(root, d, True)
        for f in files:
            if not (f.startswith('.') or f.endswith('~') or f.endswith('.pyc')):
                self.insertElement(root, f)
    
    def unpostMenu(self):
        self.menu.unpost()
        
    def getSelection(self):
        selection = self.tree.selection()
        item = fm = None
        if len(selection) == 1:
            item = selection[0]
            ext = os.path.splitext(item)[1]
            if self.extSet.has_key(ext):
                fm = self.extSet[ext]
        return (item, fm)
        
    #Functions for handling events
    def onDoubleClick(self, e=None):
        item, fm = self.getSelection()
        if fm:
            fm.onDoubleClick(item, self)
    
    def onRightClick(self, e):
        item, fm = self.getSelection()
        if fm:
            self.menu.delete(0, END)
            if fm.fillMenu(item, self):
               self.menu.post(e.x_root, e.y_root)

    def onClick(self, e):
        self.unpostMenu()
        item, fm = self.getSelection()
        if fm:
            msg = ""
            if os.path.exists(item):
                self.stat = os.stat(item)
                msg = fileInfo(self)
            msg += fm.onClick(item, self)  
            self.text.clear()
            self.text.addText(msg)      
        
    def onKeyPress(self, e):
        self.unpostMenu()
    
    def createPreviewCanvas(self):
        import matplotlib
        matplotlib.use('TkAgg')
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure
        dim = 128
        Z = zeros((dim, dim), float)
        self.figure = Figure(figsize=(dim, dim), dpi=1)
        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.details)
        self.canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
        self.figureimg = self.figure.figimage(Z, cmap=cm.gray)#, origin='lower')
        #self.canvas.show()    
    
    def updatePreview(self, filename):
        if not self.canvas:
            self.createPreviewCanvas()
        
        if not filename.endswith('.png'):
            #Read image data through Xmipp
            Z = self.getImageData(filename)        
        else:
            import matplotlib.image as mpimg
            Z = mpimg.imread(filename)
        self.figureimg.set_array(Z)
        self.figureimg.autoscale()
        
        self.canvas.show()   
        print self.figure.get_children()      

    def getImageData(self, filename):
        #TODO: improve by avoiding use of getPixel
        dim = 128
        Z = zeros((dim, dim), float)
        img = xmipp.Image()
        img.readPreview(filename, dim)
        xdim, ydim, z, n = img.getDimensions()
        for x in range(xdim):
            for y in range(ydim):
                Z[x, y] = img.getPixel(x, y)
        return Z

if __name__ == '__main__':
    if len(sys.argv) > 1:
        dir = sys.argv[1]
    else:
        dir = '.'
    XmippBrowser(dir)