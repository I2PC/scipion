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
In this module a simple ObjectBrowser is implemented.
This class can be subclasses to extend its functionality.
A concrete use of ObjectBrowser is FileBrowser, where the
elements to inspect and preview are files.
"""

import os
        
import Tkinter as tk
import ttk

import gui
import tooltip
import widgets
import matplotlib_image
from pyworkflow.utils.path import findResource
from tree import BoundTree
from text import TaggedText

class ObjectBrowser(tk.Frame):
    """ This class will implement a simple object browser.
    Basically, it will display a list of elements at the left
    panel and can display a preview and description on the
    right panel for the selected element.
    An ObjectView will be used to grab information for
    each element such as: icon, preview and description.
    A TreeProvider will be used to populate the list (Tree).
    """
    def __init__(self, parent, treeProvider, 
                 showPreview=True, **args):
        tk.Frame.__init__(self, parent, **args)
        self.treeProvider = treeProvider
        gui.configureWeigths(self)
        # The main layout will be two panes, 
        # At the left containing the elements list
        # and the right containing the preview and description
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        p.grid(row=0, column=0, sticky='news')
        
        leftPanel = tk.Frame(p)
        gui.configureWeigths(leftPanel)
        self._fillLeftPanel(leftPanel)
        p.add(leftPanel, padx=5, pady=5)
        p.paneconfig(leftPanel, minsize=300)
        
        if showPreview:
            rightPanel = tk.Frame(p, bg='blue')
            gui.configureWeigths(rightPanel)
            self._fillRightPanel(rightPanel)
            p.add(rightPanel, padx=5, pady=5)    
            p.paneconfig(rightPanel, minsize=200)    
        
            # Register a callback when the item is clicked
            self.tree.itemClick = self._itemClicked
        
    def _fillLeftPanel(self, frame):
        self.tree = BoundTree(frame, self.treeProvider)
        self.tree.grid(row=0, column=0, sticky='news')
        self.itemConfig = self.tree.itemConfig
        self.getImage = self.tree.getImage
    
    def _fillRightPanel(self, frame):
        top = tk.Frame(frame)
        top.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(top)
        top.rowconfigure(0, minsize=200)
        self._fillRightTop(top)
        
        bottom = tk.Frame(frame)
        bottom.grid(row=1, column=0, sticky='news')
        gui.configureWeigths(bottom)
        bottom.rowconfigure(1, weight=1)
        self._fillRightBottom(bottom)
        
    def _fillRightTop(self, top):
        self.noImage = self.getImage('no-image128.png')
        self.label = tk.Label(top, image=self.noImage)
        self.label.grid(row=0, column=0, sticky='news')
        
    def _fillRightBottom(self, bottom):
        self.text = TaggedText(bottom, width=40, height=15, bg='white')
        self.text.grid(row=0, column=0, sticky='news')
        
    def _itemClicked(self, obj):
        img, desc = self.treeProvider.getObjectPreview(obj)
        self.text.clear()
        img = self.getImage(img)
        if img is None:
            img = self.noImage
        self.label.config(image=img)        
        if desc is not None:
            self.text.addText(desc)
                        
"""#############################################################################"""

class FileBrowser():
    def __init__(self, initialDir='.', parent=None, root=None, seltype="both", selmode="browse", allowFilter=True, filter=None, previewDim=144):
        ''' seltype is the selection type, it could be:
              - file -> only allow files selection
              - folder -> only allow folder selection
              - both -> allow any selection
              - none -> doesn't select, only explore
            selmode is the selection mode, it could be:
              - browse -> only single file selection
              - extended -> multiple file selection
        '''
        self.seltype = seltype
        self.selmode = selmode
        self.dir = abspath(initialDir)
        self.pattern = filter
        self.allowFilter = allowFilter
        self.allowRefresh = True
        self.selectedFiles = None
        self.showPath = True
        self.dim = previewDim
        self.commonRoot = "" # this will be used to avoid display of long path names
        from protlib_filesystem import findProjectPath
        self.projectDir = findProjectPath('.')
        # Check if matplotlib is available
        
#        try: 
#            import protlib_gui_figure
#            self.matplotlibEnabled = True            
#        except ImportError:
#            self.matplotlibEnabled = False
        
#    def addFileManager(self, key, icon, extensions, 
#                       fillMenu=False, 
#                       onClick=defaultOnClick,
#                       onDoubleClick=defaultOnDoubleClick):
#        img = tk.PhotoImage(file=join(RESOURCES, icon))
#        fm = FileManager(key=key, icon=icon, image=img, 
#                         fillMenu=fillMenu,
#                         onClick=onClick,
#                         onDoubleClick=onDoubleClick)
#        self.managers[key] = fm
#        for ext in extensions:
#            self.extSet[ext] = fm
#        
#    def createFileManagers(self):
#        self.managers = {}
#        self.extSet = {}
#        addFm = self.addFileManager
#        addFm('md', 'md.gif', ['.xmd', '.sel', '.doc', '.ctfparam', '.ctfdat', '.pos', '.descr', '.param', '.hist'], 
#                            mdFillMenu, mdOnClick, mdOnDoubleClick)
#        addFm('stk', 'stack.gif', ['.stk', '.mrcs', '.st', '.pif'],
#                            stackFillMenu, imgOnClick, stackOnDoubleClick)
#        addFm('img', 'image.gif', ['.xmp', '.tif', '.tiff', '.spi', '.mrc', '.map', '.raw', '.inf', '.dm3', '.em', '.pif', '.psd', '.spe', '.ser', '.img', '.hed', '.jpeg', '.jpg', '.hdf', '.hdf5', '.h5'],
#                            imgFillMenu, imgOnClick, imgOnDoubleClick)
#        addFm('vol', 'vol.gif', ['.vol', '.mrc', '.map', '.em', '.pif'], 
#                            volFillMenu, imgOnClick, volOnDoubleClick)
#        addFm('text', 'fileopen.gif', TEXT_EXTENSIONS,
#              textFillMenu, defaultOnClick, textOnDoubleClick)
#        addFm('pyfile', 'python_file.gif', ['.py'],textFillMenu, defaultOnClick, textOnDoubleClick)
#        addFm('out', 'out.gif', ['.out'],textFillMenu, defaultOnClick, textOnDoubleClick)
#        addFm('err', 'err.gif', ['.err'],textFillMenu, defaultOnClick, textOnDoubleClick)
#        addFm('log', 'log.gif', ['.log'],textFillMenu, defaultOnClick, textOnDoubleClick)
#        addFm('folder', 'folderopen.gif', [])
#        addFm('default', 'generic_file.gif', [])
#        addFm('up', 'up.gif', [])
#        addFm('pdb', 'pdbSmall.gif', CHIMERA_EXTENSIONS, pdbFillMenu, defaultOnClick, pdbOnDoubleClick)
        
    def createDetailsTop(self, parent):
        self.detailstop = tk.Frame(parent)
        self.detailstop.rowconfigure(0, weight=1)
        self.detailstop.columnconfigure(0, weight=1)
        return self.detailstop
    
    def createDetailsBottom(self, parent):
        # Add details text
        self.text = TaggedText(parent, width=20, height=5)
        self.text.frame.grid(column=0, row=1, sticky="nwes")
        return self.text.frame
            
    def createFilterFrame(self, bottomFrame):
        #Create filter frame
        filterFrame = ttk.Frame(bottomFrame)
        filterFrame.grid(column=0, row=0, sticky='ew')
        ttk.Label(filterFrame, text="Filter").pack(side=tk.LEFT,padx=2)
        self.filterVar = tk.StringVar()
        if self.pattern:
            self.filterVar.set(self.pattern)
        filterEntry = ttk.Entry(filterFrame, width=25, textvariable=self.filterVar)
        filterEntry.pack(side=tk.LEFT,padx=2)
        self.btnFilter = ScipionButton(filterFrame, "Search", 'search.gif', command=self.filterResults)
        filterEntry.bind('<Return>', self.filterResults)
        self.btnFilter.pack(side=tk.LEFT, padx=2)
        
    def createGUI(self, root=None, title='', parent=None):
        if root:
            self.root = root
        else:
            self.root = tk.Tk()
        #Create xmipp image to show preview
        self.preview = None
        self.image = Image()
        self.lastitem = None
        
        self.parent = parent   
        self.createFileManagers()
        self.root.withdraw()
        self.root.columnconfigure(0, weight=4, minsize=120)        
        self.root.columnconfigure(1, weight=1, minsize=30)
        self.root.minsize(600, 400)
        titleStr = "Scipion Browser"
        if len(title):
            titleStr += " - %s" % title
        elif self.seltype != 'none':
            titleStr += " - Choose " + self.seltype
        self.root.title(titleStr)

        contentRow = 0
        if self.showPath:
            self.pathFrame = ttk.Frame(self.root, padding="2 2 2 2")
            self.pathFrame.grid(row=contentRow, column=0, columnspan=2, sticky='we')
            self.root.columnconfigure(contentRow, weight=1)
            self.pathVar = tk.StringVar()
            self.pathEntry = tk.Entry(self.pathFrame, textvariable=self.pathVar, 
                                      state="readonly", bd=0, fg='#7D7979')
            self.pathEntry.grid(row=0, column=0, sticky="we")
            self.pathFrame.columnconfigure(0, weight=1)
            contentRow += 1
        self.root.rowconfigure(contentRow, weight=1)
        #Create frame for tree
        frame = ttk.Frame(self.root, padding="3 3 12 12")
        frame.grid(row=contentRow, column=0, sticky="nwes")
        treeRow = 0
        if self.allowRefresh:
            btn = ScipionButton(frame, "Refresh", 'refresh.gif', command=self.refresh, tooltip='Refresh   F5')
            btn.grid(row=0, column=0, padx=(0, 5), sticky='nw')
            self.root.bind("<F5>", self.refresh)
            treeRow = 1
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(treeRow, weight=1)
        tree = Tree(frame, selectmode=self.selmode)#, columns=('items'))
        tree.grid(column=0, row=treeRow, sticky="nwes")
        self.tree = tree
        # Put scroll
        vscroll = AutoScrollbar(frame)
        vscroll.grid(column=1, row=treeRow, sticky='nes')
        vscroll.config(command=tree.yview)
        tree.config(yscrollcommand=vscroll.set)
        # List for hidden tree items, (parent, item)
        self.hidden = []
        #tree.column('items', width=60, anchor='e')
        #tree.heading('items', text='Items')

        #Create frame for details
        details = ttk.Frame(self.root, padding="3 3 3 3")
        details.grid(column=1, row=contentRow, sticky="nwes")
        details.columnconfigure(0, weight=1)
        details.rowconfigure(0, weight=1, minsize=100)  
        details.rowconfigure(1, weight=1)
        top = self.createDetailsTop(details)
        if top:
            top.grid(column=0, row=0, sticky="nwes")
        bottom = self.createDetailsBottom(details)
        if bottom:
            bottom.grid(column=0, row=1, sticky="nwes")
        self.details = details
        
        #Create bottom frame with search bar and buttons
        bottomFrame = ttk.Frame(self.root, padding="3 3 0 0")
        bottomFrame.grid(column=0, row=contentRow+1, columnspan=2, sticky='ew')
        bottomFrame.columnconfigure(0, weight=1)
        bottomFrame.columnconfigure(1, weight=1)

        if self.allowFilter:
            self.createFilterFrame(bottomFrame)
        
        #Create buttons frame
        buttonsFrame = ttk.Frame(bottomFrame)
        buttonsFrame.grid(column=1, row=0, sticky='e')
        if self.seltype != 'none':
            self.btnOk = ScipionButton(buttonsFrame, "Select", command=self.selectFiles)
            self.btnOk.pack(side=tk.LEFT,padx=(0, 5))
            cancelName = "Cancel"
        else:
            cancelName = "Close"
        self.btnCancel = ScipionButton(buttonsFrame, cancelName, command=self.close)
        self.btnCancel.pack(side=tk.LEFT, padx=(0, 10))
        
        #create menu
        # create a popup menu
        self.menu = tk.Menu(self.root, tearoff=0)
        self.menu.add_command(label="Open", command=self.onDoubleClick)
        self.menu.add_command(label="Redo", command=self.onDoubleClick) 
        
        # add bindings
        tree.bind("<Double-1>", self.onDoubleClick)
        tree.bind("<Return>", self.onDoubleClick)
        tree.bind("<Button-3>", self.onRightClick)
        self.root.bind('<<TreeviewSelect>>', self.onClick)
        self.root.bind("<Key>", self.unpostMenu)
        self.root.bind('<FocusOut>', self.unpostMenu)
        self.root.bind('<Button-1>', self.unpostMenu)
        
        #Create a dictionary with extensions and icon type
        self.insertFiles(self.dir)
        #Filter result, if pattern no provided, all files will be listed
        if self.allowFilter:
            self.filterResults()
        
    def showGUI(self, loop=True):        
        centerWindows(self.root, refWindows=self.parent)
        self.root.deiconify()
        if loop:
            self.root.mainloop()

    def insertElement(self, root, elem, isFolder=False):
        if root == self.dir: 
            parent = ''
        else: 
            parent = root
        fm = None
        if elem.endswith('..'):
            fm = self.managers['up']
        elif isFolder:
            fm = self.managers['folder']
        else:
            ext = os.path.splitext(elem)[1]
            fm = self.extSet.get(ext, self.managers['default'])
        return self.tree.insert(parent, 'end', join(root, elem), text=elem, image=fm.image)

    def insertFiles(self, path):
        files = os.listdir(path)
        files.sort()
        for f in files:
            if not (f.startswith('.') or f.endswith('~') or f.endswith('.pyc')):
                self.insertElement(path, f, os.path.isdir(join(path, f)))
    
    def refresh(self, e=None):
        self.filterResults()
        
    def changeDir(self, newDir):
        self.pathVar.set(abspath(newDir))
        self.dir = newDir
        self.tree.clear()
        self.insertElement(newDir, '..', False)
        self.insertFiles(newDir)
        
    def unpostMenu(self, e=None):
        self.menu.unpost()
        
    def getSelection(self):
        selection = self.tree.selection()
        item = fm = None
        if len(selection) == 1:
            item = join(self.commonRoot, selection[0])
            ext = os.path.splitext(item)[1]
            if self.extSet.has_key(ext):
                fm = self.extSet[ext]
            else:
                fm = self.managers['default']
        return (item, fm)
        
    #Functions for handling events
    def onDoubleClick(self, e=None):
        item, fm = self.getSelection()
        if item is not None:
            if item.endswith('..'): # Move two levels up
                self.changeDir(dirname(dirname(item)))
            elif os.path.isdir(item):
                self.changeDir(item)
            if fm:
                fm.onDoubleClick(item, self)
    
    def onRightClick(self, e):
        item, fm = self.getSelection()
        if fm:
            self.menu.delete(0, tk.END)
            if fm.fillMenu(item, self):
                self.menu.post(e.x_root, e.y_root)

    def onClick(self, e):
        item, fm = self.getSelection()    
        if fm:
            msg = ""
            if exists(item):
                self.lastitem = item
                import stat
                self.stat = os.stat(item)
                msg = fileInfo(self)
                if self.seltype != "none":
                    if stat.S_ISDIR(self.stat.st_mode):
                        correct = self.seltype in ['folder', 'both']
                    else:
                        correct = self.seltype in ['file', 'both']
                    if correct:
                        self.btnOk.config(state=tk.NORMAL)
                    else:
                        self.btnOk.config(state=tk.DISABLED)
            msg += fm.onClick(item, self)  
            if self.text:
                self.text.clear()
                self.text.addText(msg) 
        return item
       
    def updatePreview(self, filename):
        if self.matplotlibEnabled:
            if not self.preview:
                from protlib_gui_figure import ImagePreview
                self.preview = ImagePreview(self.detailstop, self.dim)
            
            if not scipionExists(filename):
                prefix, suffix = splitFilename(filename)
                fn = fixPath(suffix, self.dir, self.projectDir)
                #fn = FileName(filename)
                if not fn is None:
                    if prefix:
                        filename = "%s@%s" % (prefix, fn)
                    else:
                        filename = fn
                else:
                    filename = None
                    
            if filename is None:
                Z = getPngData(findResource('no-image.png'))
                
            if not filename.endswith('.png'):
                self.image.readPreview(filename, self.dim)
                if filename.endswith('.psd'):
                    self.image.convertPSD()
                Z = getImageData(self.image)
            else:
                Z = getPngData(filename)
                
            self.preview.updateData(Z)
    
    def filterResults(self, e=None):
        filterValue = self.filterVar.get().replace(',', '')
        self.pattern = filterValue.split()
        foundDirs = {}   
        #self.changeDir(self.dir)
        if len(self.pattern):
            self.tree.clear()
            for root, dirs, files in os.walk(self.dir, followlinks=True):
                if files:
                    files.sort()
                for f in files:
                    if self.matchPattern(f):
                        relRoot = relpath(root, self.dir)
                        rootId = foundDirs.get(relRoot, None)
                        if rootId is None: 
                            rootId = self.insertElement(self.dir, relRoot, True)
                            self.tree.item(rootId, open=tk.TRUE)
                            foundDirs[relRoot] = rootId
                        self.insertElement(rootId, f)
        else:
            self.changeDir(self.dir)
            
    def matchPattern(self, item):
        ##i = basename(item)
        from fnmatch import fnmatch
        for p in self.pattern:
            if fnmatch(item, p):
                return True
        return False
    
    def filterTreeItems(self, parent, items):
        self.tree.clear()
        for i in items:
            if os.path.isdir(i):
                self.filterTreeItems(i, self.tree.get_children(i))
            elif self.matchPattern(i):
                self.tree.see(i)
                self.filterTreeItems(i, self.tree.get_children(i))
            else:
                index = self.tree.index(i)
                self.tree.detach(i)
                self.hidden.append((parent, i, index))
        
    def close(self, e=None):
        self.root.destroy()
        
    def selectFiles(self, e=None):
        self.selectedFiles = self.tree.selection()
        self.root.destroy()

"""#############################################################################"""

class FileManager():
    ''' Class to handle different types of files '''
    def __init__(self, **attributes):
        for k, v in attributes.iteritems():
            setattr(self, k, v)
            
class ScipionButton(tk.Button):
    def __init__(self, master, text, imagePath=None, tooltip=None, **opts):
        defaults = {'activebackground': "LightSkyBlue", 'bg':"LightBlue"}
        defaults.update(opts)
        btnImage = getImage(imagePath)
        
        if btnImage:
            #height=28, width=28,
            if not opts.has_key('bg'):
                del defaults['bg']
            tk.Button.__init__(self, master, image=btnImage, bd=0,  **defaults)
            self.image = btnImage
        else:
            tk.Button.__init__(self, master, text=text, font="Verdana", **defaults)
            
        if tooltip:
            ToolTip(self, tooltip, 500)
            
    def setImage(self, imagePath):
        self.image = getImage(imagePath)
        self.config(image=self.image)

def fileInfo(browser):
    msg =  "<size:> %s\n" % pretty_size(browser.stat.st_size)
    msg += "<modified:> %s\n" % pretty_date(int(browser.stat.st_mtime))
    return msg

def pretty_date(time=False):
    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """
    from datetime import datetime
    now = datetime.now()
    if type(time) is int:
        diff = now - datetime.fromtimestamp(time)
    elif isinstance(time,datetime):
        diff = now - time 
    elif not time:
        diff = now - now
    second_diff = diff.seconds
    day_diff = diff.days

    if day_diff < 0:
        return ''

    if day_diff == 0:
        if second_diff < 10:
            return "just now"
        if second_diff < 60:
            return str(second_diff) + " seconds ago"
        if second_diff < 120:
            return  "a minute ago"
        if second_diff < 3600:
            return str( second_diff / 60 ) + " minutes ago"
        if second_diff < 7200:
            return "an hour ago"
        if second_diff < 86400:
            return str( second_diff / 3600 ) + " hours ago"
    if day_diff == 1:
        return "Yesterday"
    if day_diff < 7:
        return str(day_diff) + " days ago"
    if day_diff < 31:
        return str(day_diff/7) + " weeks ago"
    if day_diff < 365:
        return str(day_diff/30) + " months ago"
    return str(day_diff/365) + " years ago"

def pretty_size(size):
    """Human friendly file size"""
    from math import log
    unit_list = zip(['bytes', 'kB', 'MB', 'GB', 'TB', 'PB'], [0, 0, 1, 2, 2, 2])
    if size > 1:
        exponent = min(int(log(size, 1024)), len(unit_list) - 1)
        quotient = float(size) / 1024**exponent
        unit, num_decimals = unit_list[exponent]
        format_string = '{:.%sf} {}' % (num_decimals)
        return format_string.format(quotient, unit)
    if size == 0:
        return '0 bytes'
    if size == 1:
        return '1 byte'
    
def splitFilename(filename):
    ''' Split filename separating by @ 
    separating in block and filename'''
    if '@' in filename:
        block, filename = filename.split('@')
    else:
        block = None
    return block, filename

def fixPath(filename, *pathList):
    if isabs(filename):
        return filename
    for path in pathList:
        filepath = join(path, filename)
        if scipionExists(filepath):
            return filepath
    return None

def scipionExists(path):
    return FileName(path).exists()

def getImageData(img):
    ''' Function to get a matrix from an Image'''
    Z = img.getData()
    return Z


"""Show File Browser and return selected files"""
def showBrowseDialog(path='.', title='', parent=None, main=False, browser=FileBrowser, **args):
    if main:
        root = tk.Tk()
    else:
        root = tk.Toplevel()
    args['initialDir'] = path
    xb = browser(**args)
    xb.createGUI(root, title, parent)
    xb.showGUI(loop=False)
    root.wait_window(root)
    return xb.selectedFiles
    
if __name__ == '__main__':
    showBrowseDialog(path='.', main=True, seltype="none")
