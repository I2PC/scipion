'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
 *  e-mail address 'jmdelarosa@cnb.csic.es'
 ***************************************************************************/
 '''

import os
from os.path import join, exists, basename, dirname, commonprefix, relpath, abspath
import Tkinter as tk

from tkSimpleDialog import Dialog
import ttk
from config_protocols import LabelBgColor, ButtonBgColor, ButtonActiveBgColor, SectionTextColor
from protlib_filesystem import getXmippPath, xmippExists, removeFilenamePrefix, fixPath, splitFilename, getExt
from Tkinter import TclError
from protlib_gui_ext import *

      

'''**************  Implementation of Xmipp Browser *************'''
# Some helper functions for the browser
def showj(filename, mode="default"):
    from protlib_utils import runShowJ
    runShowJ(filename, extraParams="--mode %s" % mode)
    
def chimera(filename):
    runChimera(filename)
    
def vmd(filename):
    runVMD(filename)
    
def fileInfo(browser):
    from protlib_utils import pretty_date, pretty_size
    msg =  "<size:> %s\n" % pretty_size(browser.stat.st_size)
    msg += "<modified:> %s\n" % pretty_date(int(browser.stat.st_mtime))
    return msg
    
# Event handlers for each action and each type of file
def defaultOnClick(filename, browser):
    import stat
    mode = browser.stat.st_mode
    if stat.S_ISDIR(mode):
        img = 'folder.png'
        files = os.listdir(filename)
        total = 0
        for f in files:
            st = os.stat(join(filename, f))
            total += st.st_size
        from protlib_utils import pretty_size
        msg = "<%d> items: %s\n" % (len(files), pretty_size(total)) 
    else:
        img = 'file.png'
        msg = ''
    fn = join(RESOURCES, img)
    browser.updatePreview(fn)
    return msg
        
def defaultFillMenu(filename, browser):
    return False

def defaultOnDoubleClick(filename, browser):
    pass

def isChimeraSession(filename):
    ''' check if the file is a chimera .py session '''
    # Check if it is a chimera .py session
    if filename.endswith('.py'):
        f = open(filename)
        for l in f:
            if 'chimera' in l and 'import' in l:
                return True            
    return False
    
def textFillMenu(filename, browser):
    menu = browser.menu
    menu.add_command(label="Open as Text", command=lambda: showTextfileViewer(filename, [filename], browser.parent))
    if isChimeraSession(filename):
        menu.add_command(label="Open with Chimera", command=lambda:chimera(filename))
    return True

def textOnDoubleClick(filename, browser):
    filelist = [filename]
    loglist = ['.log', '.out', '.err']
    prefix, ext = os.path.splitext(filename)
    if ext in loglist:
        filelist = [prefix + ext for ext in loglist if os.path.exists(prefix + ext)]
    showTextfileViewer(filename, filelist, browser.parent)

def getMdString(filename, browser):
    from xmipp import MetaData, MDL_IMAGE, label2Str, labelIsImage
    md = MetaData()
    md.read(filename, 1)
    labels = md.getActiveLabels()
    msg =  "  <%d items>\n" % md.getParsedLines()
    msg += "  <labels:>" + ''.join(["\n   - %s" % label2Str(l) for l in labels])
    
    img = 'no-image.png'
    for label in labels:
        if labelIsImage(label):
            img = md.getValue(label, md.firstObject())
            break
    browser.updatePreview(img)
    return msg
    
def mdOnClick(filename, browser):
    if '@' not in filename:
        import xmipp
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
                
                #bnameSuffix = "@" + relpath(filename, browser.dir)
                bnameSuffix = "@" + filename
                btextSuffix = "@" + basename(filename)
                for b in blocks:
                    bname = b + bnameSuffix                    
                    btext = b + btextSuffix
                    browser.tree.insert(filename, 'end', bname, text=btext, image=fm.image)
    else:
        block, filename = splitFilename(filename)
        filename = join(browser.dir, filename)
        msg = "<Metadata Block>\n" + getMdString("%s@%s" % (block, filename), browser)
    return msg
        
def mdFillMenu(filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'metadata'))
    menu.add_command(label="Open as Images table", command=lambda:showj(filename, 'gallery'))
    menu.add_command(label="Open as ImageJ gallery", command=lambda: showj(filename, 'image'))
    menu.add_command(label="Open as Text", command=lambda: showTextfileViewer(filename, [filename], browser.parent))
    try:
        from xmipp import MetaData, MDL_MICROGRAPH
        md = MetaData(filename)
        if md.containsLabel(MDL_MICROGRAPH):
            menu.add_command(label="PSD preview", command=lambda: showCTFPreview(filename, parent=browser.parent, md=md))
    except Exception, e:
        print e
        
    return True

def mdOnDoubleClick(filename, browser):
    if '@' in filename:
        block, filename = splitFilename(filename)
        filename = join(browser.dir, filename)
        filename = '%(block)s@%(filename)s' % locals()
    showj(filename, 'metadata')

def imgOnClick(filename, browser):
    import xmipp
    x, y, z, n = xmipp.getImageSize(filename)
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
    return True

def imgOnDoubleClick(filename, browser):
    showj(filename, 'image')

def stackFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'gallery'))
    menu.add_command(label="Open in ImageJ", command=lambda:showj(filename, 'image'))
    menu.add_command(label="Open as MetaData", command=lambda: showj(filename, 'metadata'))
    return True

def stackOnDoubleClick(filename, browser):
    showj(filename, 'gallery')
        
def volFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open", command=lambda: showj(filename, 'gallery'))
    menu.add_command(label="Open as ImageJ gallery", command=lambda:showj(filename, 'image'))
    menu.add_command(label="Open with Chimera", command=lambda:chimera(filename))
    menu.add_command(label="Open mask wizard", command=lambda:showBrowseDialog(parent=browser.parent, browser=XmippBrowserMask, 
                                                                               allowFilter=False, extra={'fileList': [filename]}))
    return True

def volOnDoubleClick(filename, browser):
    showj(filename, 'gallery')
    
def pdbFillMenu( filename, browser):
    menu = browser.menu
    menu.add_command(label="Open as text", command=lambda: showTextfileViewer(filename, [filename], browser.parent))
    menu.add_command(label="Open with Chimera", command=lambda:chimera(filename))
    menu.add_command(label="Open with VMD", command=lambda:vmd(filename))
    return True
        
def pdbOnDoubleClick(filename, browser):
    chimera(filename)
    
class FileManager():
    ''' Class to handle different types of files '''
    def __init__(self, **attributes):
        for k, v in attributes.iteritems():
            setattr(self, k, v)


TEXT_EXTENSIONS = ['.txt', '.c', '.h', '.cpp', '.java', '.sh', '.star', '.emx']
CHIMERA_EXTENSIONS = ['.pdb']

class XmippBrowser():
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
        try: 
            import protlib_gui_figure
            self.matplotlibEnabled = True            
        except ImportError:
            self.matplotlibEnabled = False
        
    def addFileManager(self, key, icon, extensions, 
                       fillMenu=defaultFillMenu, 
                       onClick=defaultOnClick,
                       onDoubleClick=defaultOnDoubleClick):
        img = tk.PhotoImage(file=join(RESOURCES, icon))
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
        addFm('md', 'md.gif', ['.xmd', '.sel', '.doc', '.ctfparam', '.ctfdat', '.pos', '.descr', '.param', '.hist'], 
                            mdFillMenu, mdOnClick, mdOnDoubleClick)
        addFm('stk', 'stack.gif', ['.stk', '.mrcs', '.st', '.pif'],
                            stackFillMenu, imgOnClick, stackOnDoubleClick)
        addFm('img', 'image.gif', ['.xmp', '.tif', '.tiff', '.spi', '.mrc', '.map', '.raw', '.inf', '.dm3', '.em', '.pif', '.psd', '.spe', '.ser', '.img', '.hed', '.jpeg', '.jpg', '.hdf', '.hdf5', '.h5'],
                            imgFillMenu, imgOnClick, imgOnDoubleClick)
        addFm('vol', 'vol.gif', ['.vol', '.mrc', '.map', '.em', '.pif'], 
                            volFillMenu, imgOnClick, volOnDoubleClick)
        addFm('text', 'fileopen.gif', TEXT_EXTENSIONS,
              textFillMenu, defaultOnClick, textOnDoubleClick)
        addFm('pyfile', 'python_file.gif', ['.py'],textFillMenu, defaultOnClick, textOnDoubleClick)
        addFm('out', 'out.gif', ['.out'],textFillMenu, defaultOnClick, textOnDoubleClick)
        addFm('err', 'err.gif', ['.err'],textFillMenu, defaultOnClick, textOnDoubleClick)
        addFm('log', 'log.gif', ['.log'],textFillMenu, defaultOnClick, textOnDoubleClick)
        addFm('folder', 'folderopen.gif', [])
        addFm('default', 'generic_file.gif', [])
        addFm('up', 'up.gif', [])
        addFm('pdb', 'pdbSmall.gif', CHIMERA_EXTENSIONS, pdbFillMenu, defaultOnClick, pdbOnDoubleClick)
        
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
        self.btnFilter = XmippButton(filterFrame, "Search", 'search.gif', command=self.filterResults)
        filterEntry.bind('<Return>', self.filterResults)
        self.btnFilter.pack(side=tk.LEFT, padx=2)
        
    def createGUI(self, root=None, title='', parent=None):
        if root:
            self.root = root
        else:
            self.root = tk.Tk()
        #Create xmipp image to show preview
        import xmipp
        self.preview = None
        self.image = xmipp.Image()
        self.lastitem = None
        
        self.parent = parent   
        self.createFileManagers()
        self.root.withdraw()
        self.root.columnconfigure(0, weight=4, minsize=120)        
        self.root.columnconfigure(1, weight=1, minsize=30)
        self.root.minsize(600, 400)
        titleStr = "Xmipp Browser"
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
            btn = XmippButton(frame, "Refresh", 'refresh.gif', command=self.refresh, tooltip='Refresh   F5')
            btn.grid(row=0, column=0, padx=(0, 5), sticky='nw')
            self.root.bind("<F5>", self.refresh)
            treeRow = 1
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(treeRow, weight=1)
        tree = XmippTree(frame, selectmode=self.selmode)#, columns=('items'))
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
            self.btnOk = XmippButton(buttonsFrame, "Select", command=self.selectFiles)
            self.btnOk.pack(side=tk.LEFT,padx=(0, 5))
            cancelName = "Cancel"
        else:
            cancelName = "Close"
        self.btnCancel = XmippButton(buttonsFrame, cancelName, command=self.close)
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
            
            if not xmippExists(filename):
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
                filename = getXmippPath('resources', 'no-image.png')
                
            if not filename.endswith('.png'):
                self.image.readPreview(filename, self.dim)
                if filename.endswith('.psd'):
                    self.image.convertPSD()
                from protlib_xmipp import getImageData
                Z = getImageData(self.image)
            else:
                from protlib_gui_figure import getPngData
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

  
def initXmippBrowser(self, **args):
        if args.has_key('extra'):
            for k, v in  args['extra'].iteritems():
                setattr(self, k, v)
            del args['extra']
        XmippBrowser.__init__(self, **args)
          
class XmippBrowserPreview(XmippBrowser):
    ''' This subclass is specific preview some operations
        the extra dict will be used for personalized parameters
        that will not be passed to XmippBrowser constructor
    '''
    def __init__(self, **args):
        initXmippBrowser(self, **args)
        self.clearCallbacks = []
        self.fillCallbacks = []
        
    def createDetailsTop(self, parent):
        XmippBrowser.createDetailsTop(self, parent)
        self.root.minsize(800, 400)        
        from protlib_gui_figure import ImagePreview
        self.preview = ImagePreview(self.detailstop, self.dim, label=self.previewLabel)
        #self.frame2 = tk.Frame(self.detailstop, bg='red')
        self.detailstop.columnconfigure(1, weight=1)
#        self.frame2.grid(column=1, row=0, sticky='nsew', padx=5, pady=5)
#        self.frame2.columnconfigure(0, weight=1)
#        self.frame2.rowconfigure(0, weight=1)
        # Create result preview and clear it
        self.resultPreview = self.createResultPreview()
        self.clearResultPreview()
                
        return self.detailstop
    
    def clearResultPreview(self):
        if self.resultPreview:
            self.resultPreview.clear()
        for c in self.clearCallbacks: c()
        
    def fillResultPreview(self, e=None):
        #Read image data through Xmipp
        if self.lastitem:
            from protlib_xmipp import getImageData
            f = self.getComputeFunction()
            FlashMessage(self.root, self.computingMessage, func=f)
            #bandPassFilter(self.image, self.lastitem, 0.2, 0.4, downsampling, self.dim)
            Z = getImageData(self.image)
            self.resultPreview.updateData(Z)  
            for c in self.fillCallbacks: c() 
        else:
            showWarning("Operation failed","Select an image to preview", self.root)
    
    def onClick(self, e):
        self.lastitem = join(self.commonRoot, self.tree.selection_first())
        self.updatePreview(self.lastitem)
        self.clearResultPreview()
    
    def selectFiles(self, e=None):
        self.selectedFiles = self.getResults()
        self.root.destroy()   

    def createResultPreview(self):
        pass
    
    def getComputeFunction(self):
        pass
    
    def getResults(self):
        pass
        
    def createDetailsBottom(self, parent):
        pass
                             

'''Show Xmipp Browser and return selected files'''
def showBrowseDialog(path='.', title='', parent=None, main=False, browser=XmippBrowser, **args):
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

def showTextfileViewer(title, filelist, parent=None, main=False):
    if main:
        root = tk.Tk()
    else:
        root = tk.Toplevel()
    root.withdraw()
    root.title(title)
    from xmipp import FileName
    files = [FileName(f).removeBlockName() for f in filelist]
    l = TextfileViewer(root, files)
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    l.grid(column=0, row=0, sticky='nsew')
    centerWindows(root, refWindows=parent)
    root.deiconify()
    root.mainloop()
    
def showTiltPairsDialog(pairsList, parent=None):
    root = tk.Toplevel()
    tiltPairs = TiltPairsDialog(pairsList)
    tiltPairs.createGUI(root, parent)
    tiltPairs.showGUI()
    root.wait_window(root)
    return tiltPairs.result

def showTable(columns, rows, title='Table', root=None, width=50):
    if root is None:
        root = tk.Tk()
    else:
        root = tk.Toplevel()
    root.withdraw()
    root.title(title)
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    
    tree = ttk.Treeview(root, show='headings')
    tree['columns'] = columns
    for c in columns:
        tree.column(c, width=width)
        tree.heading(c, text=c)
    for r in rows:
        tree.insert('', 'end', text='a', values=r)
    tree.grid(row=0, column=0, sticky='news')
    
    centerWindows(root, refWindows=root)
    root.deiconify()
    root.mainloop()    
    
    
#Helper function to select Downsampling wizards
def showCTFPreview(mdPath, parent=None, md=None):  
    path = os.path.dirname(mdPath)
    from xmipp import MetaData
    if md is None:
        md = MetaData(mdPath)
    results = showBrowseDialog(path=path, browser=XmippBrowserCTF, title="PSD Preview", parent=parent,
                                    seltype="file", selmode="browse", filter=None, previewDim=256, 
                                    extra={'freqs':None, 'downsampling':1, 'previewLabel': 'Micrograph', \
                                           'computingMessage': 'Estimating PSD...', 'md':md}) # a list is returned

