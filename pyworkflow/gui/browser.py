# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
In this module a simple ObjectBrowser is implemented.
This class can be subclasses to extend its functionality.
A concrete use of ObjectBrowser is FileBrowser, where the
elements to inspect and preview are files.
"""

import os.path
import stat

import Tkinter as tk

import xmipp
from pyworkflow.utils.properties import Icon
import gui
from pyworkflow.utils import dirname, getHomePath, prettySize, getExt, dateStr
from tree import BoundTree, TreeProvider
from text import TaggedText, openTextFileEditor
from widgets import Button, HotButton


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
                 showPreview=True, showPreviewTop=True,
                 **args):
        tk.Frame.__init__(self, parent, **args)
        self.treeProvider = treeProvider
        self._lastSelected = None
        gui.configureWeigths(self)
        self.showPreviewTop = showPreviewTop
        # The main layout will be two panes, 
        # At the left containing the elements list
        # and the right containing the preview and description
        p = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        p.grid(row=0, column=0, sticky='news')
        
        leftPanel = tk.Frame(p)
        self._fillLeftPanel(leftPanel)
        p.add(leftPanel, padx=5, pady=5)
        p.paneconfig(leftPanel, minsize=300)
        
        if showPreview:
            rightPanel = tk.Frame(p)            
            self._fillRightPanel(rightPanel)
            p.add(rightPanel, padx=5, pady=5)    
            p.paneconfig(rightPanel, minsize=200)    
        
            # Register a callback when the item is clicked
            self.tree.itemClick = self._itemClicked
        
    def _fillLeftPanel(self, frame):
        gui.configureWeigths(frame)
        self.tree = BoundTree(frame, self.treeProvider)
        self.tree.grid(row=0, column=0, sticky='news')
        self.itemConfig = self.tree.itemConfig
        self.getImage = self.tree.getImage
    
    def _fillRightPanel(self, frame):
        frame.columnconfigure(0, weight=1)
        
        if self.showPreviewTop:
            top = tk.Frame(frame)
            top.grid(row=0, column=0, sticky='news')
            frame.rowconfigure(0, weight=3)
            gui.configureWeigths(top)
            top.rowconfigure(0, minsize=200)
            self._fillRightTop(top)
        
        bottom = tk.Frame(frame)
        bottom.grid(row=1, column=0, sticky='news')
        frame.rowconfigure(1, weight=1)
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
        self._lastSelected = obj
        img, desc = self.treeProvider.getObjectPreview(obj)
        # Update image preview
        if self.showPreviewTop:
            if isinstance(img, str):
                img = self.getImage(img)
            if img is None:
                img = self.noImage
            self.label.config(image=img)
        # Update text preview
        self.text.setReadOnly(False)
        self.text.clear()
        if desc is not None:
            self.text.addText(desc)
        self.text.setReadOnly(True)
            
    def getSelected(self):
        """ Return the selected object. """
        return self._lastSelected
      

#------------- Classes and Functions related to File browsing --------------

class FileInfo(object):
    """ This class will store some information about a file.
    It will serve to display files items in the Tree.
    """
    def __init__(self, path, filename):
        self._fullpath = os.path.join(path, filename)
        self._filename = filename
        if os.path.exists(self._fullpath):
            self._stat = os.stat(self._fullpath)
        else:
            self._stat = None
        
    def isDir(self):
        return stat.S_ISDIR(self._stat.st_mode) if self._stat else False
    
    def getFileName(self):
        return self._filename
    
    def getPath(self):
        return self._fullpath
    
    def getSize(self):
        return self._stat.st_size if self._stat else 0
    
    def getSizeStr(self):
        """ Return a human readable string of the file size."""
        return prettySize(self.getSize()) if self._stat else '0'

    def getDateStr(self):
        return dateStr(self.getDate()) if self._stat else '0'

    def getDate(self):
        return self._stat.st_mtime if self._stat else 0

class FileHandler(object):
    """ This class will be used to get the icon, preview and info
    from the different types of objects.
    It should be used with FileTreeProvider, where different
    types of handlers can be registered.
    """
    def getFileIcon(self, objFile):
        """ Return the icon name for a given file. """
        if objFile.isDir():
            icon = 'file_folder.gif'
        else:
            icon = 'file_generic.gif'
        
        return icon
    
    def getFilePreview(self, objFile):
        """ Return the preview image and description for the specific object."""
        if objFile.isDir():
            return 'fa-folder-open.png', None
        return None, None
    
    def getFileActions(self, objFile):
        """ Return actions that can be done with this object.
        Actions will be displayed in the context menu 
        and the first one will be the default when double-click.
        """
        return []
    
    
class TextFileHandler(FileHandler):   
    def __init__(self, textIcon):
        FileHandler.__init__(self)
        self._icon = textIcon
         
    def getFileIcon(self, objFile):
        return self._icon
    
    
class SqlFileHandler(FileHandler):
    def getFileIcon(self, objFile):
        return 'file_sqlite.gif'    
    
    
class ImageFileHandler(FileHandler):
    _image = xmipp.Image()
    _index = ''
    
    def _getImageString(self, filename):
        if isStandardImage(filename):
            return "Image file."
        x, y, z, n = xmipp.getImageSize(filename)
        objType = 'Image'
        dimMsg = "*%(objType)s file*\n  dimensions: %(x)d x %(y)d" 
        expMsg = "Columns x Rows "
        if z > 1: 
            dimMsg += " x %(z)d"
            expMsg += " x Slices"
            objType = 'Volume'
        if n > 1:
            dimMsg += " x %(n)d" 
            expMsg += " x Objects"
            objType = 'Stack'
        return (dimMsg + "\n" + expMsg) % locals()
    
    def _getImagePreview(self, filename):
        dim = 128
        
        if isStandardImage(filename):
            self.tkImg = gui.getImage(os.path.abspath(filename),
                                      tkImage=True, maxheight=dim)
        else:
            fn = self._index + filename
            self.tkImg = gui.getTkImage(self._image, fn, dim)
            
        return self.tkImg
        
    def getFilePreview(self, objFile):
        fn = objFile.getPath()
        return self._getImagePreview(fn), self._getImageString(fn)
    
    def getFileActions(self, objFile):
        from pyworkflow.em.viewer import DataView
        fn = objFile.getPath()
        return [('Open with Xmipp viewer', lambda : DataView(fn).show(),
                 Icon.ACTION_VISUALIZE)]
    
    
class ParticleFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_image.gif'
    
class VolFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_vol.gif'
    
class StackHandler(ImageFileHandler):
    _index = '1@'
    
    def getFileIcon(self, objFile):
        return 'file_stack.gif'
    

class ChimeraHandler(FileHandler):
    
    def getFileActions(self, objFile):
        from pyworkflow.em.viewer import ChimeraView
        fn = objFile.getPath()
        return [('Open with Chimera', lambda : ChimeraView(fn).show(),
                 Icon.ACTION_VISUALIZE)]
    
    def getFileIcon(self, objFile):
        return 'file_text.gif' 
    
    
class MdFileHandler(ImageFileHandler):
    def getFileIcon(self, objFile):
        return 'file_md.gif'
    
    def _getImgPath(self, mdFn, imgFn):
        """ Get ups and ups until finding the relative location to images. """
        path = dirname(mdFn)
        index, fn = xmipp.FileName(imgFn).decompose()
        
        while path and path != '/':
            newFn = os.path.join(path, fn)
            if os.path.exists(newFn):
                if index:
                    newFn = '%d@%s' % (index, newFn)
                return newFn
            path = dirname(path)
            
        return None
            
    def _getMdString(self, filename, block=None):
        md = xmipp.MetaData()
        if block:
            md.read(block + '@' + filename)
        else:
            md.read(filename, 1)
        labels = md.getActiveLabels()
        msg =  "Metadata items: *%d*\n" % md.getParsedLines()
        msg += "Metadata labels: " + ''.join(["\n   - %s" % xmipp.label2Str(l)
                                              for l in labels])
        
        imgPath = None
        for label in labels:
            if xmipp.labelIsImage(label):
                imgPath = self._getImgPath(filename,
                                           md.getValue(label, md.firstObject()))
                break
        if imgPath:
            self._imgPreview = self._getImagePreview(imgPath)
            self._imgInfo = self._getImageString(imgPath)
        return msg
    
    def getFilePreview(self, objFile):
        self._imgPreview = None
        self._imgInfo = None
        filename = objFile.getPath()
        ext = getExt(filename)
        
        if ext == '.xmd' or ext == '.ctfparam' or ext == '.pos' or ext == '.doc':
            msg = "*Metadata File* "
            blocks = xmipp.getBlocksInMetaDataFile(filename)
            nblocks = len(blocks)
            if nblocks <= 1:
                mdStr = self._getMdString(filename)
                msg += "  (single block)\n"
                if self._imgInfo:
                    msg += "\nFirst item: \n" + self._imgInfo
                msg += '\n' + mdStr
            else:
                mdStr = self._getMdString(filename, blocks[0])
                msg += "  (%d blocks) " % nblocks
                if self._imgInfo:
                    msg += "\nFirst item: \n" + self._imgInfo
                msg += "\nFirst block: \n" + mdStr
                msg += "\nAll blocks:" + ''.join(["\n  - %s" % b for b in blocks])
        elif ext == '.star':
            msg = "*Relion STAR file* \n"
            msg += self._getMdString(filename)
            
        return self._imgPreview, msg
    
    
class FileTreeProvider(TreeProvider):
    """ Populate a tree with files and folders of a given path """
    
    _FILE_HANDLERS = {}
    _DEFAULT_HANDLER = FileHandler()
    FILE_COLUMN = 'File'
    SIZE_COLUMN = 'Size'

    @classmethod
    def registerFileHandler(cls, fileHandler, *extensions):
        """ Register a FileHandler for a given file extension. 
        Params:
            fileHandler: the FileHandler that will take care of extensions.
            *extensions: the extensions list that will be associated to this
                FileHandler.
        """
        for fileExt in extensions:
            cls._FILE_HANDLERS[fileExt] = fileHandler
        
    def __init__(self, currentDir=None, showHidden=False):
        TreeProvider.__init__(self, sortingColumnName=self.FILE_COLUMN)
        self._currentDir = os.path.abspath(currentDir)
        self._showHidden = showHidden
        self.getColumns = lambda: [(self.FILE_COLUMN, 300),
                                   (self.SIZE_COLUMN, 70), ('Time', 150)]
    
    def getFileHandler(self, obj):
        filename = obj.getFileName()
        fileExt = getExt(filename)
        return self._FILE_HANDLERS.get(fileExt, self._DEFAULT_HANDLER)
        
    def getObjectInfo(self, obj):
        filename = obj.getFileName()
        fileHandler = self.getFileHandler(obj)
        icon = fileHandler.getFileIcon(obj)
        
        info = {'key': filename, 'text': filename, 
                'values': (obj.getSizeStr(), obj.getDateStr()), 'image': icon
                }
            
        return info
    
    def getObjectPreview(self, obj):
        fileHandler = self.getFileHandler(obj)
        return fileHandler.getFilePreview(obj)
    
    def getObjectActions(self, obj):
        fileHandler = self.getFileHandler(obj)
        actions = fileHandler.getFileActions(obj)
        # Always allow the option to open as text
        # specially useful for unknown formats
        fn = obj.getPath()
        actions.append(("Open external Editor", 
                        lambda: openTextFileEditor(fn), Icon.ACTION_REFERENCES))
        
        return actions
    
    def getObjects(self):
        files = os.listdir(self._currentDir)
        #files.sort()
        fileInfoList = []
        for f in files:
            if self._showHidden or not f.startswith('.'):
                fileInfoList.append(FileInfo(self._currentDir, f))

        # Sort objects
        fileInfoList.sort(key=self.fileKey, reverse=not self.isSortingAscending())

        return fileInfoList

    def fileKey(self, f):
        sortDict = {self.FILE_COLUMN: 'getFileName',
                    self.SIZE_COLUMN: 'getSize'}
        return getattr(f, sortDict.get(self._sortingColumnName, 'getDate'))()

    def getDir(self):
        return self._currentDir
    
    def setDir(self, newPath):
        self._currentDir = newPath
        
# Some constants for the type of selection
# when the file browser is opened
SELECT_NONE = 0 # No selection, just browse files                  
SELECT_FILE = 1
SELECT_FOLDER = 2
SELECT_PATH = 3 # Can be either file or folder


class FileBrowser(ObjectBrowser):
    """ The FileBrowser is a particular class of ObjectBrowser
    where the "objects" are just files and directories.
    """
    def __init__(self, parent, initialDir='.', 
                 selectionType=SELECT_FILE, 
                 selectionSingle=True, 
                 allowFilter=True, 
                 filterFunction=None, 
                 previewDim=144,
                 showHidden=False, # Show hidden files or not?
                 selectButton='Select', # Change the Select button text
                 entryLabel=None, # Display an entry for some input
                 entryValue='' # Display a value in the entry field
                 ):
        """ 
        """
        self._provider = FileTreeProvider(initialDir, showHidden)
        self.selectButton = selectButton
        self.entryLabel = entryLabel
        self.entryVar = tk.StringVar()
        self.entryVar.set(entryValue)
        
        ObjectBrowser.__init__(self, parent, self._provider)
        
        if selectionType == SELECT_NONE:
            selectButton = None
            
    
        buttonsFrame = tk.Frame(self)
        self._fillButtonsFrame(buttonsFrame)
        buttonsFrame.grid(row=1, column=0)

    def _fillLeftPanel(self, frame):
        """ Redefine this method to include a buttons toolbar and
        also include a filter bar at the bottom of the Tree.
        """
        # Tree with files
        frame.columnconfigure(0, weight=1)
        
        treeFrame = tk.Frame(frame)
        ObjectBrowser._fillLeftPanel(self, treeFrame)
        # Register the double-click event
        self.tree.itemDoubleClick = self._itemDoubleClick
        treeFrame.grid(row=1, column=0, sticky='news')
        # Toolbar frame
        toolbarFrame = tk.Frame(frame)
        self._fillToolbar(toolbarFrame)
        toolbarFrame.grid(row=0, column=0, sticky='new')
        
        # Entry frame, could be used for filter
        if self.entryLabel:  
            entryFrame = tk.Frame(frame)
            entryFrame.grid(row=2, column=0, sticky='new') 
            tk.Label(entryFrame, text=self.entryLabel).grid(row=0, column=0,
                                                            sticky='nw')
            tk.Entry(entryFrame, 
                     textvariable=self.entryVar, 
                     bg='white',
                     width=50).grid(row=0, column=1, sticky='nw')
        
        frame.rowconfigure(1, weight=1)

    def _fillToolbar(self, frame):
        """ Fill the toolbar frame with some buttons. """
        self._col = 0
        
        def addButton(text, image, command):
            btn = tk.Label(frame, text=text, image=self.getImage(image), 
                       compound=tk.LEFT, cursor='hand2')
            btn.bind('<Button-1>', command)
            btn.grid(row=0, column=self._col, sticky='nw',
                     padx=(0, 5), pady=5)
            self._col += 1
            
        addButton('Refresh', Icon.ACTION_REFRESH, self._actionRefresh)
        addButton('Home', Icon.HOME, self._actionHome)
        addButton('Back', Icon.ARROW_LEFT, self._actionUp)
        addButton('Up', Icon.ARROW_UP, self._actionUp)
        
    def _fillButtonsFrame(self, frame):
        """ Add button to the bottom frame if the selectMode
        is distinct from SELECT_NONE.
        """
        Button(frame, "Close", Icon.BUTTON_CLOSE, 
               command=self._close).grid(row=0, column=0, padx=(0,5))   
        if self.selectButton:                     
            HotButton(frame, self.selectButton, Icon.BUTTON_SELECT,
                      command=self._select).grid(row=0, column=1)
                
    def _actionRefresh(self, e=None):
        self.tree.update()
        
    def _goDir(self, newDir):
        self.treeProvider.setDir(newDir)
        self.tree.update()
        
    def _actionUp(self, e=None):
        self._goDir(dirname(self.treeProvider.getDir()))
        
    def _actionHome(self, e=None):
        self._goDir(getHomePath())
        
    def _itemDoubleClick(self, obj):
        if obj.isDir():
            self._goDir(obj.getPath())
        else:
            actions = self._provider.getObjectActions(obj)
            if actions:
                # actions[0] = first Action, [1] = the action callback
                actions[0][1]()
            
    def onClose(self):
        pass
    
    def onSelect(self, obj):
        print obj, "type: ", type(obj)
    
    def _close(self, e=None):
        self.onClose()
        
    def _select(self, e=None):
        self.onSelect(self.getSelected())
        
    def getEntryValue(self):
        return self.entryVar.get()
    
    def getCurrentDir(self):
        return self.treeProvider.getDir()
        
        
class BrowserWindow(gui.Window):
    """ Windows to hold a browser frame inside. """
    def __init__(self, title, master=None, **kwargs):
        if 'minsize' not in kwargs:
            kwargs['minsize'] = (800, 400)
        gui.Window.__init__(self, title, master, **kwargs)
        
    def setBrowser(self, browser, row=0, column=0):
        self.browser = browser
        browser.grid(row=row, column=column, sticky='news')
        self.itemConfig = browser.tree.itemConfig
        
 
STANDARD_IMAGE_EXTENSIONS = ['.png', '.jpg', '.jpeg']

def isStandardImage(filename):
    """ Check if a filename have an standard image extension. """
    fnLower = filename.lower()
    return any(fnLower.endswith(ext) for ext in STANDARD_IMAGE_EXTENSIONS)

       
class FileBrowserWindow(BrowserWindow):
    """ Windows to hold a file browser frame inside. """
    def __init__(self, title, master=None, path=None, 
                 onSelect=None, **kwargs):
        BrowserWindow.__init__(self, title, master, **kwargs)
        self.registerHandlers()
        browser = FileBrowser(self.root, path, **kwargs)
        if onSelect:
            def selected(obj):
                onSelect(obj)
                self.close()
            browser.onSelect = selected
        browser.onClose = self.close
        self.setBrowser(browser) 
        
    def getEntryValue(self):
        return self.browser.getEntryValue()
        
    def getCurrentDir(self):
        return self.browser.getCurrentDir()
        
    def registerHandlers(self):
        register = FileTreeProvider.registerFileHandler # shortcut

        register(TextFileHandler('file_text.gif'),
                 '.txt', '.log', '.out', '.err', '.stdout', '.stderr', '.emx',
                 '.json', '.xml', '.pam')
        register(TextFileHandler('file_python.gif'), '.py')
        register(TextFileHandler('file_java.gif'), '.java')
        register(MdFileHandler(), '.xmd', '.star', '.pos', '.ctfparam', '.doc')
        register(SqlFileHandler(), '.sqlite', '.db')
        register(ParticleFileHandler(),
                 '.xmp', '.tif', '.tiff', '.spi', '.mrc', '.map', '.raw',
                 '.inf', '.dm3', '.em', '.pif', '.psd', '.spe', '.ser', '.img',
                 '.hed', *STANDARD_IMAGE_EXTENSIONS)
        register(VolFileHandler(), '.vol')
        register(StackHandler(), '.stk', '.mrcs', '.st', '.pif')
        register(ChimeraHandler(), '.bild')
    
