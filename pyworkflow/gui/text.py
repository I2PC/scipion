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
Text based widgets.
"""
import os
import sys
import time
import webbrowser
import subprocess
import Tkinter as tk
import ttk

import gui
from widgets import Scrollable, IconButton
from pyworkflow.utils import (HYPER_BOLD, HYPER_ITALIC, HYPER_LINK1, HYPER_LINK2,
                              parseHyperText, renderLine, renderTextFile, colorName,
                              which, envVarOn, expandPattern)
from pyworkflow.utils.properties import Message, Color, Icon


# Define a function to open files cleanly in a system-dependent way
if sys.platform.startswith('darwin'):  # macs use the "open" command
    _open_cmd = lambda path: subprocess.Popen(['open', path])
elif os.name == 'nt':  # there is a function os.startfile for windows
    _open_cmd = lambda path: os.startfile(path)
elif os.name == 'posix':  # linux systems and so on
    def find_prog(*args):
        "Return the first argument that is a program in PATH"
        for command in args:
            if which(command):
                return command
        return None

    x_open = find_prog('xdg-open', 'gnome-open', 'kde-open', 'gvfs-open')
    editor = find_prog('pluma', 'gedit', 'kwrite', 'geany', 'kate',
                       'emacs', 'nedit', 'mousepad')

    def _open_cmd(path):
        # If it is an url, open with browser.
        if path.startswith('http://') or path.startswith('https://'):
            try:
                webbrowser.open_new_tab(path)
                return
            except:
                pass
        if x_open:  # standard way to open
            proc = subprocess.Popen([x_open, path])
            time.sleep(1)
            if proc.poll() in [None, 0]:
                return  # yay! that's the way to do it!
        if editor:  # last card: try to open it in an editor
            proc = subprocess.Popen([editor, path])
            time.sleep(1)
            if proc.poll() in [None, 0]:
                return  # hope we found your fav editor :)
        print 'WARNING: Cannot open %s' % path  # nothing worked! :(
else:
    def _open_cmd(path):
        print 'Unknown system, so cannot open %s' % path


class HyperlinkManager:
    """ Tkinter Text Widget Hyperlink Manager, taken from:
    http://effbot.org/zone/tkinter-text-hyperlink.htm """
    def __init__(self, text):
        self.text = text
        self.text.tag_config("hyper", foreground=Color.RED_COLOR, underline=1)
        self.text.tag_bind("hyper", "<Enter>", self._enter)
        self.text.tag_bind("hyper", "<Leave>", self._leave)
        self.text.tag_bind("hyper", "<Button-1>", self._click)
        self.reset()

    def reset(self):
        self.links = {}

    def add(self, action):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        return "hyper", tag

    def _enter(self, event):
        self.text.config(cursor="hand2")

    def _leave(self, event):
        self.text.config(cursor="")

    def _click(self, event):
        for tag in self.text.tag_names(tk.CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag]()
                return


class Text(tk.Text, Scrollable):    
    """ Base Text widget with some functionalities 
    that will be used for other extensions.
    """
    def __init__(self, master, **opts):
        if 'handlers' in opts:
            self.handlers = opts.pop('handlers')
        else:
            self.handlers = {}
        opts['font'] = gui.fontNormal
        defaults = self.getDefaults()
        defaults.update(opts)
        Scrollable.__init__(self, master, tk.Text, wrap=tk.WORD, **opts)
        self._createWidgets(master, **defaults)
        self.configureTags()

    def _createWidgets(self, master, **opts):
        """This is an internal function to create the Text, the Scrollbar and the Frame"""
        
        # create a popup menu
        self.menu = tk.Menu(master, tearoff=0, postcommand=self.updateMenu)
        self.menu.add_command(label="Copy to clipboard", command=self.copyToClipboard)
        self.menu.add_command(label="Open path", command=self.openFile)
        # Associate with right click
        self.bind("<Button-1>", self.onClick)
        self.bind("<Button-3>", self.onRightClick)
        
    def getDefaults(self):
        """This should be implemented in subclasses to provide defaults"""
        return {}
    
    def configureTags(self):
        """This should be implemented to create specific tags"""
        pass
    
    def addLine(self, line):
        """Should be implemented to add a line """
        self.insert(tk.END, line + '\n')
        
    def addNewline(self):
        self.insert(tk.END, '\n')
        
    def goBegin(self):
        self.see(0.0)
        
    def goEnd(self):
        self.see(tk.END)
        
    def isAtEnd(self):
        return self.scrollbar.get() == 1.0
        
    def clear(self):
        self.delete(0.0, tk.END)

    def getText(self):
        return self.get(0.0, tk.END)
    
    def setText(self, text):
        """ Replace the current text with new one. """
        self.clear()
        self.addText(text)
        
    def addText(self, text):
        """ Add some text to the current state. """
        if isinstance(text, list):
            for line in text:
                self.addLine(line)
        else:
            for line in text.splitlines():
                self.addLine(line)
        
    def onClick(self, e=None): 
        self.selection = None
        self.selection_clear()
        self.menu.unpost()
        
    def onRightClick(self, e):
        try:
            self.selection = self.selection_get().strip()
            self.menu.post(e.x_root, e.y_root)    
        except tk.TclError, e:
            pass
    
    def copyToClipboard(self, e=None):
        self.clipboard_clear()
        self.clipboard_append(self.selection)

    def openFile(self):
        # What happens when you right-click and select "Open path"
        self.openPath(self.selection)

    def openPath(self, path):
        "Try to open the selected path"
        path = expandPattern(path)

        # If the path is a dir, open it with   scipion browser dir <path>
        if os.path.isdir(path):
            dpath = (path if os.path.isabs(path)
                     else os.path.join(os.getcwd(), path))
            subprocess.Popen(['%s/scipion' % os.environ['SCIPION_HOME'],
                              'browser', 'dir', dpath])
            return

        # If it is a file, interpret it correctly and open it with DataView
        dirname = os.path.dirname(path)
        fname = os.path.basename(path)
        if '@' in fname:
            path = os.path.join(dirname, fname.split('@', 1)[-1])
        else:
            path = os.path.join(dirname, fname)

        if os.path.exists(path):
            import xmipp
            fn = xmipp.FileName(path)
            if fn.isImage() or fn.isMetaData():
                from pyworkflow.em.viewer import DataView
                DataView(path).show()
            else:
                _open_cmd(path)
        else:
            # This is probably one special reference, like sci-open:... that
            # can be interpreted with our handlers.
            tag = path.split(':', 1)[0] if ':' in path else None
            if tag in self.handlers:
                self.handlers[tag](path.split(':', 1)[-1])
            else:
                print "Can't find %s" % path

    def updateMenu(self, e=None):
        state = 'normal'
        #if not xmippExists(self.selection):
        #    state = 'disabled'#self.menu.entryconfig(1, background="green")
        self.menu.entryconfig(1, state=state)
        
    def setReadOnly(self, value):
        state = tk.NORMAL
        if value:
            state = tk.DISABLED
        self.config(state=state) 
        
    def highlight(self, pattern, tag, start="1.0", end="end", regexp=False):
        """ Apply the given tag to all text that matches the given pattern

        If 'regexp' is set to True, pattern will be treated as a regular expression
        Taken from: 
            http://stackoverflow.com/questions/3781670/tkinter-text-highlighting-in-python
        """
        start = self.index(start)
        end = self.index(end)
        self.mark_set("matchStart",start)
        self.mark_set("matchEnd",start)
        self.mark_set("searchLimit", end)

        count = tk.IntVar()
        while True:
            index = self.search(pattern, "matchEnd","searchLimit",
                                count=count, regexp=regexp)
            if index == "": break
            self.mark_set("matchStart", index)
            self.mark_set("matchEnd", "%s+%sc" % (index,count.get()))
            self.tag_add(tag, "matchStart","matchEnd")


def configureColorTags(text):
    """ Create tags in text (of type tk.Text) for all the supported colors. """
    try:
        for color in colorName.values():
            text.tag_config(color, foreground=color)
        return True
    except Exception as e:
        print "Colors still not available (%s)" % e
        return False

       
class TaggedText(Text):  
    """
    Implement a Text that will recognize some basic tags
    *some_text* will display some_text in bold
    _some_text_ will display some_text in italic
    some_link or [[some_link][some_label]] will display some_link as hiperlink or some_label as hiperlink to some_link
    also colors are recognized if set option colors=True
    """           
    def __init__(self, master, colors=True, **opts):  
        self.colors = colors
        Text.__init__(self, master, **opts)
        self.hm = HyperlinkManager(self)

    def getDefaults(self):
        return {'bg': "white", 'bd':0}
        # It used to have also 'font': gui.fontNormal  but that stops
        # this file from running. Apparently there is no fontNormal in gui.

    def configureTags(self):
        self.tag_config('normal', justify=tk.LEFT, font=gui.fontNormal)
        self.tag_config(HYPER_BOLD, justify=tk.LEFT, font=gui.fontBold)
        self.tag_config(HYPER_ITALIC, justify=tk.LEFT, font=gui.fontItalic)
        if self.colors:            
            self.colors = configureColorTags(self) # Color can be unavailable, so disable use of colors    
        
    def openLink(self, link):
        webbrowser.open_new_tab(link)  # Open in the same browser, new tab

    def matchHyperText(self, match, tag):
        """ Process when a match a found and store indexes inside string."""
        self.insert(tk.END, self.line[self.lastIndex:match.start()])
        g1 = match.group(tag)

        if tag == HYPER_BOLD or tag == HYPER_ITALIC:
            self.insert(tk.END, ' ' + g1, tag)
        elif tag == HYPER_LINK1:
            self.insert(tk.END, g1, self.hm.add(lambda: self.openLink(g1)))
        elif tag == HYPER_LINK2:
            label = match.group('link2_label')
            if g1.startswith('http:'):
                self.insert(tk.END, label, self.hm.add(lambda: self.openLink(g1)))
            else:
                self.insert(tk.END, label, self.hm.add(lambda: self.openPath(g1)))
        self.lastIndex = match.end()
        
        return g1

    def addLine(self, line):
        self.line = line
        self.lastIndex = 0
        #if protocol has been executed BUT crashee before create something
        #line=None (ROB)
        if line is not None:
            parseHyperText(line, self.matchHyperText)
            Text.addLine(self, line[self.lastIndex:])


class OutputText(Text):
    """
    Implement a Text that will show file content
    and handle console metacharacter for colored output
    """
    def __init__(self, master, filename, colors=True, t_refresh=0, maxSize=400, **opts):
        """ colors flag indicate if try to parse color meta-characters
            t_refresh is the refresh time in seconds, 0 means no refresh
        """
        self.filename = filename
        self.colors = colors
        self.t_refresh = t_refresh
        self.maxSize = maxSize

        self.refreshAlarm = None  # Identifier returned by after()
        self.lineNo = 0
        self.offset = 0
        self.lastLine = ''
        Text.__init__(self, master, **opts)
        self.hm = HyperlinkManager(self)
        self.doRefresh()

    def getDefaults(self):
        return {'bg': "black", 'fg':'white', 'bd':0,
                'height':30,  'width':100}
        # It used to have also  'font': gui.fontNormal  but that stops this
        # file from running. Apparently there is no fontNormal in gui.
        
    def configureTags(self):
        if self.colors:
            configureColorTags(self)

    def _removeLastLine(self):
        line = int(self.index(tk.END).split('.')[0])
        if line > 0:
            line -= 1
            self.delete('%d.0' % line, tk.END)
        
    def addLine(self, line):
        renderLine(line, self._addChunk, self.lineNo)

    def _addChunk(self, txt, fmt=None):
        """
        Add text txt to the widget, with format fmt.
        fmt can be a color (like 'red') or a link that looks like 'link:url'.
        """
        if self.colors and fmt is not None:
            if fmt.startswith('link:'):
                fname = fmt.split(':', 1)[-1]
                self.insert(tk.END, txt, self.hm.add(lambda: openTextFileEditor(fname)))
            else:
                self.insert(tk.END, txt, fmt)
        else:
            self.insert(tk.END, txt)

    def _notifyLine(self, line):
        if '\r' in self.lastLine and '\r' in line:
            self._removeLastLine()
            self.addNewline()
            
        self.lastLine = line        
        
    def readFile(self, clear=False):
        self.setReadOnly(False)
        
        if clear:
            self.offset = 0
            self.lineNo = 0
            self.clear()

        if os.path.exists(self.filename):
            self.offset, self.lineNo = renderTextFile(self.filename, 
                                                      self._addChunk,
                                                      offset=self.offset, 
                                                      lineNo=self.lineNo,
                                                      maxSize=self.maxSize, 
                                                      notifyLine=self._notifyLine)
        else:
            self.insert(tk.END, "File '%s' doesn't exist" % self.filename)

        self.setReadOnly(True)
        # self.goEnd()
      
    def doRefresh(self):
        # First stop pending refreshes
        if self.refreshAlarm:
            self.after_cancel(self.refreshAlarm)
            self.refreshAlarm = None

        self.readFile()

        if self.t_refresh > 0:
            self.refreshAlarm = self.after(self.t_refresh*1000, self.doRefresh)


class TextFileViewer(tk.Frame):
    """ Implementation of a simple text file viewer """
    
    LabelBgColor = "white"
    
    def __init__(self, master, fileList=[],
                 allowSearch=True, allowRefresh=True, allowOpen=False,
                 font=None, maxSize=400):
        tk.Frame.__init__(self, master)
        self.searchList = None
        self.lastSearch = None
        self.refreshAlarm = None
        self._lastTabIndex = None
        self.fileList = []  # Files being visualized
        self.taList = []  # Text areas (OutputText, a scrollable TkText)
        self.fontDict = {}
        self._allowSearch = allowSearch
        self._allowRefresh = allowRefresh
        self._allowOpen = allowOpen
        self._font = font # allow a font to be passed as argument to be used
        self.maxSize = maxSize

        self.createWidgets(fileList)
        self.master = master
        self.addBinding()
        
    def addFile(self, filename):
        self.fileList.append(filename)
        self._addFileTab(filename)
        
    def clear(self):
        """ Remove all added files. """
        self.fileList = []
        for _ in self.taList:
            self.notebook.forget(0)       
        self.taList = []
        
    def _addFileTab(self, filename):
        tab = tk.Frame(self.notebook)
        tab.rowconfigure(0, weight=1)
        tab.columnconfigure(0, weight=1)
        kwargs = {'bg': 'black',
                  'fg': 'white'}
        
        if self._font is not None:
            kwargs['font'] = self._font

        t = OutputText(tab, filename, width=100, height=30, maxSize=self.maxSize, **kwargs)
        t.frame.grid(column=0, row=0, padx=5, pady=5, sticky='nsew')
        self.taList.append(t)
        tabText = "   %s   " % os.path.basename(filename)
        self.notebook.add(tab, text=tabText)        
    
    def createWidgets(self, fileList):
        #registerCommonFonts()
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)        

        #Create toolbar frame
        toolbarFrame = tk.Frame(self)
        toolbarFrame.grid(column=0, row=0, padx=5, sticky='new')
        gui.configureWeigths(toolbarFrame)
        #Add the search box
        right = tk.Frame(toolbarFrame)
        right.grid(column=1, row=0, sticky='ne')        
        self.searchVar = tk.StringVar()
        if self._allowSearch:
            tk.Label(right, text='Search:').grid(row=0, column=3, padx=5)
            self.searchEntry = tk.Entry(right, textvariable=self.searchVar)
            self.searchEntry.grid(row=0, column=4, sticky='ew', padx=5)
            btn = IconButton(right, "Search", Icon.ACTION_SEARCH, tooltip=Message.TOOLTIP_SEARCH,
                             command=self.findText, bg=None)
            btn.grid(row=0, column=5, padx=(0, 5))
        if self._allowRefresh:
            btn = IconButton(right, "Refresh", Icon.ACTION_REFRESH, tooltip=Message.TOOLTIP_REFRESH, 
                             command=self._onRefresh, bg=None)
            btn.grid(row=0, column=6, padx=(0, 5), pady=2)
        if self._allowOpen:
            btn = IconButton(right, "Open external", Icon.ACTION_REFERENCES, tooltip=Message.TOOLTIP_EXTERNAL,
                             command=self._openExternal, bg=None)
            btn.grid(row=0, column=7, padx=(0, 5), pady=2)

        #Create tabs frame
        tabsFrame = tk.Frame(self)
        tabsFrame.grid(column=0, row=1, padx=5, pady=(0, 5), sticky="nsew")
        tabsFrame.columnconfigure(0, weight=1)
        tabsFrame.rowconfigure(0, weight=1)
        self.notebook = ttk.Notebook(tabsFrame)  
        self.notebook.rowconfigure(0, weight=1)
        self.notebook.columnconfigure(0, weight=1)      
        for f in fileList:
            self._addFileTab(f)
        self.notebook.grid(column=0, row=0, sticky='nsew', padx=5, pady=5)   
        self.notebook.bind('<<NotebookTabChanged>>', self._tabChanged)
        
    def _tabChanged(self, e=None):
        self._lastTabIndex = self.notebook.select() 

    def addBinding(self):
        self.master.bind('<Control_L><Home>', lambda e: self.changePosition(1.0))
        self.master.bind('<Control_L><End>', lambda e: self.changePosition(tk.END))
        #self.master.bind('<Alt_L><c>', lambda e: self.master.destroy())
        #self.master.bind('<Return>', lambda e: self.findText())
        self.master.bind('<Control_L><n>', lambda e: self.findText())
        self.master.bind('<Control_L><p>', lambda e: self.findText(-1))
        #self.master.bind('<Alt_L><plus>', self.changeFont)
        #self.master.bind('<Alt_L><minus>', self.changeFont)
    
    def getIndex(self):
        """ Return the index of the selected tab. """
        selected = self.notebook.select()
        if selected:
            return self.notebook.index(selected)
        return -1
    
    def setIndex(self, index):
        """ Select the tab with the given index. """
        if index != -1:
            self.notebook.select(self.notebook.tabs()[index])
    
    def selectedText(self):
        index = self.getIndex()
        if index != -1:
            return self.taList[index]
        return None
    
    def changeFont(self, event=""):
        for font in self.fontDict.values():
            gui.changeFontSize(font, event)
              
    def refreshAll(self, clear=False, goEnd=False):
        """ Refresh all output textareas. """
        for ta in self.taList:
            ta.readFile(clear)
            if goEnd:
                ta.goEnd()
        if self._lastTabIndex is not None:
            self.notebook.select(self._lastTabIndex)
            
    def _onRefresh(self, e=None):
        """ Action triggered when the 'Refresh' icon is clicked. """
        self.refreshAll(clear=False, goEnd=True)
        
    def refreshOutput(self, e=None):
        if self.refreshAlarm:
            self.after_cancel(self.refreshAlarm)
            self.refreshAlarm = None
        text = self.selectedText()
        if text:
            text.readFile()
        #self.refreshAlarm = self.after(2000, self.refreshOutput)
        
    def changePosition(self, index):
        self.selectedText().see(index)
        
    def findText(self, direction=1):
        text = self.selectedText()
        str = self.searchVar.get()
        if str is None or str != self.lastSearch:
            self.buildSearchList(text, str)
            self.lastSearch = str
        else:
            self.nextSearchIndex(text, direction)
        self.searchEntry.focus_set()
        
    def buildSearchList(self, text, str):
        text.tag_remove('found', '1.0', tk.END)
        list = []
        if str:
            idx = '1.0'
            while True:
                idx = text.search(str, idx, nocase=1, stopindex=tk.END)
                if not idx: break
                lastidx = '%s+%dc' % (idx, len(str))
                text.tag_add('found', idx, lastidx)
                list.append((idx, lastidx))
                idx = lastidx
        text.tag_config('found', foreground='white', background='blue')
        #Set class variables
        self.searchList = list
        self.currentIndex = -1
        self.nextSearchIndex(text) #select first element
    
    def nextSearchIndex(self, text, direction=1):
        #use direction=-1 to go backward
        text.tag_remove('found_current', '1.0', tk.END)
        if len(self.searchList)==0:
            return
        self.currentIndex = (self.currentIndex + direction) % len(self.searchList)
        idx, lastidx = self.searchList[self.currentIndex]
        text.tag_config('found_current', foreground='yellow', background='red')
        text.tag_add('found_current', idx, lastidx)
        text.see(idx)
        
    def _openExternal(self):
        """ Open a new window with an external viewer. """
        if envVarOn('SCIPION_EXTERNAL_VIEWER'):
            if not self.taList:
                return
            openTextFileEditor(self.taList[max(self.getIndex(), 0)].filename)
        else:
            showTextFileViewer("File viewer", self.fileList, self.windows)
  
  
def openTextFile(filename):
    """ Open a text file with an external or default viewer. """
    if envVarOn('SCIPION_EXTERNAL_VIEWER'):
        openTextFileEditor(filename)
    else:
        showTextFileViewer("File viewer", [filename])    
    
    
def openTextFileEditor(filename):
    try:
        _open_cmd(filename)
    except:
        showTextFileViewer("File viewer", [filename])
    
    
def showTextFileViewer(title, filelist, parent=None, main=False):
    w = gui.Window(title, parent, minsize=(600, 400))
    viewer = TextFileViewer(w.root, filelist, maxSize=-1)
    viewer.grid(row=0, column=0, sticky='news')
    gui.configureWeigths(w.root)
    w.show()


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.title("View files")
    l = TextFileViewer(root, fileList=sys.argv[1:])
    l.pack(side=tk.TOP, fill=tk.BOTH)
    gui.centerWindows(root)
    root.deiconify()
    root.mainloop()
