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
        
import Tkinter as tk
import ttk, os, gui
from pyworkflow.utils import *
from widgets import Scrollable, Button


class HyperlinkManager:
    """ Tkinter Text Widget Hyperlink Manager, take from:
    http://effbot.org/zone/tkinter-text-hyperlink.htm """
    def __init__(self, text):
        self.text = text
        self.text.tag_config("hyper", foreground="blue", underline=1)
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
        #self.menu.add_command(label="Open", command=self.openFile)
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
        self.config(state=tk.NORMAL)
        self.delete(0.0, tk.END)

    def getText(self):
        return self.get(0.0, tk.END)
        
    def addText(self, text):
        self.config(state=tk.NORMAL)
        if isinstance(text, list):
            for line in text:
                self.addLine(line)
        else:
            for line in text.splitlines():
                self.addLine(line)
        self.config(state=tk.DISABLED)   
        
    def onClick(self, e=None): 
        self.selection = None
        self.selection_clear()
        self.menu.unpost()
        
    def onRightClick(self, e):
        try:
            self.selection = self.selection_get()
            self.menu.post(e.x_root, e.y_root)    
        except TclError, e:
            pass
    
    def copyToClipboard(self, e=None):
        self.clipboard_clear()
        self.clipboard_append(self.selection)
        
    def openFile(self):
        openFile(self.selection)
        
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

#---------------------------------------------------------------------------
# Colors from Xmipp binding
#--------------------------------------------------------------------------- 
from xmipp import XMIPP_MAGENTA, XMIPP_BLUE, XMIPP_GREEN, XMIPP_RED, XMIPP_YELLOW, XMIPP_CYAN, colorStr

colorMap = {'red': XMIPP_RED, 'blue': XMIPP_BLUE,
                'green': XMIPP_GREEN, 'magenta': XMIPP_MAGENTA,
                'yellow': XMIPP_YELLOW, 'cyan': XMIPP_CYAN}


blueStr = lambda s: colorStr(XMIPP_BLUE, s)
greenStr = lambda s: colorStr(XMIPP_GREEN, s)
greenLowStr = lambda s: colorStr(XMIPP_GREEN, s, 0)
failStr = redStr = lambda s: colorStr(XMIPP_RED, s)
headerStr = magentaStr = lambda s: colorStr(XMIPP_MAGENTA, s)
yellowStr = lambda s: colorStr(XMIPP_YELLOW, s)
cyanStr = warnStr = cyanStr = lambda s: colorStr(XMIPP_CYAN, s)


def findColor(color):
    '''This function will search if there are color characters present
    on string and return the color and positions on string'''
    for k, v in colorMap.iteritems():
        x, y = colorStr(v, "_..._").split("_..._")
        fx = color.find(x)
        fy = color.find(y)
        if fx != -1 and fy != -1:
            color = color.replace(x, '').replace(y, '')
            return (k, fx, fy, color)
    return None

def insertColoredLine(text, line, tag=""):
    """ Check if the color codes are present in a line
    and use the corresponding tags. The colors tags should 
    be already configured on text object"""
    ctuple = findColor(line)
    if ctuple is None:
        line = line[line.rfind("\r")+1:]
        text.insert(tk.END, line, tag)  
    else:
        color,idxInitColor,idxFinishColor,cleanText=ctuple
        if idxInitColor>0:
            text.insert(tk.END, cleanText[:(idxInitColor-1)]+" ")
        text.insert(tk.END, cleanText[idxInitColor:idxFinishColor-1], "tag_" + color)
        text.insert(tk.END, cleanText[idxFinishColor:])
        
def configureColorTags(text):
    """ Function to configure tag_colorX for all supported colors.
    It is applicable to an Text text """
    try:
        for color in colorMap.keys():
            text.tag_config("tag_" + color, foreground=color)
        return True
    except Exception, e:
        print "Colors still not available"
    return False
       
class TaggedText(Text):  
    """
    Implement a Text that will recognized some basic tags
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
        return {'bg': "white", 'bd':0, 'font': 'fontNormal'}
    
    def configureTags(self):
        self.tag_config('normal', justify=tk.LEFT)
#        self.tag_config(HYPER_BOLD, justify=tk.LEFT, font=fontBold)
#        self.tag_config(HYPER_ITALIC, justify=tk.LEFT, font=fontItalic)
        self.tag_config(HYPER_BOLD, justify=tk.LEFT, font='fontBold')
        self.tag_config(HYPER_ITALIC, justify=tk.LEFT, font='fontItalic')
        if self.colors:            
            self.colors = configureColorTags(self) # Color can be unavailable, so disable use of colors    
        
    def openLink(self, link):
        from  webbrowser import open
        open(link, new=0) # Open in the same browner, new tab
        
    def matchHyperText(self, match, tag):
        """ Process when a match a found and store indexes inside string."""
        #print "match found at: ", match.start(), match.end(), "mode: ", mode
        self.insert(tk.END, self.line[self.lastIndex:match.start()])
        g1 = match.group(tag)

        if tag == HYPER_BOLD or tag == HYPER_ITALIC:
            self.insert(tk.END, ' ' + g1, tag)
        elif tag == HYPER_LINK1:
            self.insert(tk.END, g1, self.hm.add(lambda: self.openLink(g1)))
        elif tag == HYPER_LINK2:            
            self.insert(tk.END, match.group('link2_label'), self.hm.add(lambda: self.openLink(g1)))
        self.lastIndex = match.end()
        
        return g1

    def addLine(self, line):
        self.line = line
        self.lastIndex = 0
        parseHyperText(line, self.matchHyperText)
        Text.addLine(self, line[self.lastIndex:])


class OutputText(Text):
    """
    Implement a Text that will show file content
    and handle console metacharacter for colored output
    """
    def __init__(self, master, filename, colors=True, refresh=0, goEnd=True, **opts):
        """ colors flag indicate if try to parse color meta-characters
            refresh is the refresh timedeltha in seconds, 0 means no refresh
        """
        self.filename = filename
        self.colors = colors
        self.refresh = refresh
        self.refreshAlarm = None
        self.lineNo = 0
        Text.__init__(self, master, **opts)
        self.doRefresh(refresh)

    def getDefaults(self):
        return {'bg': "black", 'fg':'white', 'bd':0, 'font': gui.fontNormal, 
                'height':30,  'width':100}
        
    def configureTags(self):
        if self.colors:
            configureColorTags(self)

    def addLine(self, line): 
        self.lineNo += 1
        if self.colors:
            self.insert(tk.END, "%05d:   " % self.lineNo, "tag_cyan")  
            insertColoredLine(self, line)
        else:
            self.insert(tk.END, "%05d:   " % self.lineNo)
            self.insert(tk.END, line)   
        
    def readFile(self, clear=False):
        if clear:
            self.lineNo = 0
            self.clear()
            
        self.config(state=tk.NORMAL)
        #self.clear()
        if exists(self.filename):
            textfile = open(self.filename)
            lineNo = 1
            for line in textfile:
                if lineNo > self.lineNo:
                    self.addLine(line)
                lineNo += 1
            textfile.close()
        else:
            self.addLine("File '%s' doesn't exist" % self.filename)
        self.config(state=tk.DISABLED)
        #if self.isAtEnd():
        self.goEnd();
        #if goEnd:
        #    self.goEnd()        
      
    def doRefresh(self, seconds):
        self.stopRefresh() #Stop pending refreshes
        self.readFile()
        if seconds:
            self.refreshAlarm = self.after(seconds*1000, self.doRefresh, seconds)
    
    def stopRefresh(self):
        if self.refreshAlarm:
            self.after_cancel(self.refreshAlarm)
            self.refreshAlarm = None

  
class TextfileViewer(tk.Frame):
    """ Implementation of a simple textfile viewer """
    
    LabelBgColor = "white"
    
    def __init__(self, master, filelist):
        tk.Frame.__init__(self, master)
        self.searchList = None
        self.lastSearch = None
        self.refreshAlarm = None
        self.filelist = filelist
        self.taList = []
        self.fontDict = {}
        self.createWidgets()
        self.master = master
        self.addBinding()
        
    def addFileTab(self, filename):
        tab = tk.Frame(self.notebook)
        tab.rowconfigure(0, weight=1)
        tab.columnconfigure(0, weight=1)
        t = OutputText(tab, filename, width=100, height=30, bg='black', fg='white')
        t.frame.grid(column=0, row=0, padx=5, pady=5, sticky='nsew')
        self.taList.append(t)
        tabText = "   %s   " % os.path.basename(filename)
        self.notebook.add(tab, text=tabText)        
    
    def createWidgets(self):
        #registerCommonFonts()
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)        
        
        #Create tabs frame
        tabsFrame = tk.Frame(self)
        tabsFrame.grid(column=0, row=1, padx=5, pady=5, sticky="nsew")
        tabsFrame.columnconfigure(0, weight=1)
        tabsFrame.rowconfigure(0, weight=1)
        self.notebook = ttk.Notebook(tabsFrame)  
        self.notebook.rowconfigure(0, weight=1)
        self.notebook.columnconfigure(0, weight=1)      
        for f in self.filelist:
            self.addFileTab(f)
        self.notebook.grid(column=0, row=0, sticky='nsew', padx=5, pady=5)        

    def addBinding(self):
        self.master.bind('<Control_L><Home>', lambda e: self.changePosition(1.0))
        self.master.bind('<Control_L><End>', lambda e: self.changePosition(tk.END))
        self.master.bind('<Alt_L><c>', lambda e: self.master.destroy())
        self.master.bind('<Return>', lambda e: self.findText())
        self.master.bind('<Control_L><n>', lambda e: self.findText())
        self.master.bind('<Control_L><p>', lambda e: self.findText(-1))
        self.master.bind('<Alt_L><plus>', self.changeFont)
        self.master.bind('<Alt_L><minus>', self.changeFont)
    
    def selectedText(self):
        return self.taList[self.notebook.index(self.notebook.select())]
    
    def changeFont(self, event=""):
        for font in self.fontDict.values():
            gui.changeFontSize(font, event)
              
    def refreshAll(self, e=None):
        """ Refresh all output textareas. """
        for ta in self.taList:
            ta.readFile(clear=True)
            ta.goEnd()
        
    def refreshOutput(self, e=None):
        if self.refreshAlarm:
            self.after_cancel(self.refreshAlarm)
            self.refreshAlarm = None
        self.selectedText().readFile()
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
    gui.centerWindows(root, refWindows=parent)
    root.deiconify()
    root.mainloop()

    
if __name__ == '__main__':
    import sys
    root = tk.Tk()
    root.withdraw()
    root.title("View files")
    l = TextfileViewer(root, filelist=sys.argv[1:])
    l.pack(side=tk.TOP, fill=tk.BOTH)
    gui.centerWindows(root)
    gui.setCommonFonts()
    root.deiconify()
    root.mainloop()               
