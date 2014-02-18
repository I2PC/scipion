'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

import os
from os.path import join, exists, basename, dirname, commonprefix, relpath, abspath
import Tkinter as tk

from tkSimpleDialog import Dialog
import ttk
from config_protocols import LabelBgColor, ButtonBgColor, ButtonActiveBgColor, SectionTextColor
from protlib_filesystem import getXmippPath, xmippExists, removeFilenamePrefix, fixPath, splitFilename, getExt
from protlib_utils import runChimera, runVMD
from Tkinter import TclError
from protlib_gui_ext import *

RESOURCES = getXmippPath('resources')

Fonts = {}

def registerFont(name, **opts):
    import tkFont
    global Fonts
    Fonts[name] = tkFont.Font(**opts)

def registerCommonFonts():
    from config_protocols import FontName, FontSize
    if 'normal' not in Fonts.keys():
        registerFont('normal', family=FontName, size=FontSize)
    if 'button' not in Fonts.keys():
        registerFont('button', family=FontName, size=FontSize, weight='bold')
    if 'label' not in Fonts.keys():
        registerFont('label', family=FontName, size=FontSize+1, weight='bold')
        
def configDefaults(opts, defaults):
    for key in defaults.keys():
        if not opts.has_key(key):
            opts[key] = defaults[key]
            
def openLink(link):
    ''' Open a link in default web browser '''
    if os.path.isdir(link):
        showBrowseDialog(link, link, seltype="none", selmode="browse")
    elif xmippExists(link):
        ext = getExt(link)
        if ext in TEXT_EXTENSIONS:
            showTextfileViewer(link, [link])
        elif ext in CHIMERA_EXTENSIONS:
            chimera(link)
        else: # VALIDATE THAT showj can visualize the extesion
            showj(link)
    else:
        from  webbrowser import open
        open(link)
    
def openFile(filename):
    openLink(filename)
    
def changeFontSizeByDeltha(font, deltha, min=-999, max=999):
    size = font['size']
    new_size = size + deltha
    if new_size >= min and new_size <= max:
        font.configure(size=new_size)
            
def changeFontSize(font, event, min=-999, max=999):
    deltha = 2
    if event.char == '-':
        deltha = -2
    changeFontSizeByDeltha(font, deltha, min, max)

# **********************************************
# *     Following are some helper functions    *
# **********************************************

def getGeometry(win):
    ''' Return the geometry information of the windows
    It will be a tuple (width, height, x, y)
    '''
    return win.winfo_reqwidth(), win.winfo_reqheight(), win.winfo_x(), win.winfo_y()

def centerWindows(root, dim=None, refWindows=None):
    """Center a windows in the middle of the screen 
    or in the middle of other windows(refWindows param)"""
    root.update_idletasks()
    gw, gh, gx, gy = getGeometry(root)
    if refWindows:
        rw, rh, rx, ry = getGeometry(refWindows)
        x = rx + (rw - gw) / 2
        y = ry + (rh - gh) / 2 
    else:
        w = root.winfo_screenwidth()
        h = root.winfo_screenheight()
        if not dim is None:
            gw, gh = dim
        x = (w - gw) / 2
        y = (h - gh) / 2
        
    root.geometry("%dx%d+%d+%d" % (gw, gh, x, y))

def ProjectLabel(master, **opts):
    return tk.Label(master, font=Fonts['label'], fg=SectionTextColor, **opts)

# ******************************************************************
# *     Following are implementation of some widgets extensions   *
# ****************************************************************** 
class XmippTree(ttk.Treeview): 
    def __init__(self, master, **opts):
        ttk.Treeview.__init__(self, master, **opts)
        
    def selection_first(self):
        ''' Return first selected item or None if selection empty'''
        selection = self.selection()
        if len(selection):
            return selection[0]
        return None
    
    def _selection_move(self, moveFunc):
        item = self.selection_first()
        if item:
            item = moveFunc(item)
            if item != '':
                self.selection_set(item)
        
    def selection_up(self, e=None):
        ''' change selection to previous item '''
        self._selection_move(self.prev)
    
    def selection_down(self, e=None):
        ''' change selection to to next item '''
        self._selection_move(self.next)
        
    def item_up(self, e=None):
        '''if selected item is not the first move up one position'''
        item = self.selection_first()
        if item:
            index = self.index(item)
            if index > 0:
                self.move(item, '', index-1)
                
    def item_down(self, e=None):
        '''if selected item is not the first move up one position'''
        item = self.selection_first()
        if item:
            index = self.index(item)
            if self.next(item) != '':
                self.move(item, '', index+1)
                
    def clear(self):
        ''' remove all items '''
        childs = self.get_children('')
        for c in childs:
            self.delete(c)
            
        
class AutoScrollbar(tk.Scrollbar):
    '''A scrollbar that hides itself if it's not needed.'''
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        tk.Scrollbar.set(self, lo, hi)
        
class ScrollFrame(tk.Frame):
    ''' An scrollable Frame, that will create a Canvas with
    scrolls and with other Frame to place contents'''
    def __init__(self, master, **opts):
        tk.Frame.__init__(self, master, **opts)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        vscrollbar = AutoScrollbar(self)
        vscrollbar.grid(row=0, column=1, sticky='ns')
        hscrollbar = AutoScrollbar(self, orient=tk.HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky='ew')
        self.canvas = tk.Canvas(self, background='blue',#BgColor,
                        yscrollcommand=vscrollbar.set,
                        xscrollcommand=hscrollbar.set)
        self.canvas.grid(row=0, column=0, sticky='nsew')
        vscrollbar.config(command=self.canvas.yview)
        hscrollbar.config(command=self.canvas.xview)
        self.frame = tk.Frame(self.canvas)#background=BgColor)
        self.frame.rowconfigure(0, weight=1)
        self.frame.columnconfigure(0, weight=1)
        #self.frame.grid(row=0, column=0, sticky='nsew')
        
'''Michael Lange <klappnase (at) freakmail (dot) de>
The ToolTip class provides a flexible tooltip widget for Tkinter; it is based on IDLE's ToolTip
module which unfortunately seems to be broken (at least the version I saw).
INITIALIZATION OPTIONS:
anchor :        where the text should be positioned inside the widget, must be on of "n", "s", "e", "w", "nw" and so on;
                default is "center"
bd :            borderwidth of the widget; default is 1 (NOTE: don't use "borderwidth" here)
bg :            background color to use for the widget; default is "lightyellow" (NOTE: don't use "background")
delay :         time in ms that it takes for the widget to appear on the screen when the mouse pointer has
                entered the parent widget; default is 1500
fg :            foreground (i.e. text) color to use; default is "black" (NOTE: don't use "foreground")
follow_mouse :  if set to 1 the tooltip will follow the mouse pointer instead of being displayed
                outside of the parent widget; this may be useful if you want to use tooltips for
                large widgets like listboxes or canvases; default is 0
font :          font to use for the widget; default is system specific
justify :       how multiple lines of text will be aligned, must be "left", "right" or "center"; default is "left"
padx :          extra space added to the left and right within the widget; default is 4
pady :          extra space above and below the text; default is 2
relief :        one of "flat", "ridge", "groove", "raised", "sunken" or "solid"; default is "solid"
state :         must be "normal" or "disabled"; if set to "disabled" the tooltip will not appear; default is "normal"
text :          the text that is displayed inside the widget
textvariable :  if set to an instance of Tkinter.StringVar() the variable's value will be used as text for the widget
width :         width of the widget; the default is 0, which means that "wraplength" will be used to limit the widgets width
wraplength :    limits the number of characters in each line; default is 150

WIDGET METHODS:
configure(**opts) : change one or more of the widget's options as described above; the changes will take effect the
                    next time the tooltip shows up; NOTE: follow_mouse cannot be changed after widget initialization

Other widget methods that might be useful if you want to subclass ToolTip:
enter() :           callback when the mouse pointer enters the parent widget
leave() :           called when the mouse pointer leaves the parent widget
motion() :          is called when the mouse pointer moves inside the parent widget if follow_mouse is set to 1 and the
                    tooltip has shown up to continually update the coordinates of the tooltip window
coords() :          calculates the screen coordinates of the tooltip window
create_contents() : creates the contents of the tooltip window (by default a Tkinter.Label)
'''
# Ideas gleaned from PySol

class ToolTip:
    def __init__(self, master, text='Your text here', delay=1500, **opts):
        self.master = master
        self._opts = {'anchor':'center', 'bd':1, 'bg':'lightyellow', 'delay':delay, 'fg':'black',\
                      'follow_mouse':0, 'font':None, 'justify':'left', 'padx':4, 'pady':2,\
                      'relief':'solid', 'state':'normal', 'text':text, 'textvariable':None,\
                      'width':0, 'wraplength':150}
        self.configure(**opts)
        self._tipwindow = None
        self._id = None
        self._id1 = self.master.bind("<Enter>", self.enter, '+')
        self._id2 = self.master.bind("<Leave>", self.leave, '+')
        self._id3 = self.master.bind("<ButtonPress>", self.leave, '+')
        self._follow_mouse = 0
        if self._opts['follow_mouse']:
            self._id4 = self.master.bind("<Motion>", self.motion, '+')
            self._follow_mouse = 1
    
    def configure(self, **opts):
        for key in opts:
            if self._opts.has_key(key):
                self._opts[key] = opts[key]
            else:    
                KeyError = 'KeyError: Unknown option: "%s"' %key
                raise KeyError
    
    ##----these methods handle the callbacks on "<Enter>", "<Leave>" and "<Motion>"---------------##
    ##----events on the parent widget; override them if you want to change the widget's behavior--##
    
    def enter(self, event=None):
        self._schedule()
        
    def leave(self, event=None):
        self._unschedule()
        self._hide()
    
    def motion(self, event=None):
        if self._tipwindow and self._follow_mouse:
            x, y = self.coords()
            self._tipwindow.wm_geometry("+%d+%d" % (x, y))
    
    ##------the methods that do the work:---------------------------------------------------------##
    
    def _schedule(self):
        self._unschedule()
        if self._opts['state'] == 'disabled':
            return
        self._id = self.master.after(self._opts['delay'], self._show)

    def _unschedule(self):
        id = self._id
        self._id = None
        if id:
            self.master.after_cancel(id)

    def _show(self):
        if self._opts['state'] == 'disabled':
            self._unschedule()
            return
        if not self._tipwindow:
            self._tipwindow = tw = tk.Toplevel(self.master)
            # hide the window until we know the geometry
            tw.withdraw()
            tw.wm_overrideredirect(1)

            if tw.tk.call("tk", "windowingsystem") == 'aqua':
                tw.tk.call("::tk::unsupported::MacWindowStyle", "style", tw._w, "help", "none")

            self.create_contents()
            tw.update_idletasks()
            x, y = self.coords()
            tw.wm_geometry("+%d+%d" % (x, y))
            tw.deiconify()
    
    def _hide(self):
        tw = self._tipwindow
        self._tipwindow = None
        if tw:
            tw.destroy()
                
    ##----these methods might be overridden in derived classes:----------------------------------##
    
    def coords(self):
        # The tip window must be completely outside the master widget;
        # otherwise when the mouse enters the tip window we get
        # a leave event and it disappears, and then we get an enter
        # event and it reappears, and so on forever :-(
        # or we take care that the mouse pointer is always outside the tipwindow :-)
        tw = self._tipwindow
        twx, twy = tw.winfo_reqwidth(), tw.winfo_reqheight()
        w, h = tw.winfo_screenwidth(), tw.winfo_screenheight()
        # calculate the y coordinate:
        if self._follow_mouse:
            y = tw.winfo_pointery() + 20
            # make sure the tipwindow is never outside the screen:
            if y + twy > h:
                y = y - twy - 30
        else:
            y = self.master.winfo_rooty() + self.master.winfo_height() + 3
            if y + twy > h:
                y = self.master.winfo_rooty() - twy - 3
        # we can use the same x coord in both cases:
        x = tw.winfo_pointerx() - twx / 2
        if x < 0:
            x = 0
        elif x + twx > w:
            x = w - twx
        return x, y

    def create_contents(self):
        opts = self._opts.copy()
        for opt in ('delay', 'follow_mouse', 'state'):
            del opts[opt]
        label = tk.Label(self._tipwindow, **opts)
        label.pack()

class FlashMessage():
    def __init__(self, master, msg, delay=5, relief='solid', func=None):
        self.root = tk.Toplevel(master=master)
        #hides until know geometry
        self.root.withdraw()
        self.root.wm_overrideredirect(1)
        tk.Label(self.root, text="   %s   " % msg,
                 bd=1, bg='DodgerBlue4', fg='white').pack()
        centerWindows(self.root, refWindows=master)
        self.root.deiconify()
        self.root.grab_set()
        self.msg = msg

        if func:
            self.root.update_idletasks()
            self.root.after(10, self.proccess, func)
        else:
            self.root.after(int(delay*1000), self.close)
        self.root.wait_window(self.root)
        
    def proccess(self, func):
        func()
        self.root.destroy()
        
    def close(self):
        self.root.destroy()
        
##---------demo code-----------------------------------##
def getXmippImage(imagePath):
    ''' Return a Tk image from the 'resources' folder.
    None is returned if an error occurs '''
    image = None
    if imagePath:
        try:
            imgPath = join(RESOURCES, imagePath)
            image = tk.PhotoImage(file=imgPath)
        except tk.TclError:
            pass
    return image
    
class XmippButton(tk.Button):
    def __init__(self, master, text, imagePath=None, tooltip=None, **opts):
        defaults = {'activebackground': ButtonActiveBgColor, 'bg':ButtonBgColor}
        defaults.update(opts)
        btnImage = getXmippImage(imagePath)
        
        if btnImage:
            #height=28, width=28,
            if not opts.has_key('bg'):
                del defaults['bg']
            tk.Button.__init__(self, master, image=btnImage, bd=0,  **defaults)
            self.image = btnImage
        else:
            tk.Button.__init__(self, master, text=text, font=Fonts['button'], **defaults)
            
        if tooltip:
            ToolTip(self, tooltip, 500)
            
    def setImage(self, imagePath):
        self.image = getXmippImage(imagePath)
        self.config(image=self.image)
        
'''Implement a Listbox Dialog, it will return
the index selected in the lisbox or -1 on Cancel'''
class ListboxDialog(Dialog):
    def __init__(self, master, itemList, **kargs):
        self.list = itemList
        self.kargs = kargs
        Dialog.__init__(self, master)        
        
    def body(self, master):
        self.result = []
        self.lb = tk.Listbox(master, selectmode=tk.EXTENDED, bg="white")
        self.lb.config(**self.kargs)
        self.lb.pack(fill=tk.BOTH)
        self.lb.bind('<Double-Button-1>', self.ok)
        maxLength=0
        for item in self.list:
            self.lb.insert(tk.END, item)
            maxLength=max(maxLength,len(item))
        self.lb.config(width=maxLength+3)
        if len(self.list) > 0:
            self.lb.selection_set(0)
        return self.lb # initial focus

    def buttonbox(self):
        box = tk.Frame(self)
        w = XmippButton(box, text="OK", width=7, command=self.ok)
        w.pack(side=tk.RIGHT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        box.pack()
        
    def apply(self):
        self.result = map(int, self.lb.curselection())

    
class HyperlinkManager:
    '''
    Tkinter Text Widget Hyperlink Manager
    taken from:
    http://effbot.org/zone/tkinter-text-hyperlink.htm
    '''
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
            
class XmippText(tk.Text):    
    '''
    Base Text widget with some functionalities that will be used
    for other extensions
    ''' 
    def __init__(self, master, **options):  
        registerCommonFonts()    
        defaults = self.getDefaults()
        defaults.update(options)
        self._createWidgets(master, defaults)
        self.configureTags()        

    def _createWidgets(self, master, options):
        '''This is an internal function to create the Text, the Scrollbar and the Frame'''
        frame = tk.Frame(master)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)        
        scrollbar = AutoScrollbar(frame)
        scrollbar.grid(row=0, column=1, sticky='ns')
        options['yscrollcommand'] = scrollbar.set
        options['wrap'] = tk.WORD
        tk.Text.__init__(self, frame, options)
        scrollbar.config(command=self.yview)
        self.grid(row=0, column=0, sticky='nsew')
        self.frame = frame
        self.scrollbar = scrollbar 
        # create a popup menu
        self.menu = tk.Menu(master, tearoff=0, postcommand=self.updateMenu)
        self.menu.add_command(label="Copy to clipboard", command=self.copyToClipboard)
        self.menu.add_command(label="Open", command=self.openFile)
        # Associate with right click
        self.bind("<Button-1>", self.onClick)
        self.bind("<Button-3>", self.onRightClick)
        
    def getDefaults(self):
        '''This should be implemented in subclasses to provide defaults'''
        return {}
    
    def configureTags(self):
        '''This should be implemented to create specific tags'''
        pass
    
    def addLine(self, line):
        '''Should be implemented to add a line '''
        pass
        
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

    def addText(self, text):
        self.config(state=tk.NORMAL)
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
        if not xmippExists(self.selection):
            state = 'disabled'#self.menu.entryconfig(1, background="green")
        self.menu.entryconfig(1, state=state)
        
    def setReadOnly(self, value):
        state = tk.NORMAL
        if value:
            state = tk.DISABLED
        self.config(state=state) 

def configureColorTags(text):
    ''' Function to configure tag_colorX for all supported colors.
    It is applicable to an XmippText text '''
    try:
        from protlib_xmipp import colorMap
        for color in colorMap.keys():
            text.tag_config("tag_" + color, foreground=color)
        return True
    except Exception, e:
        print "Colors still not available"
    return False
        
        
def insertColoredLine(text, line, tag=""):
    ''' Check if the color codes are present in a line
    and use the corresponding tags. The colors tags should 
    be already configured on text object'''
    from protlib_xmipp import findColor
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
        
class TaggedText(XmippText):  
    '''
    Implement a Text that will recognized some basic tags
    <some_text> will display some_text in bold
    [some_link] will display some_link as hiperlinnk
    also colors are recognized if set option colors=True
    '''           
    def __init__(self, master, colors=True, **options):  
        self.colors = colors
        XmippText.__init__(self, master, **options)
        # Create regex for tags parsing
        import re
        self.regex = re.compile('((?:\[[^]]+\])|(?:<[^>]+>))')
        self.hm = HyperlinkManager(self)

    def getDefaults(self):
        return {'bg': "white", 'bd':0, 'font':Fonts['normal']}
    
    def configureTags(self):
        self.tag_config('normal', justify=tk.LEFT)
        self.tag_config('bold', justify=tk.LEFT, font=Fonts['button'])
        if self.colors:            
            self.colors = configureColorTags(self) # Color can be unavailable, so disable use of colors    
        
    def getTaggedParts(self, parts):
        ''' Detect [] as links text and <> as bold text'''
        tagsDict = {'[': 'link', '<': 'bold'}
        tagged_parts = []
        for p in parts:
            if len(p) > 0:
                if p[0] in tagsDict.keys():
                    tagged_parts.append((p[1:-1], tagsDict[p[0]]))
                else:
                    tagged_parts.append((p, 'normal'))
        return tagged_parts
        
    def addLine(self, line):
        parts = self.getTaggedParts(self.regex.split(line))
        for p, t in parts:
            if t == 'link':
                def insertLink(link):
                    self.insert(tk.INSERT, link, self.hm.add(lambda: openLink(link)))
                insertLink(p)
            else:
                if self.colors:
                    insertColoredLine(self, p, t)
                else:
                    self.insert(tk.END, p, t)
        self.addNewline()       

class OutputText(XmippText):
    '''
    Implement a Text that will show file content
    and handle console metacharacter for colored output
    '''
    def __init__(self, master, filename, colors=True, refresh=0, goEnd=True, **opts):
        ''' colors flag indicate if try to parse color meta-characters
            refresh is the refresh timedeltha in seconds, 0 means no refresh
        '''
        self.filename = filename
        self.colors = colors
        self.refresh = refresh
        self.refreshAlarm = None
        self.lineNo = 0
        XmippText.__init__(self, master, **opts)
        self.doRefresh(refresh)

    def getDefaults(self):
        return {'bg': "black", 'fg':'white', 'bd':0, 'font':Fonts['normal'], 
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

  
''' Implementation of a simple textfile viewer '''
class TextfileViewer(tk.Frame):
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
        t = OutputText(tab, filename, width=100, height=30)
        t.frame.grid(column=0, row=0, padx=5, pady=5, sticky='nsew')
        self.taList.append(t)
        tabText = "   %s   " % basename(filename)
        self.notebook.add(tab, text=tabText)        
    
    def createWidgets(self):
        registerCommonFonts()
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)        
        #Create toolbar frame
        toolbarFrame = tk.Frame(self)
        toolbarFrame.grid(column=0, row=0, padx=5, sticky='new')
        toolbarFrame.columnconfigure(0, weight=1)
        toolbarFrame.columnconfigure(1, weight=1)
        #Add the search box
        left = tk.Frame(toolbarFrame)
        left.grid(column=0, row=0, sticky='nw')
        XmippButton(left, "Open", 'folderopen.gif').grid(row=0, column=0, padx=(0, 5))
        XmippButton(left, "Save", 'save.gif').grid(row=0, column=1, padx=(0, 5))
        XmippButton(left, "Refresh", 'refresh.gif', command=self.refreshAll).grid(row=0, column=2, padx=(0, 5))
        right = tk.Frame(toolbarFrame)
        right.grid(column=1, row=0, sticky='ne')        
        self.searchVar = tk.StringVar()
        tk.Label(right, text='Search:').grid(row=0, column=3, padx=5, pady=5, )
        self.searchEntry = tk.Entry(right, textvariable=self.searchVar, bg=LabelBgColor)
        self.searchEntry.grid(row=0, column=4, sticky='ew', padx=5, pady=5)
        XmippButton(right, "Search", 'search.gif').grid(row=0, column=5, padx=(0, 5))
        
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
        
        self.searchEntry.focus_set()

    def addBinding(self):
        self.master.bind('<Control_L><Home>', lambda e: self.changePosition(1.0))
        self.master.bind('<Control_L><End>', lambda e: self.changePosition(tk.END))
        self.master.bind('<Alt_L><c>', lambda e: self.master.destroy())
        self.master.bind('<Control_L><f>', lambda e: self.searchEntry.focus_set())
        self.master.bind('<Return>', lambda e: self.findText())
        self.master.bind('<Control_L><n>', lambda e: self.findText())
        self.master.bind('<Control_L><p>', lambda e: self.findText(-1))
        self.master.bind('<Alt_L><plus>', self.changeFont)
        self.master.bind('<Alt_L><minus>', self.changeFont)
    
    def selectedText(self):
        return self.taList[self.notebook.index(self.notebook.select())]
    
    def changeFont(self, event=""):
        for font in self.fontDict.values():
            changeFontSize(font, event)
              
    def refreshAll(self, e=None):
        ''' Refresh all output textareas '''
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
  
if __name__ == '__main__':
    import sys
    root = tk.Tk()
    root.withdraw()
    root.title("View files")
    l = TextfileViewer(root, filelist=sys.argv[1:])
    l.pack(side=tk.TOP, fill=tk.BOTH)
    centerWindows(root)
    root.deiconify()
    root.mainloop()               
            
'''Implementation of our own dialog to display messages'''
class ShowDialog(Dialog):
    def __init__(self, master, title, msg, type, defaultResult=[]):
        self.msg = msg
        self.type = type
        self.result = defaultResult
        
        Dialog.__init__(self, master, title) 
        
    def body(self, master):
        try:
            from protlib_filesystem import getXmippPath
            self.image = tk.PhotoImage(file=join(getXmippPath('resources'), self.type + '.gif'))
        except tk.TclError:
            self.image = None
        
        self.frame = tk.Frame(master, bg="white", bd=2)
        self.text = TaggedText(self.frame)
        # Insert image
        if self.image:
            self.label = tk.Label(self.frame, image=self.image, bg='white', bd=0)
            self.label.pack(side=tk.LEFT, anchor=tk.N)
        # Insert lines of text
        mylines = self.msg.splitlines()
        m = 0
        for l in mylines:
            m = max(m, len(l))
            self.text.addLine(l)
        m = min(m + 5, 80)
        h = min(len(mylines)+3, 30)
        self.text.config(height=h, width=m)
        self.text.addNewline()
        self.text.config(state=tk.DISABLED)
        self.text.frame.pack()
        self.frame.pack()
        
    def buttonbox(self):
        box = tk.Frame(self)    
        self.btnNo = XmippButton(box, text="Ok", width=7, command=self.cancel)
        self.btnNo.pack(side=tk.LEFT, padx=5, pady=5, anchor='e')
        box.pack()
        self.bind("<Return>", self.cancel)
        self.bind("<Escape>", self.cancel)

'''Implementation of our own version of Yes/No dialog'''
class YesNoDialog(ShowDialog):
    def __init__(self, master, title, msg, DefaultNo=True):
        self.root = master  
        self.defaultNo = DefaultNo
        ShowDialog.__init__(self, master, title, msg, 'warning', False)        
        
        
    def buttonbox(self):
        box = tk.Frame(self)    
        self.btnYes = XmippButton(box, text="Yes", width=7, command=self.ok)
        self.btnYes.pack(side=tk.LEFT, padx=5, pady=5, anchor='e')
        self.btnNo = XmippButton(box, text="No", width=7, command=self.cancel)
        self.btnNo.pack(side=tk.LEFT, padx=5, pady=5, anchor='e')
        box.pack()        
        if self.defaultNo:
            self.btnNo.focus_set()
        else: 
            self.btnYes.focus_set()
        self.bind("<Return>", self.handleReturn)
        self.bind("<Escape>", self.cancel)
        self.bind("<Left>", lambda e: self.btnYes.focus_set())
        self.bind("<Right>", lambda e: self.btnNo.focus_set())

    def ok(self, e=None):
        self.result = True
        self.destroy()
        #self.root.destroy()
        
    def cancel(self, e=None):
        self.result = False
        self.destroy()
        
    def handleReturn(self, event=None):
        w = self.focus_get()
        if w is self.btnYes:
            self.ok()
        else:
            self.cancel()

class TiltPairsDialog(Dialog):
    ''' Dialog to create/edit micrographs tilt pairs
    It will receive two lists and will return the same
    '''
    def __init__(self, pairsList, **kargs):
        self.uList, self.tList = pairsList
        self.kargs = kargs
        self.result = None 
        
    def fillTreeFromList(self, tree, list):
        for item in list:
            self.insertItem(tree, item)
        
    def getListFromTree(self, tree):
        return [tree.item(i, 'text') for i in tree.get_children()]
            
    def up(self, e=None):
        self.lastTree.item_up()
        
    def down(self, e=None):
        self.lastTree.item_down()
        
    def insertItem(self, tree, text):
        tree.insert('', 'end', text=text)
        
    def moveItems(self, fromTree, toTree):
        for item in fromTree.selection():
            self.insertItem(toTree, fromTree.item(item, 'text'))
            fromTree.delete(item)
    
    def right(self, e=None):
        if self.lastTree == self.uTree:
            self.moveItems(self.uTree, self.tTree)

    def left(self, e=None):
        if self.lastTree == self.tTree:
            self.moveItems(self.tTree, self.uTree)
        
    def setLast(self, tree):
        self.lastTree = tree
        
    def remove(self, e=None):
        for item in self.lastTree.selection():
            self.lastTree.delete(item)
    
    def cancel(self, e=None):
        self.result = None
        self.root.destroy()
               
    def ok(self, e=None):
        self.result = (self.getListFromTree(self.uTree), self.getListFromTree(self.tTree))
        self.root.destroy()
           
    def createGUI(self, root=None, parent=None):
        if root:
            self.root = root
        else:
            self.root = tk.Tk()

        # Main setup            
        self.parent = parent   
        self.root.withdraw()
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.minsize(300, 300)
        self.root.title('Choose Tilt Pairs')

        frame = tk.Frame(self.root)
        # Create list for tilted micrographs
        self.tTree = XmippTree(frame)
        self.lastTree = self.tTree
        self.tTree.bind('<Button-1>', lambda e: self.setLast(self.tTree))
        #self.tTree.column('Tilted', width=100, anchor='e')
        self.tTree.heading('#0', text='Tilted')
        self.fillTreeFromList(self.tTree, self.tList)
        self.tTree.grid(row=0, column=2, sticky='nsew')
        # Create list for untilted micrographs
        self.uTree = XmippTree(frame)
        self.uTree.bind('<Button-1>', lambda e: self.setLast(self.uTree))
        #self.uTree.column('Untilted', width=100, anchor='e')
        self.uTree.heading('#0', text='Untilted')
        self.fillTreeFromList(self.uTree, self.uList)
        self.uTree.grid(column=0, row=0, sticky='nsew')
        frame.grid(column=0, row=0, sticky='nsew')
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(2, weight=1)               
        
        btnFrame = tk.Frame(frame)
        XmippButton(btnFrame, text="Up", imagePath='up.gif', command=self.up).grid(column=0, row=0, columnspan=2)#, padx=5, pady=5)
        XmippButton(btnFrame, text="Left", imagePath='left.gif', command=self.left).grid(column=0, row=1)#, padx=5, pady=5)
        XmippButton(btnFrame, text="Right", imagePath='right.gif', command=self.right).grid(column=1, row=1)#, padx=5, pady=5)
        XmippButton(btnFrame, text="Down", imagePath='down.gif', command=self.down).grid(column=0, row=2, columnspan=2)#, padx=5, pady=5)
        XmippButton(btnFrame, text="Remove", imagePath='delete.gif', command=self.remove).grid(column=0, row=3, columnspan=2, pady=5)#, padx=5, pady=5)
        btnFrame.grid(column=1, row=0)
        
        box = tk.Frame(frame)
        XmippButton(box, text="Cancel", width=7, command=self.cancel).pack(side=tk.RIGHT, padx=5, pady=5)
        XmippButton(box, text="OK", width=7, command=self.ok).pack(side=tk.RIGHT, padx=5, pady=5)
        box.grid(column=0, row=1, columnspan=3, sticky='sew')
        
    def showGUI(self):        
        centerWindows(self.root, refWindows=self.parent)
        self.root.deiconify()
        #self.root.mainloop()
    
    
''' Functions to display dialogs '''
def askYesNo(title, msg, parent):
    d = YesNoDialog(parent, title, msg)
    return d.result

def showInfo(title, msg, parent):
    ShowDialog(parent, title, msg, 'info')

def showWarning(title, msg, parent):
    ShowDialog(parent, title, msg, 'warning')
    
def showError(title, msg, parent):
    ShowDialog(parent, title, msg, 'error')

        
'''Implement an Entry that support autocomplete on <Tab>
it will have a list of choices, that will be updated
through an update function everytime the entry get the focus'''
class AutoCompleteEntry(tk.Entry):
    def __init__(self, master, **opts):
        tk.Entry.__init__(self, master, **opts)
        
    def setBuildListFunction(self, listFunc, refreshOnTab=False):
        self.listFunc = listFunc
        self.refreshOnTab = refreshOnTab
        self.bind('<FocusIn>', self.refreshList)
        self.bind('<Tab>', self.handle_tab)   
        
    def refreshList(self, event=None):
        self.choices = self.listFunc()
        
    def handle_tab(self, event):
        if self.refreshOnTab:
            self.choices = self.listFunc()
        hits = []
        for e in self.choices:
            if e.startswith(self.get()):
                hits.append(e)
        prefix = commonprefix(hits)
        if len(prefix):
            self.delete(0, tk.END)
            self.insert(0, prefix)
            self.select_clear()
            self.icursor(tk.END)
        return "break"
    

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
                files.sort()
                dirs.sort()
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

class XmippSlider(ttk.Frame):
    ''' Create a personalized frame that contains label, slider and label value
        it also keeps a variable with the value
    '''
    def __init__(self, master, label, from_=0, to=100, value=50, callback=None, step=0.01):
        self.var = tk.DoubleVar()
        self.var.set(float(value))
        ttk.Frame.__init__(self, master)
        ttk.Label(self, text=label).pack(side=tk.LEFT, padx=2, pady=2, anchor='s')
        self.slider = tk.Scale(self, from_=from_, to=to, variable=self.var, 
                                bigincrement=step, resolution=step, orient=tk.HORIZONTAL)
        if callback:
            self.var.trace('w', callback)
        self.slider.pack(side=tk.LEFT, padx=2)
        
    def getValue(self):
        return self.var.get()       
  
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
        
                              
class XmippBrowserCTF(XmippBrowserPreview):
    ''' This subclass is specific preview some operations
        the extra dict will be used for personalized parameters
        that will not be passed to XmippBrowser constructor
    '''
    def __init__(self, **args):
        XmippBrowserPreview.__init__(self, **args)
        # If metadata was provided, populate file list from there
        self.allowFilter = False
        if hasattr(self, 'md') and self.md: 
            self.insertFiles =  self.insertFilesFromMd
        self.maxFreq = 0.5
        self.unitLabel = '1/px'
        self.unit = getattr(self, 'unit', 'pixel')
        
        if hasattr(self, 'freqs') and self.unit == 'angstrom':
            self.unitLabel = '1/A'
#            self.freqs = [float(f) * self.sampling for f in self.freqs]
            self.maxFreq /= self.sampling

    def insertFilesFromMd(self, path):
        from xmipp import MDL_MICROGRAPH
        files = [self.md.getValue(MDL_MICROGRAPH, objId) for objId in self.md]
        prefix = dirname(commonprefix(files))
        self.commonRoot = prefix
        prefix = join(prefix, '')
        for fn in files:
            self.insertElement('', fn.replace(prefix, ''))
            
    def createResultPreview(self):
        self.lf, self.hf = 0, 0
        if self.freqs:
            self.lf = self.freqs[0]
            self.hf = self.freqs[1]

        from protlib_gui_figure import PsdPreview
        return PsdPreview(self.detailstop, self.dim, self.lf, self.hf, dpi=64, col=1)
    
    def getComputeFunction(self):
        from xmipp import fastEstimateEnhancedPSD
        downsampling = float(self.downsamplingVar.get())
        f = lambda :fastEstimateEnhancedPSD(self.image, self.lastitem, downsampling, self.dim, 2)
        return f
    
    def getResults(self):
        results = [self.downsamplingVar.get()]
        if self.freqs:
            results += [self.lfSlider.getValue(), self.hfSlider.getValue()]
        return results        
        
    def createPreviewButton(self, frame):
        self.btnTest = XmippButton(frame, "Preview", command=self.fillResultPreview)
        self.btnTest.pack(side=tk.LEFT, padx=2)
        
    def createLabeledEntry(self, frame, label, key, doUpdate=True):
        var = tk.StringVar()
        var.set(getattr(self, key))
        entry = ttk.Entry(frame, width=10, textvariable=var)
        entry.pack(side=tk.LEFT,padx=2)
        if doUpdate:
            entry.bind('<Return>', self.fillResultPreview)
        setattr(self, key+'Var', var) # register var
        
    def createDetailsBottom(self, parent):
        self.text = None
        frame = ttk.Frame(parent)
        ttk.Label(frame, text="Downsampling").pack(side=tk.LEFT,padx=2)
        self.createLabeledEntry(frame, 'Downsampling', 'downsampling')
        self.createPreviewButton(frame)
        if self.freqs:
            self.addFrequenciesBox(frame)
            self.clearCallbacks.append(lambda: self.updateSliderState(tk.DISABLED))
            self.fillCallbacks.append(lambda: self.updateSliderState(tk.NORMAL))
        return frame
    
    def addFrequenciesBox(self, parent):
        self.freqFrame = ttk.LabelFrame(parent, text="Frequencies", padding="5 5 5 5")
        self.freqFrame.pack(side=tk.LEFT)
        self.lfSlider = self.addFreqSlider('Low freq', self.lf)
        self.hfSlider = self.addFreqSlider('High freq', self.hf)
        
    def addFreqSlider(self, label, value):
            #def __init__(self, master, label, from_=0, to=100, value=50, callback=None, step=0.01):

        slider = XmippSlider(self.freqFrame, 
                             '%s (%s)' % (label, self.unitLabel), 
                             from_=0, 
                             to=self.maxFreq, 
                             value=value, 
                             step=self.maxFreq/100.,
                             callback=lambda a, 
                             b, 
                             c:self.updateFreqRing())
        slider.pack(side=tk.LEFT, padx=3, pady=3)
        return slider
        
    def updateFreqRing(self):
        self.resultPreview.updateFreq(self.lfSlider.getValue(), self.hfSlider.getValue())
        
    def updateSliderState(self, state=tk.NORMAL):
        self.lfSlider.slider.config(state=state)
        self.hfSlider.slider.config(state=state)
        
        
class XmippBrowserBandpassFilter(XmippBrowserCTF):
    ''' This subclass is specific preview some operations
        the extra dict will be used for personalized parameters
        that will not be passed to XmippBrowser constructor
    '''
    def __init__(self, **args):
        XmippBrowserCTF.__init__(self, **args)
        # This is need because of bad inheritance of Non-Frequency filters
        if hasattr(self, 'freqs'):
            self.lf, self.hf, self.decay = self.freqs
        for varName in ['showLowFreq', 'showHighFreq', 'showDecay']:
            setattr(self, varName, getattr(self, varName, True))
            
    def insertFiles(self, path):
        self.commonRoot = dirname(self.fileList[0])
        for f in self.fileList:
            if '@' in f:
                self.insertElement('', f)
                self.commonRoot = ""
            else:
                self.insertElement('', basename(f))
    
    def createResultPreview(self):
        from protlib_xmipp import getFirstImage
        from protlib_gui_figure import ImagePreview
        #mdFn = join(self.dir, self.pattern)
        self.lastitem = getFirstImage(self.fileList[0])
        self.updatePreview(self.lastitem)
        return ImagePreview(self.detailstop, self.dim, dpi=64, label='Filtered', col=1)
    
    def getComputeFunction(self):
        from xmipp import bandPassFilter
        lf, hf, decay = self.getResults()
        return lambda: bandPassFilter(self.image, self.lastitem, lf*self.sampling, hf*self.sampling, decay, self.dim)

    def getResults(self):
        return [self.lfSlider.getValue(), self.hfSlider.getValue(), self.decaySlider.getValue()]
        
    def createDetailsBottom(self, parent):
        frame = ttk.Frame(parent)
        self.addFrequenciesBox(frame)
        self.updateFreqRing = self.fillResultPreview
        self.fillResultPreview()
        return frame
    
    def addFrequenciesBox(self, parent):
        XmippBrowserCTF.addFrequenciesBox(self, parent)
        self.decaySlider = self.addFreqSlider('Decay', self.decay)
        if not self.showDecay:
            self.decaySlider.pack_forget()
        if not self.showLowFreq:
            self.lfSlider.pack_forget()
        if not self.showHighFreq:
            self.hfSlider.pack_forget()
        
    def onClick(self, e):
        XmippBrowserCTF.onClick(self, e)
        self.updateFreqRing()
        

class XmippBrowserGaussianFilter(XmippBrowserBandpassFilter):
    ''' This subclass is specific preview some operations
        the extra dict will be used for personalized parameters
        that will not be passed to XmippBrowser constructor
    '''
    def __init__(self, **args):
        XmippBrowserBandpassFilter.__init__(self, **args)
        self.label = 'Frequency Sigma'
        self.key = 'freqSigma'
    
    def getComputeFunction(self):
        from xmipp import gaussianFilter
        freqSigma = float(self.getResults())
        return lambda: gaussianFilter(self.image, self.lastitem, freqSigma, self.dim)

    def getResults(self):
        return getattr(self, self.key+'Var').get()
        
    def createDetailsBottom(self, parent):
        frame = ttk.Frame(parent)
        self.createLabeledEntry(frame, self.label, self.key)
        self.createPreviewButton(frame)
        self.fillResultPreview()
        return frame

class XmippBrowserRealGaussianFilter(XmippBrowserGaussianFilter):
    def __init__(self, **args):
        XmippBrowserGaussianFilter.__init__(self, **args)
        self.label = 'Real Space Sigma'
        self.key = 'realSigma'
    
    def getComputeFunction(self):
        from xmipp import realGaussianFilter
        realSigma = float(self.getResults())
        return lambda: realGaussianFilter(self.image, self.lastitem, realSigma, self.dim)

"""
class XmippBrowserCropSizeFilter(XmippBrowserCropSizeFilter):
    ''' This subclass is specific preview some operations
        the extra dict will be used for personalized parameters
        that will not be passed to XmippBrowser constructor
    '''
    def __init__(self, **args):
        XmippBrowserPreview.__init__(self, **args)
        self.label = 'CropSize'
        self.key = 'cropSize'
    
    def getComputeFunction(self):
        from xmipp import cropsize
        cropSize = float(self.getResults())
        return lambda: cropsize(self.image, self.lastitem, cropSize, self.dim)
"""

class XmippBrowserBadpixelFilter(XmippBrowserGaussianFilter):
    ''' This subclass is specific preview some operations
        the extra dict will be used for personalized parameters
        that will not be passed to XmippBrowser constructor
    '''
    def __init__(self, **args):
        XmippBrowserGaussianFilter.__init__(self, **args)
        self.label = 'DustRemovalThreshold'
        self.key = 'dustRemovalThreshold'
    
    def getComputeFunction(self):
        from xmipp import badPixelFilter
        dustRemovalThreshold = float(self.getResults())
        return lambda: badPixelFilter(self.image, self.lastitem, dustRemovalThreshold, self.dim)

class XmippBrowserMask(XmippBrowser):
    ''' Preview the middle slice and select a circular mask '''
    def __init__(self, **args):
        initXmippBrowser(self, **args)
        
    def createDetailsTop(self, parent):
        from xmipp import Image, HEADER
        from protlib_xmipp import getImageData
        from protlib_gui_figure import MaskPreview
        
        XmippBrowser.createDetailsTop(self, parent)
        # Read real dimension
        self.dim = 196
        img = Image()
        self.imgFn = self.fileList[0]
        img.read(self.imgFn, HEADER)
        self.real_dim = float(img.getDimensions()[0])
        self.rate = self.dim / self.real_dim
        self.image = Image()
        self.image.readPreview(self.imgFn, self.dim)
        #self.root.minsize(600, 400)
        self.outerRadius = getattr(self, 'outerRadius', int(self.real_dim / 2))
        self.innerRadius = getattr(self, 'innerRadius', 0)
        self.showInner = getattr(self, 'showInner', False)
        self.preview = MaskPreview(self.detailstop, self.dim, label="Central slice", 
                                   outerRadius=self.outerRadius * self.rate, innerRadius=self.innerRadius * self.rate)
        self.detailstop.columnconfigure(1, weight=1)
        ## Read volume preview and update
        Z = getImageData(self.image)
        self.preview.updateData(Z)
        return self.detailstop
    
    def createDetailsBottom(self, parent):
        self.text = None
        frame = ttk.Frame(parent)
        xdim = self.real_dim 
        step = 1
        if self.unit == 'angstrom':
            xdim *= self.sampling
            step *= self.sampling
            self.innerRadius *= self.sampling
            self.outerRadius *= self.sampling
            self.rate /= self.sampling
            
        self.innerRadiusSlider = XmippSlider(frame, "Inner radius", 
                                             from_=0, 
                                             to=xdim/2, 
                                             value=self.innerRadius, step=step,
                                             callback=lambda a, b, c:self.updateMaskRadius())
        self.outerRadiusSlider = XmippSlider(frame, "Outer radius", 
                                             from_=1, 
                                             to=xdim/2, 
                                             value=self.outerRadius, step=step,
                                             callback=lambda a, b, c:self.updateMaskRadius())
#        if self.unit == 'angstrom':
#            self.innerRadiusSlider /= self.sampling
#            self.outerRadiusSlider /= self.sampling

        if self.showInner:
            self.innerRadiusSlider.grid(row=0, column=0)#, padx=3, pady=3)
        self.outerRadiusSlider.grid(row=1, column=0)#, padx=3, pady=3)
        return frame
    
    def insertFiles(self, path):
        self.commonRoot = dirname(self.fileList[0])
        for f in self.fileList:
            if '@' in f:
                self.insertElement('', f)
                self.commonRoot = ""
            else:
                self.insertElement('', basename(f))
        
    def updateMaskRadius(self):
        innerRadius = self.innerRadiusSlider.getValue() * self.rate
        outerRadius = self.outerRadiusSlider.getValue() * self.rate
        self.preview.updateMask(outerRadius, innerRadius)

    def selectFiles(self, e=None):
        self.selectedFiles = (int(self.innerRadiusSlider.getValue()), int(self.outerRadiusSlider.getValue()))
        self.root.destroy()  

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

