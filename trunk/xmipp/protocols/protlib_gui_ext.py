import Tkinter as tk

'''
Taken from following forum:
http://code.activestate.com/recipes/52266-multilistbox-tkinter-widget/
This is a compound widget that gangs multiple Tk Listboxes to a single scrollbar to achieve a simple 
multi-column scrolled listbox. Most of the Listbox API is mirrored to make it act like the normal 
Listbox but with multiple values per row.
'''
class MultiListbox(tk.PanedWindow):
    """MultiListbox class for creating a Grid widget"""
    def __init__(self,master,lists):
        tk.PanedWindow.__init__(self,master,borderwidth=1,showhandle=False,sashwidth=2,sashpad=0,relief=tk.SUNKEN)
        self.lists = []
        self.columns=[]
        for l, w in lists:
            self.columns.append(l)
            frame = tk.Frame(self); frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
            tl=tk.Label(frame, text=l, borderwidth=2, relief=tk.GROOVE)
            tl.pack(fill=tk.X)
            tl.bind('<Button-1>',self.clickon)
            lb = tk.Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,relief=tk.FLAT, exportselection=tk.FALSE, selectmode=tk.MULTIPLE ,bg='white')
            lb.pack(expand=tk.YES, fill=tk.BOTH)
            self.lists.append(lb)
            lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
            lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
            lb.bind('<Double-Button-1>', self._doubleClick)
            lb.bind('<Leave>', lambda e: 'break')
            lb.bind('<B2-Motion>', lambda e, s=self: s._b2motion(e.x, e.y))
            lb.bind('<Button-2>', lambda e, s=self: s._button2(e.x, e.y))
            lb.bind('<Button-3>', lambda e, s=self: s._button3(e.x, e.y))
            lb.bind('<Button-4>', lambda e, s=self: s._scroll(tk.SCROLL, 1, tk.PAGES))
            lb.bind('<Button-5>', lambda e, s=self: s._scroll(tk.SCROLL, -1, tk.PAGES))
            self.add(frame)
        sb = tk.Scrollbar(master, orient=tk.VERTICAL, command=self._scroll,borderwidth=1)
        sb.pack(fill=tk.Y,side=tk.RIGHT,expand=tk.NO)
        for l in self.lists:
            l['yscrollcommand']=sb.set
        self.add(frame)
        self.pack(expand=tk.YES,fill=tk.BOTH, side=tk.TOP)
        self.sortedBy = -1
        self.previousWheel = 0
        self.SelectCallback = None
        self.DoubleClickCallback = None
        self.AllowSort = True

    def _doubleClick(self, event=''):
        if self.DoubleClickCallback:
            self.DoubleClickCallback()

    def _select(self, y,state=16):
        row = self.lists[0].nearest(y)
        if state==16:self.selection_clear(0, tk.END)
        self.selection_set(row)
        return 'break'

    def _button2(self, x, y):
        for l in self.lists: l.scan_mark(x, y)
        return 'break'

    def _button3(self, x, y):
        print "right click"
        #self._select(y, state)
        return 'break'
    

    def _b2motion(self, x, y):
        for l in self.lists: l.scan_dragto(x, y)        
        return 'break'


    def _scroll(self, *args):
        for l in self.lists:
            apply(l.yview, args)
        return 'break'


    def clickon(self,e):
        if self.AllowSort:
            self._sortBy(self.columns.index(e.widget['text']))


    def _sortBy(self, column):
        """ Sort by a given column. """
        if column == self.sortedBy:
            direction = -1 * self.direction
        else:
            direction = 1

        elements = self.get(0, tk.END)
        self.delete(0, tk.END)
        elements.sort(lambda x, y: self._sortAssist(column, direction, x, y))
        self.insert(tk.END, *elements)

        self.sortedBy = column
        self.direction = direction


    def _sortAssist(self, column, direction, x, y):
        c = cmp(x[column], y[column])
        if c:
            return direction * c
        else:
            return direction * cmp(x, y)

    def curselection(self):
        return self.lists[0].curselection()


    def delete(self, first, last=None):
        for l in self.lists:
            l.delete(first, last)


    def get(self, first, last=None):
        result = []
        for l in self.lists:
            result.append(l.get(first,last))
        if last: return apply(map, [None] + result)
        return result


    def index(self, index):
        self.lists[0].index(index)


    def insert(self, index, *elements):
        for e in elements:
            i = 0
            for l in self.lists:
                l.insert(index, e[i])
                i = i + 1
    
    def changetext(self, index, *elements):
        for e in elements:
            for l in self.lists:
                l.itemconfig(index, text=e)        

    def size(self):
        return self.lists[0].size()

    def see(self, index):
        for l in self.lists:
            l.see(index)

    def selection_anchor(self, index):
        for l in self.lists:
            l.selection_anchor(index)

    def selection_clear(self, first, last=None):
        for l in self.lists:
            l.selection_clear(first, last)


    def selection_includes(self, index):
        return self.lists[0].selection_includes(index)


    def selectedIndex(self):
        if len(self.curselection()) > 0:
            return int(self.curselection()[0])
        return -1
    
    def selection_set(self, first, last=None):
        for l in self.lists:
            l.selection_set(first, last)
        #print self.curselection()
        index = self.selectedIndex()
        if self.SelectCallback:
            self.SelectCallback(index)
            
    def selection_move_up(self):
        index = self.selectedIndex() 
        if index > 0:
            self.selection_clear(0, tk.END)
            #self.selection_clear(0, END)
            self.selection_set(index - 1)
            
    def selection_move_down(self):
        index = self.selectedIndex()
        if index != -1 and index < self.size() - 1:
            self.selection_clear(0, tk.END)
            #self.selection_clear(0, END)
            self.selection_set(index + 1)
            

def getGeometry(win):
    return win.winfo_reqwidth(), win.winfo_reqheight(), win.winfo_x(), win.winfo_y()

def centerWindows(root, dim=None, refWindows=None):
    root.update_idletasks()
    """Center a windows in the middle of the screen"""
    gw, gh, gx, gy = getGeometry(root)
    if refWindows:
        rw, rh, rx, ry = getGeometry(refWindows)
        x = rx + rw / 2 - gw / 2
        y = ry + rh / 2 - gh / 2 
    else:
        w = root.winfo_screenwidth()
        h = root.winfo_screenheight()
        if not dim is None:
            gw, gh = dim
        x = (w - gw) / 2
        y = (h - gh) / 2
    root.geometry("%dx%d+%d+%d" % (gw, gh, x, y))

__tipwindow = None

# Creates a tooptip box for a widget.
def createToolTip( widget, text ):
    def enter( event ):
        global __tipwindow
        if __tipwindow or not text:
            return
        x, y, cx, cy = widget.bbox( "insert" )
        x += widget.winfo_rootx() + 27
        y += widget.winfo_rooty() + 27
        # Creates a toplevel window
        __tipwindow = tw = tk.Toplevel( widget )
        # Leaves only the label and removes the app window
        tw.wm_overrideredirect( 1 )
        tw.wm_geometry( "+%d+%d" % ( x, y ) )
        label = tk.Label( tw, text = text, justify = tk.LEFT,
                       background = "#ffffe0", relief = tk.SOLID, borderwidth = 1,
                       font = ( "tahoma", "8", "normal" ) )
        label.pack( ipadx = 1 )
        
    def close( event ):
        global __tipwindow
        tw = __tipwindow
        __tipwindow = None
        if tw:
            tw.destroy()
            
    widget.bind( "<Enter>", enter )
    widget.bind( "<Leave>", close )


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

##---------demo code-----------------------------------##

from tkSimpleDialog import Dialog
ButtonBgColor = "LightBlue"
ButtonActiveBgColor = "LightSkyBlue"
ButtonSelectedColor = "DeepSkyBlue2"

Fonts = {}

def registerFont(name, **opts):
    global Fonts
    Fonts[name] = tkFont.Font(**opts)

def registerCommonFonts():
    if 'normal' not in Fonts.keys():
        registerFont('normal', family=FontName, size=FontSize)
    if 'button' not in Fonts.keys():
        registerFont('button', family=FontName, size=FontSize, weight=tkFont.BOLD)
    if 'label' not in Fonts.keys():
        registerFont('label', family=FontName, size=FontSize+1, weight=tkFont.BOLD)
        
def configDefaults(opts, defaults):
    for key in defaults.keys():
        if not opts.has_key(key):
            opts[key] = defaults[key]
            
def MyButton(master, text, imagePath=None, **opts):
    configDefaults(opts, {'activebackground': ButtonActiveBgColor})
    btnImage = None
    if imagePath:
        try:
            from protlib_filesystem import getXmippPath
            imgPath = os.path.join(getXmippPath('resources'), imagePath)
            btnImage = tk.PhotoImage(file=imgPath)
        except tk.TclError:
            pass
    
    if btnImage:
        btn = tk.Button(master, image=btnImage, bd=0, height=28, width=28, **opts)
        btn.image = btnImage
    else:
        btn = tk.Button(master, text=text, font=Fonts['button'], bg=ButtonBgColor, **opts)
    return btn

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
        w = MyButton(box, text="OK", width=7, command=self.ok)
        w.pack(side=tk.RIGHT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        box.pack()
        
    def apply(self):
        self.result = map(int, self.lb.curselection())
        
def openLink(link):
    ''' Open a link in default web browser '''
    from  webbrowser import open
    open(link)      

'''Implementation of our own version of Yes/No dialog'''
class ShowDialog(Dialog):
    def __init__(self, master, title, msg, type):
        self.msg = msg
        self.type = type
        Dialog.__init__(self, master, title)        
        
    def getTaggedParts(self, parts):
        tagsDict = {'[': 'link', '<': 'bold'}
        tagged_parts = []
        for p in parts:
            if len(p) > 0:
                if p[0] in tagsDict.keys():
                    tagged_parts.append((p[1:-1], tagsDict[p[0]]))
                else:
                    tagged_parts.append((p, 'normal'))
        return tagged_parts
        
    def body(self, master):
        try:
            from protlib_filesystem import getXmippPath
            self.image = tk.PhotoImage(file=os.path.join(getXmippPath('resources'), self.type + '.gif'))
        except tk.TclError:
            self.image = None
        self.result = []
        self.frame = tk.Frame(master, bg="white", bd=2)
        self.text = tk.Text(self.frame, bg="white", bd=0, font=Fonts['normal'])
        self.text.tag_config('normal', justify=tk.LEFT)
        self.text.tag_config('bold', justify=tk.LEFT, font=Fonts['button'])
        
        # Insert image
        if self.image:
            self.label = tk.Label(self.frame, image=self.image, bg='white', bd=0)
            self.label.pack(side=tk.LEFT, anchor=tk.N)

        # Insert lines of text
        # Find some tags to pretty text presentation
        import re
        r = re.compile('((?:\[[^]]+\])|(?:<[^>]+>))')
        mylines = self.msg.splitlines()
        m = 0
        hm = HyperlinkManager(self.text)
                
        for l in mylines:
            m = max(m, len(l))
            parts = self.getTaggedParts(r.split(l))
            for p, t in parts:
                if t == 'link':
                    def insertLink(link):
                        self.text.insert(tk.INSERT, link, hm.add(lambda: openLink(link)))
                    insertLink(p)
                else:
                    self.text.insert(tk.END, p, t)
            self.text.insert(tk.END, '\n')
        m = min(m, 80)
        self.text.config(height=len(mylines)+3, width=m-7)
        self.text.insert(tk.END, '\n')
        self.text.pack()
        self.frame.pack()
        
    def buttonbox(self):
        box = tk.Frame(self)    
        self.btnNo = MyButton(box, text="Ok", width=7, command=self.cancel)
        self.btnNo.pack(side=tk.LEFT, padx=5, pady=5, anchor='e')
        box.pack()
        self.bind("<Return>", self.cancel)
        self.bind("<Escape>", self.cancel)

'''Implementation of our own version of Yes/No dialog'''
class YesNoDialog(ShowDialog):
    def __init__(self, master, title, msg, DefaultNo=True):
        self.defaultNo = DefaultNo
        self.result = False
        ShowDialog.__init__(self, master, title, msg, 'warning')        
        
    def buttonbox(self):
        box = tk.Frame(self)    
        self.btnYes = MyButton(box, text="Yes", width=7, command=self.ok)
        self.btnYes.pack(side=tk.LEFT, padx=5, pady=5, anchor='e')
        self.btnNo = MyButton(box, text="No", width=7, command=self.cancel)
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

    def apply(self):
        self.result = True
        
    def handleReturn(self, event=None):
        w = self.focus_get()
        if w is self.btnYes:
            self.ok()
        else:
            self.cancel()

''' Functions to display dialogs '''
def askYesNo(title, msg, parent):
    d = YesNoDialog(parent, title, msg)
    return d.result

def showInfo(title, msg, parent):
    ShowDialog(parent, title, msg, 'warning')

def showWarning(title, msg, parent):
    ShowDialog(parent, title, msg, 'warning')
    
def showError(title, msg, parent):
    ShowDialog(parent, title, msg, 'error')
    

import os

'''Implement a Text that will show file content'''
class FilePollTextArea(tk.Frame):
    def __init__(self, master, filename, height=30, width=100, colorOn=True):
        tk.Frame.__init__(self, master)
        self.filename = filename
        #Allow not to use color, avoiding dependency with xmipp lib
        self.colorOn = colorOn
        self.createWidgets(height, width)
        self.fillTextArea()
        self.refreshAlarm = None

    def createWidgets(self, h, w):        
        #define a new frame and put a text area in it
        textfr = tk.Frame(self)
        self.text = tk.Text(textfr, height=h, width=w, background='black', fg="white")
        # put a scroll bar in the frame
        yscroll = tk.Scrollbar(textfr)
        self.text.configure(yscrollcommand=yscroll.set)
        yscroll.config(command=self.text.yview)
        #pack everything
        self.text.pack(side=tk.LEFT)
        yscroll.pack(side=tk.RIGHT,fill=tk.Y)
        textfr.pack(side=tk.TOP)
        if self.colorOn:
            from protlib_xmipp import colorMap
            for color in colorMap.keys():
                self.text.tag_config("tag_" + color, foreground=color)

    def fillTextArea(self, goEnd=False):
        self.text.config(state=tk.NORMAL)
        self.text.delete(1.0, tk.END)
        if os.path.exists(self.filename):
            textfile = open(self.filename)
            lineNo = 1
            for line in textfile:
                if self.colorOn:
                    from protlib_xmipp import findColor
                    ctuple = findColor(line)
                    self.text.insert(tk.END, "%05d:   " % lineNo,"tag_cyan")  
                    if ctuple is None:
                        self.text.insert(tk.END, line[line.rfind("\r")+1:])  
                    else:
                        color,idxInitColor,idxFinishColor,cleanText=ctuple
                        if idxInitColor>0:
                            self.text.insert(tk.END, cleanText[:(idxInitColor-1)]+" ")
                        self.text.insert(tk.END, cleanText[idxInitColor:idxFinishColor-1], "tag_" + color)
                        self.text.insert(tk.END, cleanText[idxFinishColor:])
                else:
                    self.text.insert(tk.END, "%05d:   " % lineNo)
                    self.text.insert(tk.END, line)   
                lineNo += 1
            textfile.close()
        else:
            self.text.insert(tk.END, "File '%s' doesn't exist" % self.filename) 
        self.text.config(state=tk.DISABLED)
        if goEnd:
            self.goEnd()
          
    def goEnd(self):
        self.text.see(tk.END)
          
    def doRefresh(self, seconds):
        self.stopRefresh() #Stop pending refreshes
        self.fillTextArea(goEnd=True)
        self.refreshAlarm = self.after(seconds*1000, self.doRefresh, seconds)
    
    def stopRefresh(self):
        if self.refreshAlarm:
            self.after_cancel(self.refreshAlarm)
            self.refreshAlarm = None
        
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
    
from config_protocols import LabelBgColor, FontName, FontSize
import tkFont

class OutputTextArea(tk.Frame):
    def __init__(self, master, fileprefix):
        tk.Frame.__init__(self, master)
        self.exts = ['.log', '.out', '.err']
        self.files = [fileprefix + ext for ext in self.exts]
        self.searchList = None
        self.lastSearch = None
        self.refreshAlarm = None
        self.fontDict = {}
        self.createWidgets()
        self.master = master
        self.addBinding()
        
    def createWidgets(self):
        normal = tkFont.Font(family=FontName, size=FontSize)
        bold = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        self.indexVar = tk.IntVar()
        self.indexVar.set(0)
        self.taList = []
        self.lastIndex = 0
        options = ['Log', 'Standard output', 'Standard error']
        #Label to show the output type
        frame = tk.Frame(self, bd=2, relief=tk.RIDGE)
        frame.grid(row=0, column=0, padx=5, pady=5, sticky='ew')
        tk.Label(frame, text='Select output file:', 
                 font=bold).grid(row=0, column=0, padx=5, pady=5)
        #Add the search box
        self.searchVar = tk.StringVar()
        tk.Label(frame, text='Search:', 
                 font=bold).grid(row=0, column=3, padx=5, pady=5)
        self.searchEntry = tk.Entry(frame, textvariable=self.searchVar, bg=LabelBgColor)
        self.searchEntry.grid(row=0, column=4, sticky='ew', padx=5, pady=5)
        #Add radiobuttons and textareas
        for i in range(3):
            tk.Radiobutton(frame, text='%s (%s)' % (options[i], self.exts[i]), 
                           variable=self.indexVar, command=lambda:self.changeSelection(self.indexVar.get()),
                           value=i, font=normal).grid(row=1, column=i, padx=5, pady=(5,0))
            ta = FilePollTextArea(self, self.files[i])
            ta.text.config(font=normal)
            ta.grid(row=1, column=0, columnspan=5, padx=5, pady=5)
            if i != self.indexVar.get():
                ta.grid_remove()
            self.taList.append(ta)
        self.searchEntry.focus_set()
        self.fontDict['normal'] = normal
        self.fontDict['bold'] = bold

    def addBinding(self):
        self.master.bind('<Control_L><Home>', lambda e: self.changePosition(1.0))
        self.master.bind('<Control_L><End>', lambda e: self.changePosition(tk.END))
        self.master.bind('<Alt_L><c>', lambda e: self.master.destroy())
        self.master.bind('<Alt_L>1' , lambda e: self.changeSelection(0))
        self.master.bind('<Alt_L>2' , lambda e: self.changeSelection(1))
        self.master.bind('<Alt_L>3' , lambda e: self.changeSelection(2))
        self.master.bind('<Control_L><f>', lambda e: self.searchEntry.focus_set())
        self.master.bind('<Return>', lambda e: self.findText())
        self.master.bind('<Control_L><n>', lambda e: self.findText())
        self.master.bind('<Control_L><p>', lambda e: self.findText(-1))
        self.master.bind('<Alt_L><plus>', self.changeFont)
        self.master.bind('<Alt_L><minus>', self.changeFont)
    
    def changeFont(self, event=""):
        for font in self.fontDict.values():
            changeFontSize(font, event)
              
    def refreshOutput(self):
        if self.refreshAlarm:
            self.after_cancel(self.refreshAlarm)
            self.refreshAlarm = None
        self.taList[self.lastIndex].fillTextArea()
        #self.refreshAlarm = self.after(2000, self.refreshOutput)
        
    def changeSelection(self, newValue):
        if self.lastIndex != newValue:
            self.indexVar.set(newValue)
            self.taList[self.lastIndex].grid_remove()
            self.lastIndex = newValue
            self.refreshOutput()
            self.taList[self.lastIndex].grid()
            # Reset search stuff
            self.searchList = None
            self.lastSearch = None
            
    def changePosition(self, index):
        self.taList[self.lastIndex].text.see(index)
        
    def findText(self, dir=1):
        text = self.taList[self.lastIndex].text
        str = self.searchVar.get()
        if str is None or str != self.lastSearch:
            self.buildSearchList(text, str)
            self.lastSearch = str
        else:
            self.nextSearchIndex(text, dir)
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
    
    def nextSearchIndex(self, text, dir=1):
        #use dir=-1 to go backward
        text.tag_remove('found_current', '1.0', tk.END)
        if len(self.searchList)==0:
            return
        self.currentIndex = (self.currentIndex + dir) % len(self.searchList)
        idx, lastidx = self.searchList[self.currentIndex]
        text.tag_config('found_current', foreground='yellow', background='red')
        text.tag_add('found_current', idx, lastidx)
        text.see(idx)
         
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
        prefix = os.path.commonprefix(hits)
        if len(prefix):
            self.delete(0, tk.END)
            self.insert(0, prefix)
            #self.position = self.index(Tkinter.END)
            #self.update_idletasks()
            #self.select_range(Tkinter.END,Tkinter.END)
            self.select_clear()
            self.icursor(tk.END)
        return "break"
    
# Tkinter Text Widget Hyperlink Manager
# taken from:
# http://effbot.org/zone/tkinter-text-hyperlink.htm
class HyperlinkManager:
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
            
def demoAutoComplete():
        def buildList():
            import glob
            test_list = [os.path.basename(p) for p in glob.glob(os.path.join(os.getcwd(), '*'))]
            return test_list
        """Run a mini application to test the AutocompleteEntry Widget."""
        root = tk.Tk(className=' AutocompleteEntry demo')
        entry = AutoCompleteEntry(root, width=50, bg='white')
        entry.setBuildListFunction(buildList)
        entry.pack()
        entry.focus_set()
        root.mainloop()
    
def demo():
    root = tk.Tk(className='ToolTip-demo')
    l = tk.Listbox(root)
    l.insert('end', "I'm a listbox")
    l.pack(side='top')
    ToolTip(l, follow_mouse=1, text="I'm a tooltip with follow_mouse set to 1, so I won't be placed outside my parent")
    b = tk.Button(root, text='Quit', command=root.quit)
    b.pack(side='bottom')
    ToolTip(b, text='Enough of this')
    root.mainloop()
    
def demo2():
    import sys
    root = tk.Tk()
    root.title("LogTextArea demo")
    l = OutputTextArea(root, sys.argv[1])
    l.pack(side=tk.TOP)
    #b = tk.Button(root, text='Quit', command=root.quit)
    #b.pack(side='bottom')
    root.mainloop()    

if __name__ == '__main__':
    demoAutoComplete() 

