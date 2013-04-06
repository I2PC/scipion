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
Module to handling Dialogs
some code was taken from tkSimpleDialog
"""
        
import Tkinter as tk

import gui
from widgets import Button
from text import TaggedText


class Dialog(tk.Toplevel):
    '''Implementation of our own dialog to display messages
    It will have by default a three buttons: YES, NO and CANCEL
    Subclasses can rename the labels of the buttons like: OK, CLOSE or others
    The buttons(and theirs order) can be changed.
    An image name can be passed to display left to the message.
    '''
    RESULT_YES = 0
    RESULT_NO = 1
    RESULT_CANCEL = 2
    
    def __init__(self, parent, title, **args):
        """Initialize a dialog.
        Arguments:
            parent -- a parent window (the application window)
            title -- the dialog title
        **args accepts:
            buttons -- list of buttons tuples containing which buttons to display
        """
        tk.Toplevel.__init__(self, parent)

        self.withdraw() # remain invisible for now
        # If the master is not viewable, don't
        # make the child transient, or else it
        # would be opened withdrawn
        if parent.winfo_viewable():
            self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        bodyFrame = tk.Frame(self)
        self.initial_focus = self.body(bodyFrame)
        bodyFrame.grid(row=0, column=0, sticky='news',
                       padx=5, pady=5)

        self.buttons = args.get('buttons', [('OK', Dialog.RESULT_YES),
                                            ('Cancel', Dialog.RESULT_CANCEL)])
        self.defaultButton = args.get('default', 'OK')
        btnFrame = tk.Frame(self)
        self.buttonbox(btnFrame)
        btnFrame.grid(row=1, column=0, sticky='sew',
                      padx=5, pady=(0, 5))
        
        gui.configureWeigths(self)


        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        if self.parent is not None:
            self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                      parent.winfo_rooty()+50))

        self.deiconify() # become visibile now

        self.initial_focus.focus_set()

        # wait for window to appear on screen before calling grab_set
        self.wait_visibility()
        self.grab_set()
        self.wait_window(self)
        

    def destroy(self):
        '''Destroy the window'''
        self.initial_focus = None
        tk.Toplevel.destroy(self)

    #
    # construction hooks

    def body(self, master):
        '''create dialog body.
        return widget that should have initial focus.
        This method should be overridden, and is called
        by the __init__ method.
        '''
        pass

        
    def _createButton(self, frame, text, result):
        return  Button(frame, text=text, width=7, 
                         command=lambda: self._handleResult(result))
        
    def buttonbox(self, btnFrame):
        frame = tk.Frame(btnFrame)
        btnFrame.columnconfigure(0, weight=1)
        frame.grid(row=0, column=0)
        col = 0
        for btnLabel, btnResult in self.buttons:
            btn = self._createButton(frame, btnLabel, btnResult)
            btn.grid(row=0, column=col, padx=5, pady=5)
            if btnLabel == self.defaultButton:
                self.initial_focus = btn
            col += 1
        self.bind("<Return>", self._handleReturn)
        self.bind("<Escape>", lambda: self._handleResult(Dialog.RESULT_CANCEL))


    def _handleResult(self, resultValue):
        """This method will be called when any button is pressed.
        It will set the resultValue associated with the button
        and close the Dialog"""
        self.result = resultValue
        noCancel = self.result != Dialog.RESULT_CANCEL
        
        if noCancel and not self.validate():
            self.initial_focus.focus_set() # put focus back
            return
        

        self.withdraw()
        self.update_idletasks()

        try:
            if noCancel:
                self.apply()
        finally:
            self.cancel()
            
    def _handleReturn(self, e=None):
        """Handle press return key"""
        self._handleResult(Dialog.RESULT_CANCEL)

    def cancel(self, event=None):
        # put focus back to the parent window
        if self.parent is not None:
            self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):
        '''validate the data
        This method is called automatically to validate the data before the
        dialog is destroyed. By default, it always validates OK.
        '''
        return 1 # override

    def apply(self):
        '''process the data
        This method is called automatically to process the data, *after*
        the dialog is destroyed. By default, it does nothing.
        '''
        pass # override

        
class MessageDialog(Dialog):
    """Dialog subclasses to show message info, questions or errors.
    It can display an icon with the message"""
    def __init__(self, parent, title, msg, iconPath, **args):
        self.msg = msg
        self.iconPath = iconPath
        if not 'buttons' in args:
            args['buttons'] = [('OK', Dialog.RESULT_YES)]
            args['default'] = 'OK'
        Dialog.__init__(self, parent, title, **args)
        
    def body(self, bodyFrame):
        self.image = gui.getImage(self.iconPath)
        bodyFrame.config(bg='white', bd=2)
        text = TaggedText(bodyFrame, bg='white', bd=0)
        # Insert image
        if self.image:
            label = tk.Label(bodyFrame, image=self.image, bg='white', bd=0)
            label.grid(row=0, column=0, sticky='nw')
        # Insert lines of text
        mylines = self.msg.splitlines()
        m = 0
        for l in mylines:
            m = max(m, len(l))
            text.addLine(l)
        m = min(m + 5, 80)
        h = min(len(mylines)+3, 30)
        text.config(height=h, width=m)
        text.addNewline()
        text.config(state=tk.DISABLED)
        text.frame.grid(row=0, column=1, sticky='news', padx=5, pady=5)
        
        bodyFrame.rowconfigure(0, weight=1)
        bodyFrame.columnconfigure(1, weight=1)


class YesNoDialog(MessageDialog):
    '''Ask a question with YES/NO answer'''
    def __init__(self, master, title, msg):
        MessageDialog.__init__(self, master, title, msg, 'warning.gif', default='No',
                               buttons=[('Yes', Dialog.RESULT_YES), ('No', Dialog.RESULT_NO)])        


class EntryDialog(Dialog):
    """Dialog to ask some entry"""
    def __init__(self, parent, title, entryLabel, entryWidth=20):
        self.entryLabel = entryLabel
        self.entryWidth = entryWidth
        self.value = None
        Dialog.__init__(self, parent, title)
        
    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        frame = tk.Frame(bodyFrame, bg='white')
        frame.grid(row=0, column=0, padx=20, pady=20)
        label = tk.Label(bodyFrame, text=self.entryLabel, bg='white', bd=0)
        label.grid(row=0, column=0, sticky='nw', padx=(15, 10), pady=15)
        self.entry = tk.Entry(bodyFrame, bg=gui.cfgEntryBgColor, width=self.entryWidth)
        self.entry.grid(row=0, column=1, sticky='new', padx=(0,15), pady=15)
        
    def apply(self):
        self.value = self.entry.get()
        
    def validate(self):
        if len(self.entry.get().strip()) == 0:
            showError("Validation error", "Value is empty", self)
            return False
        return True
        
''' Functions to display dialogs '''
def askYesNo(title, msg, parent):
    d = YesNoDialog(parent, title, msg)
    return d.result == Dialog.RESULT_YES

def showInfo(title, msg, parent):
    MessageDialog(parent, title, msg, 'info.gif')

def showWarning(title, msg, parent):
    MessageDialog(parent, title, msg, 'warning.gif')
    
def showError(title, msg, parent):
    MessageDialog(parent, title, msg, 'error.gif')
    
def askString(title, label, parent, entryWidth=20):
    d = EntryDialog(parent, title, label, entryWidth)
    return d.value
    
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
        w = Button(box, text="OK", width=7, command=self.ok)
        w.pack(side=tk.RIGHT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        box.pack()
        
    def apply(self):
        self.result = map(int, self.lb.curselection())
        
        
if __name__ == '__main__':
    import sys
    root = tk.Tk()
    root.withdraw()
    gui.setCommonFonts()
    #result = askYesNo("Confirm DELETE", "Are you sure to delete this?", root)
    #print result
    #showInfo('Testing Info', "this is a [really] important infor", root)
    
    #showError('Testing error', "Fatal Error due to <problems>", root)
    print askString("Enter project name", "Project name:", root)
    #root.mainloop() 
