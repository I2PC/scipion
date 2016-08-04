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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Module to handling Dialogs
some code was taken from tkSimpleDialog
"""
        
import Tkinter as tk

import gui
from tree import BoundTree
from text import Text, TaggedText
from pyworkflow.utils.properties import Message, Icon
from tkColorChooser import askcolor as _askColor

# Possible result values for a Dialog
RESULT_YES = 0
RESULT_NO = 1
RESULT_CANCEL = 2


class Dialog(tk.Toplevel):
    _images = {} #Images cache
    """Implementation of our own dialog to display messages
    It will have by default a three buttons: YES, NO and CANCEL
    Subclasses can rename the labels of the buttons like: OK, CLOSE or others
    The buttons(and theirs order) can be changed.
    An image name can be passed to display left to the message.
    """
    
    def __init__(self, parent, title, **args):
        """Initialize a dialog.
        Arguments:
            parent -- a parent window (the application window)
            title -- the dialog title
        **args accepts:
            buttons -- list of buttons tuples containing which buttons to display
        """
        
        if parent is None:
            parent = tk.Tk()
            parent.withdraw()
            gui.setCommonFonts()
    
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
        self.initial_focus = None

        bodyFrame = tk.Frame(self)
        # Call subclass method body to create that region
        self.body(bodyFrame)
        bodyFrame.grid(row=0, column=0, sticky='news',
                       padx=5, pady=5)

        self.icons = {RESULT_YES: Icon.BUTTON_SELECT, 
                      RESULT_NO: Icon.BUTTON_CLOSE,
                      RESULT_CANCEL: Icon.BUTTON_CANCEL}
        
        self.buttons = args.get('buttons', [('OK', RESULT_YES),
                                            ('Cancel', RESULT_CANCEL)])
        self.defaultButton = args.get('default', 'OK')
        btnFrame = tk.Frame(self)
        # Create buttons 
        self.buttonbox(btnFrame)
        btnFrame.grid(row=1, column=0, sticky='sew',
                      padx=5, pady=(0, 5))
        
        gui.configureWeigths(self)


        if self.initial_focus is None:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        if self.parent is not None:
            self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                      parent.winfo_rooty()+50))

        self.deiconify() # become visible now

        self.initial_focus.focus_set()

        # wait for window to appear on screen before calling grab_set
        self.wait_visibility()
        self.grab_set()
        self.wait_window(self)
        

    def destroy(self):
        """Destroy the window"""
        self.initial_focus = None
        tk.Toplevel.destroy(self)

    #
    # construction hooks

    def body(self, master):
        """create dialog body.
        return widget that should have initial focus.
        This method should be overridden, and is called
        by the __init__ method.
        """
        pass

        
    def _createButton(self, frame, text, result):
        icon = self.icons[result]
        return  tk.Button(frame, text=text, image=self.getImage(icon), compound=tk.LEFT,
                         command=lambda: self._handleResult(result))
        
    def buttonbox(self, btnFrame):
        frame = tk.Frame(btnFrame)
        btnFrame.columnconfigure(0, weight=1)
        frame.grid(row=0, column=0)
        col = 0
        for btnLabel, btnResult in self.buttons:
            btn = self._createButton(frame, btnLabel, btnResult)
            btn.grid(row=0, column=col, padx=5, pady=5)
            if (btnLabel == self.defaultButton and 
                self.initial_focus is None):
                self.initial_focus = btn
            col += 1
        self.bind("<Return>", self._handleReturn)
        self.bind("<Escape>", lambda e: self._handleResult(RESULT_CANCEL))


    def _handleResult(self, resultValue):
        """This method will be called when any button is pressed.
        It will set the resultValue associated with the button
        and close the Dialog"""
        self.result = resultValue
        noCancel = self.result != RESULT_CANCEL
        
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
        # Check which of the buttons is the default
        for button, result in self.buttons:
            if self.defaultButton == button:
                self._handleResult(result)

    def cancel(self, event=None):
        # put focus back to the parent window
        if self.parent is not None:
            self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):
        """validate the data
        This method is called automatically to validate the data before the
        dialog is destroyed. By default, it always validates OK.
        """
        return 1 # override

    def apply(self):
        """process the data
        This method is called automatically to process the data, *after*
        the dialog is destroyed. By default, it does nothing.
        """
        pass # override
    
    def getImage(self, imgName):
        """A shortcut to get an image from its name"""
        return gui.getImage(imgName, self._images)
    
    def getResult(self):
        return self.result
    
    def resultYes(self):
        return self.result == RESULT_YES
    
    def resultNo(self):
        return self.result == RESULT_NO
    
    def resultCancel(self):
        return self.result == RESULT_CANCEL

 
def fillMessageText(text, message):
    # Insert lines of text
    if isinstance(message, list):
        lines = message
    else:
        lines = message.splitlines()
    text.setReadOnly(False)
    text.clear()
    m = 0
    for l in lines:
        m = max(m, len(l))
        text.addLine(l)
    m = min(m + 5, 80)
    h = min(len(lines)+3, 30)
    text.config(height=h, width=m)
    text.addNewline()
    text.setReadOnly(True)
    
    
def createMessageBody(bodyFrame, message, image, 
                      frameBg='white',
                      textBg='white',
                      textPad=5):
    """ Create a Text containing the message.
    Params:
        bodyFrame: tk.Frame to be filled.
        msg: a str or list with the lines.
    """
    bodyFrame.config(bg=frameBg, bd=0)
    text = TaggedText(bodyFrame, bg=textBg, bd=0, highlightthickness=0)
    # Insert image
    if image:
        label = tk.Label(bodyFrame, image=image, bg=textBg, bd=0)
        label.grid(row=0, column=0, sticky='nw')
        
    text.frame.grid(row=0, column=1, sticky='news', 
                    padx=textPad, pady=textPad)
    fillMessageText(text, message)
    bodyFrame.rowconfigure(0, weight=1)
    bodyFrame.columnconfigure(1, weight=1)  
    
    return text  
        
           
class MessageDialog(Dialog):
    """Dialog subclasses to show message info, questions or errors.
    It can display an icon with the message"""
    def __init__(self, parent, title, msg, iconPath, **args):
        self.msg = msg
        self.iconPath = iconPath
        if not 'buttons' in args:
            args['buttons'] = [('OK', RESULT_YES)]
            args['default'] = 'OK'
        Dialog.__init__(self, parent, title, **args)
        
    def body(self, bodyFrame):
        self.image = gui.getImage(self.iconPath)
        createMessageBody(bodyFrame, self.msg, self.image)


class YesNoDialog(MessageDialog):
    """Ask a question with YES/NO answer"""
    def __init__(self, master, title, msg, **kwargs):
        buttonList =  [('Yes', RESULT_YES), ('No', RESULT_NO)]
        
        if kwargs.get('showCancel', False):
            buttonList.append(('Cancel', RESULT_CANCEL))
            
        MessageDialog.__init__(self, master, title, msg, 
                               'fa-exclamation-triangle_alert.png', default='No',
                               buttons=buttonList)        


class EntryDialog(Dialog):
    """Dialog to ask some entry"""
    def __init__(self, parent, title, entryLabel, entryWidth=20, defaultValue='', headerLabel=None):
        self.entryLabel = entryLabel
        self.entryWidth = entryWidth
        self.headerLabel = headerLabel
        self.tkvalue = tk.StringVar()
        self.tkvalue.set(defaultValue)
        self.value = None
        Dialog.__init__(self, parent, title)
        
    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        frame = tk.Frame(bodyFrame, bg='white')
        frame.grid(row=0, column=0, padx=20, pady=20)
        row = 0
        if self.headerLabel:
            label = tk.Label(bodyFrame, text=self.headerLabel, bg='white', bd=0)
            label.grid(row=row, column=0, columnspan=2, sticky='nw', padx=(15, 10), pady=15)
            row += 1
        label = tk.Label(bodyFrame, text=self.entryLabel, bg='white', bd=0)
        label.grid(row=row, column=0, sticky='nw', padx=(15, 10), pady=15)
        self.entry = tk.Entry(bodyFrame, bg=gui.cfgEntryBgColor, width=self.entryWidth, textvariable=self.tkvalue)
        self.entry.grid(row=row, column=1, sticky='new', padx=(0,15), pady=15)
        self.initial_focus = self.entry
        
    def apply(self):
        self.value = self.entry.get()
        
    def validate(self):
        if len(self.entry.get().strip()) == 0:
            showError("Validation error", "Value is empty", self)
            return False
        return True


class EditObjectDialog(Dialog):
    """Dialog to edit some text"""
    def __init__(self, parent, title, obj, mapper, **kwargs):
        
        self.obj = obj
        self.mapper = mapper
        
        self.textWidth = 5
        self.textHeight = 1
        self.labelText = kwargs.get('labelText', Message.TITLE_LABEL)
        self.valueText= self.obj.getObjLabel()
        
        self.commentLabel = Message.TITLE_COMMENT
        self.commentWidth = 50
        self.commentHeight = 15
        self.valueComment = self.obj.getObjComment()
        
        Dialog.__init__(self, parent, title)
        
    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        frame = tk.Frame(bodyFrame, bg='white')
        frame.grid(row=0, column=0, padx=20, pady=20)
        
        # Label
        label_text = tk.Label(bodyFrame, text=self.labelText, bg='white', bd=0)
        label_text.grid(row=0, column=0, sticky='nw', padx=(15, 10), pady=15)
        # Label box
        var = tk.StringVar()
        var.set(self.valueText)
        self.textLabel = tk.Entry(bodyFrame, width=self.textWidth, textvariable=var)
        self.textLabel.grid(row=0, column=1, sticky='news', padx=5, pady=5)
        
        # Comment
        label_comment = tk.Label(bodyFrame, text=self.commentLabel, bg='white', bd=0)
        label_comment.grid(row=1, column=0, sticky='nw', padx=(15, 10), pady=15)
        # Comment box
        self.textComment = Text(bodyFrame, height=self.commentHeight, 
                         width=self.commentWidth)
        self.textComment.setReadOnly(False)
        self.textComment.setText(self.valueComment)
        self.textComment.grid(row=1, column=1, sticky='news', padx=5, pady=5)
        self.initial_focus = self.textComment
        
    def getLabel(self):
        return self.textLabel.get()
    
    def getComment(self):
        return self.textComment.getText()
    
    def apply(self):
        self.obj.setObjLabel(self.getLabel())
        self.obj.setObjComment(self.getComment())
        
        if self.obj.hasObjId():
            self.mapper.store(self.obj)
            self.mapper.commit()
        
    def buttonbox(self, btnFrame):
        # Cancel the binding of <Return> key
        Dialog.buttonbox(self, btnFrame)
        #self.bind("<Return>", self._noReturn)
        self.unbind("<Return>")
        
    def _noReturn(self, e):
        pass
        

""" Functions to display dialogs """
def askYesNo(title, msg, parent):
    d = YesNoDialog(parent, title, msg)
    return d.resultYes()

def askYesNoCancel(title, msg, parent):
    d = YesNoDialog(parent, title, msg, showCancel=True)
    return d.result

def showInfo(title, msg, parent):
    MessageDialog(parent, title, msg, 'fa-info-circle_alert.png')

def showWarning(title, msg, parent):
    MessageDialog(parent, title, msg, 'fa-exclamation-triangle_alert.png')
    
def showError(title, msg, parent):
    MessageDialog(parent, title, msg, 'fa-times-circle_alert.png')


def askString(title, label, parent, entryWidth=20, defaultValue='', headerLabel=None):
    d = EntryDialog(parent, title, label, entryWidth, defaultValue, headerLabel)
    return d.value


def askColor(defaultColor='black'):
    (rgbcolor, hexcolor) = _askColor(defaultColor)
    return hexcolor

                
class ListDialog(Dialog):
    """
    Dialog to select an element from a list.
    It is implemented using the Tree widget.
    """
    def __init__(self, parent, title, provider, message=None, **kwargs):
        """ From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback:
                    a callback function to validate selected items.
                allowSelect: if set to False, the 'Select' button will not
                    be shown.
                allowsEmptySelection: if set to True, it will not validate
                    that at least one element was selected.
        """
        self.values = []
        self.provider = provider
        self.message = message
        self.validateSelectionCallback = kwargs.get('validateSelectionCallback',
                                                    None)
        self._selectmode = kwargs.get('selectmode', 'extended')
        self._allowsEmptySelection = kwargs.get('allowsEmptySelection', False)

        buttons = []
        if kwargs.get('allowSelect', True):
            buttons.append(('Select', RESULT_YES))
        buttons.append(('Cancel', RESULT_CANCEL))

        Dialog.__init__(self, parent, title, buttons=buttons)
        
    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        gui.configureWeigths(bodyFrame)
        self._createTree(bodyFrame)
        if self.message:
            label = tk.Label(bodyFrame, text=self.message, bg='white',
                     image=self.getImage(Icon.LIGHTBULB), compound=tk.LEFT)
            label.grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        self.initial_focus = self.tree
        
    def _createTree(self, parent):
        self.tree = BoundTree(parent, self.provider, selectmode=self._selectmode)
        
    def apply(self):
        self.values = self.tree.getSelectedObjects()
    
    def validate(self):
        self.apply() # load self.values with selected items
        err = ''
        
        if self.values:
            if self.validateSelectionCallback:
                err = self.validateSelectionCallback(self.values)
        else:
            if not self._allowsEmptySelection:
                err = "Please select an element"
            
        if err:
            showError("Validation error", err, self)
            return False
        
        return True


class ToolbarButton():
    """
    Store information about the buttons that will be added to the toolbar.
    """
    def __init__(self, text, command, icon=None, tooltip=None):
        self.text = text
        self.command = command
        self.icon = icon
        self.tooltip = tooltip


class ToolbarListDialog(ListDialog):
    """
    This class extend from ListDialog to allow an
    extra toolbar to handle operations over the elements
    in the list (e.g. Edit, New, Delete).
    """
    def __init__(self, parent, title, provider,
                 message=None, toolbarButtons=None, **kwargs):
        """ From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback:
                    a callback function to validate selected items.
                allowSelect: if set to False, the 'Select' button will not
                    be shown.
        """
        self.toolbarButtons = toolbarButtons
        self._itemDoubleClick = kwargs.get('itemDoubleClick', None)
        ListDialog.__init__(self, parent, title, provider, message, **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        gui.configureWeigths(bodyFrame, 1, 0)

        # Add an extra frame to insert the Toolbar
        # and another one for the ListDialog's body
        self.toolbarFrame = tk.Frame(bodyFrame, bg='white')
        self.toolbarFrame.grid(row=0, column=0, sticky='new')
        if self.toolbarButtons:
            for i, b in enumerate(self.toolbarButtons):
                self._addButton(b, i)

        subBody = tk.Frame(bodyFrame)
        subBody.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        ListDialog.body(self, subBody)
        if self._itemDoubleClick:
            self.tree.itemDoubleClick = self._itemDoubleClick

    def _addButton(self, button, col):
        btn = tk.Label(self.toolbarFrame, text=button.text,
                       image=self.getImage(button.icon),
                       compound=tk.LEFT, cursor='hand2', bg='white')
        btn.grid(row=0, column=col, sticky='nw', padx=(5, 0), pady=(5, 0))
        btn.bind('<Button-1>', button.command)


class FlashMessage():
    def __init__(self, master, msg, delay=5, relief='solid', func=None):
        self.root = tk.Toplevel(master=master)
        #hides until know geometry
        self.root.withdraw()
        self.root.wm_overrideredirect(1)
        tk.Label(self.root, text="   %s   " % msg,
                 bd=1, bg='DodgerBlue4', fg='white').pack()
        gui.centerWindows(self.root, refWindows=master)
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
        

class FileBrowseDialog(Dialog):
    """Dialog to select files from the filesystem."""
    def __init__(self, parent, title, provider, message=None, **args):
        """ From args:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
        """
        self.value = None
        self.provider = provider
        self.message = args.get('message', None)
        Dialog.__init__(self, parent, title,
                        buttons=[('Select', RESULT_YES), ('Cancel', RESULT_CANCEL)])
        
    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        gui.configureWeigths(bodyFrame)
        self._createTree(bodyFrame)
        if self.message:
            label = tk.Label(bodyFrame, text=self.message, bg='white',
                     image=self.getImage('fa-lightbulb-o.png'), compound=tk.LEFT)
            label.grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        self.initial_focus = self.tree
        
    def _createTree(self, parent):
        self.tree = BoundTree(parent, self.provider)
        
    def apply(self):
        index = self.tree.index(self.tree.getFirst())
        self.value = self.tree._objects[index]
    
    def validate(self):
        if self.tree.getFirst() is None:
            showError("Validation error", "Please select an element", self)
            return False
        return True                
        
