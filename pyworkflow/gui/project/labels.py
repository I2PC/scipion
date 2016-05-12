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
Tree widget implementation.
"""
        
import os
import Tkinter as tk
import ttk

from pyworkflow.config import Label
from pyworkflow.gui import Icon, configureWeigths
from pyworkflow.gui.tree import TreeProvider
import pyworkflow.gui.dialog as dialog


class LabelsTreeProvider(TreeProvider):
    """ Populate Tree from Labels. """
    def __init__(self, objList=None):
        self.objList = objList
        self._parentDict = {}

    def getColumns(self):
        return [('name', 300), ('color', 150)]

    def getObjectInfo(self, label):

        return {'key': label.getId(), 'parent': None,
                'text': label.getName(), 'values': (label.getColor()),
                'tags': label.getColor()}

    def getObjectPreview(self, obj):
        return (None, None)

    def getObjectActions(self, obj):
        return []

    def _getObjectList(self):
        """Retrieve the object list"""
        return self.objList

    def getObjects(self):
        objList = self._getObjectList()
        return objList

    def configureTags(self, tree):
        for label in self.getObjects():
            tree.tag_configure(label.getColor(), background=label.getColor())


class LabelsDialog(dialog.ToolbarListDialog):
    """
    This class extend from ListDialog to allow an
    extra toolbar to handle operations over the elements
    in the list (e.g. Edit, New, Delete).
    """
    def __init__(self, parent, labels, **kwargs):
        """ From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback:
                    a callback function to validate selected items.
                allowSelect: if set to False, the 'Select' button will not
                    be shown.
        """
        self.labels = labels
        toolbarButtons = [
            dialog.ToolbarButton('Add', self._addLabel, Icon.ACTION_NEW),
            dialog.ToolbarButton('Edit', self._editLabel, Icon.ACTION_EDIT),
            dialog.ToolbarButton('Delete', self._deleteLabel, Icon.ACTION_DELETE)
        ]

        dialog.ToolbarListDialog.__init__(self, parent,
                                   "Manage labels",
                                   LabelsTreeProvider(labels),
                                   "Select the label to edit or delete",
                                   toolbarButtons,
                                   allowsEmptySelection=True,
                                   itemDoubleClick=self._editLabel,
                                   **kwargs)

    def _newColor(self):
        """ Pick a color by default for a given label from a predefined list.
         Check that the color have not been used by other label.
        """
        colors = ["#e57373", "#4fc3f7", "#81c784", "#ff8a65", "#9575cd",
                  "#a1887f", "#ffd54f", "#dce775", "#4db6ac"]

        for c in colors:
            if all(l.getColor().lower() != c for l in self.labels):
                return c

        return 'red'


    def _addLabel(self, e=None):
        label = Label(color=self._newColor())
        dlg = EditLabelDialog(self, "Add label", label)
        if dlg.resultYes():
            self.labels.addLabel(label)
            self.tree.update()

    def _editLabel(self, e=None):
        selection = self.tree.getSelectedObjects()
        if selection:
            label = selection[0]
            dlg = EditLabelDialog(self, "Edit label", label)
            if dlg.resultYes():
                self.tree.update()

    def _deleteLabel(self, e=None):
        selection = self.tree.getSelectedObjects()
        if selection:
            labelsStr = '\n'.join('- %s' % l.getName() for l in selection)
            if dialog.askYesNo("Delete a label",
                               "Are you sure to delete the "
                               "following label(s)?\n %s" % labelsStr, self):
                for label in selection:
                    self.labels.deleteLabel(label)
                self.tree.update()


class EditLabelDialog(dialog.Dialog):
    """ Dialog to edit a label (name, color) """
    def __init__(self, parent, title, label, **kwargs):
        self.label = label
        dialog.Dialog.__init__(self, parent, title)

    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        configureWeigths(bodyFrame, 1, 1)

        # Label
        label_text = tk.Label(bodyFrame, text="Name", bg='white', bd=0)
        label_text.grid(row=0, column=0, sticky='nw', padx=(15, 10), pady=15)
        # Label box
        var = tk.StringVar()
        var.set(self.label.getName())
        self.textVar = var
        self.textLabel = tk.Entry(bodyFrame, width=20, textvariable=var)
        self.textLabel.grid(row=0, column=1, sticky='news', padx=5, pady=5)

        # Comment
        colorLabel = tk.Label(bodyFrame, text='Color \n(Click to change)',
                              bg='white', bd=0)
        colorLabel.grid(row=1, column=0, sticky='nw', padx=(15, 10), pady=15)
        self.colorVar = tk.StringVar()
        self.colorVar.set(self.label.getColor())
        self.colorBox = tk.Frame(bodyFrame, bg=self.colorVar.get())
        self.colorBox.grid(row=1, column=1, sticky='news', padx=5, pady=5)
        colorLabel.bind('<Button-1>', self._changeColor)
        self.colorBox.bind('<Button-1>', self._changeColor)

    def apply(self):
        self.label.setName(self.textVar.get())
        self.label.setColor(self.colorVar.get())

    def _changeColor(self, e=None):
        hexColor = dialog.askColor(self.colorVar.get())
        if hexColor is not None:
            self.colorBox.config(bg=hexColor)
            self.colorVar.set(hexColor)

# class TableDialogButtonDefinition:
#     """
#     TableDialog button definition: text, value, handler and icon
#     """
#
#     def __init__(self, text, value, icon=None, handler= None):
#         self.text = text
#         self.value = value
#         self.icon = icon
#         self.handler = handler
#
#
# class TableDialogConfiguration:
#     """
#     Configuration to be passed to a Table Dialog
#      It holds button definitions (text, icons,..) and handlers
#     """
#     def __init__(self, unloadButtonsList=None, onOkHandler=None, selectmode='extended', title='Table', message='', onDoubleClickHandler = None):
#
#         self.unloadButtons = unloadButtonsList or []
#         self.toolBarButtons = []
#         self.onOkHandler = onOkHandler
#         self.onDoubleClickHandler = onDoubleClickHandler
#         self.selectmode = selectmode
#         self.title = title
#         self.message = message
#
#     def addUnloadButton(self, tableButton):
#
#         self.unloadButtons.append(tableButton)
#
#     def addToolBarButton(self, tableButton):
#         self.toolBarButtons.append(tableButton)
#
#
# def createDefaultTableDialogConfiguration(okHandler, doubleClickHandler=None):
#
#     conf = TableDialogConfiguration(onDoubleClickHandler=doubleClickHandler)
#
#     btnDefYes = TableDialogButtonDefinition('Select', RESULT_YES)
#     conf.addUnloadButton(btnDefYes)
#
#     btnDefCancel = TableDialogButtonDefinition('Cancel', RESULT_CANCEL)
#     conf.addUnloadButton(btnDefCancel)
#
#     conf.onOkHandler = okHandler or defaultOkHandler
#
#
#     return conf
#
#
# def defaultOkHandler(values):
#
#     if not values:
#         return "Please select an element"
#     else:
#         return ''
#
#
# def emptyOkHandler(values):
#     return ''
#
#
# class TableDialog(Dialog):
#     """
#     Dialog to show a table (It is implemented using a Tree widget.)
#     it has 2 modes:
#         'select': It will show the table and will allow you to select elements from the list (1 or more)
#         'edit': In this case it will show buttons to perform GUI CUD actions (Create, Update and Delete).
#     """
#     EVENT_ON_DOUBLE_CLICK = 'onDoubleClick'
#     IS_SELECTED = 'isSelected'
#
#     def __init__(self, parent, provider=None, conf=None):
#         """ From kwargs:
#                 message: message tooltip to show when browsing.
#                 selected: the item that should be selected.
#                 validateSelectionCallback: a callback function to validate selected items.
#         """
#         self.values = []
#         self.provider = provider
#
#         if conf is None:
#             conf = createDefaultTableDialogConfiguration()
#
#         self.conf = conf
#         self._selectmode = self.conf.selectmode
#         self.message = conf.message
#         Dialog.__init__(self, parent, conf.title,
#                         buttons=self._getUnloadButtons())
#
#
#     def _getUnloadButtons(self):
#
#         buttons = []
#
#         for buttonDef in self.conf.unloadButtons:
#
#             buttons.append((buttonDef.text, buttonDef.value))
#
#         return buttons
#
#     def _createButton(self, frame, text, buttonId):
#         icon = self.icons[buttonId]
#
#         # Get the handler
#         handler = self._getHandler(buttonId)
#
#         return tk.Button(frame, text=text, image=self.getImage(icon), compound=tk.LEFT,
#                          command=handler)
#
#     def _getHandler(self, buttonId):
#
#         return lambda: self._handleResult(buttonId)
#
#     def body(self, bodyFrame):
#         bodyFrame.config(bg='white')
#         gui.configureWeigths(bodyFrame)
#
#         # Create the table (Tree widget)
#         self._createTree(bodyFrame)
#
#         # Add the message
#         self.addMessage(bodyFrame)
#
#         #Add a tool bar
#         self.addToolbar(bodyFrame)
#
#         # Set the focus
#         self.initial_focus = self.tree
#
#     def addMessage(self, bodyFrame):
#         if self.message:
#             label = tk.Label(bodyFrame, text=self.message, bg='white',
#                              image=self.getImage(Icon.LIGHTBULB), compound=tk.LEFT)
#             label.grid(row=2, column=0, sticky='nw', padx=5, pady=5)
#
#     def addToolbar(self, bodyFrame):
#
#         if self.conf.toolBarButtons:
#
#             # Create the Action Buttons TOOLBAR
#             toolbar = tk.Frame(self, bg='white')
#             toolbar.grid(row=2, column=0, sticky='nw')
#             gui.configureWeigths(toolbar)
#             innerToolbar = tk.Frame(toolbar, bg='white')
#             innerToolbar.grid(row=0, column=0, sticky='sw')
#
#             def addToolbarButton(tableButton):
#                 btn = tk.Label(innerToolbar, text=tableButton.text,
#                                image=self.getImage(tableButton.icon or None),
#                                compound=tk.LEFT, cursor='hand2', bg='white')
#                 btn.bind('<Button-1>', lambda e: self._toolbarButtonHandler(tableButton.value))
#                 return btn
#
#             index = 1
#             # For each toolbar button
#             for tbButton in self.conf.toolBarButtons:
#
#
#                 button = addToolbarButton(tbButton)
#                 button.grid(row=0, column=index, sticky='ns')
#                 index += 1
#
#     def _toolbarButtonHandler(self, value):
#
#         handler = self._getToolBarHandlerByValue(value)
#
#         if handler is not None:
#             handler(self)
#
#     def _getToolBarHandlerByValue(self, value):
#
#         for button in self.conf.toolBarButtons:
#             if button.value == value:
#                 return button.handler
#
#     def _createTree(self, parent):
#         self.tree = BoundTree(parent, self.provider, selectmode=self._selectmode)
#
#     def apply(self):
#         self.values = self.tree.getSelectedObjects()
#
#     def validate(self):
#         self.apply()  # load self.values with selected items
#         err = ''
#
#         if self.conf.onOkHandler:
#             err = self.conf.onOkHandler(self.values)
#
#         if err:
#             showError("Validation error", err, self)
#             return False
#
#         return True


