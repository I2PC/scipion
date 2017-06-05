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
# *  e-mail address 'scipion@cnb.csic.es'
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
        TreeProvider.__init__(self)
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

    def validate(self):

        validationMsg = None

        if len(self.textVar.get().strip()) == 0:
            validationMsg = "Label name can't be empty.\n"

        if validationMsg is not None:
            dialog.showError("Validation error", validationMsg, self)
            return False

        return True



