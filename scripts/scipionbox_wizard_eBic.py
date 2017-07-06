#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
import os
import re
import Tkinter as tk
import ttk
import tkFont
from collections import OrderedDict
from ConfigParser import ConfigParser

import pyworkflow as pw
import pyworkflow.utils as pwutils
from pyworkflow.gui import Message, Icon
from pyworkflow.config import ProjectSettings

import pyworkflow.gui as pwgui
from pyworkflow.gui.project.base import ProjectBaseWindow
from pyworkflow.gui.widgets import HotButton, Button

VIEW_WIZARD = 'wizardview'

# Session id
PROJECT_NAME = 'PROJECT_NAME'
# Where the json template is
SCIPION_WORKFLOW_TEMPLATE = 'SCIPION_WORKFLOW_TEMPLATE'
# Where to save the new custimized workflow
JSON_DESTINATION = 'JSON_DESTINATION'

# Variables to enter and write in the template
DOSE = 'DOSE'
PIXEL_SIZE = 'PIXEL_SIZE'

vars2Use = [DOSE, PIXEL_SIZE]

# Define some string constants
LABELS = {
    PROJECT_NAME: "Session id",
    DOSE: "Dose",
    PIXEL_SIZE: "Pixel size",
}

MICROSCOPE = "Microscope"
MESSAGE = 'Message'

# Project regex to validate the session id name
PROJECT_PATTERN = "em\d{5}-\d$"
PROJECT_REGEX = re.compile(PROJECT_PATTERN)


class BoxWizardWindow(ProjectBaseWindow):
    """ Windows to manage all projects. """

    def __init__(self, config, **kwargs):
        try:
            title = '%s (%s on %s)' % ('Workflow template customizer',
                                       pwutils.getLocalUserName(),
                                       pwutils.getLocalHostName())
        except Exception:
            title = Message.LABEL_PROJECTS

        settings = ProjectSettings()
        self.generalCfg = settings.getConfig()

        self.config = config
        ProjectBaseWindow.__init__(self, title, minsize=(400, 550), **kwargs)
        self.viewFuncs = {VIEW_WIZARD: BoxWizardView}
        self.switchView(VIEW_WIZARD)


class BoxWizardView(tk.Frame):
    def __init__(self, parent, windows, **kwargs):
        tk.Frame.__init__(self, parent, bg='white', **kwargs)
        self.windows = windows
        self.root = windows.root
        self.vars = {}
        self.checkvars = []
        self.microscope = None
        # Regular expression to validate username and sample name
        self.re = re.compile('\A[a-zA-Z][a-zA-Z0-9_-]+\Z')

        # tkFont.Font(size=12, family='verdana', weight='bold')
        bigSize = pwgui.cfgFontSize + 2
        smallSize = pwgui.cfgFontSize - 2
        fontName = pwgui.cfgFontName

        self.bigFont = tkFont.Font(size=bigSize, family=fontName)
        self.bigFontBold = tkFont.Font(size=bigSize, family=fontName,
                                       weight='bold')

        self.projDateFont = tkFont.Font(size=smallSize, family=fontName)
        self.projDelFont = tkFont.Font(size=smallSize, family=fontName,
                                       weight='bold')
        # Header section
        headerFrame = tk.Frame(self, bg='white')
        headerFrame.grid(row=0, column=0, sticky='new')
        headerText = "Enter your desired values"

        label = tk.Label(headerFrame, text=headerText,
                         font=self.bigFontBold,
                         borderwidth=0, anchor='nw', bg='white',
                         fg=pwgui.Color.DARK_GREY_COLOR)
        label.grid(row=0, column=0, sticky='nw', padx=(20, 5), pady=10)

        # Body section
        bodyFrame = tk.Frame(self, bg='white')
        bodyFrame.grid(row=1, column=0, sticky='news')
        self._fillContent(bodyFrame)

        # Add the create project button
        btnFrame = tk.Frame(self, bg='white')
        btn = HotButton(btnFrame, text="Create template",
                        font=self.bigFontBold,
                        command=self._onAction)
        btn.grid(row=0, column=1, sticky='ne', padx=10, pady=10)

        # Add the Import project button
        btn = Button(btnFrame, Message.LABEL_BUTTON_CANCEL,
                     Icon.ACTION_CLOSE,
                     font=self.bigFontBold,
                     command=self.windows.close)
        btn.grid(row=0, column=0, sticky='ne', padx=10, pady=10)

        btnFrame.grid(row=2, column=0, sticky='sew')
        btnFrame.columnconfigure(0, weight=1)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

    def _fillContent(self, frame):
        labelFrame = tk.LabelFrame(frame, text=' General ', bg='white',
                                   font=self.bigFontBold)
        labelFrame.grid(row=0, column=0, sticky='nw', padx=20)

        def _addPair(key, r, lf, widget='entry', traceCallback=None,
                     mouseBind=False):
            t = LABELS.get(key, key)
            label = tk.Label(lf, text=t, bg='white',
                             font=self.bigFont)
            label.grid(row=r, column=0, padx=(10, 5), pady=2, sticky='ne')

            if not widget:
                return

            var = tk.StringVar()

            if widget == 'entry':
                widget = tk.Entry(lf, width=30, font=self.bigFont,
                                  textvariable=var)
                if traceCallback:
                    if mouseBind:  # call callback on click
                        widget.bind("<Button-1>", traceCallback, "eee")
                    else:  # call callback on type
                        var.trace('w', traceCallback)
            elif widget == 'label':
                widget = tk.Label(lf, font=self.bigFont, textvariable=var)

            self.vars[key] = var
            widget.grid(row=r, column=1, sticky='nw', padx=(5, 10), pady=2)

        def _addComboPair(key, r, lf, traceCallback=None):
            t = LABELS.get(key, key)
            label = tk.Label(lf, text=t, bg='white',
                             font=self.bigFont)
            label.grid(row=r, column=0, padx=(10, 5), pady=2, sticky='ne')

            var = tk.StringVar()
            combo = ttk.Combobox(lf, textvariable=var, state='readonly')
            combo['values'] = self._getConfig().keys()

            if traceCallback:
                combo.bind('<<ComboboxSelected>>', traceCallback)
            self.vars[key] = var
            combo.grid(row=r, column=1, sticky='nw', padx=(5, 10), pady=2)

            return combo

        self.micCombo = _addComboPair(MICROSCOPE, 0, labelFrame,
                                      traceCallback=self._onMicroscopeChanged)
        _addPair(PROJECT_NAME, 1, labelFrame, traceCallback=self._onInputChange)
        _addPair(MESSAGE, 4, labelFrame, widget='label')

        labelFrame.columnconfigure(0, weight=1)
        labelFrame.columnconfigure(0, minsize=120)
        labelFrame.columnconfigure(1, weight=1)

        labelFrame2 = tk.LabelFrame(frame, text=' Pre-processing ', bg='white',
                                    font=self.bigFontBold)

        labelFrame2.grid(row=1, column=0, sticky='nw', padx=20, pady=10)
        labelFrame2.columnconfigure(0, minsize=120)

        _addPair(DOSE, 1, labelFrame2)
        _addPair(PIXEL_SIZE, 2, labelFrame2)

        frame.columnconfigure(0, weight=1)

        self._onInputChange()

    def _getConfValue(self, key, default=''):
        return self.windows.config[self.microscope].get(key, default)

    def _getConfig(self):
        return self.windows.config

    def _getVar(self, varKey):
        return self.vars[varKey]

    def _getValue(self, varKey):
        return self.vars[varKey].get()

    def _setValue(self, varKey, value):
        return self.vars[varKey].set(value)

    # noinspection PyUnusedLocal
    def _onAction(self, e=None):
        errors = []
        # Check the entered data
        # Do more checks only if there are not previous errors

        if not self._validProjectName():
            errors.append("Session ID should start with em "
                          "followed by five digits, a - and a final number."
                          "Examples: \n   em12345-1\n   em55555-1\n")

        if errors:
            errors.insert(0, "*Errors*:")
            self.windows.showError("\n  - ".join(errors))
        else:
            if self._createTemplate():
                self.windows.close()
                return

    def _createTemplate(self):

        try:
            # Compose the destination using the session id (
            destination = self._getConfValue(JSON_DESTINATION)
            sessionId = self._getValue(PROJECT_NAME)
            destination = destination.format(sessionId)

            if not self._validateDestination(destination): return False

            # Get a dictionary with only the values to write in the template
            replacementDict = self._vars2Dict()

            # Get the template and open it
            template = self._getConfValue(SCIPION_WORKFLOW_TEMPLATE)
            templateStr = open(template, 'r').read()

            # Replace values
            workflowStr = templateStr.format(**replacementDict)

            text_file = open(destination, "w")
            text_file.write(workflowStr)
            text_file.close()

            self.windows.showInfo("New workflow saved at " + destination)

        except Exception as e:
            self.windows.showError(
                "Couldn't create the template.\n" + e.message)
            return False

        return True

    def _vars2Dict(self):

        replacementDict = {}
        for varName in vars2Use:
            self._var2Dict(varName, replacementDict)

        return replacementDict

    def _var2Dict(self, varName, replacementDict):

        value = self._getValue(varName)
        replacementDict[varName] = value

    def _getProjectName(self):
        return self._getValue(PROJECT_NAME)

    def _validProjectName(self):
        return PROJECT_REGEX.match(self._getProjectName())

    def _isInputReady(self):
        return self.microscope is not None and self._validProjectName()

    def _onInputChange(self, *args):
        if self.microscope is None:
            msg = 'Select microscope'
        elif not self._validProjectName():
            msg = ('Invalid session id. Not matching ' + PROJECT_PATTERN)
        else:
            msg = ''

        self._setValue(MESSAGE, msg)

    def _onMicroscopeChanged(self, *args):
        self.microscope = self._getValue(MICROSCOPE)
        self.micCombo.selection_clear()
        # Check/uncheck different options
        for key in self.checkvars:
            self.vars[key].set(int(self._getConfValue(key, 0)))

        self._onInputChange()

    def _validateDestination(self, destination):

        destinationDir = os.path.dirname(destination)

        print destinationDir

        if os.path.exists(destination):
            answer = self.windows.askYesNo("Replace?",
                                           "There is already a customized workflow at " +
                                           destination +"\n Do you want to replace it?")
            if answer:
                return True
            else:
                return False

        elif not os.path.exists(destinationDir):
            self.windows.showError("The destination folder does not exists:\n" +
                                   destinationDir, "Folder missing")
            return False

        return True


def createDictFromConfig():
    """ Read the configuration from scipion/config/scipionbox.conf.
     A dictionary will be created where each key will be a section starting
     by MICROSCOPE:, all variables that are in the GLOBAL section will be
     inherited by default.
    """
    # Read from users' config file.
    confGlobalDict = OrderedDict()
    confDict = OrderedDict()

    cp = ConfigParser()
    cp.optionxform = str  # keep case (stackoverflow.com/questions/1611799)

    confFile = pw.getConfigPath("scipionbox.conf")

    print "Reading conf file: ", confFile
    cp.read(confFile)

    GLOBAL = "GLOBAL"
    MICROSCOPE = "MICROSCOPE:"

    if not GLOBAL in cp.sections():
        raise Exception("Missing section %s in %s" % (GLOBAL, confFile))
    else:
        for opt in cp.options(GLOBAL):
            confGlobalDict[opt] = cp.get(GLOBAL, opt)

    for section in cp.sections():
        if section != GLOBAL and section.startswith(MICROSCOPE):
            sectionDict = OrderedDict(confGlobalDict)
            for opt in cp.options(section):
                sectionDict[opt] = cp.get(section, opt)
            confDict[section.replace(MICROSCOPE, '')] = sectionDict

    return confDict


if __name__ == "__main__":
    confDict = createDictFromConfig()

    wizWindow = BoxWizardWindow(confDict)
    wizWindow.show()
