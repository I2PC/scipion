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
from __future__ import print_function

"""
Creates a scipion workflow file (json formatted) base on a template.
The template may have some ~placeholders~ that will be overwritten with values
Template may look like this, separator is "~" and within it you can define:
~title|value|type~
Template string sits at the end of the file ready for a running streaming demo.
"""
import os
import re
import Tkinter as tk
import tempfile
import tkFont

import datetime
import pyworkflow.utils as pwutils
from pyworkflow.gui import Message, Icon
from pyworkflow.config import ProjectSettings

import pyworkflow.gui as pwgui
from pyworkflow.gui.project.base import ProjectBaseWindow
from pyworkflow.gui.widgets import HotButton, Button

import traceback

FIELD_SEP = '~'
VIEW_WIZARD = 'wizardview'

# Session id
PROJECT_NAME = 'Project name'
MESSAGE = 'Message'

# Project regex to validate the session id name
PROJECT_PATTERN = "^\w{2}\d{4,6}-\d+$"
PROJECT_REGEX = re.compile(PROJECT_PATTERN)


class BoxWizardWindow(ProjectBaseWindow):
    """ Windows to manage all projects. """

    def __init__(self, **kwargs):
        try:
            title = '%s (%s on %s)' % ('Workflow template customizer',
                                       pwutils.getLocalUserName(),
                                       pwutils.getLocalHostName())
        except Exception:
            title = Message.LABEL_PROJECTS

        settings = ProjectSettings()
        self.generalCfg = settings.getConfig()

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
        btn = HotButton(btnFrame, text="Start demo",
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

        defaultProjectName = "demo_" + datetime.datetime.now().strftime("%I%M%S")
        self._addPair(PROJECT_NAME, 1, labelFrame, value=defaultProjectName)
        self._addPair(MESSAGE, 4, labelFrame, widget='label')

        labelFrame.columnconfigure(0, weight=1)
        labelFrame.columnconfigure(0, minsize=120)
        labelFrame.columnconfigure(1, weight=1)

        labelFrame2 = tk.LabelFrame(frame, text=' Acquisition values ', bg='white',
                                    font=self.bigFontBold)

        labelFrame2.grid(row=1, column=0, sticky='nw', padx=20, pady=10)
        labelFrame2.columnconfigure(0, minsize=120)

        self.addFieldsFromTemplate(labelFrame2)

        frame.columnconfigure(0, weight=1)

    def _addPair(self, text, r, lf, widget='entry', traceCallback=None,
                 mouseBind=False, value=None):
        label = tk.Label(lf, text=text, bg='white',
                         font=self.bigFont)
        label.grid(row=r, column=0, padx=(10, 5), pady=2, sticky='ne')

        if not widget:
            return

        var = tk.StringVar()

        if value is not None:
            var.set(value)

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

        self.vars[text] = var
        widget.grid(row=r, column=1, sticky='nw', padx=(5, 10), pady=2)

    def _addCheckPair(self, label, r, lf, col=1):

        var = tk.IntVar()

        cb = tk.Checkbutton(lf, text=label, font=self.bigFont, bg='white',
                            variable=var)
        self.vars[label] = var
        self.checkvars.append(label)
        cb.grid(row=r, column=col, padx=5, sticky='nw')

    def addFieldsFromTemplate(self, labelFrame2):

        self._template = getTemplateSplit()

        self._fields = getFields(self._template)

        row = 2
        for field in self._fields.values():
            self._addPair(field.getTitle(), row, labelFrame2,
                          value=field.getValue())
            row += 1

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
        for field in self._fields.values():
            newValue = self._getValue(field.getTitle())
            field.setValue(newValue)
            if not field.validate():
                errors.append("%s value does not validate. Value: %s, Type: %s"
                              % (field.getTitle(), field.getValue(),
                                 field.getType()))

        # Do more checks only if there are not previous errors
        if errors:
            errors.insert(0, "*Errors*:")
            self.windows.showError("\n  - ".join(errors))
        else:

            workflow = self._createTemplate()
            if workflow is not None:
                # Create the project
                self.createProjectFromWorkflow(workflow)

                self.windows.close()
                return

    def createProjectFromWorkflow(self, workflow):

        projectName = self._getValue(PROJECT_NAME)

        scipion = os.path.join(os.environ.get('SCIPION_HOME'), 'scipion')
        scriptsPath = os.path.join(os.environ.get('SCIPION_HOME'), 'scripts')

        # Download the required data
        # pwutils.runCommand(scipion +
        #                     " testdata --download jmbFalconMovies")

        # Create the project
        createProjectScript = os.path.join(scriptsPath, 'create_project.py')
        pwutils.runCommand(scipion + " python " + createProjectScript + " " +
                           projectName + " " + workflow)

        # Schedule the project
        scheduleProjectScript = os.path.join(scriptsPath, 'schedule_project.py')
        pwutils.runCommand(scipion + " python " + scheduleProjectScript + " " +
                           projectName)

        # Launch scipion
        pwutils.runCommand(scipion + " project " + projectName)

    def _createTemplate(self):

        try:
            # Where to write the json file.
            (fileHandle, path) = tempfile.mkstemp()

            replaceFields(self._fields.values(), self._template)

            finalJson = "".join(self._template)

            os.write(fileHandle, finalJson)
            os.close(fileHandle)

            print("New workflow saved at " + path)

        except Exception as e:
            self.windows.showError(
                "Couldn't create the template.\n" + e.message)
            traceback.print_exc()
            return None

        return path


class FormField(object):
    def __init__(self, index, title, value=None, varType=None):
        self._index = index
        self._title = title
        self._value = value
        self._type = varType

    def getTitle(self):
        return self._title

    def getIndex(self):
        return self._index

    def getType(self):
        return self._type

    def getValue(self):
        return self._value

    def setValue(self, value):
        self._value = value

    def validate(self):
        return validate(self._value, self._type)


""" FIELDS VALIDATION """
""" FIELDS TYPES"""
FIELD_TYPE_STR = "0"
FIELD_TYPE_BOOLEAN = "1"
FIELD_TYPE_PATH = "2"
FIELD_TYPE_INTEGER = "3"
FIELD_TYPE_DECIMAL = "4"


""" VALIDATION METHODS"""
def validate(value, fieldType):
    if fieldType == FIELD_TYPE_BOOLEAN:
        return validBoolean(value)
    elif fieldType == FIELD_TYPE_DECIMAL:
        return validDecimal(value)
    elif fieldType == FIELD_TYPE_INTEGER:
        return validInteger(value)
    elif fieldType == FIELD_TYPE_PATH:
        return validPath(value)
    elif fieldType == FIELD_TYPE_STR:
        return validString(value)

    else:
        print("Type %s for %snot recognized. Review the template." \
              % (type, value))
        return


def validString(value):
    return value is not None


def validInteger(value):
    return value.isdigit()


def validPath(value):
    return os.path.exists(value)


def validDecimal(value):

    try:
        float(value)
        return True
    except Exception as e:
        return False


def validBoolean(value):
    return value is True or value is False


def getFields(template):

    def fieldStr2Field(fieldIndex, fieldString):
        fieldLst = fieldString.split('|')

        title = fieldLst[0]
        defaultValue = fieldLst[1] if len(fieldLst) >= 2 else None
        varType = fieldLst[2] if len(fieldLst) >= 3 else None

        return FormField(fieldIndex, title, defaultValue, varType)

    fields = {}
    # For each field found in the template
    for index in xrange(1, len(template), 2):
        field = fieldStr2Field(index, template[index])
        fields[field.getTitle()] = field

    return fields


def replaceFields(fields, template):

    for field in fields:
        template[field.getIndex()] = field.getValue()


def getTemplateSplit():
    # Get the fields definition from the template
    templateStr = getTemplate()

    # Split the template by the field separator
    return templateStr.split(FIELD_SEP)


def getTemplate():

    template = """[
    {
        "object.className": "ProtImportMovies",
        "object.id": "2",
        "object.label": "Import movies",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "importFrom": 0,
        "filesPath": "~Movie's folder|%(SCIPION_HOME)s/data/tests/jmbFalconMovies/|2~",
        "filesPattern": "*.mrcs",
        "copyFiles": false,
        "haveDataBeenPhaseFlipped": false,
        "acquisitionWizard": null,
        "voltage": ~Voltage|300|4~,
        "sphericalAberration": ~Spherical aberration|2.7|4~,
        "amplitudeContrast": ~Amplitude contrast|0.1|4~,
        "magnification": 59000,
        "samplingRateMode": 0,
        "samplingRate": ~Sampling rate|1.34|4~,
        "scannedPixelSize": 7.0,
        "doseInitial": 0.0,
        "dosePerFrame": ~Dose per frame|0.0|4~,
        "gainFile": null,
        "darkFile": null,
        "dataStreaming": true,
        "timeout": 7200,
        "fileTimeout": 30,
        "inputIndividualFrames": false,
        "numberOfIndividualFrames": null,
        "stackFrames": false,
        "writeMoviesInProject": false,
        "movieSuffix": "_frames.mrcs",
        "deleteFrames": false,
        "streamingSocket": false,
        "socketPort": 5000
    },
    {
        "object.className": "XmippProtMovieGain",
        "object.id": "313",
        "object.label": "xmipp3 - movie gain",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "frameStep": 30,
        "movieStep": ~Gain estimation movie step|1|4~,
        "useExistingGainImage": false,
        "hostName": "localhost",
        "numberOfThreads": 1,
        "numberOfMpi": 1,
        "inputMovies": "2.outputMovies"
    },
    {
        "object.className": "ProtMotionCorr",
        "object.id": "56",
        "object.label": "Motioncorr",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "gpuMsg": "True",
        "GPUIDs": "0",
        "alignFrame0": 1,
        "alignFrameN": 0,
        "useAlignToSum": true,
        "sumFrame0": 1,
        "sumFrameN": 0,
        "binFactor": 1.0,
        "cropOffsetX": 0,
        "cropOffsetY": 0,
        "cropDimX": 0,
        "cropDimY": 0,
        "doSaveAveMic": true,
        "doSaveMovie": false,
        "doComputePSD": false,
        "doComputeMicThumbnail": false,
        "computeAllFramesAvg": false,
        "extraParams": "",
        "useMotioncor2": true,
        "doApplyDoseFilter": false,
        "patchX": 5,
        "patchY": 5,
        "group": 1,
        "tol": 0.5,
        "doMagCor": false,
        "useEst": true,
        "scaleMaj": 1.0,
        "scaleMin": 1.0,
        "angDist": 0.0,
        "defectFile": null,
        "extraParams2": "",
        "doSaveUnweightedMic": true,
        "hostName": "localhost",
        "numberOfThreads": 1,
        "numberOfMpi": 1,
        "inputMovies": "2.outputMovies"
    },
    {
        "object.className": "ProtCTFFind",
        "object.id": "118",
        "object.label": "Ctffind",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "recalculate": false,
        "sqliteFile": null,
        "ctfDownFactor": 1.0,
        "useCtffind4": true,
        "astigmatism": 100.0,
        "findPhaseShift": false,
        "minPhaseShift": 0.0,
        "maxPhaseShift": 3.15,
        "stepPhaseShift": 0.2,
        "resamplePix": true,
        "lowRes": 0.05,
        "highRes": 0.35,
        "minDefocus": 0.25,
        "maxDefocus": 4.0,
        "windowSize": 256,
        "hostName": "localhost",
        "numberOfThreads": 1,
        "numberOfMpi": 1,
        "inputMicrographs": "56.outputMicrographs"
    },
    {
        "object.className": "ProtMonitorSummary",
        "object.id": "164",
        "object.label": "Summary Monitor",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "inputProtocols": ["2","56","118","313"],
        "samplingInterval": 60,
        "stddevValue": 0.04,
        "ratio1Value": 1.15,
        "ratio2Value": 4.5,
        "maxDefocus": 40000.0,
        "minDefocus": 1000.0,
        "astigmatism": 2000.0,
        "monitorTime": 30000.0,
        "cpuAlert": 101.0,
        "memAlert": 101.0,
        "swapAlert": 101.0,
        "doGpu": false,
        "gpusToUse": "0",
        "doNetwork": false,
        "netInterfaces": 1,
        "doDiskIO": false,
        "doMail": true,
        "emailFrom": "noreply-biocomp@cnb.csic.es",
        "emailTo": "user@domain",
        "smtp": "localhost",
        "publishCmd": ""
    },
    {
        "object.className": "SparxGaussianProtPicking",
        "object.id": "284",
        "object.label": "eman2 - sparx gaussian picker",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "boxSize": 200,
        "lowerThreshold": 1.0,
        "higherThreshold": 1.2,
        "gaussWidth": 1.0,
        "inputMicrographs": "56.outputMicrographs"
    },
    {
        "object.className": "XmippProtExtractParticles",
        "object.id": "323",
        "object.label": "xmipp3 - extract particles",
        "object.comment": "",
        "runName": null,
        "runMode": 0,
        "downsampleType": 0,
        "downFactor": 1.0,
        "boxSize": 200,
        "doBorders": true,
        "doRemoveDust": true,
        "thresholdDust": 5.0,
        "doInvert": false,
        "doFlip": true,
        "doNormalize": true,
        "normType": 2,
        "backRadius": -1,
        "hostName": "localhost",
        "numberOfThreads": 4,
        "numberOfMpi": 1,
        "ctfRelations": "118.outputCTF",
        "inputCoordinates": "284.outputCoordinates"
    }
]"""
    # Replace environment variables
    template = template % os.environ

    return template


if __name__ == "__main__":
    wizWindow = BoxWizardWindow()
    wizWindow.show()
