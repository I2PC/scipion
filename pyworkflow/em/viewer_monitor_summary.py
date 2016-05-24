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

import os
import Tkinter as tk

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.gui.widgets import Button, HotButton
import pyworkflow.gui.text as text
import pyworkflow.gui as pwgui

import pyworkflow.protocol.params as params
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer

from protocol.protocol_import import ProtImportImages
from protocol.protocol_monitor_summary import ProtMonitorSummary



class ViewerMonitorSummary(Viewer):
    """ Wrapper to visualize PDF objects. """
    _environments = [DESKTOP_TKINTER]
    _targets = [ProtMonitorSummary]

    def visualize(self, obj, **kwargs):
        self.summaryWindow = self.tkWindow(SummaryWindow,
                                           title='Scipion-Box Summary',
                                           protocol=obj
                                           )
        self.summaryWindow.show()


class SummaryProvider(TreeProvider):
    """Create the tree elements for a Protocol run"""
    def __init__(self, protocol):
        self.protocol = protocol
        self.getColumns = lambda: [('Name', 300), ('Output', 150),
                                   ('Number', 100)]
        self._parentDict = {}
        self.acquisition = []
        self.refreshObjects()

    def getObjects(self):
        return self._objects

    def refreshObjects(self):
        objects = []

        def addObj(objId, name, output='', size='', parent=None):
            obj = pwobj.Object(objId=objId)
            obj.name = name
            obj.output = output
            obj.outSize = size
            obj._objParent = parent
            objects.append(obj)
            return obj

        runs = [p.get() for p in self.protocol.inputProtocols]
        g = self.protocol.getProject().getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        for n in nodes:
            prot = n.run
            pobj = addObj(prot.getObjId(),
                          '%s (id=%s)' % (prot.getRunName(), prot.strId()))

            for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                outSet.load()
                outSet.loadAllProperties()
                addObj(outSet.getObjId(), '', outName, outSet.getSize(), pobj)
                outSet.close()
                # Store acquisition parameters in case of the import protocol
                if isinstance(prot, ProtImportImages):
                    self.acquisition = [("Microscope Voltage: ",
                                         prot.voltage.get()),
                                        ("Spherical aberration: ",
                                         prot.sphericalAberration.get()),
                                        ("Magnification: ",
                                         prot.magnification.get()),
                                        ("Pixel Size (A/px): ",
                                         outSet.getSamplingRate())
                                        ]

        self._objects = objects

    def getObjectInfo(self, obj):
        info = {'key': obj.strId(),
                'parent': obj._objParent,
                'text': obj.name,
                'values': (obj.output, obj.outSize),
                'open': True
               }

        return info


class SummaryWindow(pwgui.Window):
    def __init__(self, **kwargs):
        pwgui.Window.__init__(self, **kwargs)

        self.protocol = kwargs.get('protocol')
        self.refresh = self.protocol.samplingInterval.get()
        self.provider = SummaryProvider(self.protocol)

        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)

    def _createContent(self, content):
        topFrame = tk.Frame(content)
        content.columnconfigure(0, weight=1)
        topFrame.grid(row=0, column=0, sticky='new', padx=5, pady=5)

        treeFrame = tk.Frame(content)
        content.rowconfigure(1, weight=1)
        treeFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)

        buttonsFrame = tk.Frame(content)
        buttonsFrame.grid(row=2, column=0, sticky='ne', padx=5, pady=5)

        self._fillTreeFrame(treeFrame)
        # JMRT: We fill the top frame after the tree, to make sure
        # the provider has updated the Acquisition info
        self._fillTopFrame(topFrame)
        self._fillButtonsFrame(buttonsFrame)

    def _fillTopFrame(self, frame):
        p1 = tk.Label(frame, text='Project: ')
        p1.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        projName = self.protocol.getProject().getShortName()
        p2 = tk.Label(frame, text=projName, font=self.fontBold)
        p2.grid(row=0, column=1, sticky='nw', padx=5, pady=0)

        lf = tk.LabelFrame(frame, text='Acquisition')
        lf.grid(row=1, column=0, columnspan=2, sticky='new')
        lf.columnconfigure(0, weight=1)
        lf.columnconfigure(1, weight=1)
        self.r = 0

        def add(t1, t2):
            tk.Label(lf, text=t1).grid(row=self.r, column=0, sticky='ne',
                                       padx=(10, 5), pady=(5, 0))
            tk.Label(lf, text=t2, font=self.fontBold).grid(row=self.r, column=1,
                                                           sticky='nw',
                                                           padx=(5, 25), pady=0)
            self.r += 1

        for k, v in self.provider.acquisition:
            add(k, v)

    def _fillTreeFrame(self, frame):
        self.tree = BoundTree(frame, self.provider)
        self.tree.grid(row=0, column=0)
        self.updateVar = tk.StringVar()
        updateLabel = tk.Label(frame, textvariable=self.updateVar)
        updateLabel.grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        self._updateLabel()

    def _fillButtonsFrame(self, frame):
        pdfBtn = HotButton(frame, 'Generate PDF Report',
                           command=self._generatePDF)
        pdfBtn.grid(row=0, column=0, sticky='nw', padx=10, pady=2)

        closeBtn = self.createCloseButton(frame)
        closeBtn.grid(row=0, column=1, sticky='nw')

    def _updateLabel(self):
        self.updateVar.set('Updated: %s' % pwutils.prettyTime(secs=True))
        # Schedule a refresh in some seconds
        self.tree.after(self.refresh * 1000, self._updateData)

    def _updateData(self):
        self.provider.refreshObjects()
        self.tree.update()
        self._updateLabel()

    def _generatePDF(self, e=None):
        reportName = 'report_%s.tex' % pwutils.prettyTimestamp()
        reportPath = self.protocol._getExtraPath(reportName)

        acquisitionLines = ''
        for item in self.provider.acquisition:
            acquisitionLines += '%s & %s \\\\\n' % item

        runLines = ''
        for obj in self.provider.getObjects():
            if obj.name:
                runLines += ' %s &  &  \\\\\n' % obj.name
            else:
                runLines += ' & %s & %s \\\\\n' % (obj.output, obj.outSize)

        args = {'acquisitionLines': acquisitionLines,
                'runLines': runLines,
                'dateStr': pwutils.prettyTime(secs=True),
                'projectName': self.protocol.getProject().getShortName(),
                'scipionVersion': os.environ['SCIPION_VERSION']
                }

        reportTemplate = """
\documentclass{article}

\usepackage[table]{xcolor}
\setlength{\\tabcolsep}{15pt}
\\renewcommand{\\arraystretch}{1.5}
\\begin{document}

%% Heading
\hfil{\Huge\\bf Scipion Report}\hfil
\hfill \\break
\hrule

\\begin{itemize}
\item Date: \\textbf{%(dateStr)s}
\item Project: \\textbf{%(projectName)s}
\item Scipion version: \\textbf{%(scipionVersion)s}
\end{itemize}

\hfill \\break
\hfill \\break

\\begin{tabular}{ |p{5cm} p{5cm}|  }
\hline
\multicolumn{2}{|c|}{\\textbf{Acquisition}} \\\\
\hline
%(acquisitionLines)s
\hline
\end{tabular}

\hfill \\break
\hfill \\break
\hfill \\break

\\begin{tabular}{ |p{5cm} p{3cm} p{3cm}|  }
\hline
\multicolumn{3}{|c|}{\\textbf{Runs Summary}} \\\\
\hline
%(runLines)s
\hline
\end{tabular}

\end{document}

        """
        reportFile = open(reportPath, 'w')
        reportFile.write(reportTemplate % args)
        reportFile.close()
        pwutils.runJob(None, 'pdflatex',
                       '-interaction=nonstopmode ' + reportName,
                       cwd=self.protocol._getExtraPath())
        text._open_cmd(reportPath.replace('.tex', '.pdf'))