# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import sys
import Tkinter as tk
from matplotlib import animation

from pyworkflow.gui.plotter import plt
from pyworkflow.gui.tree import BoundTree
from pyworkflow.gui.widgets import Button, HotButton
import pyworkflow.gui as pwgui
import pyworkflow.gui.text as text
import pyworkflow.utils as pwutils

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, Viewer
from pyworkflow.em.protocol import (
    ProtMonitorCTF, ProtMonitorMovieGain, ProtMonitorSystem,
    ProtMonitorSummary, SummaryProvider)

from .plotter import EmPlotter


class ProtMonitorCTFViewer(Viewer):
    _environments = [DESKTOP_TKINTER]
    _label = 'ctf monitor'
    _targets = [ProtMonitorCTF]

    def _visualize(self, obj, **kwargs):
        return [CtfMonitorPlotter(obj.createMonitor())]


class CtfMonitorPlotter(EmPlotter):
    def __init__(self, monitor):
        EmPlotter.__init__(self, windowTitle="CTF Monitor")
        self.monitor = monitor
        self.y2 = 0.
        self.y1 = 100.
        self.win = 250  # number of samples to be ploted
        self.step = 50  # self.win  will be modified in steps of this size

        self.createSubPlot(self._getTitle(), "Micrographs", "Defocus (A)")
        self.fig = self.getFigure()
        self.ax = self.getLastSubPlot()
        self.ax.margins(0.05)
        self.ax.grid(True)
        self.oldWin = self.win

        self.lines = {}
        self.init = True
        self.stop = False

    def _getTitle(self):
        return ("Use scrool wheel to change view window (win=%d)\n "
                "S stops, C continues plotting" % self.win)

    def onscroll(self, event):

        if event.button == 'up':
            self.win += self.step
        else:
            self.win -= self.step
            if self.win < self.step:
               self.win = self.step

        if self.oldWin != self.win:
            self.ax.set_title(self._getTitle())
            self.oldWin= self.win
        self.animate()
        EmPlotter.show(self)

    def press(self,event):

        sys.stdout.flush()
        if event.key == 'S':
            self.stop = True
            self.ax.set_title('Plot is Stopped. Press C to continue plotting')
        elif event.key == 'C':
            self.ax.set_title(self._getTitle())
            self.stop = False
            self.animate()
        EmPlotter.show(self)

    def has_been_closed(self,ax):
        fig = ax.figure.canvas.manager
        active_fig_managers = plt._pylab_helpers.Gcf.figs.values()
        return fig not in active_fig_managers

    def animate(self, i=0): #do NOT remove i
                           #FuncAnimation adds it as argument
        if self.stop:
            return

        data = self.monitor.getData()
        self.x = data['idValues']
        for k,v in self.lines.iteritems():
            self.y = data[k]

            lenght = len(self.x)
            imin = max(0,len(self.x) - self.win)
            xdata = self.x[imin:lenght]
            ydata = self.y[imin:lenght]
            v.set_data(xdata,ydata)

        self.ax.relim()
        self.ax.autoscale()
        self.ax.grid(True)
        self.ax.legend(loc=2).get_frame().set_alpha(0.5)
        self.ax.axhline(y=self.monitor.minDefocus, c="red",
                        linewidth=0.5, linestyle='dashed', zorder=0)
        self.ax.axhline(y=self.monitor.maxDefocus, c="red",
                        linewidth=0.5, linestyle='dashed', zorder=0)

    def show(self):
        self.paint(['defocusU','defocusV'])

    def paint(self, labels):
        for label in labels:
            if (label == 'defocusU'):
                self.lines[label], = self.ax.plot([], [], '-o',
                                                  label=label, color='b')
            else:
                self.lines[label], = self.ax.plot([], [], '-o',
                                                  label=label, color='r')

        anim = animation.FuncAnimation(self.fig, self.animate,
                                       interval=self.monitor.samplingInterval*1000)  # miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)


class ProtMonitorMovieGainViewer(Viewer):
    _environments = [DESKTOP_TKINTER]
    _label = 'movie gain monitor'
    _targets = [ProtMonitorMovieGain]

    def _visualize(self, obj, **kwargs):
        return [MovieGainMonitorPlotter(obj.createMonitor())]


class MovieGainMonitorPlotter(EmPlotter):
    def __init__(self, monitor):
        EmPlotter.__init__(self, windowTitle="Movie Gain Monitor")
        self.monitor = monitor
        self.y2 = 0.
        self.y1 = 100.
        self.win = 250 # number of samples to be plotted
        self.step = 50 # self.win will be modified in steps of this size

        self.createSubPlot(self._getTitle(), "", "")
        self.fig = self.getFigure()
        self.ax = self.getLastSubPlot()
        self.ax2 = self.ax.twinx()
        self.ax2.margins(0.05)
        self.oldWin = self.win

        self.lines = {}
        self.init = True
        self.stop = False

    def _getTitle(self):
        return ("Use scrool wheel to change view window (win=%d)\n "
                "S stops, C continues plotting" % self.win)

    def onscroll(self, event):

        if event.button == 'up':
            self.win += self.step
        else:
            self.win -= self.step
            if self.win < self.step:
               self.win = self.step

        if self.oldWin != self.win:
            self.ax.set_title(self._getTitle())
            self.oldWin= self.win
        self.animate()
        EmPlotter.show(self)

    def press(self,event):

        sys.stdout.flush()
        if event.key == 'S':
            self.stop = True
            self.ax.set_title('Plot is Stopped. Press C to continue plotting')
        elif event.key == 'C':
            self.ax.set_title(self._getTitle())
            self.stop = False
            self.animate()
        EmPlotter.show(self)

    def has_been_closed(self,ax):
        fig = ax.figure.canvas.manager
        active_fig_managers = plt._pylab_helpers.Gcf.figs.values()
        return fig not in active_fig_managers

    def animate(self, i=0): #do NOT remove i
                           #FuncAnimation adds it as argument
        if self.stop:
            return

        data = self.monitor.getData()
        self.x = data['idValues']
        for k,v in self.lines.iteritems():
            self.y = data[k]

            lenght = len(self.x)
            imin = max(0,len(self.x) - self.win)
            xdata = self.x[imin:lenght]
            ydata = self.y[imin:lenght]
            v.set_data(xdata,ydata)

        self.ax.set_ylabel('Ratios between specified percentiles',
                           color='b', size=10)
        self.ax.spines['left'].set_color('blue')
        self.ax.spines['right'].set_color('red')
        self.ax.tick_params(axis='y', labelsize=8, colors='blue')
        self.ax.relim()
        self.ax.autoscale()
        self.ax.grid(False)
        self.ax.get_xaxis().set_visible(False)
        self.ax2.set_ylabel('Residual gain standard deviation',
                            color='r', size=10)
        self.ax2.tick_params(axis='y', labelsize=8, colors='red')
        self.ax2.relim()
        self.ax2.autoscale()
        lines, labels = self.ax.get_legend_handles_labels()
        lines2, labels2 = self.ax2.get_legend_handles_labels()
        self.ax2.legend(lines + lines2, labels + labels2,
                        loc=2, prop={'size':10}).get_frame().set_alpha(0.5)

    def show(self):
        self.paint(['ratio1','ratio2','standard_deviation'])

    def paint(self, labels):
        for label in labels:
            if (label == 'standard_deviation'):
                self.lines[label], = \
                    self.ax2.plot([], [], '-o',
                                  label='Standard deviation',
                                  color='r')
            if (label == 'ratio1'):
                self.lines[label], = \
                    self.ax.plot([], [], '-o',
                                 label='97.5/2.5 percentile',
                                 color='b')
            if (label == 'ratio2'):
                self.lines[label], = \
                    self.ax.plot([], [], '-*',
                                 label='max/97.5 percentile',
                                 color='b')

        anim = animation.FuncAnimation(self.fig, self.animate,
                                       interval=self.monitor.samplingInterval*1000)#miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

    def empty(self):
        return self.init



class ProtMonitorSystemViewer(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'system monitor'
    _targets = [ProtMonitorSystem]

    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def _visualize(self, obj, **kwargs):
        return [SystemMonitorPlotter(obj.createMonitor(),
                                     nifName=self.protocol.nifsNameList[
                                         self.protocol.netInterfaces.get()])]


class SystemMonitorPlotter(EmPlotter):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'System Monitor'
    _targets = [ProtMonitorSystem]

    def __init__(self, monitor, nifName=None):
        EmPlotter.__init__(self, windowTitle="system Monitor")
        self.monitor = monitor
        self.y2 = 0.
        self.y1 = 100.
        self.win = 250  # number of samples to be ploted
        self.step = 50  # self.win  will be modified in steps of this size
        self.createSubPlot(self._getTitle(),
                           "time (hours)",
                           "percentage (or MB for IO or NetWork)")
        self.fig = self.getFigure()
        self.ax = self.getLastSubPlot()
        self.ax.margins(0.05)
        self.ax.grid(True)
        self.oldWin = self.win

        self.lines = {}
        self.init = True
        self.stop = False

        self.nifName = nifName

    def _getTitle(self):
        return ("Use scrool wheel to change view window (win=%d)\n "
                "S stops, C continues plotting. Toggle ON/OFF GPU_X "
                "by pressing X\n"
                "c/n/d toggle ON-OFF cpu/network/disk usage\n" % self.win)

    def onscroll(self, event):

        if event.button == 'up':
            self.win += self.step
        else:
            self.win -= self.step
            if self.win < self.step:
                self.win = self.step

        if self.oldWin != self.win:
            self.ax.set_title(self._getTitle())
            self.oldWin = self.win
        self.animate()
        EmPlotter.show(self)

    def press(self, event):
        def numericKey(key):
            self.colorChanged = True
            number = int(key)
            index = 3+number*3
            if (index + 3) > self.lenPlots:
                return
            if self.color['gpuMem_%d' % number] != 'w':
                self.oldColor['gpuMem_%d' % number] = \
                    self.color['gpuMem_%d' % number]
                self.oldColor['gpuUse_%d' % number] = \
                    self.color['gpuUse_%d' % number]
                self.oldColor['gpuTem_%d' % number] = \
                    self.color['gpuTem_%d' % number]
                self.color['gpuMem_%d' % number] = "w"
                self.color['gpuUse_%d' % number] = "w"
                self.color['gpuTem_%d' % number] = "w"
            else:
                self.color['gpuMem_%d' % number] = \
                    self.oldColor['gpuMem_%d' % number]
                self.color['gpuUse_%d' % number] = \
                    self.oldColor['gpuUse_%d' % number]
                self.color['gpuTem_%d' % number] = \
                    self.oldColor['gpuTem_%d' % number]

        def cpuKey(key):
            self.colorChanged = True
            if self.color['cpu'] != 'w':
                self.oldColor['cpu'] = self.color['cpu']
                self.oldColor['mem'] = self.color['mem']
                self.oldColor['swap'] = self.color['swap']
                self.color['cpu'] = "w"
                self.color['mem'] = "w"
                self.color['swap'] = "w"
            else:
                self.color['cpu'] = self.oldColor['cpu']
                self.color['swap'] = self.oldColor['swap']
                self.color['mem'] = self.oldColor['mem']

        def netKey(key):
            self.colorChanged = True
            if self.color['%s_send' % self.nifName] != 'w':
                self.oldColor['%s_send' % self.nifName] = \
                    self.color['%s_send' % self.nifName]
                self.oldColor['%s_recv' % self.nifName] = \
                    self.color['%s_recv' % self.nifName]
                self.color['%s_send' % self.nifName] = "w"
                self.color['%s_recv' % self.nifName] = "w"
            else:
                self.color['%s_send' % self.nifName] = \
                    self.oldColor['%s_send' % self.nifName]
                self.color['%s_recv' % self.nifName] = \
                    self.oldColor['%s_recv' % self.nifName]

        def diskKey(key):
            self.colorChanged = True
            if self.color['disk_read'] != 'w':
                self.oldColor['disk_read'] = self.color['disk_read']
                self.oldColor['disk_write'] = self.color['disk_write']
                self.color['disk_read'] = "w"
                self.color['disk_write'] = "w"
            else:
                self.color['disk_read'] = self.oldColor['disk_read']
                self.color['disk_write'] = self.oldColor['disk_write']

        sys.stdout.flush()
        if event.key == 'S':
            self.stop = True
            self.ax.set_title('Plot has been Stopped. '
                              'Press C to continue plotting')
        elif event.key == 'C':
            self.ax.set_title(self._getTitle())
            self.stop = False
            self.animate()
        elif event.key.isdigit():
            numericKey(event.key)
            self.animate()
        elif event.key == 'c':
            cpuKey(event.key)
            self.animate()
        elif event.key == 'n':
            netKey(event.key)
            self.animate()
        elif event.key == 'd':
            diskKey(event.key)
            self.animate()
        EmPlotter.show(self)

    def has_been_closed(self, ax):
        fig = ax.figure.canvas.manager
        active_fig_managers = plt._pylab_helpers.Gcf.figs.values()
        return fig not in active_fig_managers

    def animate(self, i=0):  # do NOT remove i

        if self.stop:
            return

        data = self.monitor.getData()
        self.x = data['idValues']
        for k, v in self.lines.iteritems():
            self.y = data[k]

            lenght = len(self.x)
            imin = max(0, len(self.x) - self.win)
            xdata = self.x[imin:lenght]
            ydata = self.y[imin:lenght]
            v.set_data(xdata, ydata)
            if self.colorChanged:
                v.set_color(self.color[k])
        self.colorChanged = False
        self.ax.relim()
        self.ax.autoscale()
        self.ax.grid(True)
        self.ax.legend(loc=2).get_frame().set_alpha(0.5)

    def paint(self, labels):
        for label in labels:
            self.lines[label], = self.ax.plot([], [], '-', label=label)

        anim = animation.FuncAnimation(
                self.fig, self.animate,
                interval=self.monitor.samplingInterval * 1000)  # miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

    def show(self):
        colortypes = ["k", "b", "r", "g", "y", "c", "m"]
        lenColortypes = len(colortypes)
        self.colorChanged = True
        self.color = {}
        self.oldColor = {}
        counter = 0
        for key in self.monitor.getLabels():
            self.color[key] = colortypes[counter % lenColortypes]
            self.oldColor[key] = colortypes[counter % lenColortypes]
            counter += 1
        self.lenPlots = len(self.color)
        self.paint(self.monitor.getLabels())


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
        buttonsFrame.grid(row=2, column=0, sticky='new', padx=5, pady=5)

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
            tk.Label(lf, text=t2, font=self.fontBold).grid(row=self.r,
                                                           column=1,
                                                           sticky='nw',
                                                           padx=(5, 25),
                                                           pady=0)
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
        subframe = tk.Frame(frame)
        subframe.grid(row=0, column=0, sticky='nw')
        frame.columnconfigure(1, weight=1)

        ctfBtn = Button(subframe, "CTF Monitor", command=self._monitorCTF)
        ctfBtn.grid(row=0, column=0, sticky='nw', padx=(0, 5))
        if self.protocol.createCtfMonitor() is None:
            ctfBtn['state'] = 'disabled'

        movieGainBtn = Button(subframe, "Movie Gain Monitor",
                              command=self._monitorMovieGain)
        movieGainBtn.grid(row=0, column=1, sticky='nw', padx=(0, 5))
        if self.protocol.createMovieGainMonitor() is None:
            movieGainBtn['state'] = 'disabled'

        sysBtn = Button(subframe, "System Monitor",
                        command=self._monitorSystem)
        sysBtn.grid(row=0, column=2, sticky='nw', padx=(0, 5))
        if self.protocol.createSystemMonitor() is None:
            sysBtn['state'] = 'disabled'

        htmlBtn = HotButton(subframe, 'Open HTML Report',
                            command=self._openHTML)
        htmlBtn.grid(row=0, column=3, sticky='nw', padx=(0, 5))

        closeBtn = self.createCloseButton(frame)
        closeBtn.grid(row=0, column=1, sticky='ne')

    def _monitorCTF(self, e=None):
        CtfMonitorPlotter(self.protocol.createCtfMonitor()).show()

    def _monitorMovieGain(self, e=None):
        MovieGainMonitorPlotter(self.protocol.createMovieGainMonitor()).show()

    def _monitorSystem(self, e=None):
        nifName = self.protocol.nifsNameList[self.protocol.netInterfaces.get()]
        SystemMonitorPlotter(self.protocol.createSystemMonitor(),
                             nifName).show()

    def _updateLabel(self):
        self.updateVar.set('Updated: %s' % pwutils.prettyTime(secs=True))
        # Schedule a refresh in some seconds
        self.tree.after(self.refresh * 1000, self._updateData)

    def _updateData(self):
        self.provider.refreshObjects()
        self.tree.update()
        self._updateLabel()

    def _openHTML(self, e=None):
        reportHtml = self.protocol.createHtmlReport()
        reportPath = reportHtml.reportPath
        if pwutils.exists(reportPath):
            text._open_cmd(reportPath)
        else:
            self.showInfo('Your html file is not ready yet. Please try again in a minute.')

