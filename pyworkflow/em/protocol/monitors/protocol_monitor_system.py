# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor, Monitor
import sqlite3 as lite
import psutil
import time, sys

from pyworkflow import VERSION_1_1
from pyworkflow.gui.plotter import plt
from pyworkflow.protocol.constants import STATUS_RUNNING, STATUS_FINISHED
from pyworkflow.protocol import getProtocolFromDb
import sys
from pyworkflow.em.plotter import EmPlotter

SYSTEM_LOG_SQLITE = 'system_log.sqlite'


class ProtMonitorSystem(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'system_monitor'
    _lastUpdateVersion = VERSION_1_1

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self.dataBase = 'log.sqlite'
        self.tableName = 'log'

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        ProtMonitor._defineParams(self, form)
        form.addParam('cpuAlert', params.FloatParam,default=101,
              label="Raise Alarm if CPU > XX%",
              help="Raise alarm if memory allocated is greater "
                   "than given percentage")
        form.addParam('memAlert', params.FloatParam,default=101,
              label="Raise Alarm if Memory > XX%",
              help="Raise alarm if cpu allocated is greater "
                   "than given percentage")
        form.addParam('swapAlert', params.FloatParam,default=101,
              label="Raise Alarm if Swap > XX%",
              help="Raise alarm if swap allocated is greater "
                   "than given percentage")

        form.addParam('monitorTime', params.FloatParam,default=300,
              label="Total Logging time (min)",
              help="Log during this interval")

        ProtMonitor._sendMailParams(self, form)

    #--------------------------- STEPS functions ------------------------------

    def monitorStep(self):
        self.createMonitor().loop()

    def createMonitor(self):
        protocols = self.getInputProtocols()
        sysMonitor = MonitorSystem(protocols,
                                   workingDir=self.workingDir.get(),
                                   samplingInterval=self.samplingInterval.get(),
                                   monitorTime=self.monitorTime.get(),
                                   email=self.createEmailNotifier(),
                                   stdout=True,
                                   cpuAlert=self.cpuAlert.get(),
                                   memAlert=self.memAlert.get(),
                                   swapAlert=self.swapAlert.get())
        return sysMonitor

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        return ['Stores CPU, memory and swap ussage in percentage']

    def _methods(self):
        return []


class MonitorSystem(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocols, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocols = protocols
        self.cpuAlert = kwargs['cpuAlert']
        self.memAlert = kwargs['memAlert']
        self.swapAlert = kwargs['swapAlert']
        self._dataBase = kwargs.get('dbName', SYSTEM_LOG_SQLITE)
        self._tableName = kwargs.get('tableName', 'log')

        self.conn = lite.connect(os.path.join(self.workingDir, self._dataBase),
                                 isolation_level=None)
        self.cur = self.conn.cursor()

    def warning(self, msg):
        self.notify("Scipion System Monitor WARNING", msg)

    def initLoop(self):
        self._createTable()
        psutil.cpu_percent(True)
        psutil.virtual_memory()

    def step(self):
        cpu = psutil.cpu_percent(interval=0)
        mem = psutil.virtual_memory()
        swap = psutil.swap_memory()

        if self.cpuAlert < 100 and cpu > self.cpuAlert:
            self.warning("CPU allocation =%f." % cpu.percent)
            self.cpuAlert = cpu

        if self.memAlert < 100 and mem.percent > self.memAlert:
            self.warning("Memory allocation =%f." % mem.percent)
            self.memAlert = mem.percent

        if self.swapAlert < 100 and swap.percent > self.swapAlert:
            self.warning("SWAP allocation =%f." % swap.percent)
            self.swapAlert = swap.percent

        sql = """INSERT INTO %s(mem,cpu,swap) VALUES(%f,%f,%f);""" % (
            self._tableName, mem.percent, cpu, swap.percent)

        try:
            self.cur.execute(sql)
        except Exception as e:
            print("ERROR: saving one data point (monitor). I continue")

        # Return finished = True if all protocols have finished
        return all(self.getUpdatedProtocol(prot).getStatus() != STATUS_RUNNING
                   for prot in self.protocols)

    def _createTable(self):
        self.cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now','localtime')),
                                cpu FLOAT,
                                mem FLOAT,
                                swap FLOAT)
                                """ % self._tableName)

    def getData(self):
        cur = self.cur
        # I guess this can be done in a single call
        # I am storing the first meassurement
        cur.execute("select julianday(timestamp)  from %s where id=1" %
                    self._tableName )
        initTime = cur.fetchone()[0]
        cur.execute("select timestamp  from %s where id=1" % self._tableName)
        initTimeTitle = cur.fetchone()[0]
        cur.execute("select (julianday(timestamp) - %f)*24  from %s" %
                    (initTime,self._tableName) )
        idValues = [r[0] for r in cur.fetchall()]

        def get(name):
            try:
                cur.execute("select %s from %s" % (name, self._tableName))
            except Exception as e:
                print("ERROR readind data (plotter). I continue")
            return [r[0] for r in cur.fetchall()]

        data = {'initTime': initTime,
                'initTimeTitle': initTimeTitle,
                'idValues': idValues,
                'cpu': get('cpu'),
                'mem': get('mem'),
                'swap': get('swap'),
                }
        #conn.close()
        return data


from pyworkflow.viewer import ( DESKTOP_TKINTER, WEB_DJANGO, Viewer)
from matplotlib import animation


class ProtMonitorSystemViewer(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'system monitor'
    _targets = [ProtMonitorSystem]

    def _visualize(self, obj, **kwargs):
        return [SystemMonitorPlotter(obj.createMonitor())]


class SystemMonitorPlotter(EmPlotter):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'System Monitor'
    _targets = [ProtMonitorSystem]

    def __init__(self, monitor):
        EmPlotter.__init__(self, windowTitle="system Monitor")
        self.monitor = monitor
        self.y2 = 0.
        self.y1 = 100.
        self.win = 250 # number of samples to be ploted
        self.step = 50 # self.win  will be modified in steps of this size
        self.createSubPlot(self._getTitle(),
                           "time (hours)", "percentage (or Mb for IO)")
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

    def animate(self,i=0): #do NOT remove i
                                         #FuncAnimation add it as argument

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

    def paint(self,labels):
        for label in labels:
            self.lines[label],=self.ax.plot([], [], '-',label=label)

        anim = animation.FuncAnimation(self.fig, self.animate,
                                       interval=self.monitor.samplingInterval * 1000)#miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

    def show(self):
        self.paint(['mem','cpu', 'swap'])

