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
from pyworkflow.em.plotter import EmPlotter
from tkMessageBox import showerror

from pynvml import nvmlInit, nvmlDeviceGetCount, nvmlDeviceGetHandleByIndex,\
    nvmlDeviceGetName, nvmlDeviceGetMemoryInfo, nvmlDeviceGetUtilizationRates,\
    NVMLError, nvmlDeviceGetTemperature, NVML_TEMPERATURE_GPU,\
    nvmlDeviceGetComputeRunningProcesses


SYSTEM_LOG_SQLITE = 'system_log.sqlite'


def errorWindow(tkParent, msg):
    try:
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)

def initGPU():
    nvmlInit()

class ProtMonitorSystem(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'system_monitor'
    _version = VERSION_1_1

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
        group = form.addGroup('GPU')
        group.addParam('doGpu', params.BooleanParam, default=False,
                       label="Check GPU",
                       help="Set to true if you want to monitor the GPU")
        group.addParam('gpusToUse', params.StringParam, default='0',
                          label='Which GPUs to use:', condition='doGpu',
                          help='Providing a list of which GPUs '
                               '(0,1,2,3, etc). Default is monitor GPU 0 only')


    #--------------------------- STEPS functions ------------------------------

    def monitorStep(self):
        self.createMonitor().loop()

    def createMonitor(self):
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)

        sysMonitor = MonitorSystem(protocols,
                                   workingDir=self.workingDir.get(),
                                   samplingInterval=self.samplingInterval.get(),
                                   monitorTime=self.monitorTime.get(),
                                   email=self.createEmailNotifier(),
                                   stdout=True,
                                   cpuAlert=self.cpuAlert.get(),
                                   memAlert=self.memAlert.get(),
                                   swapAlert=self.swapAlert.get(),
                                   doGpu=self.doGpu.get(),
                                   gpusToUse = self.gpusToUse.get(),
                                    )
        return sysMonitor

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        summary = []
        summary.append("GPU running Processes:")
        initGPU()
        try:
            gpusToUse = [int(n) for n in (self.gpusToUse.get()).split()]
            for i in gpusToUse:
                handle = nvmlDeviceGetHandleByIndex(i)
                cps = nvmlDeviceGetComputeRunningProcesses(handle)
                for ps in cps:
                    #p_tags['pid'] = ps.pid
                    msg = " %d) "%i + psutil.Process(ps.pid).name()
                    msg += " (mem =%.2f MB)"%(float(ps.usedGpuMemory)/1048576.)
                    summary.append(msg)
        except NVMLError as err:
                summary.append(str(err))

        return summary

    def _methods(self):
        return []


class MonitorSystem(Monitor):
    """ This will will be monitoring a System  protocol.
    It will internally handle a database to store produced
    system values.
    """
    def __init__(self, protocols, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocols = protocols
        self.cpuAlert = kwargs['cpuAlert']
        self.memAlert = kwargs['memAlert']
        self.swapAlert = kwargs['swapAlert']
        self._dataBase = kwargs.get('dbName', SYSTEM_LOG_SQLITE)
        self._tableName = kwargs.get('tableName', 'log')
        self.doGpu = kwargs['doGpu']
        self.labelList=["cpu","mem","swap"]
        if self.doGpu:
            self.gpuLabelList=[]
            #get Gpus to monitor
            self.gpusToUse = [int(n) for n in (kwargs['gpusToUse']).split()]
            for i in self.gpusToUse:
                self.gpuLabelList.append("gpuMem_%d"%i)
                self.gpuLabelList.append("gpuUse_%d"%i)
                self.gpuLabelList.append("gpuTem_%d"%i)
            #init GPUs
            nvmlInit()
            self.labelList += self.gpuLabelList
        else:
            self.gpusToUse = None
        self.conn = lite.connect(os.path.join(self.workingDir, self._dataBase),
                                 isolation_level=None)
        self.cur = self.conn.cursor()

    def _getUpdatedProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def warning(self, msg):
        self.notify("Scipion System Monitor WARNING", msg)

    def initLoop(self):
        self._createTable()
        psutil.cpu_percent(True)
        psutil.virtual_memory()

    def step(self):
        valuesDict={}
        valuesDict['table'] = self._tableName
        cpu  = valuesDict['cpu'] = psutil.cpu_percent(interval=0)
        mem  = valuesDict['mem'] = psutil.virtual_memory().percent
        swap = valuesDict['swap'] = psutil.swap_memory().percent
        #some code examples: https://github.com/ngi644/datadog_nvml/blob/master/nvml.py
        if self.doGpu:
            for i in self.gpusToUse:
                try:
                    handle = nvmlDeviceGetHandleByIndex(i)
                    memInfo = nvmlDeviceGetMemoryInfo(handle)
                    valuesDict["gpuMem_%d"%i] = float(memInfo.used)*100./float(memInfo.total)
                    util = nvmlDeviceGetUtilizationRates(handle)
                    valuesDict["gpuUse_%d"%i] = util.gpu
                    temp = nvmlDeviceGetTemperature(handle,NVML_TEMPERATURE_GPU)
                    valuesDict["gpuTem_%d"%i] = temp
                except NVMLError as err:
                    handle = nvmlDeviceGetHandleByIndex(i)
                    msg = "Device %d -> %s not suported\n Remove device %d from FORM"%\
                          (i,nvmlDeviceGetName(handle),i)
                    errorWindow(None, msg)

        if self.cpuAlert < 100 and cpu > self.cpuAlert:
            self.warning("CPU allocation =%f." % cpu)
            self.cpuAlert = cpu

        if self.memAlert < 100 and mem.percent > self.memAlert:
            self.warning("Memory allocation =%f." % mem)
            self.memAlert = mem

        if self.swapAlert < 100 and swap.percent > self.swapAlert:
            self.warning("SWAP allocation =%f." % swap)
            self.swapAlert = swap

        sqlCommand = "INSERT INTO %(table)s ("
        for label in self.labelList:
            sqlCommand += "%s, "%label
        #remove last comma
        sqlCommand = sqlCommand[:-2]
        sqlCommand += ") VALUES("
        for label in self.labelList:
            sqlCommand += "%"+"(%s)f, "%label
        #remove last comma
        sqlCommand = sqlCommand[:-2]
        sqlCommand += ");"
        print ("sqlCommand", sqlCommand,valuesDict)
        sql = sqlCommand%valuesDict

        try:
            self.cur.execute(sql)
        except Exception as e:
            print("ERROR: saving one data point (monitor). I continue")

        # Return finished = True if all protocols have finished
        return all(self._getUpdatedProtocol(prot).getStatus() != STATUS_RUNNING
                   for prot in self.protocols)

    def _createTable(self):
        sqlCommand = """CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now','localtime')),
                                """ % self._tableName
        for label in self.labelList:
            sqlCommand += "%s FLOAT,\n"%label
        #remove last comma and new line
        sqlCommand = sqlCommand[:-3]
        sqlCommand +=")"
        print("sqlCommand",sqlCommand)
        self.cur.execute(sqlCommand)

    def getLabels(self):
        return self.labelList

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
                'idValues': idValues}
        for label in self.labelList:
            data[label] = get(label)

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
        self.paint(self.monitor.getLabels())

