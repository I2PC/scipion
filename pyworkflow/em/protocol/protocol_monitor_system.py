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
# *  e-mail address 'roberto@cnb.csic.es'
# *
# **************************************************************************


import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor
import sqlite3 as lite
import psutil
import time, sys
from matplotlib import pyplot
from pyworkflow.protocol.constants import STATUS_RUNNING, STATUS_FINISHED
from pyworkflow.protocol import getProtocolFromDb
import sys
from pyworkflow.em.plotter import EmPlotter

class ProtMonitorSystem(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'system_monitor'

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self.dataBase = 'log.sqlite'
        self.tableName = 'log'

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        ProtMonitor._defineParams(self, form)
#        form.addParam('doDisk', params.BooleanParam, default=False,
#              label='Disk IO',
#              help='By default only CPU and Mem usage are stored. '
#                   'Set to yes to store disk related info.' )
        form.addParam('cpuAlert', params.FloatParam,default=101,
              label="Raise Alarm if CPU > XX%",
              help="Raise alarm if memory allocated is greater than given percentage")
        form.addParam('memAlert', params.FloatParam,default=101,
              label="Raise Alarm if Memory > XX%",
              help="Raise alarm if cpu allocated is greater than given percentage")
        form.addParam('swapAlert', params.FloatParam,default=101,
              label="Raise Alarm if Swap > XX%",
              help="Raise alarm if swap allocated is greater than given percentage")
        form.addParam('interval', params.FloatParam,default=300,
              label="Total Logging time (min)",
              help="Log during this interval")
        ProtMonitor._sendMailParams(self, form)

    #--------------------------- STEPS functions --------------------------------------------
    def monitorStep(self):
        baseFn = self._getExtraPath(self.dataBase)
        conn = lite.connect(baseFn, isolation_level=None)
        cur = conn.cursor()
        self.createTable(conn,self.tableName)
        #TODO: interval should be protocol
        interval = self.interval.get() # logging during these minutes
        sleepSec = self.samplingInterval.get() # wait this seconds between logs
        doDisk = False#self.doDisk
        self.loopPsUtils(cur,self.tableName,interval,sleepSec,doDisk)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        return ['Stores CPU, memory and swap ussage in percentage']

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(self.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def loopPsUtils(self, cur, tableName, interval,sleepSec,doDisk):
         timeout = time.time() + 60.*interval   # interval minutes from now
         #ignore first meassure because is very unrealiable
         psutil.cpu_percent(True)
         psutil.virtual_memory()
         disks_before = psutil.disk_io_counters(perdisk=False)

         while True:
             finished = True
             for protPointer in self.inputProtocols:
                 prot = self._updateProtocol(protPointer.get())
                 finished = finished and (prot.getStatus()!=STATUS_RUNNING)
             if (time.time() > timeout) or finished:
                break

             cpu = psutil.cpu_percent(interval=0)
             mem = psutil.virtual_memory()
             swap = psutil.swap_memory()


             if self.cpuAlert < 100 and cpu > self.cpuAlert :
                 print("Error Message\n CPU allocation =%f."%cpu)
                 sys.stdout.flush()

                 self.cpuAlert = cpu
                 if self.doMail:
                     self.sendEMail("scipion system monitor warning", "CPU allocation =%f."%cpu.percent)


             if self.memAlert < 100 and mem.percent > self.memAlert :
                 print("Error Message", "Memory allocation =%f."%mem.percent)
                 sys.stdout.flush()

                 self.memAlert = mem.percent
                 if self.doMail:
                     self.sendEMail("scipion system monitor warning", "Memory allocation =%f."%mem.percent)

             if self.swapAlert < 100 and swap.percent > self.swapAlert :
                 print("Error Message", "SWAP allocation =%f."%swap.percent)
                 sys.stdout.flush()

                 self.swapAlert = swap.percent
                 if self.doMail:
                     self.sendEMail("scipion system monitor warning", "SWAP allocation =%f."%swap.percent)

             if doDisk:
                disks_after = psutil.disk_io_counters(perdisk=False)
                disks_read_per_sec   = (disks_after.read_bytes  - disks_before.read_bytes)/(sleepSec * 1024.*1024.)
                disks_write_per_sec  = (disks_after.write_bytes - disks_before.write_bytes)/(sleepSec * 1024.*1024.)
                disks_read_time_sec  = (disks_after.read_time   - disks_before.read_time)/(sleepSec*1000.)
                disks_write_time_sec = (disks_after.write_time  - disks_before.write_time)/(sleepSec*1000.)

                disks_before = disks_after

                sql = """INSERT INTO %s(mem,cpu,swap,
                                        disks_read_per_sec,
                                        disks_write_per_sec,
                                        disks_read_time_sec,
                                        disks_write_time_sec) VALUES(%f,%f,%f,%f,%f,%f,%f);""" % (tableName, mem.percent, cpu, swap.percent,disks_read_per_sec,disks_write_per_sec,disks_read_time_sec,disks_write_time_sec )
             else:
                sql = """INSERT INTO %s(mem,cpu,swap) VALUES(%f,%f,%f);""" % (tableName, mem.percent, cpu, swap.percent)
             try:
                cur.execute(sql)
             except Exception as e:
                 print("ERROR: saving one data point (monitor). I continue")
             time.sleep(sleepSec)
         self.setStatus(STATUS_FINISHED)

    def _methods(self):
        return []

    def createTable(self,cur,tableName):

        #cur.execute("DROP TABLE IF EXISTS %s"%tableName)
        sql = """CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now','localtime')),
                                mem FLOAT,
                                cpu FLOAT,
                                disks_read_per_sec FLOAT,
                                disks_write_per_sec FLOAT,
                                disks_read_time_sec FLOAT,
                                disks_write_time_sec FLOAT,
                                openfile INT,
                                swap FLOAT)""" % tableName
        cur.execute(sql)


from pyworkflow.viewer import ( DESKTOP_TKINTER, WEB_DJANGO, Viewer)
from matplotlib import animation

class ProtMonitorSystemViewer(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'system monitor'
    _targets = [ProtMonitorSystem]

    def _visualize(self, obj, **kwargs):
        return [systemMonitorPlotter(obj)]

class systemMonitorPlotter(EmPlotter):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'System Monitor'
    _targets = [ProtMonitorSystem]

    #add param time refresh
    #boolean time refresh

    def __init__(self, protocol):
        EmPlotter.__init__(self, windowTitle="system Monitor")
        self.protocol = protocol
        self.y2 = 0.; self.y1 = 100.
        self.win = 250 # number of samples to be ploted
        self.step = 50 # self.win  will be modified in steps of this size
        self.createSubPlot("Use scrool wheel to change view window (win=%d)\n S stops, C continues plotting" % self.win, "time (hours)", "percentage (or Mb for IO)")
        self.fig = self.getFigure()
        self.ax = self.getLastSubPlot()
        self.ax.margins(0.05)
        self.ax.grid(True)
        self.oldWin = self.win

        self.lines = {}
        self.init = True
        self.stop = False

        baseFn = self.protocol._getExtraPath(self.protocol.dataBase)
        self.tableName = self.protocol.tableName
        self.samplingInterval = self.protocol.samplingInterval
        self.conn  = lite.connect(baseFn, isolation_level=None)
        self.cur   = self.conn.cursor()

    def onscroll(self, event):

        if event.button == 'up':
            self.win += self.step
        else:
            self.win -= self.step
            if self.win < self.step:
               self.win = self.step

        if self.oldWin != self.win:
            self.ax.set_title('use scroll wheel to change view window (win=%d)\n S stops, C continues plotting'%self.win)
            self.oldWin= self.win
        self.animate()
        EmPlotter.show(self)

    def press(self,event):

        sys.stdout.flush()
        if event.key == 'S':
            self.stop = True
            self.ax.set_title('Plot is Stopped. Press C to continue plotting')
        elif event.key == 'C':
	    self.ax.set_title('use scroll wheel to change view window (win=%d)\n S stops, C continues plotting'%self.win)
            self.stop = False
        self.animate()
        EmPlotter.show(self)

    def has_been_closed(self,ax):
        fig = ax.figure.canvas.manager
        active_fig_managers = pyplot._pylab_helpers.Gcf.figs.values()
        return fig not in active_fig_managers

    def animate(self,i=0): #do NOT remove i
                                         #FuncAnimation add it as argument

        if self.stop:
            return

        data = self.getData()
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
        self.doDisk = False

#    def _disk(self, e=None):
#        self.initAnimate('disks_read_per_sec')
#        self.initAnimate('disks_write_per_sec')
#        self.initAnimate('disks_read_time_sec')
#        self.initAnimate('disks_write_time_sec')
    def show(self):
        self.paint(['cpu','mem', 'swap'])
        #EmPlotter.show(self)
        #if EmPlotter.show() is here instead of in paint it does not work



    def paint(self,labels):
        for label in labels:
            self.lines[label],=self.ax.plot([], [], '-',label=label)

        anim = animation.FuncAnimation(self.fig, self.animate, interval=self.samplingInterval.get() * 1000)#miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

    def show(self):
        self.paint(['mem','cpu', 'swap'])
        #EmPlotter.show(self)

    def getData(self):

        cur = self.cur
        #I guess this can be done in a single call
        #I am storing the first meassurement
        cur.execute("select julianday(timestamp)  from %s where id=1"%self.tableName )
        initTime = cur.fetchone()[0]
        cur.execute("select timestamp  from %s where id=1"%self.tableName )
        initTimeTitle = cur.fetchone()[0]
        cur.execute("select (julianday(timestamp) - %f)*24  from %s"%(initTime,self.tableName) )
        idValues = [r[0] for r in cur.fetchall()]


        def get(name):
            try:
                cur.execute("select %s from %s" % (name, self.tableName))
            except Exception as e:
                print("ERROR readind data (plotter). I continue")
            return [r[0] for r in cur.fetchall()]

        data = {'initTime': initTime,
                'initTimeTitle': initTimeTitle,
                'idValues': idValues,
                'cpu': get('cpu'),
                'mem': get('mem'),
                'swap': get('swap'),
                'disks_read_per_sec': get('disks_read_per_sec'),
                'disks_write_per_sec': get('disks_write_per_sec'),
                'disks_read_time_sec': get('disks_read_time_sec'),
                'disks_write_time_sec': get('disks_write_time_sec'),
                }
        #conn.close()
        return data
