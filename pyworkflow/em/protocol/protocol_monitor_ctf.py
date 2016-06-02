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


import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor
import sqlite3 as lite
import time, sys
from matplotlib import pyplot
import tkMessageBox
from pyworkflow.protocol.constants import STATUS_RUNNING, STATUS_FINISHED
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.plotter import EmPlotter



class ProtMonitorCTF(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'ctf_monitor'

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self.dataBase = 'log.sqlite'
        self.tableName = 'log'
        self.readCTFs = set()

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        #ProtMonitor._defineParams(self, form)
        form.addSection(label='Input')

        form.addParam('inputProtocol', params.PointerParam,
                      label="Input protocols", important=True,
                      pointerClass='ProtCTFFind, XmippProtCTFMicrographs',
                      help="this protocol will be monitorized")
        form.addParam('samplingInterval', params.IntParam,default=60,
                      label="Sampling Interval (sec)",
                      pointerClass='EMProtocol',
                      help="Take one sample each SamplinInteval seconds")
        form.addParam('maxDefocus', params.FloatParam,default=40000,
              label="Raise Alarm if maximum defocus (A) >",
              help="Raise alarm if defocus is greater than given value")
        form.addParam('minDefocus', params.FloatParam,default=1000,
              label="Raise Alarm if minimum defocus (A) <",
              help="Raise alarm if defocus is smaller than given value")
        form.addParam('astigmatism', params.FloatParam,default=0.2,
              label="Raise Alarm if astigmatism >",
              help="Raise alarm if astigmatism is greater than given value")
        form.addParam('interval', params.FloatParam,default=300,
              label="Total Logging time (min)",
              help="Log during this interval")
        ProtMonitor._sendMailParams(self, form)

    #--------------------------- STEPS functions --------------------------------------------
    def monitorStep(self):
        baseFn = self._getPath(self.dataBase)
        conn = lite.connect(baseFn, isolation_level=None)
        cur = conn.cursor()
        self.createTable(conn,self.tableName)
        #TODO: interval should be protocol
        interval = self.interval.get() # logging during these minutes
        sleepSec = self.samplingInterval.get() # wait this seconds between logs
        self.loopCTFCollect(cur,self.tableName,interval,sleepSec)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        return ['Monitor CTF defocus']

    def _updateProtocol(self, prot):
        prot2 = getProtocolFromDb(self.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def loopCTFCollect(self, cur, tableName, interval,sleepSec):
        timeout = time.time() + 60.*interval   # interval minutes from now
        #ignore first meassure because is very unrealiable
        while True:
            prot = self._updateProtocol(self.inputProtocol.get())
            finished = prot.getStatus()!=STATUS_RUNNING
            #create set of processed CTF from CTF protocol
            CTFset=prot.outputCTF.getIdSet()#set([ctf.getObjId() for ctf in prot.outputCTF])
            #find difference
            import sys
            sys.stdout.flush()
            diffSet = CTFset -  self.readCTFs
            setOfCTFs = prot.outputCTF#.get()
            astigmatism = self.astigmatism.get()
            for ctfID in diffSet:
                defocusU = setOfCTFs[ctfID].getDefocusU()
                defocusV = setOfCTFs[ctfID].getDefocusV()
                defocusAngle = setOfCTFs[ctfID].getDefocusAngle()

                if defocusU < defocusV:
                    aux = defocusV
                    defocusV = defocusU
                    defocusU = aux
                    #TODO: check if this is always true
                    defocusAngle = 180. - defocusAngle
                    print("ERROR: defocusU should be greater than defocusV")

                #get CTFs with this ids a fill table
                #do not forget to compute astigmatism
                sql = """INSERT INTO %s(defocusU,defocusV,astigmatism,ratio) VALUES(%f,%f,%f,%f);""" \
                      %(tableName, defocusU, defocusV, defocusAngle, defocusU/defocusV)
                try:
                    cur.execute(sql)
                except Exception as e:
                    print("ERROR: saving one data point (CTF monitor). I continue")

                if  (defocusU/defocusV) > (1. + astigmatism):
                     print("Error Message", "defocusU/defocusV =%f."%defocusU/defocusV)
                     sys.stdout.flush()

                     self.swapAlert = swap.percent
                     if self.doMail:
                          self.sendEMail("scipion system monitor warning", "defocus ratio  =%f."%defocusU/defocusV)

                if  defocusU > self.maxDefocus:
                     print("Error Message", "defocusU =%f."%defocusU)
                     sys.stdout.flush()

                     self.maxDefocus = defocusU

                     if self.doMail:
                         self.sendEMail("scipion system monitor warning", "defocusU  =%f."%defocusU)

                if  defocusV < self.minDefocus:
                     print("Error Message", "defocusV =%f."%defocusV)
                     sys.stdout.flush()

                     self.maxDefocus = defocusV

                     if self.doMail:
                         self.sendEMail("scipion system monitor warning", "defocusV  =%f."%defocusV)

            self.readCTFs.update(diffSet)
            if (time.time() > timeout) or finished:
               break
            time.sleep(sleepSec)

         #self.setStatus(STATUS_FINISHED)

    def _methods(self):
        return []

    def createTable(self,cur,tableName):

        sql = """CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now','localtime')),
                                defocusU FLOAT,
                                defocusV FLOAT,
                                defocus FLOAT,
                                astigmatism FLOAT,
                                ratio FLOAT)
                                """ % tableName
        cur.execute(sql)


from pyworkflow.viewer import ( DESKTOP_TKINTER, WEB_DJANGO, Viewer)
from pyworkflow.protocol.params import (LabelParam, NumericRangeParam,
                                        EnumParam, FloatParam, IntParam)
from matplotlib import animation


class ProtMonitorCTFViewer(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'ctf monitor'
    _targets = [ProtMonitorCTF]

    def _visualize(self, obj, **kwargs):
        return [CtfMonitorPlotter(obj)]


class CtfMonitorPlotter(EmPlotter):
    def __init__(self, protocol):
        EmPlotter.__init__(self, windowTitle="CTF Monitor")
        self.protocol = protocol
        self.y2 = 0.; self.y1 = 100.
        self.win = 250 # number of samples to be ploted
        self.step = 50 # self.win  will be modified in steps of this size

        self.createSubPlot("Use scrool wheel to change view window (win=%d)\n S stops, C continues plotting" % self.win, "Micrographs", "Defocus (A)")
        self.fig = self.getFigure()
        self.ax = self.getLastSubPlot()
        self.ax.margins(0.05)
        self.ax.grid(True)
        self.oldWin = self.win

        self.lines = {}
        self.init = True
        self.stop = False

        baseFn = self.protocol._getPath(self.protocol.dataBase)
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
                           #FuncAnimation adds it as argument
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
	self.ax.axhline(y=self.protocol.minDefocus, c="red", linewidth=0.5, linestyle='dashed', zorder=0)
	self.ax.axhline(y=self.protocol.maxDefocus, c="red", linewidth=0.5, linestyle='dashed', zorder=0)
    def show(self):
        self.paint(['defocusU','defocusV'])


    def paint(self, labels):
        for label in labels:
            self.lines[label],=self.ax.plot([], [], '-',label=label)

        anim = animation.FuncAnimation(self.fig, self.animate, interval = self.samplingInterval.get() * 1000)#miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

    def getData(self):

        cur = self.cur

        def get(name):
            try:
                cur.execute("select %s from %s" % (name, self.tableName))
            except Exception as e:
                print("ERROR readind data (plotter). I continue")
            return [r[0] for r in cur.fetchall()]

        data = {
                'defocusU': get('defocusU'),
                'defocusV': get('defocusV'),
                'astigmatism': get('astigmatism'),
                'ratio': get('ratio'),
                'idValues': get('id')
                }
        #conn.close()
        return data
