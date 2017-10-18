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
import sys

import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor, Monitor, EmailNotifier
from report_html import PSD_PATH
import sqlite3 as lite
import time, sys

from pyworkflow import VERSION_1_1
from pyworkflow.gui.plotter import plt
import tkMessageBox
from pyworkflow.protocol.constants import STATUS_RUNNING
from pyworkflow.protocol import getUpdatedProtocol

from pyworkflow.em.plotter import EmPlotter

CTF_LOG_SQLITE = 'ctf_log.sqlite'


class ProtMonitorCTF(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'ctf monitor'
    _lastUpdateVersion = VERSION_1_1
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

        form.addParam('monitorTime', params.FloatParam, default=300,
              label="Total Logging time (min)",
              help="Log during this interval")

        ProtMonitor._sendMailParams(self, form)

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):

        self.createMonitor().loop()

    def createMonitor(self):

        ctfProt = self.inputProtocol.get()
        ctfProt.setProject(self.getProject())

        ctfMonitor = MonitorCTF(ctfProt,
                                workingDir=self.workingDir.get(),
                                samplingInterval=self.samplingInterval.get(),
                                monitorTime=self.monitorTime.get(),
                                email=self.createEmailNotifier(),
                                stdout=True,
                                minDefocus=self.minDefocus.get(),
                                maxDefocus=self.maxDefocus.get(),
                                astigmatism=self.astigmatism.get())
        return ctfMonitor

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        return ['Monitor CTF defocus']


class MonitorCTF(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocol, **kwargs):
        Monitor.__init__(self, **kwargs)

        # The CTF protocol to monitor
        self.protocol = protocol

        self.maxDefocus = kwargs['maxDefocus']
        self.minDefocus = kwargs['minDefocus']
        self.astigmatism = kwargs['astigmatism']
        self._dataBase = kwargs.get('dbName', CTF_LOG_SQLITE)
        self._tableName = kwargs.get('tableName', 'log')
        self.readCTFs = set()

        self.conn = lite.connect(os.path.join(self.workingDir, self._dataBase),
                                 isolation_level=None)
        self.cur = self.conn.cursor()

    def warning(self, msg):
        self.notify("Scipion CTF Monitor WARNING", msg)

    def initLoop(self):
        self._createTable()

    def step(self):
        prot = getUpdatedProtocol(self.protocol)
        # Create set of processed CTF from CTF protocol
        if hasattr(prot, 'outputCTF'):
            CTFset = prot.outputCTF.getIdSet()
        else:
            return False
        # find difference
        sys.stdout.flush()
        diffSet = CTFset - self.readCTFs
        setOfCTFs = prot.outputCTF
        astigmatism = self.astigmatism

        for ctfID in diffSet:
            ctf = setOfCTFs[ctfID]
            defocusU = ctf.getDefocusU()
            defocusV = ctf.getDefocusV()
            defocusAngle = ctf.getDefocusAngle()
            psdPath = os.path.abspath(ctf.getPsdFile())
            micPath = os.path.abspath(ctf.getMicrograph().getFileName())
            shiftPlot = (getattr(ctf.getMicrograph(), 'plotCart', None)
                         or getattr(ctf.getMicrograph(), 'plotGlobal', None))
            if shiftPlot is not None:
                shiftPlotPath = os.path.abspath(shiftPlot.getFileName())
            else:
                shiftPlotPath = ""

            if defocusU < defocusV:
                aux = defocusV
                defocusV = defocusU
                defocusU = aux
                # TODO: check if this is always true
                defocusAngle = 180. - defocusAngle
                print("ERROR: defocusU should be greater than defocusV")

            # get CTFs with this ids a fill table
            # do not forget to compute astigmatism
            sql = """INSERT INTO %s(defocusU,defocusV,astigmatism,ratio,psdPath)
                     VALUES(%f,%f,%f,%f,"%s");""" % (self._tableName, defocusU,
                     defocusV, defocusAngle, defocusU / defocusV, psdPath)
            try:
                self.cur.execute(sql)
            except Exception as e:
                print("ERROR: saving one data point (CTF monitor). I continue")
                print e

            if (defocusU / defocusV) > (1. + astigmatism):
                self.warning("Defocus ratio (defocusU / defocusV)  = %f."
                             % (defocusU / defocusV))

            if defocusU > self.maxDefocus:
                self.warning("DefocusU (%f) is larger than defocus "
                             "maximum (%f)" % (defocusU, self.maxDefocus))
                self.maxDefocus = defocusU

            if defocusV < self.minDefocus:
                self.warning("DefocusV (%f) is smaller than defocus "
                             "minumum (%f)" % (defocusV, self.maxDefocus))
                self.minDefocus = defocusV

        self.readCTFs.update(diffSet)
        # Finish when protocol is not longer running
        return prot.getStatus() != STATUS_RUNNING

    def _createTable(self):
        self.cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now','localtime')),
                                defocusU FLOAT,
                                defocusV FLOAT,
                                defocus FLOAT,
                                astigmatism FLOAT,
                                ratio FLOAT, 
                                psdPath STRING)
                                """ % self._tableName)

    def getData(self):
        def get(name):
            try:
                self.cur.execute("select %s from %s" % (name, self._tableName))
            except Exception as e:
                print("MonitorCTF, ERROR reading data from db: %s" %
                      os.path.join(self.workingDir, self._dataBase))
            return [r[0] for r in self.cur.fetchall()]

        data = {
            'defocusU': get('defocusU'),
            'defocusV': get('defocusV'),
            'astigmatism': get('astigmatism'),
            'ratio': get('ratio'),
            'idValues': get('id'),
            PSD_PATH: get('psdPath'),
        }
        # conn.close()
        return data


from pyworkflow.viewer import ( DESKTOP_TKINTER, Viewer)
from pyworkflow.protocol.params import (LabelParam, NumericRangeParam,
                                        EnumParam, FloatParam, IntParam)
from matplotlib import animation


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
        self.y2 = 0.; self.y1 = 100.
        self.win = 250 # number of samples to be ploted
        self.step = 50 # self.win  will be modified in steps of this size

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
                                       interval=self.monitor.samplingInterval*1000)#miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

