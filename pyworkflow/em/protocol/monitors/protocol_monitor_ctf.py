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
from math import isinf

import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor, Monitor
import sqlite3 as lite

from pyworkflow import VERSION_1_1

from pyworkflow.protocol.constants import STATUS_RUNNING
from pyworkflow.protocol import getUpdatedProtocol

PHASE_SHIFT = 'phaseShift'
TIME_STAMP = 'timeStamp'
DEFOCUS_U = 'defocusU'
RESOLUTION = 'resolution'
CTF_LOG_SQLITE = 'ctf_log.sqlite'


class ProtMonitorCTF(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'ctf monitor'
    _lastUpdateVersion = VERSION_1_1

    # -------------------------- DEFINE param functions ---------------------
    def _defineParams(self, form):    
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
        form.addParam('astigmatism', params.FloatParam,default=2000,
              label="Raise Alarm if astigmatism (A) >",
              help="Raise alarm if astigmatism is greater than given value")

        form.addParam('monitorTime', params.FloatParam, default=300,
              label="Total Logging time (min)",
              help="Log during this interval")

        ProtMonitor._sendMailParams(self, form)

    # -------------------------- STEPS functions ------------------------------
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

    # -------------------------- INFO functions -------------------------------
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

            # Defocus angle
            defocusAngle = ctf.getDefocusAngle()
            if defocusAngle > 360 or defocusAngle< -360:
                defocusAngle = 0

            # Astigmatism
            astig = abs(defocusU - defocusV)

            # Resolution
            resolution = ctf.getResolution()
            if isinf(resolution):
                resolution = 0.
            
            # Fit quality
            fitQuality = ctf.getFitQuality()
            if fitQuality is None or isinf(fitQuality): 
                fitQuality = 0.

            # PhaseShift
            phaseShift = ctf.getPhaseShift() if ctf.hasPhaseShift() else 0.

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

            ctfCreationTime = ctf.getObjCreation()

            # get CTFs with this ids a fill table
            # do not forget to compute astigmatism
            sql = """INSERT INTO %s(timestamp, ctfID, defocusU,defocusV,astigmatism,ratio, resolution, fitQuality, phaseShift, micPath,psdPath,shiftPlotPath )
                     VALUES("%s",%d,%f,%f,%f,%f,%f,%f,%f,"%s","%s","%s");""" % (self._tableName, ctfCreationTime, ctfID, defocusU,
                     defocusV, astig, defocusU / defocusV, resolution, fitQuality, phaseShift,  micPath, psdPath, shiftPlotPath)
            try:
                self.cur.execute(sql)
            except Exception as e:
                print("ERROR: saving one data point (CTF monitor). I continue")
                print(e)
                print(sql)

            if abs(defocusU - defocusV) > astigmatism:
                self.warning("Astigmatism (defocusU - defocusV)  = %f."
                             % abs(defocusU - defocusV))

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
                                ctfID INTEGER,
                                defocusU FLOAT,
                                defocusV FLOAT,
                                defocus FLOAT,
                                astigmatism FLOAT,
                                ratio FLOAT,
                                resolution FLOAT,
				                fitQuality FLOAT,
				                phaseShift FLOAT,
                                micPath STRING,
                                psdPath STRING,
                                shiftPlotPath STRING)
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
            DEFOCUS_U: get('defocusU'),
            'defocusV': get('defocusV'),
            'astigmatism': get('astigmatism'),
            'ratio': get('ratio'),
            'idValues': get('ctfID'),
            RESOLUTION: get('resolution'),
            'fitQuality': get('fitQuality'),
            PHASE_SHIFT: get('phaseShift'),
            'imgMicPath': get('micPath'),
            'imgPsdPath': get('psdPath'),
            'imgShiftPath': get('shiftPlotPath'),
            TIME_STAMP: get("strftime('%s', timestamp) * 1000")
         }
        # conn.close()
        return data


