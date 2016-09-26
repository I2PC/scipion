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

import sys
import os
from os.path import abspath, dirname

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtMonitor, Monitor
from pyworkflow.em.protocol import ProtImportMovies


class ProtMonitorISPyB(ProtMonitor):
    """ Provide some summary of the basic steps of the Scipion-Box:
    - Import movies
    - Align movies (global and/or local)
    - CTF estimation.
    """
    _label = 'monitor to ISPyB'

    def _defineParams(self, form):
        ProtMonitor._defineParams(self, form)

        group = form.addGroup('Experiment')
        group.addParam('groupid', params.StringParam,
                      label="Group Id",
                      help="Group Id")
        group.addParam('visit', params.StringParam,
                      label="Visit",
                      help="Visit")
        group.addParam('sampleid', params.StringParam,
                      label="Sample Id",
                      help="Sample Id")
        group.addParam('detectorid', params.StringParam,
                      label="Detector Id",
                      help="Detector Id")

        form.addParam('db', params.EnumParam,
                      choices=["production", "devel", "test"],
                      label="Database",
                      help="Select which ISPyB database you want to use.")

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        monitor = MonitorISPyB(self, workingDir=self._getPath(),
                               samplingInterval=self.samplingInterval.get(),
                               monitorTime=100,
                               stdout=True,
                               )
        #monitor.initLoop()
        print "just created the monitor and loop"
        monitor.loop()


class MonitorISPyB(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocol, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocol = protocol

    def step(self):
        print "MonitorISPyB: only one step"

        prot = self.protocol
        db = ISPyBdb(prot.db.get(), prot.groupid.get(), prot.visit.get(),
                     prot.sampleid.get(), prot.detectorid.get())

        try:
            for p in prot.inputProtocols:
                obj = p.get()
                print "protocol: ", obj.getRunName()

                if isinstance(obj, ProtImportMovies):
                    outSet = obj.outputMovies
                    outSet.load()
                    for movie in outSet:
                        print "movieId: ", movie.getObjId()
                        db.put_movie(movie)
                    outSet.close()
        except Exception as ex:
            print "ERROR: ", ex

        return False


class ISPyBdb():
    def __init__(self, db, parentid, visit, sampleid, detectorid):
        self.params =  {}
        self.params['parentid'] = parentid
        self.params['visitid'] = visit
        self.params['sampleid'] = sampleid
        self.params['detectorid'] = detectorid

    def put_movie(self, movie):
        movieFn = movie.getFileName()
        self.params['imgdir'] = dirname(movieFn)
        self.params['imgprefix'] = pwutils.removeBaseExt(movieFn)
        self.params['imgsuffix'] = pwutils.getExt(movieFn)
        self.params['file_template'] = movieFn

        # TODO: Use ispyb-api some day
        cmd = 'module load python/ana; python --version' # em_put_movie.py'
        for k, v in self.params.iteritems():
            cmd += ' --%s %s' % (k, v)
        print cmd
        #os.system(cmd)
