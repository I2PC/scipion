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

import pyworkflow as pw
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtMonitor, Monitor, PrintNotifier
from pyworkflow.em.protocol import ProtImportMovies

from ispyb_proxy import ISPyBProxy


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
                               monitorTime=100)

        monitor.addNotifier(PrintNotifier())
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
        self.info("MonitorISPyB: only one step")

        prot = self.protocol

        proxy = ISPyBProxy(["prod", "dev", "test"][prot.db.get()],
                           experimentParams={'parentid': prot.groupid.get(),
                                             'visit': prot.visit.get(),
                                             'sampleid': prot.sampleid.get(),
                                             'detectorid': prot.detectorid.get()
                                             }
                           )
        try:
            for p in prot.inputProtocols:
                obj = p.get()
                self.info("protocol: %s" % obj.getRunName())

                if isinstance(obj, ProtImportMovies):
                    outSet = obj.outputMovies
                    outSet.load()
                    for movie in outSet:
                        self.info("movieId: %s" % movie.getObjId())
                        self.put_movie(proxy, movie)
                    outSet.close()
        except Exception as ex:
            print "ERROR: ", ex

        self.info("Closing proxy")
        proxy.close()

        return False

    def put_movie(self, proxy, movie):
        movieFn = movie.getFileName()
        proxy.sendMovieParams({'imgdir': dirname(movieFn),
                               'imgprefix': pwutils.removeBaseExt(movieFn),
                               'imgsuffix': pwutils.getExt(movieFn),
                               'file_template': movieFn
                               })


class FileNotifier():
    def __init__(self, filename):
        self.f = open(filename, 'w')

    def notify(self, title, message):
        print >> self.f, title, message
        self.f.flush()