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

from os.path import abspath, dirname
from collections import OrderedDict

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtMonitor, Monitor, PrintNotifier
from pyworkflow.em.protocol import ProtImportMovies, ProtAlignMovies, ProtCTFMicrographs

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
        group.addParam('visit', params.StringParam,
                      label="Visit",
                      help="Visit")

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
        self.allIds = OrderedDict()
        self.numberOfFrames = None

    def step(self):
        self.info("MonitorISPyB: only one step")

        prot = self.protocol

        proxy = ISPyBProxy(["prod", "dev", "test"][prot.db.get()],
                           experimentParams={'visit': prot.visit.get()})


        runs = [p.get() for p in self.protocol.inputProtocols]
        g = self.protocol.getProject().getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        allParams = OrderedDict()

        for n in nodes:
            prot = n.run
            self.info("protocol: %s" % prot.getRunName())

            if isinstance(prot, ProtImportMovies):
                self.create_movie_params(prot, allParams)
            elif isinstance(prot, ProtAlignMovies) and hasattr(prot, 'outputMicrographs'):
                self.update_align_params(prot, allParams)
            elif isinstance(prot, ProtCTFMicrographs):
                self.update_ctf_params(prot, allParams)

        for itemId, params in allParams.iteritems():
            ispybId = proxy.sendMovieParams(params)
            self.allIds[itemId] = ispybId

        self.info("Closing proxy")
        proxy.close()

        return False

    def iter_updated_set(self, objSet):
        objSet.load()
        objSet.loadAllProperties()
        for obj in objSet:
            yield obj
        objSet.close()

    def create_movie_params(self, prot, allParams):

        for movie in self.iter_updated_set(prot.outputMovies):
            if self.numberOfFrames is None:
                self.numberOfFrames = movie.getNumberOfFrames()
            movieFn = movie.getFileName()
            movieId = movie.getObjId()

            allParams[movieId] = {
                'id': self.allIds.get(movieId, None),
                'imgdir': dirname(movieFn),
                'imgprefix': pwutils.removeBaseExt(movieFn),
                'imgsuffix': pwutils.getExt(movieFn),
                'file_template': movieFn,
                'n_images': self.numberOfFrames
             }

    def update_align_params(self, prot, allParams):
        for mic in self.iter_updated_set(prot.outputMicrographs):
            allParams[mic.getObjId()].update({
                'comments': 'aligned'
            })

    def update_ctf_params(self, prot, allParams):
        for ctf in self.iter_updated_set(prot.outputCTF):
            allParams[ctf.getObjId()].update({
            'min_defocus': ctf.getDefocusU(),
            'max_defocus': ctf.getDefocusV(),
            'amount_astigmatism': ctf.getDefocusRatio()
            })


class FileNotifier():
    def __init__(self, filename):
        self.f = open(filename, 'w')

    def notify(self, title, message):
        print >> self.f, title, message
        self.f.flush()