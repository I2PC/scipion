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

import time

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params

from protocol_monitor import ProtMonitor



class ProtMonitorSummary(ProtMonitor):
    """ Provide some summary of the basic steps of the Scipion-Box:
    - Import movies
    - Align movies (global and/or local)
    - CTF estimation.
    """
    _label = 'monitor summary'

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        finished = False
        completedDict = {}
        f = open(self._getPath('log.txt'), 'w')

        def log(msg):
            print >> f, msg
            f.flush()

        log("Inside monitorStep......")

        while not finished:
            log("Iterating over protocols...")
            for protPointer in self.inputProtocols:
                prot = protPointer.get()
                log(" prot.getObjId(): %s" % prot.getObjId())
                for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                    outSet.load()
                    outSet.loadAllProperties()
                    completedDict[(prot.getObjId(), outName)] = outSet.getSize()
                    outSet.close()
            log("=" * 80)
            log(pwutils.prettyTime())
            log(completedDict)
            # Wait some seconds before checking for new data
            time.sleep(self.samplingInterval.get())

        f.close()


