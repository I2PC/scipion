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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************


import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from protocol import EMProtocol



class ProtMonitor(EMProtocol):
    """ This is the base class for implementing 'Monitors', a special type
    of protocols intended to be used in the context of the Scipion-Box,
    where several protocols will be run in streaming and the monitor will
    be 'observing' their progress.
    """
    _label = 'monitor'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputProtocols', params.MultiPointerParam,
                      label="Input protocols", important=True,
                      pointerClass='EMProtocol',
                      help="this protocol/s will be monitorized")

        form.addParam('samplingInterval', params.IntParam,default=60,
                      label="Sampling Interval (sec)",
                      pointerClass='EMProtocol',
                      help="Take one sample each SAmplinInteval seconds")

    #--------------------------- INSERT steps functions --------------------------------------------   
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions --------------------------------------------
    def monitorStep(self):
        finished = False
        completedDict = {}

        while not finished:
            for protPointer in self.inputProtocols:
                prot = protPointer.get()
                for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                    outSet.load()
                    completedDict[(prot.getObjId(), outName)] = outSet.getSize()
                    outSet.close()
            print "=" * 80
            print pwutils.prettyTime()
            print completedDict

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []


