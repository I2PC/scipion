# **************************************************************************
# *
# * Authors:     R. Marabini (roberto@cnb.csic.es)
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

from collections import OrderedDict

import pyworkflow.protocol.params as params
from protocol import EMProtocol
import time
import xmipp
import random
from pyworkflow.em.data import SetOfMicrographs, Micrograph, Acquisition
from pyworkflow.protocol.constants import STEPS_PARALLEL
from os.path import basename, exists


class ProtCreateStreamData(EMProtocol):
    """ create  setofXXXX in streaming mode.
    """
    _label="create stream data"
    _singleImageFn = "singleImage.xmp"
    _magnification = 500000
    _voltage = 200
    _sphericalAberration = 2.0
    _amplitudeContrast = 0.07

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.xmippMic = OrderedDict()
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        #TODO: add movie creation in streaming
        form.addParam('setof', params.EnumParam, default=0,
                      choices=['setOfMicrographs'],
                      label='set Of',
                      help='create set of')
        form.addParam('xDim', params.IntParam,default=1024,
                      label="xdim",
                      help="X dim ")
        form.addParam('yDim', params.IntParam,default=1024,
                      label="ydim",
                      help="Y dim ")
        form.addParam('nDim', params.IntParam,default=10,
                      label="ndim",
                      help="setofXX size")
        form.addParam('samplingRate', params.FloatParam,default=4,
                      label="samplingRate",
                      help="Sampling rate")
        form.addParam('creationInteval', params.IntParam,default=60,
              label="Create Object each (sec)",
              pointerClass='EMProtocol',
              help="create one object each creationInteval seconds")
        form.addParam('delay', params.IntParam,default=0,
                      label="delay (sec)",
                      help="wait this seconds before creating stram data")

        form.addParallelSection(threads=2, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        deps = []

        if self.delay != 0:
            delayId = self._insertFunctionStep('delayStep')
            deps.append(delayId)

        for mic in range (1, self.nDim.get() +1):
            self._insertFunctionStep('createStep', mic, prerequisites=deps)

    #--------------------------- STEPS functions --------------------------------------------
    def delayStep(self):
        time.sleep(self.delay)

    def _checkNewMics(self, micSet):
        """ Check for already computed CTF and update the output set. """
        micDict = {}
        newMic = False
        for mic in micSet:
            micDict[mic.getFileName()] = True

        if micDict:
            if micSet.getSize():
                micSet.enableAppend()
                micSet.loadAllProperties()
        else:
            micSet.setStreamState(micSet.STREAM_OPEN)
            acquisition = Acquisition()
            acquisition.setMagnification(self._magnification)
            acquisition.setVoltage(self._voltage)
            acquisition.setSphericalAberration(self._sphericalAberration)
            acquisition.setAmplitudeContrast(self._amplitudeContrast)
            micSet.setAcquisition(acquisition)
            micSet.setSamplingRate(self.samplingRate.get())


        mic = Micrograph()

        counter = 0
        for k,v in self.xmippMic.iteritems():
            counter += 1
            if (k not in micDict):
                mic.setFileName(k)
                mic.setMicName(basename(k))
                mic.setObjId(counter)
                micSet.append(mic)
                newMic= True

        return micSet, newMic#why a dictionary, a boolean may be enought

    def _updateOutput(self, micSet):
        #firstTime = not self.hasAttribute('outputMicrographs')
        self._defineOutputs(outputMicrographs=micSet)


    def _stepsCheck(self):
        setOfMicFn = self._getPath('micrographs.sqlite')
        micSet = SetOfMicrographs(filename=setOfMicFn)
        micSet, newMic = self._checkNewMics(micSet)

        #check if end ....
        endMics = micSet.getSize() == self.nDim.get()

        if newMic:
            if endMics:
                micSet.setStreamState(micSet.STREAM_CLOSED)
            self._updateOutput(micSet)
        micSet.close()

    def createStep(self, mic):
        from pyworkflow.em.packages.xmipp3 import getEnviron

        #create image
        img = xmipp.Image()
        img.setDataType(xmipp.DT_FLOAT)
        img.resize(self.xDim, self.yDim)
        img.initRandom(0., 1., xmipp.XMIPP_RND_UNIFORM)
        baseFn = self._getExtraPath(self._singleImageFn)
        img.write(baseFn)

        md1 = xmipp.MetaData()
        md1.setColumnFormat(False)
        idctf = md1.addObject()

        baseFnCtf = self._getTmpPath("ctf_%d.param"%mic)
        baseFnImageCTF = self._getExtraPath("imageCTF_%d.xmp"%mic)

        md1.setValue(xmipp.MDL_CTF_SAMPLING_RATE, 1., idctf)
        md1.setValue(xmipp.MDL_CTF_VOLTAGE, 200., idctf);
        defocus = 20000 + 10000 * random.random()
        udefocus = defocus + 1000 * random.random()
        vdefocus = defocus + 1000 * random.random()
        if udefocus < vdefocus:
            aux = vdefocus
            vdefocus = udefocus
            udefocus = aux
        md1.setValue(xmipp.MDL_CTF_DEFOCUSU, udefocus, idctf);
        md1.setValue(xmipp.MDL_CTF_DEFOCUSV, vdefocus, idctf);
        md1.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, 180.0 * random.random(), idctf);
        md1.setValue(xmipp.MDL_CTF_CS, 2., idctf);
        md1.setValue(xmipp.MDL_CTF_Q0, 0.07, idctf);
        md1.setValue(xmipp.MDL_CTF_K, 1., idctf);

        md1.write(baseFnCtf)

        #apply ctf
        args  = " -i %s"%baseFn
        args += " -o %s"%baseFnImageCTF
        args += " -f ctf %s"%baseFnCtf
        args += " --sampling %f"%self.samplingRate
        self.runJob("xmipp_transform_filter", args, env=getEnviron())
        self.xmippMic[baseFnImageCTF] = True
        time.sleep(self.creationInteval.get())

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []


