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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from collections import OrderedDict

import pyworkflow.protocol.params as params
from protocol import EMProtocol
import time
import xmipp
import random

from pyworkflow import VERSION_1_2
from pyworkflow.em.data import SetOfMicrographs, Micrograph, Acquisition, Movie, SetOfMovies
from pyworkflow.protocol.constants import STEPS_PARALLEL
from os.path import basename
from pyworkflow.em.convert import ImageHandler

SET_OF_MOVIES=0
SET_OF_MICROGRAPHS=1
SET_OF_RANDOM_MICROGRAPHS=2

class ProtCreateStreamData(EMProtocol):
    """ create  setofXXXX in streaming mode.
        micrograph -> read a micrograph in memory and writes it nDim times
        movie      -> read a movie in memory and writes it nDim times
        randomMicrographs -> creates a micrograph with random values and aplies a reandom CTF
    """
    _label="create stream data"
    _lastUpdateVersion = VERSION_1_2
    _singleImageFn = "singleImage.xmp"
    _magnification = 500000
    _voltage = 200
    _sphericalAberration = 2.0
    _amplitudeContrast = 0.07
    object = None

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.dictObj = OrderedDict()
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        #TODO: add movie creation in streaming
        form.addParam('setof', params.EnumParam, default=0,
                      choices=['Movies', 'Micrographs', 'RandomMicrographs'],
                      label='set Of',
                      help='create set of')
        form.addParam('inputMovie', params.PointerParam, pointerClass='Movie',
                      condition="setof==%d"%SET_OF_MOVIES,
                      label="movie",
                      help='This movie will be copied "number of items" times')
        form.addParam('inputMic', params.PointerParam, pointerClass='Micrograph',
                      condition="setof==%d"%SET_OF_MICROGRAPHS,
                      label="micrograph",
                      help='This micrograph will be copied "number of items" times')
        form.addParam('xDim', params.IntParam,default=1024,
                      condition="setof==%d"%SET_OF_RANDOM_MICROGRAPHS,
                      label="xdim",
                      help="X dim ")
        form.addParam('yDim', params.IntParam,default=1024,
                      condition="setof==%d"%SET_OF_RANDOM_MICROGRAPHS,
                      label="ydim",
                      help="Y dim ")
        form.addParam('nDim', params.IntParam,default=10,
                      label="number of items",
                      help="setofXX size")
        form.addParam('samplingRate', params.FloatParam,default=4,
                      label="samplingRate",
                      help="Sampling rate")
        form.addParam('creationInterval', params.IntParam,default=60,
              label="Create Object each (sec)",
              pointerClass='EMProtocol',
              help="create one object each creationInterval seconds")
        form.addParam('delay', params.IntParam,default=0,
                      label="delay (sec)",
                      help="wait this seconds before creating stream data")

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        deps = []

        if self.delay != 0:
            delayId = self._insertFunctionStep('delayStep')
            deps.append(delayId)

        step = None
        for mic in range (1, self.nDim.get() +1):
            if self.setof == SET_OF_MOVIES:
                step = 'createStep'
            elif self.setof == SET_OF_MICROGRAPHS:
                step = 'createStep'
            elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
                step = 'createStep_random_mic'
            else:
                raise Exception('Unknown data type')
            self._insertFunctionStep(step, mic, prerequisites=deps)

    #--------------------------- STEPS functions --------------------------------------------
    def delayStep(self):
        time.sleep(self.delay)

    def _checkNewItems(self, objSet):
        """ Check for already computed micrograph/movie and update the output set. """
        objDict = {}
        newObj = False
        for obj in objSet:
            objDict[obj.getFileName()] = True

        if objDict:
            if objSet.getSize():
                objSet.enableAppend()
                objSet.loadAllProperties()
        else:
            objSet.setStreamState(objSet.STREAM_OPEN)
            acquisition = Acquisition()
            acquisition.setMagnification(self._magnification)
            acquisition.setVoltage(self._voltage)
            acquisition.setSphericalAberration(self._sphericalAberration)
            acquisition.setAmplitudeContrast(self._amplitudeContrast)
            objSet.setAcquisition(acquisition)
            objSet.setSamplingRate(self.samplingRate.get())


        if self.setof == SET_OF_MOVIES:
            obj = Movie()
        elif self.setof == SET_OF_MICROGRAPHS:
            obj = Micrograph()
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            obj = Micrograph()
        else:
            raise Exception('Unknown data type')

        counter = 0
        for k,v in self.dictObj.iteritems():
            counter += 1
            if (k not in objDict):
                obj.setFileName(k)
                obj.setMicName(basename(k))
                obj.setObjId(counter)
                objSet.append(obj)
                newObj= True

        return objSet, newObj#why a dictionary, a boolean may be enought

    def _updateOutput(self, objSet):
        if self.setof == SET_OF_MOVIES:
            self._defineOutputs(outputMovies=objSet)
        elif self.setof == SET_OF_MICROGRAPHS:
            self._defineOutputs(outputMicrographs=objSet)
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            self._defineOutputs(outputMicrographs=objSet)

    def _stepsCheck(self):
        if self.setof == SET_OF_MOVIES:
            objSet = SetOfMovies(filename=self._getPath('movies.sqlite'))
        elif self.setof == SET_OF_MICROGRAPHS:
            objSet = SetOfMicrographs(filename=self._getPath('micrographs.sqlite'))
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            objSet = SetOfMicrographs(filename=self._getPath('micrographs.sqlite'))
        else:
            raise Exception('Unknown data type')

        newObjSet, newObj = self._checkNewItems(objSet)

        #check if end ....
        endObjs = newObjSet.getSize() == self.nDim.get()

        if newObj:
            if endObjs:
                newObjSet.setStreamState(newObjSet.STREAM_CLOSED)
            self._updateOutput(newObjSet)
        newObjSet.close()

    def createStep(self, counter):

        if not ProtCreateStreamData.object:
            print("read object")
            if self.setof == SET_OF_MOVIES:
                ProtCreateStreamData.object = ImageHandler().read(self.inputMovie.get().getLocation())
                self.name = "movie"
            elif self.setof == SET_OF_MICROGRAPHS:
                ProtCreateStreamData.object = ImageHandler().read(self.inputMic.get().getLocation())
                self.name = "micro"

        #save file
        destFn = self._getExtraPath("%s_%05d" % (self.name,counter))
        ProtCreateStreamData.object.write(destFn)
        self.dictObj[destFn] = True
        print "self.dictObj", self.dictObj
        time.sleep(self.creationInterval.get())

    def createStep_random_mic(self, mic):
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
        self.dictObj[baseFnImageCTF] = True
        time.sleep(self.creationInterval.get())

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []


