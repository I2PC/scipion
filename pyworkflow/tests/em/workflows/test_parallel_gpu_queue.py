# **************************************************************************
# *
# * Authors:    David Maluenda Niubo (dmaluenda@cnb.csic.es)
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
import json  # to fill the queue form
import time

from pyworkflow.tests import *
from pyworkflow.em.protocol import *
from pyworkflow.protocol.launch import schedule
import pyworkflow.utils as pwutils

try:
    from relion.protocols import *
    from relion.convert import *
except:
    pwutils.pluginNotFound('relion')


# Set this to match with your queue system
QUEUE_PARAMS = (u'myslurmqueue', {u'JOB_TIME': u'1',  # in hours
                                  u'JOB_MEMORY': u'8192'}) # in Mb

class TestQueueBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='mda'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles, 
                                     filesPath=pattern,
                                     samplingRate=samplingRate,
                                     checkStack=checkStack)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles '
                            'is None.' % pattern)
        return protImport

    @classmethod
    def runNormalizeParticles(cls, particles):
        """ Run normalize particles protocol """
        protPreproc = cls.newProtocol(ProtRelionPreprocessParticles,
                                      doNormalize=True)
        protPreproc.inputParticles.set(particles)
        cls.launchProtocol(protPreproc)
        cls.sampling = protPreproc.outputParticles.getSamplingRate()
        return protPreproc

    def _checkAsserts(self, relionProt):
        def loadQueue():
            os.system("gpu squeue > queue.log")
            queue = open("queue.log", "r").readlines()
            os.system("rm queue.log")
            return queue[1:]  # the first line is the header

        def isJobInQueue(jobId, protId, Yes=True):
            """ returns the value of Yes if the protId and the jobId
                are found in the queue. The value of Yes must be a boolean
            """
            isQueueThere = not Yes
            queueList = loadQueue()
            for queueStr in queueList:
                if jobId in queueStr and str(protId) in queueStr:
                    isQueueThere = Yes
                    break
            return isQueueThere

        def wait_until(somepredicate, timeout, *args, **kwargs):
            mustend = time.time() + timeout
            while time.time() < mustend:
                if somepredicate(*args, **kwargs):
                    return True
                time.sleep(1)
            return False

        def checkQueue(relionProt):
            jobId = relionProt.getJobId()   # is an string
            protId = relionProt.getObjId()  # is an integer

            self.assertTrue(isJobInQueue(jobId, protId),
                            "The job %s corresponding to "
                            "the protocol %d has been not "
                            "attached to the system queue"
                             % (jobId, protId))

            print(pwutils.magentaStr(" > job %s of the protocol %d found in the "
                                     "queue, wait a sec..." % (jobId, protId)))
            isDone = wait_until(isJobInQueue, 10*60, jobId, protId, Yes=False)
            self.assertTrue(isDone, "Timeout: the job has not ended...")
            print(pwutils.magentaStr("    ...job ended!"))
            print("%s  >>  %s" % (relionProt, id(relionProt)))
            return relionProt

        if relionProt.useQueue():
            relionProt = checkQueue(relionProt)
            return  # I don't know why, but we cannot retrieve the output, permissions???

        print("%s  >>  %s" % (relionProt, id(relionProt)))
        self.assertIsNotNone(relionProt.outputClasses,
                             "There was a problem with Relion 2D classify")

        classsesPixSize = relionProt.outputClasses.getImages().getSamplingRate()
        self.assertAlmostEquals(self.sampling, classsesPixSize,
                                "There was a problem with the sampling rate "
                                "of the particles")
        for class2D in relionProt.outputClasses:
                self.assertTrue(class2D.hasAlignment2D())


    def _runRelionClassify2D(self, previousRun, label='', threads=1, MPI=1,
                             doGpu=False, GPUs='', useQueue=False):
        prot2D = self.newProtocol(ProtRelionClassify2D,
                                  doCTF=False, maskDiameterA=340,
                                  numberOfMpi=MPI, numberOfThreads=threads)
        prot2D.numberOfClasses.set(4)
        prot2D.numberOfIterations.set(3)
        prot2D.inputParticles.set(previousRun.outputParticles)
        prot2D.setObjLabel(label)

        if useQueue:
            prot2D._useQueue.set(True)
            prot2D._queueParams.set(json.dumps(QUEUE_PARAMS))

        prot2D.doGpu.set(doGpu)
        if doGpu:
            prot2D.gpusToUse.set(GPUs)

        self.launchProtocol(prot2D)
        return prot2D


class Test_Queue_ALL(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testNoGpuSerial(self):
        relionNoGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU serial",
                                                  useQueue=True)
        self._checkAsserts(relionNoGpu11)

    def testNoGpuMPI(self):
        relionNoGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI",
                                                  MPI=4, useQueue=True)
        self._checkAsserts(relionNoGpu14)

    def testNoGpuThreads(self):
        relionNoGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU Threads",
                                                  threads=4, useQueue=True)
        self._checkAsserts(relionNoGpu41)

    def testNoGpuMPIandThreads(self):
        relionNoGpu44 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI+Threads",
                                                  MPI=4, threads=4,
                                                  useQueue=True)
        self._checkAsserts(relionNoGpu44)

    def testGpuSerial(self):
        relionGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU serial",
                                                doGpu=True, useQueue=True)
        self._checkAsserts(relionGpu11)

    def testGpuMPI(self):
        relionGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=4,
                                                useQueue=True)
        self._checkAsserts(relionGpu14)

    def testGpuThreads(self):
        relionGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU Threads",
                                                doGpu=True, threads=4,
                                                useQueue=True)
        self._checkAsserts(relionGpu41)

    def testGpuMPIandThreads(self):
        relionGpu44 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI+Threads",
                                                doGpu=True, useQueue=True,
                                                MPI=4, threads=4)
        self._checkAsserts(relionGpu44)


class Test_Queue_Small(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testGpuMPI(self):
        relionGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=4,
                                                useQueue=True)
        self._checkAsserts(relionGpu14)

class Test_noQueue_ALL(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testNoGpuSerial(self):
        relionNoGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU serial")
        self._checkAsserts(relionNoGpu11)

    def testNoGpuMPI(self):
        relionNoGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI",
                                                  MPI=4)
        self._checkAsserts(relionNoGpu14)

    def testNoGpuThreads(self):
        relionNoGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU Threads",
                                                  threads=4)
        self._checkAsserts(relionNoGpu41)

    def testNoGpuMPIandThreads(self):
        relionNoGpu44 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI+Threads",
                                                  MPI=4, threads=4)
        self._checkAsserts(relionNoGpu44)

    def testGpuSerial(self):
        relionGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU serial",
                                                doGpu=True)
        self._checkAsserts(relionGpu11)

    def testGpuMPI(self):
        relionGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=4)
        self._checkAsserts(relionGpu14)

    def testGpuThreads(self):
        relionGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU Threads",
                                                doGpu=True, threads=4)
        self._checkAsserts(relionGpu41)

    def testGpuMPIandThreads(self):
        relionGpu44 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI+Threads",
                                                doGpu=True,
                                                MPI=4, threads=4)
        self._checkAsserts(relionGpu44)


class Test_noQueue_Small(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testGpuMPI(self):
        relionGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=4)
        self._checkAsserts(relionGpu14)
