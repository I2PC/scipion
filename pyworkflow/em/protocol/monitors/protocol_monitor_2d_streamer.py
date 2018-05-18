# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import os.path
import time
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils

import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor, Monitor


class ProtMonitor2dStreamer(ProtMonitor):
    """ This protocol will monitor an input set of particles
    (usually in streaming) and will run/schedule many copies
     of a given 2D classification protocol but using subsets
     of the input particles as the 2D classification input.
    """
    _label = '2d streamer'

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self._runIds = pwobj.CsvList(pType=int)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('input2dProtocol', params.PointerParam,
                      label="Input 2D protocol", important=True,
                      pointerClass='ProtClassify2D',
                      help="This protocol will serve as the template run"
                           "that will be repeated with subsets of the "
                           "input particles. ")

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Input particles",
                      help='Input particles that will be used to execute'
                           'many 2D classification runs based on the 2D '
                           'protocol template selected. ')

        form.addParam('batchSize', params.IntParam,
                      label="Batch size",
                      help="How many particles (approximately) you want to "
                           "group to make the new batch and launch a new 2d"
                           "classification job. ")

        form.addParam('startingNumber', params.IntParam, default=0,
                      label="Starting number",
                      help="Specify a value greater than 0 if you want to skip "
                           "this amount of particles from the classification "
                           "batches (e.g, if you have classified them for the "
                           "initial 2D classification template. ")

        group = form.addGroup('Monitoring')
        group.addParam('samplingInterval', params.IntParam, default=10,
                       label="Update interval (min)",
                       help="After how many minutes the protocol should look "
                            "for new input data and schedule more 2D classification"
                            "jobs if necessary. ")

    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    # --------------------------- STEPS functions ----------------------------
    def monitorStep(self):
        interval = self.samplingInterval.get() * 60
        # list of particles that will be inserted in the new set
        self._counter = 0
        self._lastMicId = None
        self._lastPartId = 0
        self._subset = self._createSubset()
        self._runPrerequisites = []
        self._streamClosed = False
        # list of runs that has been (or will) be scheduled/run

        finished = False

        while not finished:
            self._checkNewInput()
            time.sleep(interval)
            finished = self._streamClosed

    # -------------------------- UTILS functions ------------------------------
    def _createSubset(self):
        """ Create a new empty set of particles with a given suffix. """
        self._counter += 1
        subset = self._createSetOfParticles(suffix="_%03d" % self._counter)
        subset.copyInfo(self.inputParticles.get())

        return subset

    def _writeSubset(self, subset):
        """ Generated the output of this subset. """
        newSubsetName = 'outputParticles_%03d' % self._counter
        self.info("Creating new subset: %s" % newSubsetName)
        subset.write()
        self._defineOutputs(**{newSubsetName: subset})
        self._defineTransformRelation(self.inputParticles, subset)
        # The following is required to commit the changes to the database
        self._store(subset)
        subset.close()

        from pyworkflow.manager import Manager
        manager = Manager()
        project = manager.loadProject(self.getProject().getName())
        input2D = self.input2dProtocol.get()
        copyProt = project.copyProtocol(project.getProtocol(input2D.getObjId()))
        copyProt.inputParticles.set(project.getProtocol(self.getObjId()))
        copyProt.inputParticles.setExtended(newSubsetName)
        project.scheduleProtocol(copyProt, self._runPrerequisites)
        # Next schedule will be after this one
        self._runPrerequisites.append(copyProt.getObjId())

    def _checkNewInput(self):
        """ Check if there are new particles and generate a new set
        and its corresponding 2D classification. """
        self.info("Checking new input...")
        subset = self._subset

        for particle in self._iterParticles():
            micId = particle.getMicId()
            partId = particle.getObjId()
            subset.append(particle)
            self.info("micId: %03d, particle: %05s, size: %s"
                      % (micId, partId, subset.getSize()))

            # Check the following after finding particles of a new micrograph
            if micId != self._lastMicId:
                if self._lastMicId is not None and subset.getSize() > self.batchSize:
                    self._writeSubset(subset)
                    subset = self._createSubset()

                self._lastMicId = micId

            self._lastPartId = partId

        # Write last group of particles if input stream is closed
        if self._streamClosed:
            self._writeSubset(subset)

        self._subset = subset

    def _iterParticles(self):
        inputParts = self.inputParticles.get()
        inputParts.load()
        inputParts.loadAllProperties()
        self._streamClosed = inputParts.isStreamClosed()

        for p in inputParts.iterItems(orderBy=['_micId', 'id'],
                                      direction='ASC',
                                      where='id > %d' % self._lastPartId):
            yield p

        inputParts.close()
