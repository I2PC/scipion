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

import os

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em.data import Coordinate
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md



class ProtEthanPicker(ProtParticlePicking):
    """  """

    _label = 'ethan picker'

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):

        ProtParticlePicking._defineParams(self, form)
        form.addParam('radius', params.IntParam,
                      label='Radius of particle (px)')

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        deps = []

        for mic in self.getInputMicrographs():
            stepId = self._insertFunctionStep('pickMicrographStep',
                                              self.radius.get(),
                                              mic.getFileName())
            deps.append(stepId)

        self._insertFunctionStep('createOutputStep', prerequisites=deps)

    #--------------------------- STEPS functions -------------------------------
    def pickMicrographStep(self, radius, micFn):
        micDir = self._getMicDir(micFn)
        # Convert micrographs to mrc (uint8) as required by ETHAN program
        fnMicBase = pwutils.replaceBaseExt(micFn, 'mrc')
        fnMicFull = os.path.join(micDir, fnMicBase)
        fnPosBase = self._getMicPosFn(micFn)

        # FIXME: Right now we can not define the depth (8, 16 or 32 bits)
        # from the ImageHandler API
        from pyworkflow.em.packages.xmipp3 import runProgram
        args = '-i %s -o %s --depth uint8' % (micFn, fnMicFull)
        runProgram('xmipp_image_convert', args)

        # Run ethan program with the required arguments
        program = os.environ.get('ETHAN_BIN')
        args = "%d %s %s" % (radius, fnMicBase, fnPosBase)
        self.runJob(program, args, cwd=micDir)

        # Clean temporary micrograph
        pwutils.cleanPath(fnMicFull)

    def createOutputStep(self):
        coordSet = self._createSetOfCoordinates(self.getInputMicrographs())
        self.readSetOfCoordinates(self._getExtraPath(), coordSet)
        coordSet.setBoxSize(self.radius.get() * 2)
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    #--------------------------- UTILS functions -------------------------------
    def _getMicDir(self, micFn):
        return self._getExtraPath()

    def _getMicPosFn(self, micFn):
        return pwutils.replaceBaseExt(micFn, 'txt')

    def readSetOfCoordinates(self, workingDir, coordSet):

        for mic in self.getInputMicrographs():
            micFn = mic.getFileName()
            micDir = self._getMicDir(micFn)
            coordFile = os.path.join(micDir, self._getMicPosFn(micFn))
            if os.path.exists(coordFile):
                coordMd = md.MetaData()
                coordMd.readPlain(coordFile, 'xcoor ycoor')
                for objId in coordMd:
                    x = coordMd.getValue(md.MDL_XCOOR, objId)
                    y = coordMd.getValue(md.MDL_YCOOR, objId)
                    coord = Coordinate()
                    coord.setPosition(x, y)
                    coord.setMicrograph(mic)
                    coordSet.append(coord)
            else:
                print "Coordinate file '%s' not found. " % coordFile

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return []

    def _methods(self):
        methodsMsgs = []
        return methodsMsgs

