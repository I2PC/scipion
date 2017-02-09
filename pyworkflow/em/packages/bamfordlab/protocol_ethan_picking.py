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

        form.addParam('height', params.FloatParam, default=0.5,
                      label="Min. height",
                      help='The minimum height of the virus peak compared to '
                           'the average peak.')

        form.addParam('squareWidth', params.FloatParam, default=1.5,
                      label="Square width",
                      help='Minimum distance between two viruses is:\n'
                           '2 * SQUARE_WIDTH_PARAM * RADIUS.')

        form.addParam('ringWidth', params.FloatParam, default=1.2,
                      label="Ring width",
                      help='Width of the ring in ring filter is: \n'
                           'RING_WIDTH_PARAM * RADIUS - RADIUS.')

        form.addParam('dist', params.FloatParam, default=0.1,
                      label="Distance",
                      help='Distance which the peak can move in x or y '
                           'direction during center refinement is: \n'
                           'DIST_PARAM * RADIUS.')

        form.addParam('reduction', params.IntParam, default=0,
                      label="Reduction",
                      help='REDUCTION_PARAM * REDUCTION_PARAM pixels are '
                           'averaged before filtering. The value has to be '
                           'integer. Special value 0 means that the program '
                           'determines the parameter.')

        form.addParam('doRefine', params.BooleanParam, default=True,
                      label="Refine particle centers?",
                      help='')
        form.addParam('doSectorTest', params.BooleanParam, default=True,
                      label="Perform sector test?",
                      help='')
        form.addParam('doHeightTest', params.BooleanParam, default=True,
                      label="Perform height test?",
                      help='')
        form.addParam('doDistanceTest', params.BooleanParam, default=True,
                      label="Perform peak pair distance test?",
                      help='')

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
        fnMicCfg = pwutils.replaceBaseExt(micFn, 'cfg')
        fnMicFull = os.path.join(micDir, fnMicBase)
        fnPosBase = self._getMicPosFn(micFn)

        # FIXME: Right now we can not define the depth (8, 16 or 32 bits)
        # from the ImageHandler API
        from pyworkflow.em.packages.xmipp3 import runProgram
        args = '-i %s -o %s --depth uint8' % (micFn, fnMicFull)
        runProgram('xmipp_image_convert', args)

        # Create a configuration file to be used by ETHAN with the parameters
        # selected by the user
        self.writeConfigFile(os.path.join(micDir, fnMicCfg))
        # Run ethan program with the required arguments
        program = os.environ.get('ETHAN_BIN')
        args = "%d %s %s %s" % (radius, fnMicBase, fnPosBase, fnMicCfg)
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

    def writeConfigFile(self, configFn):
        f = open(configFn, 'w')
        argsDict = {
            'height': self.height.get(),
            'squareWidth': self.squareWidth.get(),
            'ringWidth': self.ringWidth.get(),
            'dist': self.dist.get(),
            'reduction': self.reduction.get(),
            'refinement': 1 if self.doRefine else 0,
            'sectorTest': 1 if self.doSectorTest else 0,
            'heightTest': 1 if self.doHeightTest else 0,
            'distanceTest': 1 if self.doDistanceTest else 0
        }

        f.write("""
# The minimum height of the virus peak compared to the average peak.
HEIGHT_PARAM %(height)s

# Minimum distance between two viruses is 2 * SQUARE_WIDTH_PARAM * RADIUS.
SQUARE_WIDTH_PARAM %(squareWidth)s

# Width of the ring in ring filter is RING_WIDTH_PARAM * RADIUS - RADIUS.
RING_WIDTH_PARAM %(ringWidth)s

# Distance which the peak can move in x or y direction during center
# refinement is DIST_PARAM * RADIUS.
DIST_PARAM %(dist)s

# REDUCTION_PARAM * REDUCTION_PARAM pixels are averaged before filtering.
# The value has to be integer
# Special value 0 means that the program determines the parameter.
REDUCTION_PARAM %(reduction)s

# Is refinement of particle centers performed? Put 1 for yes, 0 for no.
# Producing virus files is not possible when refinement is off
REFINEMENT_PARAM %(refinement)s

# Is sector test performed? Put 1 for yes, 0 for no.
SECTOR_TEST %(sectorTest)s

# Is height test performed? Put 1 for yes, 0 for no.
HEIGHT_TEST %(heightTest)s

# Is peak pair distance test performed? Put 1 for yes, 0 for no.
DISTANCE_TEST %(distanceTest)s
        """ % argsDict)
        f.close()

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return []

    def _methods(self):
        methodsMsgs = []
        return methodsMsgs

