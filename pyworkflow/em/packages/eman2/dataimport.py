# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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

import pyworkflow.utils as pwutils
from pyworkflow.em.data import Coordinate, CTFModel
from pyworkflow.em.data_tiltpairs import Angles
from pyworkflow.em.metadata import (MetaData, MDL_XCOOR, MDL_YCOOR,
                                    MDL_PICKING_PARTICLE_SIZE)
from convert import loadJson, readCTFModel, readSetOfParticles


class EmanImport():

    def __init__(self, protocol, lstFile):
        self.protocol = protocol
        self._lstFile = lstFile
        self.copyOrLink = protocol.getCopyOrLink()

    def importAngles(self, fileName, addAngles):
        if pwutils.exists(fileName):
            ext = pwutils.getExt(fileName)

            if ext == ".json":
                fnBase = pwutils.replaceBaseExt(fileName, 'hdf')
                keyName = 'tiltparams_micrographs/' + fnBase.replace('_info', '')
                jsonAngDict = loadJson(fileName)
                if jsonAngDict.has_key(keyName):
                    angles = jsonAngDict[keyName]
                    tilt, y2, y = angles[:3]  # y2=tilted, y=gamma(untilted)
                    ang = Angles()
                    # TODO: check this conversion
                    ang.setAngles(y, y2, tilt)
                    addAngles(ang)
            else:
                raise Exception('Unknown extension "%s" to import Eman tilt pair angles' % ext)

    def importCoordinates(self, fileName, addCoordinate):
        if pwutils.exists(fileName):
            ext = pwutils.getExt(fileName)

            if ext == ".json":
                jsonPosDict = loadJson(fileName)
                boxes = []

                if jsonPosDict.has_key("boxes"):
                    boxes = jsonPosDict["boxes"]
                elif jsonPosDict.has_key("boxes_rct"):
                    boxes = jsonPosDict["boxes_rct"]
                if boxes:
                    for box in boxes:
                        x, y = box[:2]
                        coord = Coordinate()
                        coord.setPosition(x, y)
                        addCoordinate(coord)

            elif ext == ".box":
                md = MetaData()
                md.readPlain(fileName, "xcoor ycoor particleSize")
                size = md.getValue(MDL_PICKING_PARTICLE_SIZE, md.firstObject())
                if size is None:
                    print ">>> WARNING: Error parsing coordinate file: %s" % fileName
                    print "             Skipping this file."
                else:
                    half = size / 2
                    for objId in md:
                        x = md.getValue(MDL_XCOOR, objId)
                        y = md.getValue(MDL_YCOOR, objId)
                        coord = Coordinate()
                        coord.setPosition(x + half, y + half)
                        addCoordinate(coord)
            else:
                raise Exception('Unknown extension "%s" to import Eman coordinates' % ext)

    def getBoxSize(self, coordFile):
        """ Try to infer the box size from the given coordinate file.
        In the case of .box files, the size is the 3rd column
        In the case of .json files, we will look for file e2boxercache/base.json
        """
        if coordFile.endswith('.box'):
            md = MetaData()
            md.readPlain(coordFile, "xcoor ycoor particleSize")
            return md.getValue(MDL_PICKING_PARTICLE_SIZE, md.firstObject())

        elif coordFile.endswith('.json'):
            infoDir = pwutils.dirname(coordFile)
            # Still go one level up of info dir
            jsonBase = pwutils.join(pwutils.dirname(infoDir), 'e2boxercache', 'base.json')
            jsonBase2 = pwutils.join(infoDir, 'project.json')
            if pwutils.exists(jsonBase):
                jsonDict = loadJson(jsonBase)
                if jsonDict.has_key('box_size'):
                    return int(jsonDict["box_size"])
            elif pwutils.exists(jsonBase2):
                jsonDict = loadJson(jsonBase2)
                if jsonDict.has_key('global.boxsize'):
                    return int(jsonDict["global.boxsize"])

        return None

    def importCTF(self, mic, fileName):
        ctf = CTFModel()
        ctf.setMicrograph(mic)
        readCTFModel(ctf, fileName)
        return ctf

    def importParticles(self):
        """ Import particles from 'imageSet.lst' file. """
        partSet = self.protocol._createSetOfParticles()
        partSet.setObjComment('Particles imported from EMAN lst file:\n%s' % self._lstFile)
        self.protocol.setSamplingRate(partSet)
        partSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
        self.protocol.fillAcquisition(partSet.getAcquisition())

        # Now read the alignment/ctf/acquisition parameters
        direc = self.protocol._getExtraPath()
        readSetOfParticles(self._lstFile, partSet, direc)

        # Register the output set of particles
        self.protocol._defineOutputs(outputParticles=partSet)

    def validateParticles(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        errors = []
        return errors
