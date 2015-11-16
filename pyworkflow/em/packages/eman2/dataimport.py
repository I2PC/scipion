# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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

from os.path import  exists, join, dirname

from pyworkflow.utils.path import getExt
from pyworkflow.em.data import Coordinate
from pyworkflow.em.packages.eman2 import loadJson
from pyworkflow.em.metadata import MetaData, MDL_XCOOR, MDL_YCOOR, MDL_PICKING_PARTICLE_SIZE


class EmanImport():

    def __init__(self, protocol):
        self.protocol = protocol

    def importCoordinates(self, fileName, addCoordinate):
        if exists(fileName):
            ext = getExt(fileName)
            
            if ext == ".json":
                jsonPosDict = loadJson(fileName)
                
                if jsonPosDict.has_key("boxes"):
                    boxes = jsonPosDict["boxes"]

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
                        coord.setPosition(x+half, y+half)
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
            infoDir = dirname(coordFile)
            # Still go one level up of info dir
            jsonBase = join(dirname(infoDir), 'e2boxercache', 'base.json')
            if exists(jsonBase):
                jsonDict = loadJson(jsonBase)
                if jsonDict.has_key('box_size'):
                    return int(jsonDict["box_size"])
                
        return None
                