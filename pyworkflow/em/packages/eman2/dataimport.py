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



from pyworkflow.em.packages.eman2.convert import readCoordinates

from os.path import  exists
from pyworkflow.em.data import Coordinate
from pyworkflow.em.packages.eman2 import loadJson

class EmanImport():

    def __init__(self, protocol):
        self.protocol = protocol

    def importCoordinates(self, fileName, addCoordinate):

        print 'import from eman ' + fileName
        if exists(fileName):
            jsonPosDict = loadJson(fileName)

            if jsonPosDict.has_key("boxes"):
                boxes = jsonPosDict["boxes"]

                for box in boxes:
                    x, y = box[:2]
                    coord = Coordinate()
                    coord.setPosition(x, y)
                    addCoordinate(coord)