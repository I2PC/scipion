# *********************************************************************
# * Authors:  Marta Martinez (mmmtnez@cnb.csic.es), May 2013
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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.packages.chimera.convert import runChimeraProgram, getProgram
from pyworkflow.em.protocol.protocol_import.volumes import ProtImportVolumes
import os
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer

class viewerProtImportVolumes(XmippViewer):
    """ Visualize the output of protocol volume strain """
    _label = 'viewer input volume'
    _targets = [ProtImportVolumes]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **kwargs):
        XmippViewer.__init__(self, **kwargs)


    # ROB: I know that there is a nice chimera interface
    # but it does not work in this case since I am interested
    # in reading the MRC header. So I will use chimera as
    # an external program

    def _visualize(self, obj, **args):
        # save temporal file
        outputVolume = self.protocol.outputVolume
        sampling = outputVolume.getSamplingRate()
        tmpFileNameBILD = self.protocol._getTmpPath("axis.bild")
        dim = outputVolume.getDim()[0]

        f = open(tmpFileNameBILD, "w")
        arrowDict = {}
        arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = dim /2
        arrowDict["r1"] = dim / 20
        arrowDict["r2"] = 4 * arrowDict["r1"]
        arrowDict["rho"] = 7.5 * arrowDict["r1"]


        f.write(".color 1 0 0\n"
                ".arrow 0 0 0 %(x)d 0 0 %(r1)f  \n"
                ".color 1 1 0\n"
                ".arrow 0 0 0 0 %(y)d 0 %(r1)f \n"
                ".color 0 0 1\n"
                ".arrow 0 0 0 0 0 %(z)d %(r1)f \n"%arrowDict)
        f.close()

        tmpFileNameCMD = self.protocol._getTmpPath("chimera.cmd")
        f = open(tmpFileNameCMD, "w")
        origin1 = self.protocol.outputVolume.getOrigin(
            returnInitIfNone=True).getShifts()
        x = origin1[0] * sampling
        y = origin1[1] * sampling
        z = origin1[2] * sampling

        f.write("volume #0 style surface level .001 origin %f,%f,%f\n" % (x, y, z))
        f.close()

        inputVol = outputVolume.getFileName().replace(':mrc', '')
        args = " "
        if os.path.exists(inputVol):
            args += inputVol + " "
        if os.path.exists(tmpFileNameBILD):
            args += tmpFileNameBILD + " "
        if os.path.exists(tmpFileNameCMD):
            args += tmpFileNameCMD + " "
        runChimeraProgram(getProgram(), args)
