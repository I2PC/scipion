# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
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

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
from pyworkflow.em.packages.chimera.convert import runChimeraProgram, \
    getProgram, chimeraPdbTemplateFileName
from protocol_fit import ChimeraProtRigidFit
import os

# TODO: very likely this should inherit from ProtocolViewer
# not from XmippViewer. But then I get an empty form :-(


class viewerChimeraProtRigidFit(XmippViewer):
    """ Visualize the output of protocol volume strain """
    _label = 'viewer fit'
    _targets = [ChimeraProtRigidFit]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

# ROB: I know that there is a nice chimera interface
# but it does not work in this case since I am interested
# in reading the MRC header. So I will use chimera as
# an external program

    def _visualize(self, obj, **args):
        # TODO if input volume is not mrc this will not work.
        inputVol = self.protocol.inputVolume.get().getFileName().\
            replace(':mrc', '')
        outputPDB = self.protocol._getExtraPath(chimeraPdbTemplateFileName) % 1
        args = " "
        if os.path.exists(outputPDB):
            args += outputPDB + " "
        if os.path.exists(inputVol):
            args += inputVol + " "
        runChimeraProgram(getProgram(), args)
