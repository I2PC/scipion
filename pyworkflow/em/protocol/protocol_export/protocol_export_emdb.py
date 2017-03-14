# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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


import re
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import EMProtocol
from pyworkflow import VERSION_1_2
from pyworkflow.em.convert import ImageHandler
import os


class ProtExportEMDB(EMProtocol):
    """ stress  will  stress  test  a  computer system in various selectable
       ways. Several options require the program stress-ng.
    """
    _label = 'exportEMDB'
    _program = "" 
    _version = VERSION_1_2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.xmippMic={}

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('exportVolume', params.PointerParam, label="Volume to export", important=True,
                      pointerClass='Volume',
                      help='This volume will be exported using mrc format')
        form.addParam('exportFSC', params.PointerParam, label="FSC to export", important=True,
                      pointerClass='FSC',
                      help='This FSC will be exported using XML emdb format')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
            self._insertFunctionStep('exportVolumeStep')
            self._insertFunctionStep('exportFSCStep')

    #--------------------------- STEPS functions --------------------------------------------

    def exportVolumeStep(self):

        self.nombreVolumen, self.nombreFsc = self.getName()

        ih = ImageHandler()
        cmdFile = self._getExtraPath(self.nombreVolumen)
        ih.convert(self.exportVolume.get().getLocation(), cmdFile)

    def exportFSCStep(self):

        x,y = self.exportFSC.get().getData()
        cmdFSCFile = self._getExtraPath(self.nombreFsc)
        fo = open(cmdFSCFile, "w")
        fo.write('<fsc title="FSC(%s)" xaxis="Resolution(A-1)" yaxis="Correlation Coefficient">\n'%self.nombreVolumen)
        for i in range(len(x)):
            fo.write("<coordinate>\n")
            fo.write("<x>%f</x>\n"%x[i])
            fo.write("<y>%f</y>\n" % y[i])
            fo.write("</coordinate>\n")

        fo.write("</fsc>\n")
        fo.close()


    def _validate(self):
        message = []
        return message

    def _summary(self):
        message = "Data Available at directory: *%s*"% os.path.abspath(self._getExtraPath())

        return [message]

    def _methods(self):
        return []

#--------------------------- UTILS functions ---------------------------------------------------

    def getName(self):
        name = self.exportVolume.get().getFileName()
        match = re.search(r'[\w.]+.mrc', name)
        if match:
            nombreVolumen = match.group()[:-4]
            nombreFsc = match.group()[:-8] + '.xml'
        else:
            nombreVolumen = "volume.mrc"
            nombreFsc = "fsc.xml"

        return nombreVolumen, nombreFsc
