# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.protocol.params as params

from pyworkflow import VERSION_1_2
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import FSC
from pyworkflow.em.convert.atom_struct import AtomicStructHandler, toCIF


class ProtExportEMDB(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'export emdb'
    _program = ""
    _lastUpdateVersion = VERSION_1_2
    VOLUMENAME = 'volume.mrc'
    COORDINATEFILENAME = 'coordinates.cif'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

        #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('exportVolume', params.PointerParam, label="Volume to export", important=True,
                      pointerClass='Volume',
                      help='This volume will be exported using mrc format')
        form.addParam('exportFSC', params.PointerParam, label="FSC to export", important=True,
                      pointerClass='FSC, SetOfFSCs',
                      help='This FSC will be exported using XML emdb format')
        form.addParam('exportAtomStruct', params.PointerParam,
                      label="Atomic structure to export", allowsNull=True,
                      pointerClass='AtomStruct',
                      help='This atomic structure will be exported using mmCIF format')
        form.addParam('filesPath', params.PathParam, important=True,
                      label="Export to directory",
                      help="Directory where the files will be generated.")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('exportVolumeStep')
        self._insertFunctionStep('exportFSCStep')
        self._insertFunctionStep('exportAtomStructStep')

    #--------------------------- STEPS functions --------------------------------------------

    def exportVolumeStep(self):
        #create directory if needed
        dirName = self.filesPath.get()
        try:
            os.makedirs(dirName)
        except OSError:
            if not os.path.isdir(dirName):
                raise
        ih = ImageHandler()
        ih.convert(self.exportVolume.get().getLocation(),
                   os.path.join(dirName, self.VOLUMENAME))

    def exportFSCStep(self):
        exportFSC = self.exportFSC.get()
        if isinstance(self.exportFSC.get(), FSC):
            fscSet = self._createSetOfFSCs()
            fscSet.append(exportFSC)
        else:
            fscSet = exportFSC

        dirName = self.filesPath.get()
        for i, exportFSC in enumerate(fscSet):

            x,y = exportFSC.getData()
            fnFSC = os.path.join(dirName, "fsc_%02d.xml" % i)
            fo = open(fnFSC, "w")
            fo.write('<fsc title="FSC(%s)" xaxis="Resolution (A-1)" '
                     'yaxis="Correlation Coefficient">\n' %
                     os.path.join(dirName, self.VOLUMENAME))
            for i in range(len(x)):
                fo.write("<coordinate>\n")
                fo.write("<x>%f</x>\n"%x[i])
                fo.write("<y>%f</y>\n" % y[i])
                fo.write("</coordinate>\n")

            fo.write("</fsc>\n")
            fo.close()

    def exportAtomStructStep(self):
        exportAtomStruct = self.exportAtomStruct.get()
        originStructPath = exportAtomStruct.getFileName()
        dirName = self.filesPath.get()
        destinyStructPath = os.path.join(dirName, self.COORDINATEFILENAME)
        if originStructPath.endswith(".cif") or originStructPath.endswith(".mmcif"):
            h = AtomicStructHandler()
            h.read(originStructPath)
            h.write(destinyStructPath)
        else:
            toCIF(originStructPath, destinyStructPath)


    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        message = []
        fnPath = self.filesPath.get()
        if fnPath == "" or fnPath is None:
            message.append("You must set a path to export.")
        return message

    def _summary(self):
        message = "Data Available at : *%s*"% self.filesPath.get()
        return [message]

    def _methods(self):
        return []

#--------------------------- UTILS functions ---------------------------------------------------

    def getFnPath(self, label='volume'):
        return os.path.join(self.filesPath.get(),
                            self._getFileName(label))
