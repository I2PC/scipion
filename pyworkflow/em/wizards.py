# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
""" This module is for Actual wizards to be able to be discovered from the
Domain. wizard.py is left for wizard models and base classes."""
import os, sys

from pyworkflow.em.convert.atom_struct import AtomicStructHandler
from pyworkflow.object import String
from pyworkflow.wizard import Wizard
from pyworkflow.em.wizard import EmWizard, ListTreeProviderString
from pyworkflow.gui import dialog
from pyworkflow.em import ProtImportImages, ProtImportCoordinates, \
    ProtImportCoordinatesPairs, ProtImportVolumes, Volume, Ccp4Header, \
    ProtImportSequence
import requests

class ImportAcquisitionWizard(EmWizard):
    _targets = [(ProtImportImages, ['acquisitionWizard'])]

    def show(self, form, *params):
        acquisitionInfo = form.protocol.loadAcquisitionInfo()

        if isinstance(acquisitionInfo, dict):
            # If acquisitionInfo is None means something is wrong.
            # Now, let's try to show a meaningful error message.
            self._setAcquisition(form, acquisitionInfo)
        else:
            # If not dict, it should be an error message
            dialog.showError("Input error", acquisitionInfo, form.root)

    @classmethod
    def _setAcquisition(cls, form, acquisitionInfo):
        """ Ask whether to set the AcquisitionInfo to the protocol parameters.
        Params:
            acquisitionInfo: Should be a dictionary with acquisition values.
                If None, show an error.
        """
        msg = ''
        for k, v in acquisitionInfo.iteritems():
            msg += '%s = %s\n' % (k, v)
        msg += '\n*Do you want to use detected acquisition values?*'
        response = dialog.askYesNo("Import acquisition",
                                   msg, form.root)
        if response:
            prot = form.protocol
            comment = ''

            for k, v in acquisitionInfo.iteritems():
                if prot.hasAttribute(k):
                    form.setVar(k, v)
                else:
                    comment += "%s = %s\n" % (k, v)
            if comment:
                prot.setObjComment(comment)


class ImportCoordinatesBoxSizeWizard(Wizard):
    _targets = [(ProtImportCoordinates, ['boxSize']),
                (ProtImportCoordinatesPairs, ['boxSize'])]

    @classmethod
    def _getBoxSize(cls, protocol):

        return protocol.getDefaultBoxSize()

    @classmethod
    def show(cls, form, *params):
        form.setVar('boxSize', cls._getBoxSize(form.protocol))


class ImportOriginVolumeWizard(Wizard):

    _targets = [(ProtImportVolumes, ['x', 'y', 'z'])]

    def show(self, form, *params):
        protocol = form.protocol
        filesPath = protocol.filesPath.get()
        filesPattern = protocol.filesPattern.get()
        if filesPattern:
            fullPattern = os.path.join(filesPath, filesPattern)
        else:
            fullPattern = filesPath

        sampling = protocol.samplingRate.get()
        for fileName, fileId in protocol.iterFiles():
            inputVol = Volume()
            inputVol.setFileName(fileName)
            if ((str(fullPattern)).endswith('mrc') or
               (str(fullPattern)).endswith('map')):
                ccp4header = Ccp4Header(fileName, readHeader=True)
                x, y, z = ccp4header.getOrigin(changeSign=True)  # In Angstroms
            else:
                x, y, z = self._halfOriginCoordinates(inputVol, sampling)

            form.setVar('x', x)
            form.setVar('y', y)
            form.setVar('z', z)

    @classmethod
    def _halfOriginCoordinates(cls, volume, sampling):
        xdim, ydim, zdim = volume.getDim()
        if zdim > 1:
            zdim = zdim / 2.
        x = xdim / 2. * sampling
        y = ydim / 2. * sampling
        z = zdim * sampling
        return x, y, z


class GetStructureChainsWizard(Wizard):
    _targets = [(ProtImportSequence, ['inputStructureChain'])
                # NOTE: be careful if you change this class since
                # chimera-wizard inherits from it.
                # (ChimeraModelFromTemplate, ['inputStructureChain'])
                # (atomstructutils, ['inputStructureChain'])
                ]

    @classmethod
    def getModelsChainsStep(cls, protocol):
        structureHandler = AtomicStructHandler()
        fileName = ""
        if hasattr(protocol, 'pdbId'):
            if protocol.pdbId.get() is not None:
                pdbID = protocol.pdbId.get()
                url = "https://www.rcsb.org/structure/"
                URL = url + ("%s" % pdbID)
                try:
                    response = requests.get(URL)
                except:
                    raise Exception("Cannot connect to PDB server")
                if ((response.status_code >= 400) and (response.status_code < 500)):
                    raise Exception("%s is a wrong PDB ID" % pdbID)
                fileName = structureHandler.readFromPDBDatabase(
                    os.path.basename(pdbID), dir="/tmp/")
            else:
                fileName = protocol.pdbFile.get()
        else:
            if protocol.pdbFileToBeRefined.get() is not None:
                fileName = os.path.abspath(protocol.pdbFileToBeRefined.get(
                ).getFileName())

        structureHandler.read(fileName)
        structureHandler.getStructure()
        models = structureHandler.getModelsChains()
        return models

    def editionListOfChains(self, models):
        self.chainList = []
        for model, chainDic in models.iteritems():
            for chainID, lenResidues in chainDic.iteritems():

                self.chainList.append(
                    '{"model": %d, "chain": "%s", "residues": %d}' %
                    (model, str(chainID), lenResidues))

    def show(self, form, *params):
        protocol = form.protocol
        try:
            models = self.getModelsChainsStep(protocol)
        except Exception as e:
            print "ERROR: ", e.message
            return

        self.editionListOfChains(models)
        finalChainList = []
        for i in self.chainList:
            finalChainList.append(String(i))
        provider = ListTreeProviderString(finalChainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                                "Select one of the chains (model, chain, "
                                "number of chain residues)")
        form.setVar('inputStructureChain', dlg.values[0].get())
