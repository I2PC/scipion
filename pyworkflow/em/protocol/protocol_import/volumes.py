# -*- coding: utf-8 -*-
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
"""
In this module are protocol base classes related to EM imports of Micrographs,
Particles, Volumes...
"""

from os.path import exists, basename
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import copyFile
import pyworkflow.protocol.constants as const
import pyworkflow.protocol.params as params
from pyworkflow.em import Volume, ImageHandler, PdbFile
from pyworkflow.em.convert import downloadPdb
from pyworkflow.em.data import Transform
from base import ProtImportFiles
from images import ProtImportImages
from pyworkflow.em.utils.ccp4_utilities.convert import Ccp4Header


class ProtImportVolumes(ProtImportImages):
    """Protocol to import a set of volumes to the project"""
    _outputClassName = 'SetOfVolumes'
    _label = 'import volumes'

    def __init__(self, **args):
        ProtImportImages.__init__(self, **args)

    def _defineAcquisitionParams(self, form):
        """ Define acquisition parameters, it can be overriden
        by subclasses to change what parameters to include.
        """
        form.addParam('samplingRate', params.FloatParam,
                      label=Message.LABEL_SAMP_RATE)
        form.addParam('setDefaultOrigin', params.BooleanParam,
                      label="setDefaultOrigin",
                      help="Set origin of coordinates in the 3D map center "
                           "by default (True) or provide it. Only Modeling "
                           "related programs support this feature so far",
                      default=True)
        line = form.addLine('Offset',
                            help= "You have to provide the map center "
                            "coordinates in Angstroms (pixels x sampling). "
                            "We follow the same convention than CCP4. Chimera"
                            " considers the same magnitude and opposite sign "
                            "than CCP4.", condition='not setDefaultOrigin')
        line.addParam('x', params.FloatParam, condition='not setDefaultOrigin',
                      label="x", help="offset along x axis (A)")
        line.addParam('y', params.FloatParam, condition='not setDefaultOrigin',
                      label="y", help="offset along y axis (A)")
        line.addParam('z', params.FloatParam, condition='not setDefaultOrigin',
                      label="z", help="offset along z axis (A)")

    def _insertAllSteps(self):
        self._insertFunctionStep('importVolumesStep', self.getPattern(),
                                 self.samplingRate.get(),
                                 self.setDefaultOrigin.get())

    # --------------------------- STEPS functions -----------------------------

    def importVolumesStep(self, pattern, samplingRate, setDefaultOrigin=True):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        self.info("Using pattern: '%s'" % pattern)

        # Create a Volume template object
        vol = Volume()
        vol.setSamplingRate(samplingRate)
        copyOrLink = self.getCopyOrLink()
        imgh = ImageHandler()

        volSet = self._createSetOfVolumes()
        volSet.setSamplingRate(samplingRate)

        for fileName, fileId in self.iterFiles():
            dst = self._getExtraPath(basename(fileName))
            copyOrLink(fileName, dst)
            x, y, z, n = imgh.getDimensions(dst)
            # First case considers when reading mrc without volume flag
            # Second one considers single volumes (not in stack)
            if (z == 1 and n != 1) or (z != 1 and n == 1):
                vol.setObjId(fileId)
                if dst.endswith('.mrc'):
                    dst += ':mrc'
                vol.setLocation(dst)
                t = Transform()
                if setDefaultOrigin:
                    if (z == 1 and n != 1):
                        zDim = n
                    else:
                        zDim = z
                    t.setShifts(x/2. * samplingRate,
                                y/2. * samplingRate,
                                zDim/2. * samplingRate)
                else:
                    t.setShifts(self.x, self.y, self.z)
                vol.setOrigin(t)
                volSet.append(vol)
            else:
                for index in range(1, n+1):
                    vol.cleanObjId()
                    vol.setLocation(index, dst)
                    if setDefaultOrigin:
                        t = Transform()
                        if setDefaultOrigin:
                            t.setShifts(x / 2., y / 2., z / 2.)
                        else:
                            t.setShifts(self.x, self.y, self.z)
                        vol.setOrigin(t)
                    volSet.append(vol)
            # ##DELETE THIS
            #
            # ccp4header = Ccp4Header(vol.getFileName(), readHeader=True)
            # sampling = ccp4header.computeSampling()
            # print "origin.getShifts: ", vol.getOrigin().getShifts()
            # print "ccp4header.getStartAngstrom(sampling): ", ccp4header.getStartAngstrom(
            #     sampling)
            # print "ccp4header: ", ccp4header
            # ccp4header.setStartAngstrom(vol.getOrigin().getShifts(), sampling)
            # # ccp4header.writeHeader()
            # ccp4header.getStartAngstrom(sampling)
            # print "ccp4header.getStartAngstrom(sampling): ", ccp4header.getStartAngstrom(
            #         sampling)
            # print "ccp4header: ", ccp4header
            # ##

        if volSet.getSize() > 1:
            self._defineOutputs(outputVolumes=volSet)
        else:
            self._defineOutputs(outputVolume=vol)

    # --------------------------- INFO functions ------------------------------

    def _getVolMessage(self):
        if self.hasAttribute('outputVolume'):
            return "Volume %s" % self.getObjectTag('outputVolume')
        else:
            return "Volumes %s" % self.getObjectTag('outputVolumes')

    def _summary(self):
        summary = []
        if self.hasAttribute('outputVolume') or \
                self.hasAttribute('outputVolumes'):
            summary.append("%s imported from:\n%s" % (self._getVolMessage(),
                           self.getPattern()))

            summary.append(u"Sampling rate: *%0.2f* (Å/px)" %
                           self.samplingRate.get())
        return summary

    def _methods(self):
        methods = []
        if self.hasAttribute('outputVolume') or \
                self.hasAttribute('outputVolumes'):
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getVolMessage(), self.samplingRate.get()),)
        return methods


class ProtImportPdb(ProtImportFiles):
    """ Protocol to import a set of pdb volumes to the project"""
    _label = 'import atomic structure'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1

    def __init__(self, **args):
        ProtImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPdbData', params.EnumParam, choices=['id', 'file'],
                      label="Import atomic structure from",
                      default=self.IMPORT_FROM_ID,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import mmCIF data from online server or local file')
        form.addParam('pdbId', params.StringParam,
                      condition='inputPdbData == IMPORT_FROM_ID',
                      label="Atomic structure ID ", allowsNull=True,
                      help='Type a mmCIF ID (four alphanumeric characters).')
        form.addParam('pdbFile', params.PathParam, label="File path",
                      condition='inputPdbData == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to desired atomic structure.')
        form.addParam('inputVolume', params.PointerParam, label="Input Volume",
                      pointerClass='Volume',
                      allowsNull=True,
                      help='Associate this volume to the mmCIF file.')

    def _insertAllSteps(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            pdbPath = self._getPath('%s.cif' % self.pdbId.get())
            self._insertFunctionStep('pdbDownloadStep', pdbPath)
        else:
            pdbPath = self.pdbFile.get()
        self._insertFunctionStep('createOutputStep', pdbPath)

    def pdbDownloadStep(self, pdbPath):
        """Download all pdb files in file_list and unzip them."""
        downloadPdb(self.pdbId.get(), pdbPath, self._log)

    def createOutputStep(self, pdbPath):
        """ Copy the PDB structure and register the output object.
        """
        if not exists(pdbPath):
            raise Exception("Atomic structure not found at *%s*" % pdbPath)

        baseName = basename(pdbPath)
        localPath = self._getExtraPath(baseName)
        copyFile(pdbPath, localPath)
        pdb = PdbFile()
        volume = self.inputVolume.get()

        # if a volume exists assign it to the pdb object
        # IMPORTANT: we DO need to if volume is not None
        # because we need to persist the pdb object
        # before we can make the last source relation
        if volume is not None:
            pdb.setVolume(volume)

        pdb.setFileName(localPath)
        self._defineOutputs(outputPdb=pdb)

        if volume is not None:
            self._defineSourceRelation(volume, pdb)

    def _summary(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            summary = ['Atomic structure imported from ID: *%s*' %
                       self.pdbId]
        else:
            summary = ['Atomic structure imported from file: *%s*' %
                       self.pdbFile]

        return summary

    def _validate(self):
        errors = []
        if (self.inputPdbData == self.IMPORT_FROM_FILES and not exists(
                self.pdbFile.get())):
            errors.append("Atomic structure not found at *%s*" %
                          self.pdbPath.get())
        # TODO: maybe also validate that if exists is a valid PDB file
        return errors
