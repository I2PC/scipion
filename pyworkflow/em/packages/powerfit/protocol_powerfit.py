# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
#             Marta Martinez (mmmtnez@cnb.csic.es), Nov 2017
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

from pyworkflow.em import *
from pyworkflow.em.headers import adaptFileToCCP4, \
    ORIGIN
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils import *


class PowerfitProtRigidFit(ProtFitting3D):
    """ Protocol for fitting a PDB into a 3D volume

    This is actually a wrapper to the program Powerfit.
    See documentation at:
       http://www.bonvinlab.org/education/powerfit
    """
    _label = 'powerfit'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPDB', PointerParam, pointerClass='PdbFile',
                      label="Input PDBx/mmCIF", important=True)
        form.addParam('inputVol', PointerParam, pointerClass='Volume',
                      label="Input volume", important=True, allowsNull=True)
        form.addParam('resolution', FloatParam, default=6,
                      label="Resolution (A)", important=True,
                      help="Resolution for the fitting. "
                           "The PDBx/mmCIF is filtered to this frequency.")
        form.addParam('angleStep', FloatParam, label="Angular step",
                      default=10.0, expertLevel=LEVEL_ADVANCED,
                      help='Angular step for the alignment search')
        form.addParam('nModels', IntParam, label="Number of models",
                      default=10, expertLevel=LEVEL_ADVANCED,
                      help='Number of models to estimate')
        form.addParam('doLaplacian', BooleanParam, label="Apply Laplacian",
                      default=False, expertLevel=LEVEL_ADVANCED,
                      help='Apply a Laplacian to the volume to highlight '
                           'borders')
        form.addParam('doCoreWeight', BooleanParam, label="Apply core weight",
                      default=False, expertLevel=LEVEL_ADVANCED,
                      help='Apply core weights')
        form.addParam('otherPowerfit', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Other parameters for Powerfit',
                      help='See http://www.bonvinlab.org/education/powerfit')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('powerfitWrapper')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def powerfitWrapper(self):
        _localInputVol = "volume.mrc"
        localInputVol = self._getExtraPath(_localInputVol)
        if self.inputVol.get() is None:
            volume = self.inputPDB.get().getVolume()
            print "Volume: Volume associated to atomic structure %s\n" % volume
        else:
            volume = self.inputVol.get()
            print "Volume: Input volume %s\n" % volume
        sampling = volume.getSamplingRate()
        origin = volume.getOrigin(force=True).getShifts()

        # powerfit needs offset in origin
        adaptFileToCCP4(volume.getFileName(), localInputVol,
                        origin, sampling, ORIGIN)
        args = "%s %f %s -d %s -p %d -a %f -n %d" % (localInputVol,
                                                     self.resolution,
                                                     self.inputPDB.get().
                                                     getFileName(),
                                                     self._getExtraPath(),
                                                     self.numberOfThreads,
                                                     self.angleStep,
                                                     self.nModels)

        if self.doLaplacian:
            args += " -l"
        if self.doCoreWeight:
            args += " -cw"
        self.runJob("powerfit", args)

        # Construct the chimera viewers
        dim = volume.getDim()[0]
        bildFileName = os.path.abspath(self._getExtraPath("axis.bild"))
        createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)

        for n in range(1, self.nModels.get() + 1):
            fnPdb = self._getExtraPath("fit_%d.pdb" % n)
            if exists(fnPdb):
                fnCmd = self._getExtraPath("chimera_%d.cmd" % n)
                fhCmd = open(fnCmd, 'w')
                fhCmd.write("open %s\n" % _localInputVol)
                fhCmd.write("open lcc.mrc\n")
                fhCmd.write("open fit_%d.pdb\n" % n)
                fhCmd.write("vol #1 hide\n")
                fhCmd.write("scolor #0 volume #1 cmap rainbow\n")
                fhCmd.write("open %s\n" % bildFileName)
                fhCmd.close()

    def createOutputStep(self):
        if self.inputVol.get() is None:
            volume = self.inputPDB.get().getVolume()
        else:
            volume = self.inputVol.get()

        fnOutput = self._getExtraPath("solutions.out")
        qualifiers = {}
        if exists(fnOutput):
            lineCounter = 0
            for line in open(fnOutput, "r").readlines():
                if lineCounter > 0:
                    tokens = line.split()
                    fnPdb = self._getExtraPath("fit_%d.pdb" % int(tokens[0]))
                    qualifiers[fnPdb] = (Float(tokens[1]),
                                         Float(tokens[2]),
                                         Float(tokens[3]))
                lineCounter += 1

        setOfPDBs = self._createSetOfPDBs()

        for n in range(self.nModels.get()):
            fnPdb = self._getExtraPath("fit_%d.pdb" % (n + 1))
            if exists(fnPdb):
                pdb = PdbFile(fnPdb)
                pdb.setVolume(volume)
                pdb._powerfit_cc = qualifiers[fnPdb][0]
                pdb._powerfit_Fish_z = qualifiers[fnPdb][1]
                pdb._powerfit_rel_z = qualifiers[fnPdb][2]
                setOfPDBs.append(pdb)

        self._defineOutputs(outputPDBs=setOfPDBs)
        self._defineSourceRelation(self.inputPDB.get(), setOfPDBs)
        self._defineSourceRelation(volume, setOfPDBs)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        summary.append('Angular step: %f' % self.angleStep.get())
        if self.doLaplacian:
            summary.append("Apply Laplacian")
        if self.doLaplacian:
            summary.append("Apply core weights")
        if not self.otherPowerfit.empty():
            summary.append('Other powerfit parameters: %s' %
                           self.otherPowerfit)
        return summary

    def _methods(self):
        summary = []
        summary.append('We rigidly fitted the structure %s into the volume %s'
                       'with an angular step of %f using Powerfit '
                       '[vanZundert2015].'
                       % (self.inputPDB.get().getNameId(),
                          self.inputVol.get().getNameId(),
                          self.angleStep))
        if self.doLaplacian:
            summary.append(
                "We applied a Laplacian filter to the input volume.")
        if self.doLaplacian:
            summary.append("We used core weighted local cross-correlations.")
        return summary

    def _citations(self):
        return ['vanZundert2015']

    def _validate(self):
        errors = []
        if which('powerfit') is '':
            errors.append('You should have the program powerfit in the PATH')
        # Check that the input volume exist
        if (not self.inputPDB.get().hasVolume()) \
                and (self.inputVol.get() is None):
            errors.append("Error: You should provide a volume.\n")
        return errors
