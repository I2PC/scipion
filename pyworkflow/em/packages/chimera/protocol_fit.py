# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

from protocol_operate import ChimeraProtOperate
from pyworkflow.protocol.params import MultiPointerParam, PointerParam, \
    StringParam


class ChimeraProtRigidFit(ChimeraProtOperate):
    """Protocol to perform rigid fit using Chimera.
        Execute command *scipionwrite [model #n] [refmodel #p]
        [saverefmodel 0|1]* from command line in order to transferm fitted
        pdb to scipion. Default values are model=#0,
        refmodel =#1 and saverefmodel 0 (false).
        model refers to the pdb file. refmodel to a 3Dmap"""
    _label = 'chimera rigid fit'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Volume to process")
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="PdbFile",
                      label='PDBx/mmCIF file to be refined',
                      help="PDBx/mmCIF file to be refined. This structure "
                           "object will be saved after refinement")
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="PdbFile",
                      label='Other reference PDBx/mmCIF files',
                      help="Other PDBx/mmCIF files used as reference. "
                           "This structure objects will not be saved")
        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="Add extra commands in cmd file. Use for testing")
        form.addSection(label='Help')
        form.addLine('''Execute command *scipionwrite [model #n] [refmodel
        #p] [saverefmodel 0|1]* from command line in order to transfer fitted 
        structure to scipion. Default values are model=#2,
        refmodel =#1 and saverefmodel 0 (false).
        Model refers to the structure file. refmodel to a 3Dmap''')

        # --------------------------- INSERT steps functions --------------------

    def prerequisitesStep(self):
        """
        """
        if self.inputVolume.get() is None:
            fnVol = self.pdbFileToBeRefined.get().getVolume()
            index, fn = fnVol.getLocation()
            print "Volume: Volume associated to atomic structure %s(%d)\n" \
                  % (fn, index)
        else:
            fnVol = self.inputVolume.get()
            print "Volume: Input volume %s\n" % fnVol

    def _validate(self):
        errors = super(ChimeraProtRigidFit, self)._validate()

        # Check that the input volume exist
        if (not self.pdbFileToBeRefined.get().hasVolume()) and (
                    self.inputVolume.get() is None):
            errors.append("Error: You should provide a volume.\n")

        return errors
