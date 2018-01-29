# **************************************************************************
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
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


from pyworkflow import VERSION_1_2
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam, EnumParam, \
    IntParam
from pyworkflow.em.constants import SYM_I222, SYM_I222r, SYM_In25, SYM_In25r, \
    SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL, SYM_OCTAHEDRAL
from pyworkflow.em.packages.xmipp3 import XMIPP_SYM_NAME
from pyworkflow.em.constants import SCIPION_SYM_NAME
from pyworkflow.em import Volume
from pyworkflow.em.packages.ccp4.convert import Ccp4Header
from pyworkflow.em.data import Transform

DEBUG = True


class XmippProtExtractUnit(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'extract unit cell'
    _program = ""
    _version = VERSION_1_2

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, label="Input Volume",
                      important=True, pointerClass='Volume',
                      help='This volume will be cropped')
        form.addParam('symmetryGroup', EnumParam,
                      choices=[XMIPP_SYM_NAME[SYM_CYCLIC] +
                               " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                               XMIPP_SYM_NAME[SYM_DIHEDRAL] +
                               " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL] + ")",
                               XMIPP_SYM_NAME[SYM_TETRAHEDRAL] +
                               " (" + SCIPION_SYM_NAME[SYM_TETRAHEDRAL] + ")",
                               XMIPP_SYM_NAME[SYM_OCTAHEDRAL] +
                               " (" + SCIPION_SYM_NAME[SYM_OCTAHEDRAL] + ")",
                               XMIPP_SYM_NAME[SYM_I222] +
                               " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                               XMIPP_SYM_NAME[SYM_I222r] +
                               " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                               XMIPP_SYM_NAME[SYM_In25] +
                               " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                               XMIPP_SYM_NAME[SYM_In25r] +
                               " (" + SCIPION_SYM_NAME[SYM_In25r] + ")"],
                      default=SYM_I222r,
                      label="Symmetry",
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/"
                           "Symmetry for a description of the symmetry groups "
                           "format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )
        form.addParam('symmetryOrder', IntParam, default=1,
                      condition='symmetryGroup<=%d' % SYM_DIHEDRAL,
                      label='Symmetry Order',
                      help='Order of cyclic symmetry.')
        form.addParam('offset', FloatParam, default=0.,
                      condition='symmetryGroup<=%d' % SYM_DIHEDRAL,
                      label="offset",
                      help="rotate unit cell around z-axis by offset degrees")
        form.addParam('innerRadius', FloatParam, default=-1,
                      label="Inner Radius (px)",
                      help="inner Mask radius, if -1, the radius will be 0")
        form.addParam('outerRadius', FloatParam, default=-1,
                      label="Outer Radius (px)",
                      help="outer Mask radius, if -1, the radius will be "
                           "volume_size/2")
        form.addParam('expandFactor', FloatParam, default=0.,
                      label="Expand Factor",
                      help="Increment cropped region by this factor")

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('extractUnit')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------

    def extractUnit(self):
        sym = self.symmetryGroup.get()
        if sym == SYM_CYCLIC:
            sym = "%s%d" % (XMIPP_SYM_NAME[SYM_CYCLIC][:1], self.symmetryOrder)
        elif sym == SYM_DIHEDRAL:
            sym = "%s%d" %\
                  (XMIPP_SYM_NAME[SYM_DIHEDRAL][:1], self.symmetryOrder)
        elif sym == SYM_TETRAHEDRAL:
            sym = "%s" % (XMIPP_SYM_NAME[SYM_TETRAHEDRAL])
        elif sym == SYM_OCTAHEDRAL:
            sym = "%s" % (XMIPP_SYM_NAME[SYM_OCTAHEDRAL])
        elif sym >= SYM_I222 and sym <= SYM_In25r:
            sym = XMIPP_SYM_NAME[self.symmetryGroup.get()]
        args = "-i %s -o %s" % \
               (self.inputVolumes.get().getFileName(), self._getOutputVol())
        args += " --unitcell %s " % sym
        args += " %f " % self.innerRadius.get()
        args += " %f " % self.outerRadius.get()
        args += " %f " % self.expandFactor.get()
        args += " %f " % self.offset.get()
        args += " %f " % self.inputVolumes.get().getSamplingRate()
        origin = self.inputVolumes.get().getOrigin().getShifts()
        # x origin coordinate
        args += " %f " % origin[0]
        # y origin coordinate
        args += " %f " % origin[1]
        # z origin coordinate
        args += " %f " % origin[2]

        self.runJob("xmipp_transform_window", args)

    def createOutputStep(self):
        vol = Volume()
        vol.setLocation(self._getOutputVol())
        sampling = self.inputVolumes.get().getSamplingRate()
        vol.setSamplingRate(sampling)
        #
        ccp4header = Ccp4Header(self._getOutputVol(), readHeader=True)
        t = Transform()
        x, y, z = ccp4header.getOffset()  # origin output vol coordinates
        #_inputVol = self.inputVolumes.get()
        origin = self.inputVolumes.get().getOrigin().getShifts()
        # x, y, z origin input vol coordinates
        x_origin = origin[0]
        y_origin = origin[1]
        z_origin = origin[2]
        # x, y, z origin output vol coordinates
        dim = self.inputVolumes.get().getDim()
        x += dim[0] / 2. - x_origin
        y += dim[1] / 2. - y_origin
        z += dim[2] / 2. - z_origin
        t.setShifts(-x, -y, -z)  # we follow chimera convention no MRC
        vol.setOrigin(t)
        #
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputVolumes, self.outputVolume)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        message = []
        return message

    def _summary(self):
        # message = "Data Available at : *%s*"% self.filesPath.get()
        message = ""
        return [message]

    def _methods(self):
        return []

    # --------------------------- UTILS functions -----------------------------

    def _getOutputVol(self):
        return self._getExtraPath("output_volume.mrc")

    def replace_at_index(self, tup, ix, val):
        return tup[:ix] + (val,) + tup[ix+1:]
